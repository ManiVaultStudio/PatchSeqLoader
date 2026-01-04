#include "PatchSeqDataLoader.h"

#include "InputDialog.h"

#include "MatrixDataLoader.h"
#include "MatrixData.h"

#include "EphysData/Experiment.h"
#include "Electrophysiology/NWBLoader.h"

#include "Taxonomy.h"
#include "Morphology/SWCLoader.h"
#include "FeatureNames.h"

#include <util/Timer.h>

#include "csv.h"

#include <CellMorphologyData/CellMorphology.h>

#include <Set.h>

#include <QtCore>
#include <QtDebug>
#include <QFileDialog>
#include <QDir>

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <map>

Q_PLUGIN_METADATA(IID "studio.manivault.PatchSeqDataLoader")

//#define DALLEYLEE
#define KALMBACH
//#define WALEBOER

#ifdef DALLEYLEE

#define CELL_ID_TAG "cell_id"
#define CELL_NAME_TAG "cell_name"

#define METADATA_CLUSTER_LABEL "cluster_label_Hierarchical"
#define METADATA_SUBCLASS_LABEL "subclass_label_Hierarchical"

#define TX_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_Original/IDs_w_tc_data.csv"
#define EPHYS_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_November_2025/20251107_ephys_data.csv"
#define MORPHO_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_July_Cortex/Morphology/dendrite_morphometric_features.csv"
#define META_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_July_Cortex/Dalley_simplified_metadata_4_4_2025.csv"
#define MORPHO_DIR "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_August_2025/SWC_LayerAligned"
#define TRACES_DIR "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_July_Cortex/Electrophysiology/dalley_cytosplore_nwb"

#define EPHYS_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_November_2025/e_umap.csv"
#define MORPHO_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_November_2025/m_umap.csv"
#define TX_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_July_Cortex/Transcriptomics/tx_umap.csv"
#define ME_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_November_2025/me_umap.csv"

#endif

#ifdef KALMBACH

#define CELL_ID_TAG "cell_id"
#define CELL_NAME_TAG "cell_name"

#define METADATA_CLUSTER_LABEL "Group_name"
#define METADATA_SUBCLASS_LABEL "Subclass_name"

#define TX_PATH ""
#define EPHYS_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_June_Macaque/NHP_ephys_features_20250520.csv"
#define MORPHO_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_June_Macaque/RawFeatureWide_dend_20250616.csv"
#define META_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_December_Macaque/HANN_filtered_Cytosplore_20251212.csv"
#define MORPHO_DIR "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_December_Macaque/SWC_files"
#define TRACES_DIR "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_December_Macaque/NWB_files"

#define EPHYS_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_June_Macaque/ephys_umap_s.csv"
#define MORPHO_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_June_Macaque/morpho_umap_s.csv"
#define TX_UMAP_PATH ""

#endif

#ifdef WALEBOER

#define CELL_ID_TAG "cell_id"
#define CELL_NAME_TAG "cell_name"

#define METADATA_CLUSTER_LABEL "T-type"
#define METADATA_SUBCLASS_LABEL "Subclass"

#define TX_PATH ""
#define EPHYS_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_VU/Data Julian/Processed_data/ephys_data.csv"
#define MORPHO_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_VU/Data Julian/Processed_data/morph_data.csv"
#define META_PATH "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_VU/Data Julian/Processed_data/meta_data.csv"
#define MORPHO_DIR "D:/Dropbox/Julian/Patchseq/Supplied_Data/Data_VU/Data Julian/Processed_data/SWC"
#define TRACES_DIR ""

#define EPHYS_UMAP_PATH ""
#define MORPHO_UMAP_PATH ""
#define TX_UMAP_PATH ""
#define ME_UMAP_PATH ""

#endif

using namespace mv;
using namespace mv::gui;

namespace
{
    void addPointsToTextDataset(Dataset<Points>& points, Dataset<Text>& text, std::vector<uint32_t>& indexMapping)
    {
        for (int d = 0; d < points->getNumDimensions(); d++)
        {
            std::vector<QString> column(text->getNumRows(), "?");
            QString dimName = points->getDimensionNames()[d];

            // Get a proper name if possible
            QString properHeaderName = dimName;
            if (properFeatureNames.find(dimName) != properFeatureNames.end())
                properHeaderName = properFeatureNames[dimName];

            for (size_t i = 0; i < indexMapping.size(); i++)
            {
                uint32_t row = indexMapping[i];
                column[row] = QString::number(points->getValueAt(i * points->getNumDimensions() + d));
            }
            text->addColumn(properHeaderName, column);
        }
    }

    void removeRowsNotInMetadata(DataFrame& df, QString columnToCheck, DataFrame& metadata, MatrixData& matrix)
    {
        std::vector<QString> metaColumn = metadata[columnToCheck];
        // Make a unique set out of meta column
        std::unordered_set<QString> uniqueValues;
        for (int i = 0; i < metaColumn.size(); i++)
        {
            uniqueValues.insert(metaColumn[i]);
        }

        // For every value in df column, check if its in metadata
        std::vector<QString> dfColumn = df[columnToCheck];
        std::vector<int> rowsToRemove;
        for (int i = 0; i < dfColumn.size(); i++)
        {
            QString row = dfColumn[i];
            if (uniqueValues.find(row) == uniqueValues.end())
            {
                // Remove row because its not in the metadata
                rowsToRemove.push_back(i);
            }
        }
        qDebug() << "Removing " << rowsToRemove.size() << " rows because they are not found in metadata";

        df.removeRows(rowsToRemove);
        matrix.removeRows(rowsToRemove);
    }

    void removeFeatures(MatrixData& matrix, QStringList features)
    {
        std::vector<int> colsToDelete;
        for (int i = 0; i < matrix.headers.size(); i++)
        {
            for (int j = 0; j < features.size(); j++)
            {
                if (matrix.headers[i] == features[j])
                    colsToDelete.push_back(i);
            }
        }

        matrix.removeCols(colsToDelete);
    }

    std::map<QString, std::vector<unsigned int>> makeClustersFromList(std::vector<QString> list)
    {
        std::map<QString, std::vector<unsigned int>> clusterData;

        for (int i = 0; i < list.size(); i++)
        {
            const QString& str = list[i];
            clusterData[str].push_back(i);
        }
        return clusterData;
    }

    void removeRowsWithAllDataMissing(DataFrame& df, MatrixData& matrix)
    {
        int numRows = matrix.numRows;

        // Identify rows with all missing values
        std::vector<int> badRowIndices;
        for (int row = 0; row < numRows; row++)
        {
            bool badRow = true;
            for (int col = 0; col < matrix.numCols; col++)
            {
                if (matrix.data[row * matrix.numCols + col] != MISSING_VALUE)
                {
                    badRow = false;
                }
            }

            if (badRow)
            {
                badRowIndices.push_back(row);
            }
        }

        df.removeRows(badRowIndices);
        matrix.removeRows(badRowIndices);
    }

    void removeDuplicateRows(DataFrame& df, QString columnToCheck, MatrixData& matrix)
    {
        std::vector<int> duplicateRows = df.findDuplicateRows(columnToCheck);
        qDebug() << "Removing duplicate rows: " << duplicateRows.size();
        for (int i = 0; i < duplicateRows.size(); i++)
        {
            qDebug() << df[columnToCheck][duplicateRows[i]];
        }

        df.removeRows(duplicateRows);
        matrix.removeRows(duplicateRows);
    }

    QColor hexToQColor(const QString& hex)
    {
        return QColor(hex.left(7));
    }

    void buildMapOfCellTypesColors(const DataFrame& taxonomyDf, QHash<QString, QColor>& cellTypeColors)
    {
        std::vector<QString> superTypeNames = taxonomyDf["labels"];
        std::vector<QString> superTypeColors = taxonomyDf["colors"];

        for (int i = 0; i < superTypeNames.size(); i++)
        {
            // Colors are in RGBA hex, but Qt wants ARGB, so we need to swizzle, actually just remove alpha component, make it RGB
            cellTypeColors[superTypeNames[i]] = hexToQColor(superTypeColors[i]); // #RRGGBB
        }
    }

    void addColorizedClusters(Dataset<Points> umapDataset, DataFrame& umapDf, QString columnName)
    {
        // Create a list of clusters and their indices from the list of cluster names
        Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", columnName, umapDataset);

        const std::vector<QString>& clusterAsList = umapDf[columnName];
        std::map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(clusterAsList);

        for (auto& kv : clusterMap)
        {
            Cluster cluster;

            cluster.setName(kv.first);
            cluster.setIndices(kv.second);

            clusterData->addCluster(cluster);
        }
        Cluster::colorizeClusters(clusterData->getClusters(), 0);

        mv::events().notifyDatasetDataChanged(clusterData);
        mv::events().notifyDatasetDataDimensionsChanged(clusterData);
    }

    void addColorizedClustersFromMetadata(Dataset<Text>& metadata, Dataset<Points>& umap, KeyBasedSelectionGroup& selectionGroup, QString columnName, QHash<QString, QColor> cellTypeColors = QHash<QString, QColor>())
    {
        const BiMap& umapBiMap = selectionGroup.getBiMap(umap);
        const BiMap& metadataBiMap = selectionGroup.getBiMap(metadata);

        std::vector<uint32_t> allVals(umap->getNumPoints());
        std::iota(allVals.begin(), allVals.end(), 0);
        std::vector<QString> keys = umapBiMap.getKeysByValues(allVals);
        std::vector<int> metaIndices = metadataBiMap.getValuesByKeysWithMissingValue(keys, -1);

        // Create a list of clusters and their indices from the list of cluster names
        Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", columnName, umap);

        const std::vector<QString>& column = metadata->getColumn(columnName);
        std::vector<QString> clusterKeys;
        for (int idx : metaIndices)
        {
            clusterKeys.push_back(column[idx]);
        }

        std::map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(clusterKeys);

        for (auto& kv : clusterMap)
        {
            //if (kv.first == "")
            //    continue;

            Cluster cluster;

            cluster.setName(kv.first);
            cluster.setIndices(kv.second);

            clusterData->addCluster(cluster);
        }
        Cluster::colorizeClusters(clusterData->getClusters(), 0);

        // Assign unassigned indices to unknown cluster
        std::unordered_set<uint32_t> unassignedIndicesSet(allVals.begin(), allVals.end());
        for (const Cluster& cluster : clusterData->getClusters())
        {
            for (const auto& value : cluster.getIndices()) {
                unassignedIndicesSet.erase(value);
            }
        }
        //if (columnName != "lobe")
        {
            Cluster cluster;

            if (columnName != "lobe")
                cluster.setName("Unlabeled");
            else
                cluster.setName("TemL (ref)");

            std::vector<uint32_t> unassignedIndices(unassignedIndicesSet.begin(), unassignedIndicesSet.end());
            cluster.setIndices(unassignedIndices);
            cluster.setColor(hexToQColor("#DDDDDD"));

            clusterData->addCluster(cluster);
        }
        //}
        //else
        //{
        //    for (Cluster& cluster : clusterData->getClusters())
        //    {
        //        if (cluster.getName() == "TemL")
        //        {
        //            std::vector<uint32_t> unassignedIndices(unassignedIndicesSet.begin(), unassignedIndicesSet.end());
        //            std::vector<uint32_t> combinedIndices;
        //            combinedIndices.insert(combinedIndices.end(), cluster.getIndices().begin(), cluster.getIndices().end());
        //            combinedIndices.insert(combinedIndices.end(), unassignedIndices.begin(), unassignedIndices.end());
        //            cluster.setIndices(combinedIndices);
        //        }
        //    }
        //}

        if (columnName == "lobe")
        {
            for (Cluster& cluster : clusterData->getClusters())
            {
                if (cluster.getName() == "CeG") cluster.setColor(hexToQColor("#E6B6B6"));
                if (cluster.getName() == "FroL") cluster.setColor(hexToQColor("#E8CD59"));
                if (cluster.getName() == "ParL") cluster.setColor(hexToQColor("#EC813B"));
                if (cluster.getName() == "TemL") cluster.setColor(hexToQColor("#D55C92"));
                if (cluster.getName() == "OccL") cluster.setColor(hexToQColor("#D24D48"));
                if (cluster.getName() == "InL") cluster.setColor(hexToQColor("#B94A39"));
                if (cluster.getName() == "LimL") cluster.setColor(hexToQColor("#8D7BDB"));
                if (cluster.getName() == "") { cluster.setName("Unknown"); cluster.setColor(hexToQColor("#4C72B0")); }
                if (cluster.getName() == "TemL (ref)") cluster.setColor(QColor(220, 220, 220, 128));
            }
        }
        if (columnName == "paradigm")
        {
            for (Cluster& cluster : clusterData->getClusters())
            {
                if (cluster.getName() == "acute") cluster.setColor(hexToQColor("#DBCD70"));
                if (cluster.getName() == "culture") cluster.setColor(hexToQColor("#70BEDB"));
            }
        }
#ifdef DALLEYLEE
        if (columnName == "Supertype" || columnName == "Subclass")
        {
            for (Cluster& cluster : clusterData->getClusters())
            {
                if (cellTypeColors.contains(cluster.getName()))
                    cluster.setColor(cellTypeColors[cluster.getName()]);
            }
        }
#endif

        mv::events().notifyDatasetDataChanged(clusterData);
        mv::events().notifyDatasetDataDimensionsChanged(clusterData);
    }
}

// =============================================================================
// View
// =============================================================================

PatchSeqDataLoader::~PatchSeqDataLoader(void)
{

}

void PatchSeqDataLoader::init()
{

}

int extractLayerNumber(const QString& str)
{
    QRegularExpression re("L(\\d+)");
    QRegularExpressionMatch match = re.match(str);
    if (match.hasMatch()) {
        return match.captured(1).toInt();
    }
    return -1; // Return -1 if no match is found, for safety
}

void PatchSeqDataLoader::addTaxonomyClustersForDf(DataFrame& df, DataFrame& metadata, TaxonomyLevel level, QString name, mv::Dataset<mv::DatasetImpl> parent, QString metaLabel)
{
    // Get the available cluster labels from the metadata df
    const std::vector<QString>& treeCluster = metadata[metaLabel];

    Dataset<Clusters> treeClusterData = mv::data().createDataset<Clusters>("Cluster", properFeatureNames[metaLabel], parent);

    //std::unordered_map<QString, std::vector<unsigned int>> clusterData;

    const std::vector<QString>& metadataCellIds = metadata[CELL_ID_TAG];
    const std::vector<QString>& dataCellIds = df[CELL_ID_TAG];
    std::vector<QString> clusterNames(dataCellIds.size());

    // FIXME Do this more efficiently
    for (int i = 0; i < dataCellIds.size(); i++)
    {
        const QString& cellId = dataCellIds[i];

        bool foundMatch = false;
        for (int j = 0; j < metadataCellIds.size(); j++)
        {
            const QString& metaCellId = metadataCellIds[j];

            if (cellId == metaCellId)
            {
                clusterNames[i] = treeCluster.at(j);
                foundMatch = true;
            }
        }
        if (!foundMatch)
        {
            clusterNames[i] = "Undefined";
        }
    }

    // Create a list of clusters and their indices from the list of cluster names
    std::map<QString, std::vector<unsigned int>> clusterData = makeClustersFromList(clusterNames);

//#ifdef DALLEYLEE
//    // Try to sort cluster names according to layer
//    std::sort(clusterNames.begin(), clusterNames.end(), [](const QString& a, const QString& b) {
//        return extractLayerNumber(a) < extractLayerNumber(b);
//        });
//    qDebug() << clusterNames;
//#endif

    for (auto& kv : clusterData)
    {
        if (kv.first.isEmpty())
        {
            qWarning() << "Skipping cluster with no name..";
            continue;
        }

        Cluster cluster;

        cluster.setName(kv.first);

#if defined(DALLEYLEE) || defined(WALEBOER)
        if (_cellTypeColors.contains(cluster.getName()))
            cluster.setColor(_cellTypeColors[cluster.getName()]);
#else
        if (level == TaxonomyLevel::GROUP)
        {
            const QColor& color = _colorTaxonomy.getGroupColor(cluster.getName());
            cluster.setColor(color);
        }
        if (level == TaxonomyLevel::SUBCLASS)
        {
            const QColor& color = _colorTaxonomy.getSubclassColor(cluster.getName());
            cluster.setColor(color);
        }
#endif

        cluster.setIndices(kv.second);

        treeClusterData->addCluster(cluster);
    }
    
    mv::events().notifyDatasetDataChanged(treeClusterData);
    mv::events().notifyDatasetDataDimensionsChanged(treeClusterData);
}

void PatchSeqDataLoader::createClusterData(std::vector<QString> stringList, QString dataName, mv::Dataset<mv::DatasetImpl> parent)
{
    Dataset<Clusters> clusterData = mv::data().createDataset<Points>("Cluster", dataName, parent);

    std::map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(stringList);

    for (auto& kv : clusterMap)
    {
        if (kv.first.isEmpty())
        {
            qWarning() << "[PatchSeqDataLoader::createClusterData] Skipping cluster with empty name..";
            continue;
        }

        Cluster cluster;

        cluster.setName(kv.first);

        cluster.setIndices(kv.second);

        clusterData->addCluster(cluster);
    }

    mv::events().notifyDatasetDataChanged(clusterData);
    mv::events().notifyDatasetDataDimensionsChanged(clusterData);
}

void PatchSeqDataLoader::loadData()
{
    Q_INIT_RESOURCE(met_loader_resources);

    InputDialog inputDialog(nullptr, *this);
    inputDialog.setModal(true);

    // open dialog and wait for user input
    int ok = inputDialog.exec();

    PatchSeqFilePaths filePaths;
    QDir morphologiesDir;
    QDir ephysTracesDir;
    if (ok == QDialog::Accepted)
    {
        //filePaths.gexprFilePath = "D:/Dropbox/Julian/Patchseq/ProvidedData/IDs_w_tc_data.csv";// inputDialog.getTranscriptomicsFilePath();
        //filePaths.ephysFilePath = "D:/Dropbox/Julian/Patchseq/ProvidedData/allen_test_human_exc_simple_ephys.csv";// "D:/Dropbox/Julian/Patchseq/NewData/240928_human_exc_dataset_rsc369_ephys_data.csv";// inputDialog.getElectrophysiologyFilePath();
        //filePaths.morphoFilePath = "D:/Dropbox/Julian/Patchseq/ProvidedData/allen_test_human_exc_simple_morpho.csv";// inputDialog.getMorphologyFilePath();
        //filePaths.metadataFilePath = "D:/Dropbox/Julian/Patchseq/ProvidedData/allen_test_human_exc_metadata_simple.csv";// inputDialog.getMetadataFilePath();
        //morphologiesDir = QDir("D:/Dropbox/Julian/Patchseq/ProvidedData/SWC_Upright"); //inputDialog.getMorphologiesDir();

        filePaths.gexprFilePath = TX_PATH;// inputDialog.getTranscriptomicsFilePath();
        filePaths.ephysFilePath = EPHYS_PATH;// "D:/Dropbox/Julian/Patchseq/NewData/240928_human_exc_dataset_rsc369_ephys_data.csv";// inputDialog.getElectrophysiologyFilePath();
        filePaths.morphoFilePath = MORPHO_PATH;// inputDialog.getMorphologyFilePath();
        filePaths.metadataFilePath = META_PATH;// inputDialog.getMetadataFilePath();
        morphologiesDir = QDir(MORPHO_DIR); //inputDialog.getMorphologiesDir();
        ephysTracesDir = QDir(TRACES_DIR);

        filePaths.ephysUMapFilePath = EPHYS_UMAP_PATH;
        filePaths.morphoUMapFilePath = MORPHO_UMAP_PATH;
        filePaths.txUMapFilePath = TX_UMAP_PATH;
#ifdef ME_UMAP_PATH
        filePaths.meUMapFilePath = ME_UMAP_PATH;
#endif
        filePaths.ephysTracesDir = TRACES_DIR;
        filePaths.morphologiesDir = MORPHO_DIR;
    }

    // Locate all the necessary patch-seq files
    //if (!filePaths.allFilesLocated())
    //{
    //    qDebug() << "Failed to locate all of the necessary patch-seq files.";
    //    return;
    //}

    _task.setEnabled(true);
    _task.setRunning();
    QCoreApplication::processEvents();

    // Load taxonomy
    qDebug() << "Reading taxonomy annotations from file..";
    //Taxonomy taxonomy = Taxonomy::fromJsonFile();
    //taxonomy.printTree();
    //_taxonomyDf.readFromFile(":met_loader/hodge_taxonomy.csv");
#if defined(DALLEYLEE) || defined(WALEBOER)
    _taxonomyDf.readFromFile(":met_loader/SEAAD_colors.csv");
    buildMapOfCellTypesColors(_taxonomyDf, _cellTypeColors);
#else
    _taxonomyDf.readFromFile(":met_loader/HMBA_BG_consensus_colors_all_levels.csv");
    _colorTaxonomy.processFromDataFrame(_taxonomyDf);
#endif

    qDebug() << "Loading CSV file: " << filePaths.metadataFilePath;

    // Read metadata file
    _metadataDf.readFromFile(filePaths.metadataFilePath);
    _metadataDf.removeDuplicateRows(CELL_ID_TAG);

    // Gene expression data
    if (filePaths.hasGeneExpressions())
    {
        //qDebug() << "Loading gene expression data..";
        //// Read gene expression file
        //_task.setSubtaskStarted("Loading Transcriptomics");
        //loadGeneExpressionData(filePaths.gexprFilePath, _metadataDf);
        //_task.setSubtaskFinished("Loading Transcriptomics");

        //_gexprMetadata = DataFrame::subsetAndReorderByColumn(_metadataDf, _transcriptomicsDf, CELL_ID_TAG, CELL_ID_TAG);

        //// Add cluster meta data
        //addTaxonomyClustersForDf(_transcriptomicsDf, _gexprMetadata, _taxonomyDf, QFileInfo(filePaths.gexprFilePath).baseName(), _geneExpressionData, METADATA_CLUSTER_LABEL);
        //addTaxonomyClustersForDf(_transcriptomicsDf, _gexprMetadata, _taxonomyDf, QFileInfo(filePaths.gexprFilePath).baseName(), _geneExpressionData, METADATA_SUBCLASS_LABEL);
    }

    if (filePaths.hasEphysFeatures())
    {
        qDebug() << "Load electrophysiology feature data..";
        // Read electrophysiology file
        loadEphysData(filePaths.ephysFilePath, _metadataDf);
        qDebug() << "bee";
        // Subset and reorder the metadata
        _ephysMetadata = DataFrame::subsetAndReorderByColumn(_metadataDf, _ephysDf, CELL_ID_TAG, CELL_ID_TAG);
        qDebug() << "bee1";
    }

    if (filePaths.hasMorphologyFeatures())
    {
        qDebug() << "bee4";
        loadMorphologyData(filePaths.morphoFilePath, _metadataDf);
    }

    _task.setFinished();
    qDebug() << "bee5";
    //----------------------------------------------------------------------------------------------------------------------
    // Make a metadata text dataset and adds its columns, tries to assign a proper header name if possible 
    //----------------------------------------------------------------------------------------------------------------------
    _metadata = mv::data().createDataset<Text>("Text", "Cell Metadata", mv::Dataset<DatasetImpl>(), "", false);
    _metadata->setProperty("PatchSeqType", "Metadata");
    qDebug() << "bee6";
    for (int i = 0; i < _metadataDf.getHeaders().size(); i++)
    {
        const QString& header = _metadataDf.getHeaders()[i];

        // Get a proper name if possible
        QString properHeaderName = header;
        if (properFeatureNames.find(header) != properFeatureNames.end())
            properHeaderName = properFeatureNames[header];

        std::vector<QString> column = _metadataDf[header];
        _metadata->addColumn(properHeaderName, column);
    }
    qDebug() << "bee7";
    //----------------------------------------------------------------------------------------------------------------------
    // Link up all the datasets
    //----------------------------------------------------------------------------------------------------------------------
    // Take columns from ephys and morpho data and order them correctly, filling in missing data
    BiMap metadataBiMap;
    std::vector<uint32_t> metaCellIdIndices(_metadataDf.numRows());
    std::iota(metaCellIdIndices.begin(), metaCellIdIndices.end(), 0);
    metadataBiMap.addKeyValuePairs(_metadataDf[CELL_ID_TAG], metaCellIdIndices);

    qDebug() << "bee8";
    std::vector<QString> ephysCellIdColumn = _ephysMetadata[CELL_ID_TAG];
    std::vector<uint32_t> ephysToMetaIndices = metadataBiMap.getValuesByKeys(ephysCellIdColumn);

    std::vector<QString> morphoCellIdColumn = _morphoMetadata[CELL_ID_TAG];
    std::vector<uint32_t> morphoToMetaIndices = metadataBiMap.getValuesByKeys(morphoCellIdColumn);
    qDebug() << "bee9";
    // Add ephys and morpho data to metadata dataset
    addPointsToTextDataset(_ephysData, _metadata, ephysToMetaIndices); qDebug() << "bee10";
    addPointsToTextDataset(_morphoData, _metadata, morphoToMetaIndices);
    qDebug() << "bee11";
    events().notifyDatasetAdded(_metadata);
    events().notifyDatasetDataDimensionsChanged(_metadata);

    //createClusterData(_metadataDf[METADATA_CLUSTER_LABEL], "tree_cluster", _metadata);

    qDebug() << ">>>>>>>>>>>>>> Loading morphology cells";
    loadMorphologyCells(morphologiesDir);
    
    qDebug() << ">>>>>>>>>>>>>> Loading ephys cells";

    if (filePaths.hasEphysTraces())
        loadEphysTraces(ephysTracesDir);
    else
        qDebug() << "No ephys traces directory set, skipping loading..";

    qDebug() << ">>>>>>>>>>>>>> Making selection group";
    // Make selection group
    _selectionGroup.addDataset(_metadata, metadataBiMap);

    if (filePaths.hasGeneExpressions())
    {
        //// Gene expression data
        //BiMap gexprBiMap;
        //std::vector<uint32_t> gexprIndices(_geneExpressionData->getNumPoints());
        //std::iota(gexprIndices.begin(), gexprIndices.end(), 0);
        //gexprBiMap.addKeyValuePairs(_transcriptomicsDf[CELL_ID_TAG], gexprIndices);
        //qDebug() << "Gexpr: " << _transcriptomicsDf[CELL_ID_TAG].size() << gexprIndices.size();

        //_selectionGroup.addDataset(_geneExpressionData, gexprBiMap);
    }

    if (filePaths.hasEphysFeatures())
    {
        // Ephys data
        BiMap ephysBiMap;
        std::vector<uint32_t> ephysIndices(_ephysData->getNumPoints());
        std::iota(ephysIndices.begin(), ephysIndices.end(), 0);
        _ephysMetadata.removeDuplicateRows(CELL_ID_TAG);
        ephysBiMap.addKeyValuePairs(_ephysDf[CELL_ID_TAG], ephysIndices);
        qDebug() << "Ephys: " << _ephysMetadata[CELL_ID_TAG].size() << ephysIndices.size();

        // Ephys UMAP
        if (filePaths.hasEphysUMap())
            loadUMap(filePaths.ephysUMapFilePath, _ephysData, "Ephys UMAP", ephysBiMap);

        _selectionGroup.addDataset(_ephysData, ephysBiMap);

        // Add cluster meta data
        addTaxonomyClustersForDf(_ephysDf, _ephysMetadata, TaxonomyLevel::GROUP, QFileInfo(filePaths.ephysFilePath).baseName(), _ephysData, METADATA_CLUSTER_LABEL); qDebug() << "bee2";
        addTaxonomyClustersForDf(_ephysDf, _ephysMetadata, TaxonomyLevel::SUBCLASS, QFileInfo(filePaths.ephysFilePath).baseName(), _ephysData, METADATA_SUBCLASS_LABEL); qDebug() << "bee3";
#ifdef DALLEYLEE
        addColorizedClustersFromMetadata(_metadata, _ephysData, _selectionGroup, "paradigm");
        addColorizedClustersFromMetadata(_metadata, _ephysData, _selectionGroup, "lobe");
#endif
    }

    if (filePaths.hasMorphologyFeatures())
    {
        // Morphology data
        BiMap morphBiMap;
        std::vector<uint32_t> morphoIndices(_morphoData->getNumPoints());
        std::iota(morphoIndices.begin(), morphoIndices.end(), 0);
        morphBiMap.addKeyValuePairs(_morphologyDf[CELL_ID_TAG], morphoIndices);
        qDebug() << "Morph: " << _morphoMetadata[CELL_ID_TAG].size() << morphoIndices.size();

        // Morphology UMAP
        if (filePaths.hasMorphoUMap())
            loadUMap(filePaths.morphoUMapFilePath, _morphoData, "Morpho UMAP", morphBiMap);

        _selectionGroup.addDataset(_morphoData, morphBiMap);

        // Add cluster meta data
        addTaxonomyClustersForDf(_morphologyDf, _morphoMetadata, TaxonomyLevel::GROUP, QFileInfo(filePaths.morphoFilePath).baseName(), _morphoData, METADATA_CLUSTER_LABEL);
        addTaxonomyClustersForDf(_morphologyDf, _morphoMetadata, TaxonomyLevel::SUBCLASS, QFileInfo(filePaths.morphoFilePath).baseName(), _morphoData, METADATA_SUBCLASS_LABEL);
#ifdef DALLEYLEE
        addColorizedClustersFromMetadata(_metadata, _morphoData, _selectionGroup, "paradigm");
        addColorizedClustersFromMetadata(_metadata, _morphoData, _selectionGroup, "lobe");
#endif
    }

    if (filePaths.hasEphysTraces())
    {
        // Ephys traces mapping
        BiMap ephysTraceBiMap;
        std::vector<uint32_t> ephysTraceIndices(_ephysTraceCellIds.size());
        std::iota(ephysTraceIndices.begin(), ephysTraceIndices.end(), 0);
        ephysTraceBiMap.addKeyValuePairs(_ephysTraceCellIds, ephysTraceIndices);

        _selectionGroup.addDataset(_ephysTraces, ephysTraceBiMap);
    }

    // Cell ID data
    BiMap cellIdBiMap;
    std::vector<uint32_t> cellIdIndices(_metadata->getNumRows());
    std::iota(cellIdIndices.begin(), cellIdIndices.end(), 0);
    cellIdBiMap.addKeyValuePairs(_metadataDf[CELL_ID_TAG], cellIdIndices);
    qDebug() << "Metadata: " << _metadataDf[CELL_ID_TAG].size() << cellIdIndices.size();

    // Morphology mapping
    BiMap cellMorphologyBiMap;
    std::vector<uint32_t> cmIndices(_cellMorphoData->getData().size());
    std::iota(cmIndices.begin(), cmIndices.end(), 0);
    QStringList cmCellIds = _cellMorphoData->getCellIdentifiers();
    std::vector<QString> cmCellIdVector(cmCellIds.size());
    for (int i = 0; i < cmCellIdVector.size(); i++)
    {
        cmCellIdVector[i] = cmCellIds[i];
    }
    cellMorphologyBiMap.addKeyValuePairs(cmCellIdVector, cmIndices);

    // ME UMAP
    if (filePaths.hasMEUMap())
    {
        MatrixData umapData;
        MatrixDataLoader matrixDataLoader(true);
        DataFrame umapDf;
        matrixDataLoader.LoadMatrixData(filePaths.meUMapFilePath, umapDf, umapData, 1);

        // Create dataset
        Dataset<Points> umapDataset = mv::data().createDataset("Points", "MorphoElectric UMAP", mv::Dataset<DatasetImpl>(), "", false);
        umapDataset->setData(umapData.data, umapData.numCols);
        umapDataset->setDimensionNames(umapData.headers);

        events().notifyDatasetAdded(umapDataset);
        events().notifyDatasetDataChanged(umapDataset);

        // UMAP BiMap
        BiMap umapBiMap;
        std::vector<uint32_t> indices(umapDf[CELL_ID_TAG].size());
        std::iota(indices.begin(), indices.end(), 0);
        umapBiMap.addKeyValuePairs(umapDf[CELL_ID_TAG], indices);

        _selectionGroup.addDataset(umapDataset, umapBiMap);

#ifdef DALLEYLEE
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "Supertype", _cellTypeColors);
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "Subclass", _cellTypeColors);
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "paradigm");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "lobe");
#endif
    }

    // Transcriptomics UMAP
    if (filePaths.hasTxUMap())
    {
        MatrixData umapData;
        MatrixDataLoader matrixDataLoader(true);
        DataFrame umapDf;
        matrixDataLoader.LoadMatrixData(filePaths.txUMapFilePath, umapDf, umapData, 6);

        //// TEMP Check data
        //float minV = std::numeric_limits<float>::max();
        //float maxV = -std::numeric_limits<float>::max();
        //for (int i = 0; i < umapData.data.size(); i++)
        //{
        //    if (umapData.data[i] < minV) minV = umapData.data[i];
        //    if (umapData.data[i] > maxV) maxV = umapData.data[i];

        //    if (umapData.data[i] > 10000)
        //    {
        //        qDebug() << i;
        //    }
        //}
        
        // Create dataset
        Dataset<Points> umapDataset = mv::data().createDataset("Points", "Tx UMAP", mv::Dataset<DatasetImpl>(), "", false);
        umapDataset->setData(umapData.data, umapData.numCols);
        umapDataset->setDimensionNames(umapData.headers);

        events().notifyDatasetAdded(umapDataset);
        events().notifyDatasetDataChanged(umapDataset);

        // Transcriptomics UMAP BiMap
        BiMap txUmapBiMap;
        std::vector<uint32_t> indices(umapDf[CELL_ID_TAG].size());
        std::iota(indices.begin(), indices.end(), 0);
        const std::vector<QString>& keys = umapDf[CELL_ID_TAG];
        for (int i = 0; i < keys.size(); i++)
        {
            if (keys[i].isNull() || keys[i].isEmpty())
                indices[i] = -1;
        }

        txUmapBiMap.addKeyValuePairs(keys, indices);

        _selectionGroup.addDataset(umapDataset, txUmapBiMap);

        // Annotation metadata
        {
            // Create a list of clusters and their indices from the list of cluster names
            Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", "Supertype", umapDataset);

            const std::vector<QString>& clusterAsList = umapDf["celltype"];
            std::map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(clusterAsList);

            for (auto& kv : clusterMap)
            {
                Cluster cluster;

                cluster.setName(kv.first);
                cluster.setIndices(kv.second);

#if defined(DALLEYLEE) || defined(WALEBOER)
                if (_cellTypeColors.contains(cluster.getName()))
                    cluster.setColor(_cellTypeColors[cluster.getName()]);
#else
                const QColor& color = _colorTaxonomy.getGroupColor(cluster.getName());
                cluster.setColor(color);
#endif

                clusterData->addCluster(cluster);
            }

            mv::events().notifyDatasetDataChanged(clusterData);
            mv::events().notifyDatasetDataDimensionsChanged(clusterData);
        }
        // Subclass
        {
            // Create a list of clusters and their indices from the list of cluster names
            Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", "Subclass", umapDataset);

            const std::vector<QString>& clusterAsList = umapDf["subclass"];
            std::map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(clusterAsList);

            for (auto& kv : clusterMap)
            {
                Cluster cluster;

                cluster.setName(kv.first);
                cluster.setIndices(kv.second);

#if defined(DALLEYLEE) || defined(WALEBOER)
                if (_cellTypeColors.contains(cluster.getName()))
                    cluster.setColor(_cellTypeColors[cluster.getName()]);
#else
                const QColor& color = _colorTaxonomy.getGroupColor(cluster.getName());
                cluster.setColor(color);
#endif

                clusterData->addCluster(cluster);
            }

            mv::events().notifyDatasetDataChanged(clusterData);
            mv::events().notifyDatasetDataDimensionsChanged(clusterData);
        }

#ifdef DALLEYLEE
        //addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "Supertype", _cellTypeColors);
        //addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "Subclass", _cellTypeColors);
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "paradigm");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "days_in_culture");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "patched_cell_container");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "donor_name");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "medical_conditions");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "lobe");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "age_label(yrs)");
        addColorizedClustersFromMetadata(_metadata, umapDataset, _selectionGroup, "sex");
#endif
    }

    // Set cell morphology colors
    {
        // Get std::vector of cell ids of morphologies
        QStringList ids = _cellMorphoData->getCellIdentifiers();
        std::vector<QString> cellIds(ids.size());
        for (int i = 0; i < ids.size(); i++)
        {
            cellIds[i] = ids[i];
        }

        // Map these cell ids to metadata indices
        std::vector<int> indices = cellIdBiMap.getValuesByKeysWithMissingValue(cellIds, -1);

        std::vector<QString> clusterLabels = _metadataDf[METADATA_CLUSTER_LABEL];

        for (int i = 0; i < indices.size(); i++)
        {
            if (indices[i] == -1) continue;

            QString label = clusterLabels[indices[i]];

#if defined(DALLEYLEE) || defined(WALEBOER)
            QColor color = _cellTypeColors[label];
            _cellMorphoData->getData()[i].cellTypeColor.set(color.redF(), color.greenF(), color.blueF());
#else
            const QColor& color = _colorTaxonomy.getGroupColor(label);
            _cellMorphoData->getData()[i].cellTypeColor.set(color.redF(), color.greenF(), color.blueF());
#endif
        }
    }

    _selectionGroup.addDataset(_metadata, cellIdBiMap);
    _selectionGroup.addDataset(_cellMorphoData, cellMorphologyBiMap);

    events().addSelectionGroup(_selectionGroup);
}

void PatchSeqDataLoader::loadGeneExpressionData(QString filePath, const DataFrame& metadata)
{
    qDebug() << "Loading transcriptomic data..";

    MatrixData matrixData;
    MatrixDataLoader matrixDataLoader;
    matrixDataLoader.LoadMatrixData(filePath, _transcriptomicsDf, matrixData, 3);
    //matrixData.standardize();

    removeDuplicateRows(_transcriptomicsDf, CELL_ID_TAG, matrixData);

    _transcriptomicsDf.printFirstFewDimensionsOfDataFrame();

    _geneExpressionData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName(), mv::Dataset<DatasetImpl>(), "", false);
    _geneExpressionData->setProperty("PatchSeqType", "T");
    _geneExpressionData->setData(matrixData.data, matrixData.numCols);
    _geneExpressionData->setDimensionNames(matrixData.headers);

    events().notifyDatasetAdded(_geneExpressionData);
    events().notifyDatasetDataChanged(_geneExpressionData);
    events().notifyDatasetDataDimensionsChanged(_geneExpressionData);
}

void PatchSeqDataLoader::loadEphysData(QString filePath, const DataFrame& metadata)
{
    qDebug() << "Loading electrophysiology data..";

    MatrixData matrixData;
    MatrixDataLoader matrixDataLoader(true);
    matrixDataLoader.LoadMatrixData(filePath, _ephysDf, matrixData, 2);

    removeDuplicateRows(_ephysDf, CELL_ID_TAG, matrixData);
    removeRowsNotInMetadata(_ephysDf, CELL_ID_TAG, _metadataDf, matrixData);
    removeFeatures(matrixData, featuresToDelete);
    //removeRowsWithAllDataMissing(_ephysDf, matrixData);
    matrixData.imputeMissingValues();
    //matrixData.standardize();

    _ephysDf.printFirstFewDimensionsOfDataFrame();
    qDebug() << "PreCreate";
    _ephysData = mv::data().createDataset<Points>("Points", "Ephys Feature Data", mv::Dataset<DatasetImpl>(), "", false); //QFileInfo(filePath).baseName()
    qDebug() << "PostCreate";
    _ephysData->setProperty("PatchSeqType", "E");
    _ephysData->setData(matrixData.data, matrixData.numCols);

    // Replace feature names with proper names
    for (int i = 0; i < matrixData.headers.size(); i++)
    {
        QString& header = matrixData.headers[i];
        if (properEphysFeatureNames.find(header) != properEphysFeatureNames.end())
        {
            header = properEphysFeatureNames[header];
        }
    }
    _ephysData->setDimensionNames(matrixData.headers);

    qDebug() << "PreAdd";
    events().notifyDatasetAdded(_ephysData);
    events().notifyDatasetDataChanged(_ephysData);
    events().notifyDatasetDataDimensionsChanged(_ephysData);
    qDebug() << "PostAdd";
}

void PatchSeqDataLoader::loadMorphologyData(QString filePath, const DataFrame& metadata)
{
    qDebug() << "Loading morphology data..";

    MatrixData matrixData;
    MatrixDataLoader matrixDataLoader(true);
    matrixDataLoader.LoadMatrixData(filePath, _morphologyDf, matrixData, 1);

    // Find 
    std::vector<int> badRows;
    //badRows.push_back(_morphologyDf.findRowWithColumnValue("cell_id", "770268110"));

    _morphologyDf.removeRows(badRows);
    matrixData.removeRows(badRows);

    removeDuplicateRows(_morphologyDf, CELL_ID_TAG, matrixData);
    matrixData.fillMissingValues(0);
    //matrixData.standardize();

    _morphoData = mv::data().createDataset<Points>("Points", "Morphology Feature Data", mv::Dataset<DatasetImpl>(), "", false); //QFileInfo(filePath).baseName()
    _morphoData->setProperty("PatchSeqType", "M");
    _morphoData->setData(matrixData.data, matrixData.numCols);

    // Replace feature names with proper names
    for (int i = 0; i < matrixData.headers.size(); i++)
    {
        QString& header = matrixData.headers[i];
        if (properMorphologyFeatureNames.find(header) != properMorphologyFeatureNames.end())
        {
            header = properMorphologyFeatureNames[header];
        }
    }
    _morphoData->setDimensionNames(matrixData.headers);

    events().notifyDatasetAdded(_morphoData);
    events().notifyDatasetDataChanged(_morphoData);
    events().notifyDatasetDataDimensionsChanged(_morphoData);

    // Subset and reorder the metadata
    _morphoMetadata = DataFrame::subsetAndReorderByColumn(metadata, _morphologyDf, CELL_ID_TAG, CELL_ID_TAG);

    qDebug() << "Successfully loaded" << _morphoData->getNumPoints() << "cell morphologies";
}

void PatchSeqDataLoader::loadMorphologyCells(QDir dir)
{
    Timer timer("SWC Morphology Loading");

    // Load morphology cells
    _cellMorphoData = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", "Cell Morphologies", mv::Dataset<DatasetImpl>(), "", false);
    _cellMorphoData->setProperty("PatchSeqType", "Morphologies");
#if defined(DALLEYLEE) || defined(WALEBOER)
    _cellMorphoData->setProperty("isCortical", true);
#endif

    QDir morphologyDir(dir);

    QStringList swcFiles = morphologyDir.entryList(QStringList() << "*.swc" << "*.SWC", QDir::Files);

    QStringList cellIds;
    std::vector<CellMorphology> cellMorphologies(swcFiles.size());

    SWCLoader loader;
    for (int i = 0; i < cellMorphologies.size(); i++)
    {
        QString swcFile = swcFiles[i];
        CellMorphology& cellMorphology = cellMorphologies[i];

        loader.LoadSWC(morphologyDir.filePath(swcFile), cellMorphology);

        cellMorphology.findCentroid();
        cellMorphology.findExtents();
        cellMorphology.cellTypeColor.set(0.11f, 0.79f, 0);

        cellIds.append(QFileInfo(swcFile).baseName());
    }

    _cellMorphoData->setCellIdentifiers(cellIds);
    _cellMorphoData->setData(cellMorphologies);

    events().notifyDatasetAdded(_cellMorphoData);
    events().notifyDatasetDataChanged(_cellMorphoData);
    events().notifyDatasetDataDimensionsChanged(_cellMorphoData);
}

void PatchSeqDataLoader::loadEphysTraces(QDir dir)
{
    // Create ephys dataset
    _ephysTraces = mv::data().createDataset<EphysExperiments>("Electrophysiology Data", "Ephys Traces", mv::Dataset<DatasetImpl>(), "", false);
    _ephysTraces->setProperty("PatchSeqType", "EphysTraces");

    // Find all .nwb files in given directory
    QDir ephysTracesDir(dir);

    QStringList nwbFiles = ephysTracesDir.entryList(QStringList() << "*.nwb" << "*.NWB", QDir::Files);

    // Map metadata cell names to cell_ids
    std::vector<QString> metaCellSpecimenNames = _metadataDf[CELL_NAME_TAG];
    std::vector<QString> metaCellIds = _metadataDf[CELL_ID_TAG];

    std::unordered_map<QString, QString> specimenNameToCellId;
    for (size_t i = 0; i < metaCellSpecimenNames.size(); ++i)
        specimenNameToCellId[metaCellSpecimenNames[i]] = metaCellIds[i];

    qDebug() << "Found" << nwbFiles.size() << "NWB files, attempting to load them..";

    // Load NWB files and add them to dataset
    LoadInfo loadInfo;
    NWBLoader loader;
    for (int i = 0; i < nwbFiles.size(); i++)
    {
        Experiment experiment;

        QString fileName = nwbFiles[i];
        
        qDebug() << "Loading NWB file: " << ephysTracesDir.filePath(fileName);

        loader.LoadNWB(ephysTracesDir.filePath(fileName), experiment, loadInfo);

        QString specimenName = fileName;
        specimenName.chop(4); // Cut off the .nwb part
        qDebug() << "Specimen name: " << specimenName;
        if (specimenNameToCellId.find(specimenName) == specimenNameToCellId.end())
            continue;

        _ephysTraceCellIds.push_back(specimenNameToCellId[specimenName]);
        qDebug() << "Cell ID: " << _ephysTraceCellIds[_ephysTraceCellIds.size()-1];
        _ephysTraces->addExperiment(std::move(experiment));
    }

    qDebug() << "Ignored stimsets: " << loadInfo.ignoredStimsets;
    qDebug() << "Loaded stimsets: " << loadInfo.loadedStimsets;

    events().notifyDatasetAdded(_ephysTraces);
    events().notifyDatasetDataChanged(_ephysTraces);
}

void PatchSeqDataLoader::loadUMap(QString filePath, mv::Dataset<Points> parent, QString datasetName, BiMap& bimap)
{
    MatrixData umapData;
    MatrixDataLoader matrixDataLoader(false);
    DataFrame umapDf;
    matrixDataLoader.LoadMatrixData(filePath, umapDf, umapData, 1);
    
    // Reorder UMAP points according to parent dataset order
    std::vector<int> parentOrder = bimap.getValuesByKeysWithMissingValue(umapDf[CELL_ID_TAG], -1);
    std::vector<float> reorderedData(parent->getNumPoints() * 2, 0);
    for (int i = 0; i < umapData.numRows; i++)
    {
        int parentIndex = parentOrder[i];
        if (parentIndex == -1)
            continue;

        reorderedData[parentIndex * 2 + 0] = umapData.data[i * 2 + 0];
        reorderedData[parentIndex * 2 + 1] = umapData.data[i * 2 + 1];
    }

    // Create dataset
    Dataset<Points> umapDataset = mv::data().createDerivedDataset(datasetName, parent, parent, false);
    umapDataset->setData(reorderedData, umapData.numCols);
    umapDataset->setDimensionNames(umapData.headers);

    events().notifyDatasetAdded(umapDataset);
    events().notifyDatasetDataChanged(umapDataset);
}

// =============================================================================
// Factory
// =============================================================================

LoaderPlugin* PatchSeqDataLoaderFactory::produce()
{
    return new PatchSeqDataLoader(this);
}

DataTypes PatchSeqDataLoaderFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;
    supportedTypes.append(PointType);
    return supportedTypes;
}

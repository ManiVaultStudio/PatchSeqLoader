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
#include <SelectionGroup.h>

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

Q_PLUGIN_METADATA(IID "studio.manivault.PatchSeqDataLoader")

#define CELL_ID_TAG "cell_id"
#define METADATA_CLUSTER_LABEL "HANN_cluster_label_assignment_winner"
#define METADATA_SUBCLASS_LABEL "HANN_subclass_label_assignment_winner"

#define EPHYS_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Data_February_2025/Electrophysiology/umap_2d_nn15_md0.1_ephys.csv"
#define MORPHO_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Data_February_2025/Morphology/umap_2d_nn15_md0.1_morphometric_dendrite.csv"
#define TRANS_UMAP_PATH "D:/Dropbox/Julian/Patchseq/Data_February_2025/Transcriptomics/umap_seuratQC_mapping_results_SEAAD_NEW_Nov_19th_2024.csv"

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

    std::unordered_map<QString, std::vector<unsigned int>> makeClustersFromList(std::vector<QString> list)
    {
        std::unordered_map<QString, std::vector<unsigned int>> clusterData;

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

        }
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

void PatchSeqDataLoader::addTaxonomyClustersForDf(DataFrame& df, DataFrame& metadata, DataFrame& taxonomyDf, QString name, mv::Dataset<mv::DatasetImpl> parent, QString metaLabel)
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
    std::unordered_map<QString, std::vector<unsigned int>> clusterData = makeClustersFromList(clusterNames);

    for (auto& kv : clusterData)
    {
        Cluster cluster;

        cluster.setName(kv.first);

        if (_cellTypeColors.contains(cluster.getName()))
            cluster.setColor(_cellTypeColors[cluster.getName()]);

        cluster.setIndices(kv.second);

        treeClusterData->addCluster(cluster);
    }

    mv::events().notifyDatasetDataChanged(treeClusterData);
    mv::events().notifyDatasetDataDimensionsChanged(treeClusterData);
}

void PatchSeqDataLoader::createClusterData(std::vector<QString> stringList, QString dataName, mv::Dataset<mv::DatasetImpl> parent)
{
    Dataset<Clusters> clusterData = mv::data().createDataset<Points>("Cluster", dataName, parent);

    std::unordered_map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(stringList);

    for (auto& kv : clusterMap)
    {
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

        filePaths.gexprFilePath = "D:/Dropbox/Julian/Patchseq/Data_Original/IDs_w_tc_data.csv";// inputDialog.getTranscriptomicsFilePath();
        filePaths.ephysFilePath = "D:/Dropbox/Julian/Patchseq/Data_February_2025/Electrophysiology/ephys_features.csv";// "D:/Dropbox/Julian/Patchseq/NewData/240928_human_exc_dataset_rsc369_ephys_data.csv";// inputDialog.getElectrophysiologyFilePath();
        filePaths.morphoFilePath = "D:/Dropbox/Julian/Patchseq/Data_February_2025/Morphology/dendrite_morphometric_features.csv";// inputDialog.getMorphologyFilePath();
        filePaths.metadataFilePath = "D:/Dropbox/Julian/Patchseq/Data_February_2025/metadata_simple_MMC_11_6_24.csv";// inputDialog.getMetadataFilePath();
        morphologiesDir = QDir("D:/Dropbox/Julian/Patchseq/Data_February_2025/Morphology/SWC_LayerAligned"); //inputDialog.getMorphologiesDir();
        ephysTracesDir = QDir("D:/Dropbox/Julian/Patchseq/Data_May_2025/Electrophysiology/NWB_files");
    }

    // Locate all the necessary patch-seq files
    if (!filePaths.allFilesLocated())
    {
        qDebug() << "Failed to locate all of the necessary patch-seq files.";
        return;
    }

    _task.setEnabled(true);
    _task.setRunning();
    QCoreApplication::processEvents();

    // Load taxonomy
    qDebug() << "Reading taxonomy annotations from file..";
    //Taxonomy taxonomy = Taxonomy::fromJsonFile();
    //taxonomy.printTree();
    //_taxonomyDf.readFromFile(":met_loader/hodge_taxonomy.csv");
    _taxonomyDf.readFromFile(":met_loader/SEAAD_colors.csv");
    buildMapOfCellTypesColors(_taxonomyDf, _cellTypeColors);

    qDebug() << "Loading CSV file: " << filePaths.metadataFilePath;

    // Read metadata file
    _metadataDf.readFromFile(filePaths.metadataFilePath);
    _metadataDf.removeDuplicateRows(CELL_ID_TAG);

    // Read gene expression file
    _task.setSubtaskStarted("Loading Transcriptomics");
    loadGeneExpressionData(filePaths.gexprFilePath, _metadataDf);
    _task.setSubtaskFinished("Loading Transcriptomics");

    // Read electrophysiology file
    loadEphysData(filePaths.ephysFilePath, _metadataDf);

    // Subset and reorder the metadata
    qDebug() << "<<<<<<<<<<<<<<<<<<<<<< GEXPR";
    DataFrame gexpr_metadata = DataFrame::subsetAndReorderByColumn(_metadataDf, _transcriptomicsDf, CELL_ID_TAG, CELL_ID_TAG);
    qDebug() << "<<<<<<<<<<<<<<<<<<<<<< EPHYS";
    // Subset and reorder the metadata
    DataFrame ephys_metadata = DataFrame::subsetAndReorderByColumn(_metadataDf, _ephysDf, CELL_ID_TAG, CELL_ID_TAG);

    // Add cluster meta data
    addTaxonomyClustersForDf(_transcriptomicsDf, gexpr_metadata, _taxonomyDf, QFileInfo(filePaths.gexprFilePath).baseName(), _geneExpressionData, METADATA_CLUSTER_LABEL);
    addTaxonomyClustersForDf(_transcriptomicsDf, gexpr_metadata, _taxonomyDf, QFileInfo(filePaths.gexprFilePath).baseName(), _geneExpressionData, METADATA_SUBCLASS_LABEL);

    // Add cluster meta data
    addTaxonomyClustersForDf(_ephysDf, ephys_metadata, _taxonomyDf, QFileInfo(filePaths.ephysFilePath).baseName(), _ephysData, METADATA_CLUSTER_LABEL);
    addTaxonomyClustersForDf(_ephysDf, ephys_metadata, _taxonomyDf, QFileInfo(filePaths.ephysFilePath).baseName(), _ephysData, METADATA_SUBCLASS_LABEL);

    loadMorphologyData(filePaths.morphoFilePath, _metadataDf);

    _task.setFinished();

    //----------------------------------------------------------------------------------------------------------------------
    // Make a metadata text dataset and adds its columns, tries to assign a proper header name if possible 
    //----------------------------------------------------------------------------------------------------------------------
    _metadata = mv::data().createDataset<Text>("Text", "cell_metadata", mv::Dataset<DatasetImpl>(), "", false);
    _metadata->setProperty("PatchSeqType", "Metadata");

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

    //----------------------------------------------------------------------------------------------------------------------
    // Link up all the datasets
    //----------------------------------------------------------------------------------------------------------------------
    // Take columns from ephys and morpho data and order them correctly, filling in missing data
    BiMap metadataBiMap;
    std::vector<uint32_t> metaCellIdIndices(_metadataDf.numRows());
    std::iota(metaCellIdIndices.begin(), metaCellIdIndices.end(), 0);
    metadataBiMap.addKeyValuePairs(_metadataDf[CELL_ID_TAG], metaCellIdIndices);

    std::vector<QString> ephysCellIdColumn = ephys_metadata[CELL_ID_TAG];
    std::vector<uint32_t> ephysToMetaIndices = metadataBiMap.getValuesByKeys(ephysCellIdColumn);

    std::vector<QString> morphoCellIdColumn = _morphoMetadata[CELL_ID_TAG];
    std::vector<uint32_t> morphoToMetaIndices = metadataBiMap.getValuesByKeys(morphoCellIdColumn);

    // Add ephys and morpho data to metadata dataset
    addPointsToTextDataset(_ephysData, _metadata, ephysToMetaIndices);
    addPointsToTextDataset(_morphoData, _metadata, morphoToMetaIndices);

    events().notifyDatasetAdded(_metadata);
    events().notifyDatasetDataDimensionsChanged(_metadata);

    createClusterData(_metadataDf[METADATA_CLUSTER_LABEL], "tree_cluster", _metadata);

    qDebug() << ">>>>>>>>>>>>>> Loading morphology cells";
    loadMorphologyCells(morphologiesDir);

    qDebug() << ">>>>>>>>>>>>>> Loading ephys cells";
    loadEphysTraces(ephysTracesDir);

    qDebug() << ">>>>>>>>>>>>>> Making selection group";
    // Make selection group
    KeyBasedSelectionGroup selectionGroup;

    // Gene expression data
    BiMap gexprBiMap;
    std::vector<uint32_t> gexprIndices(_geneExpressionData->getNumPoints());
    std::iota(gexprIndices.begin(), gexprIndices.end(), 0);
    gexprBiMap.addKeyValuePairs(_transcriptomicsDf[CELL_ID_TAG], gexprIndices);
    qDebug() << "Gexpr: " << _transcriptomicsDf[CELL_ID_TAG].size() << gexprIndices.size();

    // Morphology data
    BiMap morphBiMap;
    std::vector<uint32_t> morphoIndices(_morphoData->getNumPoints());
    std::iota(morphoIndices.begin(), morphoIndices.end(), 0);
    morphBiMap.addKeyValuePairs(_morphologyDf[CELL_ID_TAG], morphoIndices);
    qDebug() << "Morph: " << _morphoMetadata[CELL_ID_TAG].size() << morphoIndices.size();

    // Ephys data
    BiMap ephysBiMap;
    std::vector<uint32_t> ephysIndices(_ephysData->getNumPoints());
    std::iota(ephysIndices.begin(), ephysIndices.end(), 0);
    ephys_metadata.removeDuplicateRows(CELL_ID_TAG);
    ephysBiMap.addKeyValuePairs(_ephysDf[CELL_ID_TAG], ephysIndices);
    qDebug() << "Ephys: " << ephys_metadata[CELL_ID_TAG].size() << ephysIndices.size();

    // Ephys traces mapping
    BiMap ephysTraceBiMap;
    std::vector<uint32_t> ephysTraceIndices(_ephysTraceCellIds.size());
    std::iota(ephysTraceIndices.begin(), ephysTraceIndices.end(), 0);
    ephysTraceBiMap.addKeyValuePairs(_ephysTraceCellIds, ephysTraceIndices);

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

    // Ephys UMAP
    loadUMap(EPHYS_UMAP_PATH, _ephysData, "Ephys UMAP", ephysBiMap);

    // Morphology UMAP
    loadUMap(MORPHO_UMAP_PATH, _morphoData, "Morpho UMAP", morphBiMap);

    // Transcriptomics UMAP
    {
        MatrixData umapData;
        MatrixDataLoader matrixDataLoader(true);
        DataFrame umapDf;
        matrixDataLoader.LoadMatrixData("D:/Dropbox/Julian/Patchseq/Data_February_2025/Transcriptomics/umap_seuratQC_mapping_results_SEAAD_NEW_Nov_19th_2024.csv", umapDf, umapData, 10);

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
        txUmapBiMap.addKeyValuePairs(umapDf[CELL_ID_TAG], indices);

        selectionGroup.addDataset(umapDataset, txUmapBiMap);

        // Annotation metadata
        {
            // Create a list of clusters and their indices from the list of cluster names
            Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", "HANNtype", umapDataset);

            const std::vector<QString>& clusterAsList = umapDf["HANNtype"];
            std::unordered_map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(clusterAsList);

            for (auto& kv : clusterMap)
            {
                Cluster cluster;

                cluster.setName(kv.first);
                cluster.setIndices(kv.second);

                if (_cellTypeColors.contains(cluster.getName()))
                    cluster.setColor(_cellTypeColors[cluster.getName()]);

                clusterData->addCluster(cluster);
            }

            mv::events().notifyDatasetDataChanged(clusterData);
            mv::events().notifyDatasetDataDimensionsChanged(clusterData);
        }
        {
            // Create a list of clusters and their indices from the list of cluster names
            Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", "platform", umapDataset);

            const std::vector<QString>& clusterAsList = umapDf["platform"];
            std::unordered_map<QString, std::vector<unsigned int>> clusterMap = makeClustersFromList(clusterAsList);

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
    }

    // Set cell morphology colors
    {
        QStringList ids = _cellMorphoData->getCellIdentifiers();
        std::vector<QString> cellIds(ids.size());
        for (int i = 0; i < ids.size(); i++)
        {
            cellIds[i] = ids[i];
        }

        std::vector<uint32_t> indices = cellIdBiMap.getValuesByKeys(cellIds);

        std::vector<QString> clusterLabels = _metadataDf[METADATA_CLUSTER_LABEL];

        for (int i = 0; i < indices.size(); i++)
        {
            QString label = clusterLabels[indices[i]];

            QColor color = _cellTypeColors[label];
            _cellMorphoData->getData()[i].cellTypeColor.set(color.redF(), color.greenF(), color.blueF());
        }
    }

    selectionGroup.addDataset(_geneExpressionData, gexprBiMap);
    selectionGroup.addDataset(_morphoData, morphBiMap);
    selectionGroup.addDataset(_ephysData, ephysBiMap);
    selectionGroup.addDataset(_ephysTraces, ephysTraceBiMap);
    selectionGroup.addDataset(_metadata, cellIdBiMap);
    selectionGroup.addDataset(_cellMorphoData, cellMorphologyBiMap);

    events().addSelectionGroup(selectionGroup);
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
    //removeRowsWithAllDataMissing(_ephysDf, matrixData);
    matrixData.imputeMissingValues();
    //matrixData.standardize();

    _ephysDf.printFirstFewDimensionsOfDataFrame();

    _ephysData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName(), mv::Dataset<DatasetImpl>(), "", false);
    _ephysData->setProperty("PatchSeqType", "E");
    _ephysData->setData(matrixData.data, matrixData.numCols);
    _ephysData->setDimensionNames(matrixData.headers);

    events().notifyDatasetAdded(_ephysData);
    events().notifyDatasetDataChanged(_ephysData);
    events().notifyDatasetDataDimensionsChanged(_ephysData);
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

    _morphoData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName(), mv::Dataset<DatasetImpl>(), "", false);
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

    // Add cluster meta data
    addTaxonomyClustersForDf(_morphologyDf, _morphoMetadata, _taxonomyDf, QFileInfo(filePath).baseName(), _morphoData, METADATA_CLUSTER_LABEL);
    addTaxonomyClustersForDf(_morphologyDf, _morphoMetadata, _taxonomyDf, QFileInfo(filePath).baseName(), _morphoData, METADATA_SUBCLASS_LABEL);

    qDebug() << "Successfully loaded" << _morphoData->getNumPoints() << "cell morphologies";
}

void PatchSeqDataLoader::loadMorphologyCells(QDir dir)
{
    Timer timer("SWC Morphology Loading");

    // Load morphology cells
    _cellMorphoData = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", "cell_morphology", mv::Dataset<DatasetImpl>(), "", false);
    _cellMorphoData->setProperty("PatchSeqType", "Morphologies");

    QDir morphologyDir(dir);
    morphologyDir.cd("SWC_LayerAligned");

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
    _ephysTraces = mv::data().createDataset<EphysExperiments>("Electrophysiology Data", "EphysTraces", mv::Dataset<DatasetImpl>(), "", false);
    _ephysTraces->setProperty("PatchSeqType", "EphysTraces");

    // Find all .nwb files in given directory
    QDir ephysTracesDir(dir);
    ephysTracesDir.cd("nwb2_files");

    QStringList nwbFiles = ephysTracesDir.entryList(QStringList() << "*.nwb" << "*.NWB", QDir::Files);

    // Map metadata cell names to cell_ids
    std::vector<QString> metaCellSpecimenNames = _metadataDf["cell_specimen_name"];
    std::vector<QString> metaCellIds = _metadataDf["cell_id"];

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

void PatchSeqDataLoader::loadUMap(QString filePath, mv::Dataset<DatasetImpl> parent, QString datasetName, BiMap& bimap)
{
    MatrixData umapData;
    MatrixDataLoader matrixDataLoader(false);
    DataFrame umapDf;
    matrixDataLoader.LoadMatrixData(filePath, umapDf, umapData, 1);

    // Reorder UMAP points according to parent dataset order
    std::vector<uint32_t> parentOrder = bimap.getValuesByKeys(umapDf[CELL_ID_TAG]);
    std::vector<float> reorderedData(umapData.data.size());
    for (int i = 0; i < umapData.numRows; i++)
    {
        int parentIndex = parentOrder[i];
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

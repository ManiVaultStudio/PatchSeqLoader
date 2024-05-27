#include "PatchSeqDataLoader.h"

#include "InputDialog.h"

#include "MatrixDataLoader.h"
#include "MatrixData.h"

#include "Taxonomy.h"
#include "CellLoader.h"
#include "FeatureNames.h"

#include "csv.h"

#include <CellMorphologyData/CellMorphologyData.h>
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

void PatchSeqDataLoader::addTaxonomyClustersForDf(DataFrame& df, DataFrame& metadata, DataFrame& taxonomyDf, QString name, mv::Dataset<mv::DatasetImpl> parent)
{
    std::vector<QString> treeCluster = metadata["tree_cluster"];

    Dataset<Clusters> treeClusterData;
    treeClusterData = mv::data().createDataset<Points>("Cluster", name, parent);

    std::unordered_map<QString, std::vector<unsigned int>> clusterData = makeClustersFromList(treeCluster);

    std::vector<QString> final_cluster_names = taxonomyDf["final cluster"];
    std::vector<QString> final_cluster_colors = taxonomyDf["final cluster color"];

    for (auto& kv : clusterData)
    {
        Cluster cluster;

        cluster.setName(kv.first);

        bool found = false;
        for (int i = 0; i < final_cluster_names.size(); i++)
        {
            if (final_cluster_names[i] == kv.first)
            {
                cluster.setColor(final_cluster_colors[i]);
                found = true;
                //qDebug() << "Found attributes in taxonomy CSV!";
                break;
            }
        }
        if (found == false)
        {
            qDebug() << "ERROR: Failed to find attributes in taxonomy CSV: " << kv.first;
        }

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
    if (ok == QDialog::Accepted)
    {
        filePaths.gexprFilePath = inputDialog.getTranscriptomicsFilePath();
        filePaths.ephysFilePath = inputDialog.getElectrophysiologyFilePath();
        filePaths.morphoFilePath = inputDialog.getMorphologyFilePath();
        filePaths.metadataFilePath = inputDialog.getMetadataFilePath();
        morphologiesDir = inputDialog.getMorphologiesDir();
    }

    //// Locate all the necessary patch-seq files
    //filePaths.locateFilePaths(dir);

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
    _taxonomyDf.readFromFile(":met_loader/hodge_taxonomy.csv");

    qDebug() << "Loading CSV file: " << filePaths.metadataFilePath;

    // Read metadata file
    _metadata.readFromFile(filePaths.metadataFilePath);
    _metadata.removeDuplicateRows(CELL_ID_TAG);

    // Read gene expression file
    _task.setSubtaskStarted("Loading Transcriptomics");
    loadGeneExpressionData(filePaths.gexprFilePath, _metadata);
    _task.setSubtaskFinished("Loading Transcriptomics");

    // Read electrophysiology file
    loadEphysData(filePaths.ephysFilePath, _metadata);

    // Subset and reorder the metadata
    qDebug() << "<<<<<<<<<<<<<<<<<<<<<< GEXPR";
    DataFrame gexpr_metadata = DataFrame::subsetAndReorderByColumn(_metadata, _transcriptomicsDf, CELL_ID_TAG, CELL_ID_TAG);
    qDebug() << "<<<<<<<<<<<<<<<<<<<<<< EPHYS";
    // Subset and reorder the metadata
    DataFrame ephys_metadata = DataFrame::subsetAndReorderByColumn(_metadata, _ephysDf, CELL_ID_TAG, CELL_ID_TAG);
    qDebug() << "Ephys metadata: ";
    qDebug() << ephys_metadata[CELL_ID_TAG];

    // Add cluster meta data
    addTaxonomyClustersForDf(_transcriptomicsDf, gexpr_metadata, _taxonomyDf, QFileInfo(filePaths.gexprFilePath).baseName(), _geneExpressionData);

    // Add cluster meta data
    addTaxonomyClustersForDf(_ephysDf, ephys_metadata, _taxonomyDf, QFileInfo(filePaths.ephysFilePath).baseName(), _ephysData);

    loadMorphologyData(filePaths.morphoFilePath, _metadata);

    _task.setFinished();

    Dataset<Text> metaDataset;
    metaDataset = mv::data().createDataset<Text>("Text", "cell_metadata", mv::Dataset<DatasetImpl>(), "", false);
    metaDataset->setProperty("PatchSeqType", "Metadata");

    for (int i = 0; i < _metadata.getHeaders().size(); i++)
    {
        const QString& header = _metadata.getHeaders()[i];

        // Get a proper name if possible
        QString properHeaderName = header;
        if (properFeatureNames.find(header) != properFeatureNames.end())
            properHeaderName = properFeatureNames[header];

        std::vector<QString> column = _metadata[header];
        metaDataset->addColumn(properHeaderName, column);
    }

    // Take columns from ephys and morpho data and order them correctly, filling in missing data
    BiMap metadataBiMap;
    std::vector<uint32_t> metaCellIdIndices(_metadata.numRows());
    std::iota(metaCellIdIndices.begin(), metaCellIdIndices.end(), 0);
    metadataBiMap.addKeyValuePairs(_metadata[CELL_ID_TAG], metaCellIdIndices);

    std::vector<QString> ephysCellIdColumn = ephys_metadata[CELL_ID_TAG];
    std::vector<uint32_t> ephysToMetaIndices = metadataBiMap.getValuesByKeys(ephysCellIdColumn);

    std::vector<QString> morphoCellIdColumn = _morphoMetadata[CELL_ID_TAG];
    std::vector<uint32_t> morphoToMetaIndices = metadataBiMap.getValuesByKeys(morphoCellIdColumn);

    // Add ephys and morpho data to metadata dataset
    addPointsToTextDataset(_ephysData, metaDataset, ephysToMetaIndices);
    addPointsToTextDataset(_morphoData, metaDataset, morphoToMetaIndices);

    events().notifyDatasetAdded(metaDataset);
    events().notifyDatasetDataDimensionsChanged(metaDataset);

    createClusterData(_metadata["tree_cluster"], "tree_cluster", metaDataset);

    qDebug() << ">>>>>>>>>>>>>> Loading morphology cells";
    loadMorphologyCells(morphologiesDir);

    qDebug() << ">>>>>>>>>>>>>> Making selection group";
    // Make selection group
    KeyBasedSelectionGroup selectionGroup;

    // Gene expression data
    BiMap gexprBiMap;
    std::vector<uint32_t> gexprIndices(_geneExpressionData->getNumPoints());
    std::iota(gexprIndices.begin(), gexprIndices.end(), 0);
    gexprBiMap.addKeyValuePairs(gexpr_metadata[CELL_ID_TAG], gexprIndices);
    qDebug() << "Gexpr: " << gexpr_metadata[CELL_ID_TAG].size() << gexprIndices.size();

    // Morphology data
    BiMap morphBiMap;
    std::vector<uint32_t> morphoIndices(_morphoData->getNumPoints());
    std::iota(morphoIndices.begin(), morphoIndices.end(), 0);
    morphBiMap.addKeyValuePairs(_morphoMetadata[CELL_ID_TAG], morphoIndices);
    qDebug() << "Morph: " << _morphoMetadata[CELL_ID_TAG].size() << morphoIndices.size();

    // Ephys data
    BiMap ephysBiMap;
    std::vector<uint32_t> ephysIndices(_ephysData->getNumPoints());
    std::iota(ephysIndices.begin(), ephysIndices.end(), 0);
    ephys_metadata.removeDuplicateRows(CELL_ID_TAG);
    ephysBiMap.addKeyValuePairs(ephys_metadata[CELL_ID_TAG], ephysIndices);
    qDebug() << "Ephys: " << ephys_metadata[CELL_ID_TAG].size() << ephysIndices.size();

    // Cell ID data
    BiMap cellIdBiMap;
    std::vector<uint32_t> cellIdIndices(metaDataset->getNumRows());
    std::iota(cellIdIndices.begin(), cellIdIndices.end(), 0);
    cellIdBiMap.addKeyValuePairs(_metadata[CELL_ID_TAG], cellIdIndices);
    qDebug() << "Metadata: " << _metadata[CELL_ID_TAG].size() << cellIdIndices.size();

    selectionGroup.addDataset(_geneExpressionData, gexprBiMap);
    selectionGroup.addDataset(_morphoData, morphBiMap);
    selectionGroup.addDataset(_ephysData, ephysBiMap);
    selectionGroup.addDataset(metaDataset, cellIdBiMap);

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
    removeRowsNotInMetadata(_ephysDf, CELL_ID_TAG, _metadata, matrixData);
    removeRowsWithAllDataMissing(_ephysDf, matrixData);
    matrixData.imputeMissingValues();
    matrixData.standardize();

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

    removeDuplicateRows(_morphologyDf, CELL_ID_TAG, matrixData);
    matrixData.fillMissingValues(0);
    matrixData.standardize();

    _morphoData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName(), mv::Dataset<DatasetImpl>(), "", false);
    _morphoData->setProperty("PatchSeqType", "M");
    _morphoData->setData(matrixData.data, matrixData.numCols);
    _morphoData->setDimensionNames(matrixData.headers);

    events().notifyDatasetAdded(_morphoData);
    events().notifyDatasetDataChanged(_morphoData);
    events().notifyDatasetDataDimensionsChanged(_morphoData);

    // Subset and reorder the metadata
    _morphoMetadata = DataFrame::subsetAndReorderByColumn(metadata, _morphologyDf, CELL_ID_TAG, CELL_ID_TAG);

    // Add cluster meta data
    addTaxonomyClustersForDf(_morphologyDf, _morphoMetadata, _taxonomyDf, QFileInfo(filePath).baseName(), _morphoData);
}

void PatchSeqDataLoader::loadMorphologyCells(QDir dir)
{
    // Load morphology cells
    Dataset<CellMorphologies> cellMorphoData = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", "cell_morphology", mv::Dataset<DatasetImpl>(), "", false);
    cellMorphoData->setProperty("PatchSeqType", "Morphologies");

    QDir morphologyDir(dir);
    morphologyDir.cd("SWC_Upright");

    QStringList swcFiles = morphologyDir.entryList(QStringList() << "*.swc" << "*.SWC", QDir::Files);

    QStringList cellIds;
    std::vector<CellMorphology> cellMorphologies(swcFiles.size());

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < cellMorphologies.size(); i++)
    {
        QString swcFile = swcFiles[i];
        CellMorphology& cellMorphology = cellMorphologies[i];

        std::string fileContents;
        loadCellContentsFromFile(morphologyDir.filePath(swcFile), fileContents);
        //qDebug() << QString::fromStdString(fileContents);

        readCell(fileContents, cellMorphology);

        cellMorphology.findCentroid();
        cellMorphology.findExtents();

        cellIds.append(QFileInfo(swcFile).baseName());
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Loading Morphology elapsed time: " << elapsed.count() << " s\n";

    cellMorphoData->setCellIdentifiers(cellIds);
    cellMorphoData->setData(cellMorphologies);

    events().notifyDatasetAdded(cellMorphoData);
    events().notifyDatasetDataChanged(cellMorphoData);
    events().notifyDatasetDataDimensionsChanged(cellMorphoData);
}

QIcon PatchSeqDataLoaderFactory::getIcon(const QColor& color /*= Qt::black*/) const
{
    return Application::getIconFont("FontAwesome").getIcon("database");
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

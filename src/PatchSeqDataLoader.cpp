#include "PatchSeqDataLoader.h"

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

Q_PLUGIN_METADATA(IID "studio.manivault.PatchSeqDataLoader")

#define CELL_ID_TAG "cell_id"

using namespace mv;
using namespace mv::gui;

namespace
{
    void readDataFrame(DataFrame& df, QString fileName)
    {
        QFile file(fileName);

        bool header = true;
        if (!file.open(QIODevice::ReadOnly))
        {
            throw DataLoadException(fileName, "File was not found at location.");
        }

        QTextStream in(&file);

        while (!in.atEnd())
        {
            QString line = in.readLine();

            QStringList tokens = line.split(",");
            if (header)
            {
                df.setHeaders(tokens);
                header = false;
                continue;
            }

            std::vector<QString> row(tokens.size());
            for (int i = 0; i < tokens.size(); i++)
                row[i] = tokens[i];
                
            df.getData().push_back(row);
        }

        file.close();
    }

    void printFirstFewDimensionsOfDataFrame(const DataFrame& df)
    {
        std::cout << "Loaded file with first 20 dimensions: ";
        for (int i = 0; i < std::min(20, (int) df.getHeaders().size()); i++)
        {
            std::cout << df.getHeaders()[i].toStdString() << ", ";
        }
        std::cout << std::endl;
    }

    void readGeneExpressionDf(DataFrame& df, QString fileName, std::vector<float>& geneExpressionMatrix, mv::ModalTask& task, unsigned int& numCols)
    {
        auto start = std::chrono::high_resolution_clock::now();

        // Read text part of the CSV into a DataFrame
        QFile inputFile(fileName);
        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);

            // Read header
            if (!in.atEnd())
            {
                QString line = in.readLine();

                QStringList tokens = line.split(",");
                numCols = (tokens.size() - 3);

                for (int i = 0; i < tokens.size(); i++)
                {
                    QString& token = tokens[i];
                    token.replace("\"", "");
                }

                df.setHeaders(tokens);
            }

            inputFile.close();
        }
        else
        {
            throw DataLoadException(fileName, "File was not found at location.");
        }

        // Read the numbers part
        std::vector<QString> row(3);
        std::vector<float> geneExpressionRow(numCols);

        io::LineReader fin(fileName.toStdString());

        char* token;

        bool skippedHeader = false;

        int lineCount = 0;
        task.setSubtasks(100);

        while (char* line = fin.next_line())
        {
            if (!skippedHeader) { skippedHeader = true; continue; }
            //qDebug() << line;
            //qDebug() << "new line";
            task.setSubtaskStarted(lineCount / 18);

            int colIndex = 0;
            /* get the first three tokens */
            token = strtok(line, ",");
            row[0] = token;
            token = strtok(NULL, ",");
            row[1] = token;
            token = strtok(NULL, ",");
            row[2] = token;

            row[0].replace("\"", "");
            row[1].replace("\"", "");
            row[2].replace("\"", "");
            df.getData().push_back(row);

            token = strtok(NULL, ",");

            /* walk through other tokens */
            while (token != NULL) {
                //printf(" %s\n", token);

                geneExpressionRow[colIndex] = atof(token);
                colIndex++;

                token = strtok(NULL, ",");

                //qDebug() << colIndex;
            }
            //df.getData().push_back(row);
            geneExpressionMatrix.insert(geneExpressionMatrix.end(), geneExpressionRow.begin(), geneExpressionRow.end());

            task.setSubtaskFinished(lineCount / 18);

            QApplication::processEvents();

            lineCount++;
        }

        qDebug() << "Done";

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "[PatchSeqDataLoader] Gene expression matrix loaded in: " << elapsed.count() << " s\n";
    }

    void readPatchSeqDf(DataFrame& df, QString fileName, int numStringCols, std::vector<float>& matrix, unsigned int& numCols)
    {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Starting read" << std::endl;

        std::vector<QString> row(numStringCols);

        QFile inputFile(fileName);
        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);

            // Read header
            if (!in.atEnd())
            {
                QString line = in.readLine();

                QStringList tokens = line.split(",");
                numCols = (tokens.size() - numStringCols);

                for (int i = 0; i < tokens.size(); i++)
                {
                    QString& token = tokens[i];
                    token.replace("\"", "");
                }

                df.setHeaders(tokens);
            }

            std::vector<float> dataRow(numCols);
            while (!in.atEnd())
            {
                QString line = in.readLine();

                QStringList tokens = line.split(",");

                for (int i = 0; i < numStringCols; i++)
                {
                    QString token = tokens[i];
                    token.replace("\"", "");
                    row[i] = token;
                }

                df.getData().push_back(row);

                for (int col = 0; col < numCols; col++)
                    dataRow[col] = tokens[col + numStringCols].toFloat();

                matrix.insert(matrix.end(), dataRow.begin(), dataRow.end());
            }
            inputFile.close();
        }
        else
        {
            throw DataLoadException(fileName, "File was not found at location.");
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "[PatchSeqDataLoader] " << fileName.toStdString() << " loaded in : " << elapsed.count() << " s\n";
    }

    void readMorphologyDf(DataFrame& df, QString fileName, std::vector<float>& matrix, unsigned int& numCols)
    {
        auto start = std::chrono::high_resolution_clock::now();

        int numStringColumns = 1;

        std::vector<QString> row(numStringColumns);

        QFile inputFile(fileName);
        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);

            // Read header
            if (!in.atEnd())
            {
                QString line = in.readLine();

                QStringList tokens = line.split(",");
                numCols = (tokens.size() - numStringColumns);

                for (int i = 0; i < tokens.size(); i++)
                {
                    QString& token = tokens[i];
                    token.replace("\"", "");
                }

                df.setHeaders(tokens);
            }

            std::vector<float> dataRow(numCols);
            while (!in.atEnd())
            {
                QString line = in.readLine();

                QStringList tokens = line.split(",");

                for (int i = 0; i < numStringColumns; i++)
                {
                    QString token = tokens[i];
                    token.replace("\"", "");
                    row[i] = token;
                }

                df.getData().push_back(row);

                for (int col = 0; col < numCols; col++)
                    dataRow[col] = tokens[col + numStringColumns].toFloat();

                matrix.insert(matrix.end(), dataRow.begin(), dataRow.end());
            }
            inputFile.close();
        }
        else
        {
            throw DataLoadException(fileName, "File was not found at location.");
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "[PatchSeqDataLoader] " << fileName.toStdString() << " loaded in : " << elapsed.count() << " s\n";
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

    std::unordered_map<unsigned int, unsigned int> computeLinkedData(DataFrame& leftDf, DataFrame& rightDf, QString columnNameLeft, QString columnNameRight)
    {
        std::vector<QString> columnLeft = leftDf[columnNameLeft];
        std::vector<QString> columnRight = rightDf[columnNameRight];

        // Make a map out of meta column
        std::unordered_map<QString, int> indexMap;
        for (int i = 0; i < columnRight.size(); i++)
        {
            indexMap[columnRight[i]] = i;
        }

        // Find mapping
        std::unordered_map<unsigned int, unsigned int> mapping;
        for (int i = 0; i < columnLeft.size(); i++)
        {
            const QString& cell_id = columnLeft[i];
            if (indexMap.find(cell_id) == indexMap.end())
            {
                //qDebug() << "Failed to find cell ID: " << cell_id << " in metadata file.";
                continue;
            }
            int index = indexMap[cell_id];
            mapping[i] = index;
        }

        return mapping;
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

    QDir dir = QFileDialog::getExistingDirectory(nullptr, tr("Open Directory"),
        "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    // Locate all the necessary patch-seq files
    PatchSeqFilePaths filePaths;
    filePaths.locateFilePaths(dir);

    if (!filePaths.allFilesLocated())
    {
        qDebug() << "Failed to locate all of the necessary patch-seq files.";
        return;
    }

    // Load taxonomy
    qDebug() << "Reading taxonomy annotations from file..";
    //Taxonomy taxonomy = Taxonomy::fromJsonFile();
    //taxonomy.printTree();
    readDataFrame(_taxonomyDf, ":met_loader/hodge_taxonomy.csv");

    qDebug() << "Loading CSV file: " << filePaths.metadataFilePath;

    // Read metadata file
    DataFrame metadata;
    readDataFrame(metadata, filePaths.metadataFilePath);

    // Read gene expression file
    loadGeneExpressionData(filePaths.gexprFilePath, metadata);

    // Read electrophysiology file
    loadEphysData(filePaths.ephysFilePath, metadata);

    // Subset and reorder the metadata
    qDebug() << "<<<<<<<<<<<<<<<<<<<<<< GEXPR";
    DataFrame gexpr_metadata = DataFrame::subsetAndReorderByColumn(metadata, _geneExpressionDf, CELL_ID_TAG, CELL_ID_TAG);
    qDebug() << "<<<<<<<<<<<<<<<<<<<<<< EPHYS";
    // Subset and reorder the metadata
    DataFrame ephys_metadata = DataFrame::subsetAndReorderByColumn(metadata, _ephysDf, CELL_ID_TAG, CELL_ID_TAG);

    // Add cluster meta data
    addTaxonomyClustersForDf(_geneExpressionDf, gexpr_metadata, _taxonomyDf, QFileInfo(filePaths.gexprFilePath).baseName(), _geneExpressionData);

    // Add cluster meta data
    addTaxonomyClustersForDf(_ephysDf, ephys_metadata, _taxonomyDf, QFileInfo(filePaths.ephysFilePath).baseName(), _ephysData);

    loadMorphologyData(filePaths.morphoFilePath, metadata);

    Dataset<Text> metaDataset;
    metaDataset = mv::data().createDataset<Text>("Text", "cell_metadata", mv::Dataset<DatasetImpl>(), "", false);

    metaDataset->setProperty("PatchSeqType", "Metadata");

    for (int i = 0; i < metadata.getHeaders().size(); i++)
    {
        const QString& header = metadata.getHeaders()[i];

        // Get a proper name if possible
        QString properHeaderName = header;
        if (properFeatureNames.find(header) != properFeatureNames.end())
            properHeaderName = properFeatureNames[header];

        std::vector<QString> column = metadata[header];
        metaDataset->addColumn(properHeaderName, column);
    }

    // Take columns from ephys and morpho data and order them correctly, filling in missing data
    BiMap metadataBiMap;
    std::vector<uint32_t> metaCellIdIndices(metadata.numRows());
    std::iota(metaCellIdIndices.begin(), metaCellIdIndices.end(), 0);
    metadataBiMap.addKeyValuePairs(metadata[CELL_ID_TAG], metaCellIdIndices);

    std::vector<QString> ephysCellIdColumn = ephys_metadata[CELL_ID_TAG];
    std::vector<uint32_t> ephysToMetaIndices = metadataBiMap.getValuesByKeys(ephysCellIdColumn);

    std::vector<QString> morphoCellIdColumn = _morphoMetadata[CELL_ID_TAG];
    std::vector<uint32_t> morphoToMetaIndices = metadataBiMap.getValuesByKeys(morphoCellIdColumn);

    for (int d = 0; d < _ephysData->getNumDimensions(); d++)
    {
        std::vector<QString> column(metadata.numRows(), "?");
        QString dimName = _ephysData->getDimensionNames()[d];

        // Get a proper name if possible
        QString properHeaderName = dimName;
        if (properFeatureNames.find(dimName) != properFeatureNames.end())
            properHeaderName = properFeatureNames[dimName];

        for (size_t i = 0; i < ephysToMetaIndices.size(); i++)
        {
            uint32_t row = ephysToMetaIndices[i];
            column[row] = QString::number(_ephysData->getValueAt(i * _ephysData->getNumDimensions() + d));
        }
        metaDataset->addColumn(properHeaderName, column);
    }

    for (int d = 0; d < _morphoData->getNumDimensions(); d++)
    {
        std::vector<QString> column(metadata.numRows(), "?");
        QString dimName = _morphoData->getDimensionNames()[d];

        // Get a proper name if possible
        QString properHeaderName = dimName;
        if (properFeatureNames.find(dimName) != properFeatureNames.end())
            properHeaderName = properFeatureNames[dimName];

        for (size_t i = 0; i < morphoToMetaIndices.size(); i++)
        {
            uint32_t row = morphoToMetaIndices[i];
            column[row] = QString::number(_morphoData->getValueAt(i * _morphoData->getNumDimensions() + d));
        }
        metaDataset->addColumn(properHeaderName, column);
    }

    events().notifyDatasetAdded(metaDataset);
    events().notifyDatasetDataDimensionsChanged(metaDataset);

    createClusterData(metadata["tree_cluster"], "tree_cluster", metaDataset);

    qDebug() << ">>>>>>>>>>>>>> Loading morphology cells";
    loadMorphologyCells(dir);

    qDebug() << ">>>>>>>>>>>>>> Making selection group";
    // Make selection group
    KeyBasedSelectionGroup selectionGroup;

    // Gene expression data
    BiMap gexprBiMap;
    std::vector<uint32_t> gexprIndices(_geneExpressionData->getNumPoints());
    std::iota(gexprIndices.begin(), gexprIndices.end(), 0);
    gexprBiMap.addKeyValuePairs(gexpr_metadata[CELL_ID_TAG], gexprIndices);

    // Morphology data
    BiMap morphBiMap;
    std::vector<uint32_t> morphoIndices(_morphoData->getNumPoints());
    std::iota(morphoIndices.begin(), morphoIndices.end(), 0);
    morphBiMap.addKeyValuePairs(_morphoMetadata[CELL_ID_TAG], morphoIndices);

    // Ephys data
    BiMap ephysBiMap;
    std::vector<uint32_t> ephysIndices(_ephysData->getNumPoints());
    std::iota(ephysIndices.begin(), ephysIndices.end(), 0);
    ephysBiMap.addKeyValuePairs(ephys_metadata[CELL_ID_TAG], ephysIndices);

    // Cell ID data
    BiMap cellIdBiMap;
    std::vector<uint32_t> cellIdIndices(metaDataset->getNumRows());
    std::iota(cellIdIndices.begin(), cellIdIndices.end(), 0);
    cellIdBiMap.addKeyValuePairs(metadata[CELL_ID_TAG], cellIdIndices);

    selectionGroup.addDataset(_geneExpressionData, gexprBiMap);
    selectionGroup.addDataset(_morphoData, morphBiMap);
    selectionGroup.addDataset(_ephysData, ephysBiMap);
    selectionGroup.addDataset(metaDataset, cellIdBiMap);

    events().addSelectionGroup(selectionGroup);
}

void PatchSeqDataLoader::loadGeneExpressionData(QString filePath, const DataFrame& metadata)
{
    _task.setEnabled(true);
    _task.setRunning();
    QCoreApplication::processEvents();

    std::vector<float> geneExpressionMatrix;
    unsigned int numCols;
    readGeneExpressionDf(_geneExpressionDf, filePath, geneExpressionMatrix, _task, numCols);

    _task.setFinished();

    printFirstFewDimensionsOfDataFrame(_geneExpressionDf);

    _geneExpressionData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName());

    std::vector<QString> gexprDimNames(_geneExpressionDf.getHeaders().begin() + 3, _geneExpressionDf.getHeaders().end());

    _geneExpressionData->setData(geneExpressionMatrix, numCols);
    _geneExpressionData->setDimensionNames(gexprDimNames);

    events().notifyDatasetDataChanged(_geneExpressionData);
    events().notifyDatasetDataDimensionsChanged(_geneExpressionData);
}

void PatchSeqDataLoader::loadEphysData(QString filePath, const DataFrame& metadata)
{
    std::vector<float> ephysMatrix;
    unsigned int numEphysCols;
    readPatchSeqDf(_ephysDf, filePath, 2, ephysMatrix, numEphysCols);

    printFirstFewDimensionsOfDataFrame(_ephysDf);

    _ephysData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName());

    std::vector<QString> ephysDimNames(_ephysDf.getHeaders().begin() + 2, _ephysDf.getHeaders().end());

    _ephysData->setData(ephysMatrix, numEphysCols);
    _ephysData->setDimensionNames(ephysDimNames);

    events().notifyDatasetDataChanged(_ephysData);
    events().notifyDatasetDataDimensionsChanged(_ephysData);
}

void PatchSeqDataLoader::loadMorphologyData(QString filePath, const DataFrame& metadata)
{
    // Load morphology data
    std::vector<float> morphologyMatrix;
    unsigned int numMorphoCols;
    readMorphologyDf(_morphologyDf, filePath, morphologyMatrix, numMorphoCols);

    _morphoData = mv::data().createDataset<Points>("Points", QFileInfo(filePath).baseName());

    std::vector<QString> morphDimNames(_morphologyDf.getHeaders().begin() + 1, _morphologyDf.getHeaders().end());

    _morphoData->setData(morphologyMatrix, numMorphoCols);
    _morphoData->setDimensionNames(morphDimNames);

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
    Dataset<CellMorphologies> cellMorphoData = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", "cell_morphology");

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

InputDialog::InputDialog(QWidget* parent, PatchSeqDataLoader& binLoader, QString fileName) :
    QDialog(parent),
    _datasetNameAction(this, "Dataset name", fileName),
    _dataTypeAction(this, "Data type", { "Float", "Unsigned Byte" }),
    _numberOfDimensionsAction(this, "Number of dimensions", 1, 1000000, 1),
	_storeAsAction(this, "Store as"),
    _isDerivedAction(this, "Mark as derived", false),
    _datasetPickerAction(this, "Source dataset"),
    _loadAction(this, "Load"),
    _groupAction(this, "Settings")
{
    setWindowTitle(tr("Patch-seq Data Loader"));

    _numberOfDimensionsAction.setDefaultWidgetFlags(IntegralAction::WidgetFlag::SpinBox);

    QStringList pointDataTypes;
    for (const char* const typeName : PointData::getElementTypeNames())
    {
        pointDataTypes.append(QString::fromLatin1(typeName));
    }
    _storeAsAction.setOptions(pointDataTypes);

    // Load some settings
    _dataTypeAction.setCurrentIndex(binLoader.getSetting("DataType").toInt());
    _numberOfDimensionsAction.setValue(binLoader.getSetting("NumberOfDimensions").toInt());
    _storeAsAction.setCurrentIndex(binLoader.getSetting("StoreAs").toInt());

    _groupAction.addAction(&_datasetNameAction);
    _groupAction.addAction(&_dataTypeAction);
    _groupAction.addAction(&_numberOfDimensionsAction);
    _groupAction.addAction(&_storeAsAction);
    _groupAction.addAction(&_isDerivedAction);
    _groupAction.addAction(&_datasetPickerAction);
    _groupAction.addAction(&_loadAction);

    auto layout = new QVBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(_groupAction.createWidget(this));

    setLayout(layout);

    // Update the state of the dataset picker
    const auto updateDatasetPicker = [this]() -> void {
        if (_isDerivedAction.isChecked()) {

            // Get unique identifier and gui names from all point data sets in the core
            auto dataSets = mv::data().getAllDatasets(std::vector<mv::DataType> {PointType});

            // Assign found dataset(s)
            _datasetPickerAction.setDatasets(dataSets);
        }
        else {

            // Assign found dataset(s)
            _datasetPickerAction.setDatasets(mv::Datasets());
        }

        // Disable dataset picker when not marked as derived
        _datasetPickerAction.setEnabled(_isDerivedAction.isChecked());
    };

    // Populate source datasets once the dataset is marked as derived
    connect(&_isDerivedAction, &ToggleAction::toggled, this, updateDatasetPicker);

    // Update dataset picker at startup
    updateDatasetPicker();

    // Accept when the load action is triggered
    connect(&_loadAction, &TriggerAction::triggered, this, [this, &binLoader]() {

        // Save some settings
        binLoader.setSetting("DataType", _dataTypeAction.getCurrentIndex());
        binLoader.setSetting("NumberOfDimensions", _numberOfDimensionsAction.getValue());
        binLoader.setSetting("StoreAs", _storeAsAction.getCurrentIndex());

        accept();
    });
}

#include "PatchSeqDataLoader.h"

#include "Taxonomy.h"
#include "CellLoader.h"

#include "csv.h"

#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>
#include <TextData/TextData.h>
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

using namespace mv;
using namespace mv::gui;

namespace
{
    void readDataFrame(DataFrame& df, QString fileName)
    {
        std::ifstream file(fileName.toStdString());

        bool header = true;
        if (file)
        {
            std::string line;
            while (std::getline(file, line))
            {
                QString qline = QString::fromStdString(line);

                QStringList tokens = qline.split(",");
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
        }
        else
        {
            throw DataLoadException(fileName, "File was not found at location.");
        }
    }

    void readGeneExpressionDf(DataFrame& df, QString fileName, std::vector<float>& geneExpressionMatrix, mv::ModalTask& task, unsigned int& numCols)
    {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Starting read" << std::endl;
        std::vector<QString> row(3);

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
                token = strtok(NULL, ",");
                token = strtok(NULL, ",");
                token = strtok(NULL, ",");

                /* walk through other tokens */
                while (token != NULL) {
                    //printf(" %s\n", token);

                    geneExpressionRow[colIndex] = atof(token);
                    colIndex++;

                    token = strtok(NULL, ",");

                    //qDebug() << colIndex;
                }
                geneExpressionMatrix.insert(geneExpressionMatrix.end(), geneExpressionRow.begin(), geneExpressionRow.end());

                task.setSubtaskFinished(lineCount / 18);

                QApplication::processEvents();

                lineCount++;
            }

            qDebug() << "Done";

            //io::CSVReader in(fileName);
            //in.read_header(io::ignore_extra_column, "vendor", "size", "speed");
            //std::string vendor; int size; double speed;
            //while (in.read_row(vendor, size, speed)) {
            //    // do stuff with the data
            //}


            //int lineCount = 0;

            //task.setSubtasks(100);

            //while (!in.atEnd())
            //{
            //    QString line = in.readLine();

            //    lineCount++;

            //    if (lineCount > 1800) lineCount = 1800;

            //    //task.setSubtaskStarted()

            //    task.setSubtaskStarted(lineCount / 18);

            //    QStringList tokens = line.split(",");

            //    for (int i = 0; i < 3; i++)
            //    {
            //        QString token = tokens[i];
            //        token.replace("\"", "");
            //        row[i] = token;
            //    }

            //    df.getData().push_back(row);

            //    for (int col = 0; col < numCols; col++)
            //        geneExpressionRow[col] = tokens[col+3].toFloat();

            //    geneExpressionMatrix.insert(geneExpressionMatrix.end(), geneExpressionRow.begin(), geneExpressionRow.end());

            //    task.setSubtaskFinished(lineCount / 18);

            //    QApplication::processEvents();
            //}
            inputFile.close();
        }
        else
        {
            throw DataLoadException(fileName, "File was not found at location.");
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
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
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    }

    void readMorphologyDf(DataFrame& df, QString fileName, std::vector<float>& matrix, unsigned int& numCols)
    {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Starting read" << std::endl;
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
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
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
                qDebug() << "Found attributes in taxonomy CSV!";
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
    QDir dir = QFileDialog::getExistingDirectory(nullptr, tr("Open Directory"),
        "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    QStringList csvFiles = dir.entryList(QStringList() << "*.csv" << "*.CSV", QDir::Files);
    QString gexprFilePath;
    QString ephysFilePath;
    QString morphoFilePath;
    QString metadataFilePath;
    QString annotationFilePath;
    
    for (QString filePath : csvFiles)
    {
        if (filePath.contains("IDs"))
            gexprFilePath = dir.filePath(filePath);
        if (filePath.contains("ephys"))
            ephysFilePath = dir.filePath(filePath);
        if (filePath.contains("morpho"))
            morphoFilePath = dir.filePath(filePath);
        if (filePath.contains("metadata"))
            metadataFilePath = dir.filePath(filePath);
        if (filePath.contains("annotation"))
            annotationFilePath = dir.filePath(filePath);

        qDebug() << filePath;
    }

    qDebug() << gexprFilePath;
    qDebug() << ephysFilePath;
    qDebug() << morphoFilePath;
    qDebug() << metadataFilePath;

    // Don't try to load a file if the dialog was cancelled or the file name is empty
    if (metadataFilePath.isNull() || metadataFilePath.isEmpty())
        return;

    qDebug() << "Loading taxonomy";
    //Taxonomy taxonomy = Taxonomy::fromJsonFile();
    //taxonomy.printTree();
    
    DataFrame taxonomyDf;
    readDataFrame(taxonomyDf, annotationFilePath);

    qDebug() << "Loading CSV file: " << metadataFilePath;

    // Read metadata file
    DataFrame metadata;
    readDataFrame(metadata, metadataFilePath);

    _task.setEnabled(true);
    _task.setRunning();
    QCoreApplication::processEvents();

    // Read gene expression file
    DataFrame geneExpressionDf;
    std::vector<float> geneExpressionMatrix;
    unsigned int numCols;
    readGeneExpressionDf(geneExpressionDf, gexprFilePath, geneExpressionMatrix, _task, numCols);

    _task.setFinished();

    for (int i = 0; i < std::min(20, (int) geneExpressionDf.getHeaders().size()); i++)
    {
        qDebug() << geneExpressionDf.getHeaders()[i];
    }

    Dataset<Points> pointData;
    pointData = mv::data().createDataset<Points>("Points", QFileInfo(gexprFilePath).baseName());

    std::vector<QString> gexprDimNames(geneExpressionDf.getHeaders().begin() + 3, geneExpressionDf.getHeaders().end());

    pointData->setData(geneExpressionMatrix, numCols);
    pointData->setDimensionNames(gexprDimNames);

    qDebug() << "Notify data changed";
    events().notifyDatasetDataChanged(pointData);
    events().notifyDatasetDataDimensionsChanged(pointData);

    // Read electrophysiology file
    DataFrame ephysDf;
    std::vector<float> ephysMatrix;
    unsigned int numEphysCols;
    readPatchSeqDf(ephysDf, ephysFilePath, 2, ephysMatrix, numEphysCols);

    for (int i = 0; i < std::min(20, (int)ephysDf.getHeaders().size()); i++)
    {
        qDebug() << ephysDf.getHeaders()[i];
    }

    Dataset<Points> ephysData = mv::data().createDataset<Points>("Points", QFileInfo(ephysFilePath).baseName());

    std::vector<QString> ephysDimNames(ephysDf.getHeaders().begin() + 2, ephysDf.getHeaders().end());

    ephysData->setData(ephysMatrix, numEphysCols);
    ephysData->setDimensionNames(ephysDimNames);

    events().notifyDatasetDataChanged(ephysData);
    events().notifyDatasetDataDimensionsChanged(ephysData);

    // Subset and reorder the metadata
    DataFrame gexpr_metadata = DataFrame::subsetAndReorderByColumn(metadata, geneExpressionDf, "cell_id", "cell_id");

    // Subset and reorder the metadata
    DataFrame ephys_metadata = DataFrame::subsetAndReorderByColumn(metadata, ephysDf, "cell_id", "cell_id");

    // Add cluster meta data
    addTaxonomyClustersForDf(geneExpressionDf, gexpr_metadata, taxonomyDf, QFileInfo(gexprFilePath).baseName(), pointData);

    // Add cluster meta data
    addTaxonomyClustersForDf(ephysDf, ephys_metadata, taxonomyDf, QFileInfo(ephysFilePath).baseName(), ephysData);

    // Load morphology data
    DataFrame morphologyDf;
    std::vector<float> morphologyMatrix;
    unsigned int numMorphoCols;
    readMorphologyDf(morphologyDf, morphoFilePath, morphologyMatrix, numMorphoCols);

    Dataset<Points> morphoData;
    morphoData = mv::data().createDataset<Points>("Points", QFileInfo(morphoFilePath).baseName());

    std::vector<QString> morphDimNames(morphologyDf.getHeaders().begin() + 1, morphologyDf.getHeaders().end());

    morphoData->setData(morphologyMatrix, numMorphoCols);
    morphoData->setDimensionNames(morphDimNames);

    qDebug() << "Notify data changed";
    events().notifyDatasetDataChanged(morphoData);
    events().notifyDatasetDataDimensionsChanged(morphoData);

    // Subset and reorder the metadata
    DataFrame morpho_metadata = DataFrame::subsetAndReorderByColumn(metadata, morphologyDf, "cell_id", "cell_id");

    // Add cluster meta data
    addTaxonomyClustersForDf(morphologyDf, morpho_metadata, taxonomyDf, QFileInfo(morphoFilePath).baseName(), morphoData);

    //// Add linked selections between data
    //std::unordered_map<unsigned int, unsigned int> linkedData = computeLinkedData(geneExpressionDf, morphologyDf, "cell_id", "cell_id");

    //SelectionMap selectionMap;
    //auto& mapping = selectionMap.getMap();

    //for (auto kv = linkedData.begin(); kv != linkedData.end(); kv++)
    //{
    //    std::vector<unsigned int> a(1, kv->second);
    //    mapping[kv->first] = a;
    //}
    //pointData->getSelection()->addLinkedData(morphoData, selectionMap);
    //qDebug() << "add linked data";

    qDebug() << ">>>>>>>>>>>>>> Create metadata";
    // Create text data for metadata columns
    std::vector<QString> cellIdColumn = metadata["cell_id"];
    std::vector<QString> cellNameColumn = metadata["cell_name_label"];
    std::vector<QString> treeSubclassColumn = metadata["tree_subclass"];
    std::vector<QString> treeClusterColumn = metadata["tree_cluster"];
    std::vector<QString> paradigmColumn = metadata["paradigm"];

    Dataset<Text> metaDataset;
    metaDataset = mv::data().createDataset<Text>("Text", "cell_metadata");

    metaDataset->addColumn("cell_id", cellIdColumn);
    metaDataset->addColumn("cell_name_label", cellNameColumn);
    metaDataset->addColumn("tree_subclass", treeSubclassColumn);
    metaDataset->addColumn("tree_cluster", treeClusterColumn);
    metaDataset->addColumn("paradigm", paradigmColumn);

    // Take columns from ephys and morpho data and order them correctly, filling in missing data
    BiMap metadataBiMap;
    std::vector<uint32_t> metaCellIdIndices(metadata.numRows());
    std::iota(metaCellIdIndices.begin(), metaCellIdIndices.end(), 0);
    metadataBiMap.addKeyValuePairs(metadata["cell_id"], metaCellIdIndices);

    std::vector<QString> ephysCellIdColumn = ephys_metadata["cell_id"];
    std::vector<uint32_t> ephysToMetaIndices = metadataBiMap.getValuesByKeys(ephysCellIdColumn);

    std::vector<QString> morphoCellIdColumn = morpho_metadata["cell_id"];
    std::vector<uint32_t> morphoToMetaIndices = metadataBiMap.getValuesByKeys(morphoCellIdColumn);

    for (int d = 0; d < ephysData->getNumDimensions(); d++)
    {
        std::vector<QString> column(metadata.numRows(), "?");
        QString dimName = ephysData->getDimensionNames()[d];

        for (size_t i = 0; i < ephysToMetaIndices.size(); i++)
        {
            uint32_t row = ephysToMetaIndices[i];
            column[row] = QString::number(ephysData->getValueAt(i * ephysData->getNumDimensions() + d));
        }
        metaDataset->addColumn(dimName, column);
    }

    for (int d = 0; d < morphoData->getNumDimensions(); d++)
    {
        std::vector<QString> column(metadata.numRows(), "?");
        QString dimName = morphoData->getDimensionNames()[d];

        for (size_t i = 0; i < morphoToMetaIndices.size(); i++)
        {
            uint32_t row = morphoToMetaIndices[i];
            column[row] = QString::number(morphoData->getValueAt(i * morphoData->getNumDimensions() + d));
        }
        metaDataset->addColumn(dimName, column);
    }

    qDebug() << "Notify data changed";
    events().notifyDatasetDataChanged(metaDataset);
    events().notifyDatasetDataDimensionsChanged(metaDataset);

    createClusterData(metadata["tree_cluster"], "tree_cluster", metaDataset);

    qDebug() << ">>>>>>>>>>>>>> Loading morphology cells";
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

    qDebug() << ">>>>>>>>>>>>>> Making selection group";
    // Make selection group
    KeyBasedSelectionGroup selectionGroup;

    // Gene expression data
    BiMap gexprBiMap;
    std::vector<uint32_t> gexprIndices(pointData->getNumPoints());
    std::iota(gexprIndices.begin(), gexprIndices.end(), 0);
    gexprBiMap.addKeyValuePairs(gexpr_metadata["cell_id"], gexprIndices);

    // Morphology data
    BiMap morphBiMap;
    std::vector<uint32_t> morphoIndices(morphoData->getNumPoints());
    std::iota(morphoIndices.begin(), morphoIndices.end(), 0);
    morphBiMap.addKeyValuePairs(morpho_metadata["cell_id"], morphoIndices);

    // Ephys data
    BiMap ephysBiMap;
    std::vector<uint32_t> ephysIndices(ephysData->getNumPoints());
    std::iota(ephysIndices.begin(), ephysIndices.end(), 0);
    ephysBiMap.addKeyValuePairs(ephys_metadata["cell_id"], ephysIndices);

    // Cell ID data
    BiMap cellIdBiMap;
    std::vector<uint32_t> cellIdIndices(metaDataset->getNumRows());
    std::iota(cellIdIndices.begin(), cellIdIndices.end(), 0);
    cellIdBiMap.addKeyValuePairs(metaDataset->getColumn("cell_id"), cellIdIndices);

    selectionGroup.addDataset(pointData, gexprBiMap);
    selectionGroup.addDataset(morphoData, morphBiMap);
    selectionGroup.addDataset(ephysData, ephysBiMap);
    selectionGroup.addDataset(metaDataset, cellIdBiMap);

    events().addSelectionGroup(selectionGroup);
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

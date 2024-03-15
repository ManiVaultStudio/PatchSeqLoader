#include "PatchSeqDataLoader.h"

#include "DataFrame.h"
#include "Taxonomy.h"

#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>

#include <Set.h>

#include <QtCore>
#include <QtDebug>

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <iostream>
#include <algorithm>

Q_PLUGIN_METADATA(IID "studio.manivault.SparseDataLoader")

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

    void readGeneExpressionDf(DataFrame& df, QString fileName, std::vector<float>& geneExpressionMatrix, unsigned int& numCols)
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
            while (!in.atEnd())
            {
                QString line = in.readLine();

                QStringList tokens = line.split(",");

                for (int i = 0; i < 3; i++)
                {
                    QString token = tokens[i];
                    token.replace("\"", "");
                    row[i] = token;
                }

                df.getData().push_back(row);

                for (int col = 0; col < numCols; col++)
                    geneExpressionRow[col] = tokens[col+3].toFloat();

                geneExpressionMatrix.insert(geneExpressionMatrix.end(), geneExpressionRow.begin(), geneExpressionRow.end());
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

void PatchSeqDataLoader::loadData()
{
    const QString fileName = AskForFileName(tr("Metadata Files (*.csv)"));

    // Don't try to load a file if the dialog was cancelled or the file name is empty
    if (fileName.isNull() || fileName.isEmpty())
        return;

    qDebug() << "Loading taxonomy";
    Taxonomy taxonomy = Taxonomy::fromJsonFile();
    //taxonomy.printTree();

    DataFrame taxonomyDf;
    readDataFrame(taxonomyDf, "D:/Dropbox/Julian/Patchseq/FinalHumanMTGclusterAnnotation_update.csv");

    qDebug() << "Loading CSV file: " << fileName;

    // Read metadata file
    DataFrame metadata;
    readDataFrame(metadata, fileName);

    // Read gene expression file
    const QString gexprFileName = AskForFileName(tr("Gene Expression Files (*.csv)"));

    DataFrame geneExpressionDf;
    std::vector<float> geneExpressionMatrix;
    unsigned int numCols;
    readGeneExpressionDf(geneExpressionDf, gexprFileName, geneExpressionMatrix, numCols);

    for (int i = 0; i < std::min(20, (int) geneExpressionDf.getHeaders().size()); i++)
    {
        qDebug() << geneExpressionDf.getHeaders()[i];
    }

    Dataset<Points> pointData;
    pointData = mv::data().createDataset<Points>("Points", QFileInfo(fileName).baseName());

    pointData->setData(geneExpressionMatrix, numCols);

    qDebug() << "Notify data changed";
    events().notifyDatasetDataChanged(pointData);
    events().notifyDatasetDataDimensionsChanged(pointData);

    // Subset and reorder the metadata
    std::vector<QString> cell_id_column_gexpr = geneExpressionDf["cell_id"];
    std::vector<QString> cell_id_column_meta = metadata["cell_id"];

    // Make a map out of meta column
    std::unordered_map<QString, int> indexMap;
    for (int i = 0; i < cell_id_column_meta.size(); i++)
    {
        indexMap[cell_id_column_meta[i]] = i;
    }

    // Find ordering
    std::vector<int> ordering;
    for (const QString& cell_id : cell_id_column_gexpr)
    {
        if (indexMap.find(cell_id) == indexMap.end())
        {
            qDebug() << "Failed to find cell ID: " << cell_id << " in metadata file.";
        }
        int index = indexMap[cell_id];
        ordering.push_back(index);
    }

    // Subset and reorder metadata
    metadata.reorder(ordering);

    // Add cluster meta data
    std::vector<QString> treeCluster = metadata["tree_cluster"];
    for (int i = 0; i < treeCluster.size(); i++)
    {
        //qDebug() << treeCluster[i];
    }
    std::cout << "NUM CLUSTERS LIST: " << treeCluster.size();

    Dataset<Clusters> treeClusterData;
    treeClusterData = mv::data().createDataset<Points>("Cluster", QFileInfo(fileName).baseName(), pointData);

    std::unordered_map<QString, std::vector<unsigned int>> clusterData = makeClustersFromList(treeCluster);

    std::vector<QString> final_cluster_names = taxonomyDf["final cluster"];
    std::vector<QString> final_cluster_colors = taxonomyDf["final cluster color"];

    for (auto& kv : clusterData) {
        Cluster cluster;

        cluster.setName(kv.first);
        //cluster.setId(QString::fromStdString(dataContent.clusterIDs[i]));
        TaxonomyAttributes* attr = taxonomy.findLeafWithAttribute("label", kv.first.toStdString());
        if (attr != nullptr)
        {
            qDebug() << "Found attributes in taxonomy!";
            std::string color = attr->_map["nodePar.col"];
            cluster.setColor(QColor(QString::fromStdString(color)));
        }
        else
        {
            qDebug() << "ERROR: Failed to find attributes in taxonomy: " << kv.first;
        }

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
        //QColor color = { dataContent.clusterColors[i * 3], dataContent.clusterColors[i * 3 + 1] , dataContent.clusterColors[i * 3 + 2] };
        //cluster.setColor(color);

        //cluster.getIndices() = std::vector<std::uint32_t>(dataContent.clusterIndices.begin() + globalIndicesOffset, dataContent.clusterIndices.begin() + globalIndicesOffset + dataContent.clusterSizes[i]);
        //globalIndicesOffset += dataContent.clusterSizes[i];

        cluster.setIndices(kv.second);

        treeClusterData->addCluster(cluster);
    }

    mv::events().notifyDatasetDataChanged(treeClusterData);
    mv::events().notifyDatasetDataDimensionsChanged(treeClusterData);


    //for (std::vector<QString>& row : metadata)
    //{
    //    qDebug() << row[1];
    //}

    return;

    //// Read the 
    //std::ifstream file(fileName.toStdString(), std::ios::binary);

    //if (file)
    //{
    //    size_t numRows;
    //    file.read((char*)&numRows, sizeof(size_t));

    //    size_t numCols;
    //    file.read((char*)&numCols, sizeof(size_t));

    //    size_t numNonZero;
    //    file.read((char*)&numNonZero, sizeof(size_t));

    //    qDebug() << "Num rows:" << numRows;
    //    qDebug() << "Num cols:" << numCols;
    //    qDebug() << "Num nonzeros:" << numNonZero;

    //    std::vector<size_t> rowPointers(numRows + 1);
    //    std::vector<uint16_t> colIndices(numNonZero);
    //    std::vector<uint16_t> values(numNonZero);

    //    file.read((char*)rowPointers.data(), sizeof(size_t) * rowPointers.size());
    //    file.read((char*)colIndices.data(), sizeof(uint16_t) * colIndices.size());
    //    file.read((char*)values.data(), sizeof(uint16_t) * values.size());

    //    file.close();

    //    Dataset<Points> pointData;
    //    pointData = mv::data().createDataset<Points>("Points", QFileInfo(fileName).baseName());

    //    pointData->setSparseData<uint16_t, uint16_t>(numRows, numCols, rowPointers, colIndices, values);

    //    qDebug() << "Notify data changed";
    //    events().notifyDatasetDataChanged(pointData);
    //    events().notifyDatasetDataDimensionsChanged(pointData);
    //}
    //else
    //{
    //    throw DataLoadException(fileName, "File was not found at location.");
    //}

    //SparseMatrix<uint16_t, uint16_t> csr(numRows, numCols, numNonZero);
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

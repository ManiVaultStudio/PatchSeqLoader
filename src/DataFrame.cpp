#include "DataFrame.h"

#include "CSVReader.h"

#include <LoaderPlugin.h>

#include <QDebug>
#include <QFile>
#include <QTextStream>

#include <unordered_map>
#include <sstream>
#include <iostream>

DataFrame::DataFrame()
{

}

bool DataFrame::empty() const
{
    return _data.empty();
}

unsigned int DataFrame::numRows() const
{
    return _data.size();
}

unsigned int DataFrame::numCols() const
{
    return _data[0].size();
}

bool DataFrame::hasColumn(const QString& columnName) const
{
    return columnIndex(columnName) >= 0;
}

int DataFrame::columnIndex(const QString& columnName) const
{
    for (int i = 0; i < static_cast<int>(_headers.size()); i++)
    {
        if (_headers[i] == columnName)
            return i;
    }

    return -1;
}

QString DataFrame::getValue(int row, int col)
{
    return _data[row][col];
}

int DataFrame::findRowWithColumnValue(QString columnName, QString value)
{
    int col = getColumnIndex(columnName);

    for (int i = 0; i < numRows(); i++)
    {
        const QString& val = _data[i][col];
        if (val == value)
            return i;
    }
    qWarning() << "Failed to find value: " << value << " in column: " << columnName;
    return -1;
}

void DataFrame::readFromFile(QString fileName)
{
    QFile file(fileName);

    if (!file.open(QIODevice::ReadOnly))
    {
        throw mv::plugin::DataLoadException(fileName, "File was not found at location.");
    }

    QTextStream in(&file);
    std::stringstream csvStream;
    while (!in.atEnd()) {
        QString line = in.readLine();
        csvStream << line.toStdString() << "\n";
    }

    file.close();
    qDebug() << "Pre CSV";
    CSVReader reader;
    reader.LoadCSV(csvStream, _headers, _data);
}

std::vector<int> DataFrame::findDuplicateRows(QString columnToCheck)
{
    std::vector<int> duplicateRows;

    const int col = columnIndex(columnToCheck);

    if (col < 0)
    {
        qWarning() << "Cannot find duplicates. Missing column:" << columnToCheck;
        return duplicateRows;
    }

    QSet<QString> uniqueRows;

    for (int i = 0; i < static_cast<int>(_data.size()); i++)
    {
        if (col >= static_cast<int>(_data[i].size()))
            continue;

        const QString value = _data[i][col];

        if (!uniqueRows.contains(value))
            uniqueRows.insert(value);
        else
            duplicateRows.push_back(i);
    }

    return duplicateRows;
}

void DataFrame::removeRow(int rowIndex)
{
    _data.erase(_data.begin() + rowIndex);
}

void DataFrame::removeRows(const std::vector<int>& rowsToDelete)
{
    int rowsRemoved = 0;
    // Delete bad rows from both the dataframe and the matrix
    for (int rowToDelete : rowsToDelete)
    {
        rowToDelete -= rowsRemoved;
        removeRow(rowToDelete);
        rowsRemoved++;
    }
}

void DataFrame::removeDuplicateRows(QString columnToCheck)
{
    std::vector<int> duplicateRows = findDuplicateRows(columnToCheck);
    qDebug() << "Removing duplicate rows: " << duplicateRows.size();
    for (int i = 0; i < duplicateRows.size(); i++)
    {
        qDebug() << (*this)[columnToCheck][duplicateRows[i]];
    }
    removeRows(duplicateRows);
}

void DataFrame::addHeader(QString header)
{
    _headers.push_back(header);
}

void DataFrame::setHeaders(const QStringList& columnNames)
{
    for (const QString columnName : columnNames)
    {
        _headers.push_back(columnName);
    }
}

void DataFrame::reorder(std::vector<int> order)
{
    std::vector<std::vector<QString>> reorderedData;

    for (const int index : order)
    {
        reorderedData.push_back(_data[index]);
    }

    _data = reorderedData;
}

void DataFrame::subsetAndReorderAccordingTo(DataFrame& rightDf, QString columnNameLeft, QString columnNameRight)
{
    std::vector<QString> columnRight = rightDf[columnNameRight];
    std::vector<QString> columnLeft = (*this)[columnNameLeft];

    // Make a map out of meta column
    std::unordered_map<QString, int> indexMap;
    for (int i = 0; i < columnLeft.size(); i++)
    {
        indexMap[columnLeft[i]] = i;
    }

    // Find ordering
    std::vector<int> ordering;
    for (const QString& cell_id : columnRight)
    {
        auto it = indexMap.find(cell_id);

        if (it == indexMap.end())
        {
            qDebug() << "[subsetAndReorderAccordingTo] Failed to find cell ID:" << cell_id;
            continue;
        }

        ordering.push_back(it->second);
    }

    // Subset and reorder metadata
    reorder(ordering);
}

void DataFrame::printFirstFewDimensionsOfDataFrame()
{
    std::cout << "Loaded file with first 20 dimensions: ";
    for (int i = 0; i < std::min(20, (int) _headers.size()); i++)
    {
        std::cout << _headers[i].toStdString() << ", ";
    }
    std::cout << std::endl;
}

DataFrame DataFrame::subsetAndReorderByColumn(const DataFrame& leftDf, DataFrame& rightDf, QString columnNameLeft, QString columnNameRight)
{
    std::vector<QString> columnRight = rightDf[columnNameRight];
    std::vector<QString> columnLeft = leftDf[columnNameLeft];

    // Make a map out of meta column
    std::unordered_map<QString, int> indexMap;
    for (int i = 0; i < columnLeft.size(); i++)
    {
        indexMap[columnLeft[i]] = i;
    }

    // Find ordering
    std::vector<int> ordering;
    for (const QString& cell_id : columnRight)
    {
        if (indexMap.find(cell_id) == indexMap.end())
        {
            qDebug() << "[subsetAndReorderByColumn] Failed to find cell ID: " << cell_id << " in metadata file.";
            continue;
        }
        int index = indexMap[cell_id];
        ordering.push_back(index);
    }
    qDebug() << "Ordering: " << ordering.size();
    // Subset and reorder metadata
    std::vector<std::vector<QString>> reorderedData;

    for (const int index : ordering)
    {
        reorderedData.push_back(leftDf._data[index]);
    }

    DataFrame resultDf;
    resultDf._headers = leftDf._headers;
    resultDf._data = reorderedData;

    return resultDf;
}

std::vector<QString> DataFrame::operator[](QString columnName) const
{
    const int col = columnIndex(columnName);

    if (col < 0)
    {
        qWarning() << "Could not find column with name:" << columnName;
        return {};
    }

    std::vector<QString> column;
    column.reserve(_data.size());

    for (const auto& row : _data)
    {
        if (col < static_cast<int>(row.size()))
            column.push_back(row[col]);
        else
            column.push_back({});
    }

    return column;
}

int DataFrame::getColumnIndex(QString columnName) const
{
    for (int i = 0; i < _headers.size(); i++)
    {
        if (_headers[i] == columnName)
            return i;
    }

    qWarning() << "Could not find column with name: " << columnName;
    return -1;
}

#include "DataFrame.h"

#include <LoaderPlugin.h>

#include <QDebug>
#include <QFile>
#include <QTextStream>

#include <unordered_map>
#include <iostream>

DataFrame::DataFrame()
{

}

unsigned int DataFrame::numRows() const
{
    return _data.size();
}

unsigned int DataFrame::numCols() const
{
    return _data[0].size();
}

QString DataFrame::getValue(int row, int col)
{
    return _data[row][col];
}

std::vector<std::vector<QString>>& DataFrame::getData()
{
    return _data;
}

void DataFrame::readFromFile(QString fileName)
{
    QFile file(fileName);

    bool header = true;
    if (!file.open(QIODevice::ReadOnly))
    {
        throw mv::plugin::DataLoadException(fileName, "File was not found at location.");
    }

    QTextStream in(&file);

    while (!in.atEnd())
    {
        QString line = in.readLine();

        QStringList tokens = line.split(",");
        if (header)
        {
            setHeaders(tokens);
            header = false;
            continue;
        }

        std::vector<QString> row(tokens.size());
        for (int i = 0; i < tokens.size(); i++)
            row[i] = tokens[i];

        _data.push_back(row);
    }

    file.close();
}

std::vector<int> DataFrame::findDuplicateRows(QString columnToCheck)
{
    // Iterate over rows
    std::vector<QString> column = (*this)[columnToCheck];

    QSet<QString> uniqueRows;
    std::vector<int> duplicateRows;
    for (int i = 0; i < column.size(); i++)
    {
        if (!uniqueRows.contains(column[i]))
            uniqueRows.insert(column[i]);
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
        if (indexMap.find(cell_id) == indexMap.end())
        {
            qDebug() << "[subsetAndReorderAccordingTo] Failed to find cell ID: " << cell_id << " in metadata file.";
        }
        int index = indexMap[cell_id];
        ordering.push_back(index);
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
    int columnIndex = getColumnIndex(columnName);

    std::vector<QString> column;

    for (int row = 0; row < _data.size(); row++)
    {
        column.push_back(_data[row][columnIndex]);
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

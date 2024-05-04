#include "DataFrame.h"

#include <QDebug>

#include <unordered_map>

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

void DataFrame::removeRow(int rowIndex)
{
    _data.erase(_data.begin() + rowIndex);
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
        }
        int index = indexMap[cell_id];
        ordering.push_back(index);
    }

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

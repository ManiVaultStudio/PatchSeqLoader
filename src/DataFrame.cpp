#include "DataFrame.h"

#include <QDebug>

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

std::vector<QString> DataFrame::operator[](QString columnName)
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

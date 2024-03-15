#pragma once

#include <QString>
#include <QStringList>

#include <vector>

class DataFrame
{
public:
    DataFrame();

    unsigned int numRows() const;
    unsigned int numCols() const;
    QString getValue(int row, int col);
    const std::vector<QString>& getHeaders() { return _headers; }
    std::vector<std::vector<QString>>& getData();

    void setHeaders(const QStringList& columnNames);
    void reorder(std::vector<int> order);

    std::vector<QString> operator[](QString columnName);

private:
    int getColumnIndex(QString columnName) const;

private:
    std::vector<QString> _headers;

    std::vector<std::vector<QString>> _data;
};

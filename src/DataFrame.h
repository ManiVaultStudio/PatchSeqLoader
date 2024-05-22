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
    const std::vector<QString>& getHeaders() const { return _headers; }
    std::vector<std::vector<QString>>& getData();

    void addRow(std::vector<QString>& row) { _data.push_back(row); }
    std::vector<int> findDuplicateRows(QString columnToCheck);
    void removeRow(int rowIndex);
    void removeRows(const std::vector<int>& rowsToDelete);
    void removeDuplicateRows(QString columnToCheck);

    void addHeader(QString header);
    void setHeaders(const QStringList& columnNames);
    void reorder(std::vector<int> order);
    void subsetAndReorderAccordingTo(DataFrame& rightDf, QString columnNameLeft, QString columnNameRight);

    void printFirstFewDimensionsOfDataFrame();

    static DataFrame subsetAndReorderByColumn(const DataFrame& leftDf, DataFrame& rightDf, QString columnNameLeft, QString columnNameRight);

    std::vector<QString> operator[](QString columnName) const;

private:
    int getColumnIndex(QString columnName) const;

private:
    std::vector<QString> _headers;

    std::vector<std::vector<QString>> _data;
};

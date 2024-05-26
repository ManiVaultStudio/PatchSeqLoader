#pragma once

#include <QString>

#include <vector>

// Magic number that represents a missing value, to be imputed
constexpr float MISSING_VALUE = 1234567.0f;

class MatrixData
{
public:
    void removeRow(int row);
    void removeRows(const std::vector<int>& rowsToDelete);

    void fillMissingValues(float fillValue);
    void imputeMissingValues();
    void standardize();

    std::vector<float> operator[](QString columnName) const;

private:
    int getColumnIndex(QString columnName) const;

public:
    std::vector<QString> headers;
    std::vector<float> data;
    size_t numRows;
    size_t numCols;
};

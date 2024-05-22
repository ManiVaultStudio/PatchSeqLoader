#pragma once

#include <QString>

#include <vector>

// Magic number that represents a missing value, to be imputed
constexpr float MISSING_VALUE = 1234567.0f;

class MatrixData
{
public:
    void removeRow(int row);

    void fillMissingValues(float fillValue);
    void imputeMissingValues();

    std::vector<QString> headers;
    std::vector<float> data;
    size_t numRows;
    size_t numCols;
};

#pragma once

#include "DataFrame.h"

#include <QString>

#include <vector>

namespace mv
{
    class ModalTask;
}

// Magic number that represents a missing value, to be imputed
constexpr float MISSING_VALUE = 1234567.0f;

class MatrixData
{
public:
    void removeRow(int row)
    {
        data.erase(data.begin() + (row * numCols + 0), data.begin() + (row * numCols + numCols));
        numRows--;
    }

    std::vector<QString> headers;
    std::vector<float> data;
    size_t numRows;
    size_t numCols;
};

class MatrixDataLoader
{
public:
    MatrixDataLoader(bool handleMissingValues = false) :
        _handleMissingValues(handleMissingValues)
    {

    }

    void LoadMatrixData(QString fileName, DataFrame& df, MatrixData& matrix, int numMetaCols);

private:
    bool _handleMissingValues = false;
};

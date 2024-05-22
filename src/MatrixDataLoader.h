#pragma once

#include "DataFrame.h"

namespace mv
{
    class ModalTask;
}

class MatrixData;

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

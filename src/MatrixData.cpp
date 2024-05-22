#include "MatrixData.h"

void MatrixData::removeRow(int row)
{
    data.erase(data.begin() + (row * numCols + 0), data.begin() + (row * numCols + numCols));
    numRows--;
}

void MatrixData::fillMissingValues(float fillValue)
{
    // Compute means, ignoring missing values
    for (int col = 0; col < numCols; col++)
    {
        // Set the imputed value for the rows with missing values
        for (int row = 0; row < numRows; row++)
        {
            float v = data[row * numCols + col];
            if (v == MISSING_VALUE)
                data[row * numCols + col] = fillValue;
        }
    }
}

void MatrixData::imputeMissingValues()
{
    // Compute means, ignoring missing values
    for (int col = 0; col < numCols; col++)
    {
        // Compute the mean
        float mean = 0;
        int numContributions = 0;
        for (int row = 0; row < numCols; row++)
        {
            float v = data[row * numCols + col];
            if (v == MISSING_VALUE)
                continue;
            mean += v;
            numContributions++;
        }
        if (numContributions > 0)
            mean /= numContributions;

        // Set the imputed value for the rows with missing values
        for (int row = 0; row < numRows; row++)
        {
            float v = data[row * numCols + col];
            if (v == MISSING_VALUE)
                data[row * numCols + col] = mean;
        }
    }
}

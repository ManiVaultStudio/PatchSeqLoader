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

namespace
{
    // Function to compute the mean of each column
    std::vector<float> computeColumnMeans(const std::vector<float>& matrix, int rows, int cols) {
        std::vector<float> means(cols, 0.0f);
        for (int c = 0; c < cols; ++c) {
            float sum = 0.0f;
            for (int r = 0; r < rows; ++r) {
                sum += matrix[r * cols + c];
            }
            means[c] = sum / rows;
        }
        return means;
    }

    // Function to compute the standard deviation of each column
    std::vector<float> computeColumnStdDevs(const std::vector<float>& matrix, const std::vector<float>& means, int rows, int cols) {
        std::vector<float> stdDevs(cols, 0.0f);
        for (int c = 0; c < cols; ++c) {
            float sumSqDiffs = 0.0f;
            for (int r = 0; r < rows; ++r) {
                float diff = matrix[r * cols + c] - means[c];
                sumSqDiffs += diff * diff;
            }
            stdDevs[c] = std::sqrt(sumSqDiffs / rows);
        }
        return stdDevs;
    }
}

void MatrixData::standardize()
{
    std::vector<float> means = computeColumnMeans(data, numRows, numCols);
    std::vector<float> stdDevs = computeColumnStdDevs(data, means, numRows, numCols);

    for (int c = 0; c < numCols; ++c) {
        for (int r = 0; r < numRows; ++r) {
            int index = r * numCols + c;

            if (stdDevs[c] != 0) {
                data[index] = (data[index] - means[c]) / stdDevs[c];
            }
            else {
                data[index] = data[index] - means[c]; // Avoid division by zero
            }
        }
    }
}

#include "MatrixData.h"

#include <QDebug>

#include <unordered_set>

void MatrixData::removeRow(int row)
{
    data.erase(data.begin() + (row * numCols + 0), data.begin() + (row * numCols + numCols));
    numRows--;
}

void MatrixData::removeRows(const std::vector<int>& rowsToDelete)
{
    int rowsRemoved = 0;
    // Delete bad rows from both the dataframe and the matrix
    for (int rowToDelete : rowsToDelete)
    {
        rowToDelete -= rowsRemoved;
        removeRow(rowToDelete);
        rowsRemoved++;
    }
}

void MatrixData::removeCols(const std::vector<int>& colsToDelete)
{
    if (colsToDelete.empty())
        return;

    std::unordered_set<int> deleteSet(colsToDelete.begin(), colsToDelete.end());

    std::vector<int> colsToKeep;
    for (int i = 0; i < numCols; i++)
    {
        if (!deleteSet.count(i))
            colsToKeep.push_back(i);
    }

    // Only include cols to be kept in the new data
    std::vector<float> newData;
    for (int row = 0; row < numRows; row++)
    {
        for (int j = 0; j < colsToKeep.size(); j++)
        {
            int col = colsToKeep[j];
            newData.push_back(data[row * numCols + col]);
        }
    }

    // Only include headers to be kept
    std::vector<QString> newHeaders;
    for (int j = 0; j < colsToKeep.size(); j++)
    {
        int col = colsToKeep[j];
        newHeaders.push_back(headers[col]);
    }

    headers = newHeaders;
    data = newData;
    numCols = colsToKeep.size();
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

int MatrixData::getColumnIndex(QString columnName) const
{
    for (int i = 0; i < headers.size(); i++)
    {
        if (headers[i] == columnName)
            return i;
    }

    qWarning() << "Could not find column with name: " << columnName;
    return -1;
}

std::vector<float> MatrixData::operator[](QString columnName) const
{
    int columnIndex = getColumnIndex(columnName);

    std::vector<float> column;

    for (int row = 0; row < numRows; row++)
    {
        column.push_back(data[row * numCols + columnIndex]);
    }

    return column;
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

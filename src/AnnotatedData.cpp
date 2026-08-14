#include "AnnotatedData.h"

#include "AnnotatedDataPrinter.h"

#include <iostream>

bool AnnotationTable::HasColumn(QString columnName)
{
    return std::find(columnNames.begin(), columnNames.end(), columnName) != columnNames.end();
}

void AnnotationTable::RemoveRows(const std::vector<size_t>& rowsToDelete)
{
    if (rowsToDelete.empty())
        return;

    // AnnotationTable::values is column-major, so index.size() is the row count.
    const size_t rowCount = index.size();

    // Mark rows for deletion first, so row indices are not invalidated while erasing.
    std::vector<char> remove(rowCount, false);

    for (size_t row : rowsToDelete)
    {
        // Ignore invalid row indices defensively.
        if (row < rowCount)
            remove[row] = true;
    }

    // Erase entries whose row index was marked.
    // This works for both the index vector and each column in values.
    auto eraseMarkedRows = [&](std::vector<QString>& rows)
        {
            size_t row = 0;

            rows.erase(
                std::remove_if(
                    rows.begin(),
                    rows.end(),
                    [&](const QString&)
                    {
                        // If a column is unexpectedly longer than index, keep extra entries.
                        if (row >= remove.size())
                            return false;

                        return remove[row++] != 0;
                    }),
                rows.end());
        };

    // Remove rows from the table index.
    eraseMarkedRows(index);

    // Remove the same row positions from every metadata column.
    for (std::vector<QString>& column : values)
        eraseMarkedRows(column);
}

void NumericMatrix::RemoveRows(const std::vector<size_t>& rowsToDelete)
{
    // Create rows-length vector marking which rows should be removed
    std::vector<char> remove(rowCount, false);
    for (size_t row : rowsToDelete)
        remove[row] = true;

    // In one pass compact the value matrix
    size_t writeIndex = 0;
    for (size_t row = 0; row < rowCount; row++)
    {
        if (!remove[row])
        {
            if (writeIndex != row)
            {
                std::move(values.begin() + row * columnCount,
                          values.begin() + (row + 1) * columnCount,
                          values.begin() + writeIndex * columnCount);
            }
            writeIndex++;
        }
    }

    // Update to the new size of the matrix and row count
    values.resize(writeIndex * columnCount);
    rowCount = writeIndex;
}

void AnnotatedData::RemoveRows(const std::vector<size_t>& rowsToDelete)
{
    obs.RemoveRows(rowsToDelete);

    X.RemoveRows(rowsToDelete);
}

void AnnotatedData::Print(size_t maxRows, size_t maxColumns) const
{
    AnnotatedDataPrinter printer(std::cout);
    printer.Print(*this, maxRows, maxColumns);
}

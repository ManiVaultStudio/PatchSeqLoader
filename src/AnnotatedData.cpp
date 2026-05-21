#include "AnnotatedData.h"

#include "AnnotatedDataPrinter.h"

#include <iostream>

bool AnnotationTable::HasColumn(QString columnName)
{
    return std::find(columnNames.begin(), columnNames.end(), columnName) != columnNames.end();
}

void AnnotationTable::RemoveRows(const std::vector<size_t>& rowsToDelete)
{
    // Delete rows from index
    for (size_t row : rowsToDelete)
        index.erase(index.begin() + row);

    // Delete rows from values
    size_t rowCount = values.size();

    std::vector<char> remove(rowCount, false);

    for (std::size_t row : rowsToDelete)
        remove[row] = true;

    std::size_t row_index = 0;

    values.erase(std::remove_if(values.begin(), values.end(),
            [&](const std::vector<QString>&) {
                return remove[row_index++] != 0;
            }
        ), values.end()
    );
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

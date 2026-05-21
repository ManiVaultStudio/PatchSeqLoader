#include "AnnotatedDataPrinter.h"

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <sstream>

AnnotatedDataPrinter::AnnotatedDataPrinter(std::ostream& output)
    : output_(output)
{
}

void AnnotatedDataPrinter::Print(const AnnotatedData& data, size_t maxRows, size_t maxColumns) const
{
    PrintSummary(data);

    output_ << '\n';
    PrintAnnotationTable("obs", data.obs, maxRows, maxColumns);

    output_ << '\n';
    PrintAnnotationTable("var", data.var, maxRows, maxColumns);

    output_ << '\n';
    PrintMatrixPreview(data, maxRows, maxColumns);

    output_ << '\n';
    output_ << "* = imputed value\n";
}

void AnnotatedDataPrinter::PrintSummary(const AnnotatedData& data) const
{
    output_ << "AnnotatedData object\n";
    output_ << "------------------------------------------------------------\n";

    output_ << "Shape: "
        << data.X.rowCount << " observations x "
        << data.X.columnCount << " variables\n";

    output_ << "obs: "
        << data.obs.index.size() << " rows x "
        << data.obs.columnNames.size() << " columns\n";

    output_ << "var: "
        << data.var.index.size() << " rows x "
        << data.var.columnNames.size() << " columns\n";

    output_ << "X: "
        << data.X.rowCount << " rows x "
        << data.X.columnCount << " columns\n";
}

void AnnotatedDataPrinter::PrintAnnotationTable(const char* name, const AnnotationTable& table, size_t maxRows, size_t maxColumns) const
{
    output_ << name << " preview:\n";

    const size_t rowsToPrint = std::min(maxRows, table.index.size());
    const size_t columnsToPrint = std::min(maxColumns, table.columnNames.size());

    const QString indexHeader = table.indexName.isEmpty() ? QString("index") : table.indexName;

    output_ << std::setw(CELL_WIDTH) << FormatString(indexHeader);

    for (size_t col = 0; col < columnsToPrint; ++col)
        output_ << std::setw(CELL_WIDTH) << FormatString(table.columnNames[col]);

    if (table.columnNames.size() > columnsToPrint)
        output_ << std::setw(CELL_WIDTH) << "...";

    output_ << '\n';

    for (size_t row = 0; row < rowsToPrint; ++row)
    {
        output_ << std::setw(CELL_WIDTH) << FormatString(table.index[row]);

        for (size_t col = 0; col < columnsToPrint; ++col)
        {
            QString value;

            if (row < table.values.size() && col < table.values[row].size())
                value = table.values[row][col];

            output_ << std::setw(CELL_WIDTH) << FormatString(value);
        }

        if (table.columnNames.size() > columnsToPrint)
            output_ << std::setw(CELL_WIDTH) << "...";

        output_ << '\n';
    }

    if (table.index.size() > rowsToPrint)
    {
        output_ << std::setw(CELL_WIDTH) << "...";
        output_ << '\n';
    }
}

void AnnotatedDataPrinter::PrintMatrixPreview(const AnnotatedData& data, size_t maxRows, size_t maxColumns) const
{
    output_ << "X preview:\n";

    const NumericMatrix& X = data.X;

    const size_t rowsToPrint = std::min(maxRows, X.rowCount);
    const size_t columnsToPrint = std::min(maxColumns, X.columnCount);

    output_ << std::setw(CELL_WIDTH) << "";

    for (size_t col = 0; col < columnsToPrint; ++col)
    {
        QString columnName;

        if (col < data.var.index.size())
            columnName = data.var.index[col];
        else
            columnName = QString("var_%1").arg(col);

        output_ << std::setw(CELL_WIDTH)
            << FormatString(columnName);
    }

    if (X.columnCount > columnsToPrint)
        output_ << std::setw(CELL_WIDTH) << "...";

    output_ << '\n';

    for (size_t row = 0; row < rowsToPrint; ++row)
    {
        QString rowName;

        if (row < data.obs.index.size())
            rowName = data.obs.index[row];
        else
            rowName = QString("obs_%1").arg(row);

        output_ << std::setw(CELL_WIDTH) << FormatString(rowName);

        for (size_t col = 0; col < columnsToPrint; ++col)
        {
            const size_t flatIndex = row * X.columnCount + col;

            if (flatIndex >= X.values.size())
            {
                output_ << std::setw(CELL_WIDTH) << "<missing>";
                continue;
            }

            const bool imputed = flatIndex < X.imputed.size() && X.imputed[flatIndex] != 0;

            output_ << std::setw(CELL_WIDTH) << FormatFloat(X.values[flatIndex], imputed);
        }

        if (X.columnCount > columnsToPrint)
            output_ << std::setw(CELL_WIDTH) << "...";

        output_ << '\n';
    }

    if (X.rowCount > rowsToPrint)
    {
        output_ << std::setw(CELL_WIDTH) << "...";
        output_ << '\n';
    }
}

std::string AnnotatedDataPrinter::FormatString(const QString& value) const
{
    std::string text = value.toStdString();

    if (text.size() <= MAX_TEXT_LENGTH)
        return text;

    if (MAX_TEXT_LENGTH <= 3)
        return text.substr(0, MAX_TEXT_LENGTH);

    return text.substr(0, MAX_TEXT_LENGTH - 3) + "...";
}

std::string AnnotatedDataPrinter::FormatFloat(float value, bool imputed) const
{
    std::ostringstream stream;

    stream << std::setprecision(5) << value;

    if (imputed)
        stream << '*';

    return stream.str();
}

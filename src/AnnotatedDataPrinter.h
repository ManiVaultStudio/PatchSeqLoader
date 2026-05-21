#pragma once

#include "AnnotatedData.h"

#include <cstddef>
#include <iosfwd>
#include <string>

class AnnotatedDataPrinter
{
public:
    explicit AnnotatedDataPrinter(std::ostream& output);

    void Print(const AnnotatedData& data,
        size_t maxRows = 5,
        size_t maxColumns = 6) const;

private:
    void PrintSummary(const AnnotatedData& data) const;

    void PrintAnnotationTable(const char* name,
        const AnnotationTable& table,
        size_t maxRows,
        size_t maxColumns) const;

    void PrintMatrixPreview(const AnnotatedData& data,
        size_t maxRows,
        size_t maxColumns) const;

    std::string FormatString(const QString& value) const;
    std::string FormatFloat(float value, bool imputed) const;

private:
    std::ostream& output_;

    static constexpr size_t CELL_WIDTH = 14;
    static constexpr size_t MAX_TEXT_LENGTH = 13;
};

#include "CSVReader.h"

#include "csv.hpp"

void CSVReader::LoadCSV(QString filePath, std::vector<QString>& headers, std::vector<std::vector<QString>>& data)
{
    csv::CSVFormat format;
    format.delimiter(',').quote('"').header_row(0);

    csv::CSVReader reader(filePath.toStdString(), format);

    // Capture headers (optional, if the CSV has them)
    headers.clear();
    for (const auto& col : reader.get_col_names()) {
        headers.push_back(QString::fromStdString(col));
    }

    data.clear();
    for (csv::CSVRow& row : reader)
    {
        std::vector<QString> dataRow;
        for (csv::CSVField& field : row) {
            dataRow.push_back(QString::fromStdString(field.get<>()));
        }
        data.push_back(dataRow);
    }
}

void CSVReader::LoadCSV(std::stringstream& sstream, std::vector<QString>& headers, std::vector<std::vector<QString>>& data)
{
    csv::CSVFormat format;
    format.delimiter(',').quote('"').header_row(0);

    csv::CSVReader reader(sstream, format);

    // Capture headers (optional, if the CSV has them)
    headers.clear();
    for (const auto& col : reader.get_col_names()) {
        headers.push_back(QString::fromStdString(col));
    }

    data.clear();
    for (csv::CSVRow& row : reader)
    {
        std::vector<QString> dataRow;
        for (csv::CSVField& field : row) {
            dataRow.push_back(QString::fromStdString(field.get<>()));
        }
        data.push_back(dataRow);
    }
}

#pragma once

#include "AnnotatedData.h"

#include "Config/Config.h"

#include <QString>

#include <vector>

enum class CsvColumnKind
{
    Index,
    Metadata,
    Data
};

class CsvLoadContext
{
public:
    std::vector<CsvColumnKind> columnKinds;
    bool allColumnsAreMetadata = false;
};

//class CsvData
//{
//public:
//    std::vector<QString> metadataHeaders;
//    std::vector<QString> dataHeaders;
//
//    std::vector<CsvColumnKind> columnKinds;
//
//    // Row-major metadata storage
//    std::vector<std::vector<QString>> metadata;
//    // Row-major data storage
//    std::vector<float> data;
//
//    std::vector<uint8_t> imputed;
//
//    // Number of rows
//    size_t numRows;
//
//    // Number of data dimensions
//    size_t numDataDimensions;
//};

class CsvLoader
{
public:
    void Load(const QString& filePath, const QStringList& metadataHeaders, const QString& index, AnnotatedData& dataFile);
    void Load(const config::TableSource& config, AnnotatedData& dataFile);
private:

};

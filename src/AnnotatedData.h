#pragma once

#include <vector>

#include <QString>
#include <QStringList>
#include <QHash>

// Magic number that represents a missing value, to be imputed
constexpr float MISSING_VALUE = 1234567.0f;

class AnnotationTable
{
public:
    bool HasColumn(QString columnName);

    void RemoveRows(const std::vector<size_t>& rowsToDelete);

public:
    QString indexName;
    std::vector<QString> index;

    std::vector<QString> columnNames;
    // Column-major storage of observation annotations
    std::vector<std::vector<QString>> values;
};

class NumericMatrix
{
public:
    void RemoveRows(const std::vector<size_t>& rowsToDelete);

public:
    size_t rowCount;
    size_t columnCount;

    // Row-major [row * columnCount + col]
    std::vector<float> values;

    // Same shape as values
    // 0 = observed/original value
    // 1 = missing value that was imputed
    std::vector<uint8_t> imputed;
};

struct NamedNumericMatrix
{
    std::vector<QString> columnNames;
    NumericMatrix matrix;
};

class AnnotatedData
{
public:
    void RemoveRows(const std::vector<size_t>& rowsToDelete);
    void Print(size_t maxRows = 5, size_t maxColumns = 6) const;

public:
    AnnotationTable obs;
    QHash<QString, NamedNumericMatrix> obsm;
    AnnotationTable var;
    NumericMatrix X;
private:

};

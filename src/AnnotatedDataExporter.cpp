#include "AnnotatedDataExporter.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>

#include <stdexcept>

void AnnotatedDataExporter::Export(const AnnotatedData& data, const QString& filePath)
{
    const AnnotationTable& obs = data.obs;

    ValidateObs(obs);

    qDebug()
        << "Exporting obs:"
        << "indexName =" << obs.indexName
        << "rows =" << obs.index.size()
        << "columns =" << obs.columnNames.size()
        << "values columns =" << obs.values.size();

    QFile file(filePath);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        throw std::runtime_error(
            QString("Failed to open CSV file for writing: %1. Reason: %2")
                .arg(filePath, file.errorString())
                .toStdString()
        );
    }

    QTextStream out(&file);

    // Header
    const QString indexHeader =
        obs.indexName.isEmpty() ? QString("index") : obs.indexName;

    out << EscapeCsvField(indexHeader);

    for (const QString& columnName : obs.columnNames)
    {
        out << ',';
        out << EscapeCsvField(columnName);
    }

    out << '\n';

    // Rows
    for (size_t row = 0; row < obs.index.size(); ++row)
    {
        out << EscapeCsvField(obs.index[row]);

        for (size_t col = 0; col < obs.columnNames.size(); ++col)
        {
            out << ',';

            QString value;

            // Column-major: values[col][row]
            if (col < obs.values.size() && row < obs.values[col].size())
                value = obs.values[col][row];

            out << EscapeCsvField(value);
        }

        out << '\n';
    }
}

void AnnotatedDataExporter::ValidateObs(const AnnotationTable& obs)
{
    if (obs.values.size() != obs.columnNames.size())
    {
        throw std::runtime_error(
            QString("Cannot export obs: values column count (%1) does not match columnNames count (%2).")
                .arg(obs.values.size())
                .arg(obs.columnNames.size())
                .toStdString()
        );
    }

    for (size_t col = 0; col < obs.values.size(); ++col)
    {
        if (obs.values[col].size() != obs.index.size())
        {
            throw std::runtime_error(
                QString("Cannot export obs: column %1 has %2 values, expected %3 rows.")
                    .arg(col)
                    .arg(obs.values[col].size())
                    .arg(obs.index.size())
                    .toStdString()
            );
        }
    }
}

QString AnnotatedDataExporter::EscapeCsvField(const QString& value)
{
    QString escaped = value;

    const bool needsQuotes =
        escaped.contains(',') ||
        escaped.contains('"') ||
        escaped.contains('\n') ||
        escaped.contains('\r');

    escaped.replace("\"", "\"\"");

    if (needsQuotes)
        return "\"" + escaped + "\"";

    return escaped;
}

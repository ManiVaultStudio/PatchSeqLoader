#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "AnnotatedDataExporter.h"

#include <QSet>

#include <algorithm>
#include <stdexcept>
#include <array>

namespace
{
    bool TryParseFloat(const QString& text, float& value)
    {
        const QString trimmed = text.trimmed();

        if (trimmed.isEmpty())
            return false;

        bool ok = false;
        value = trimmed.toFloat(&ok);

        return ok;
    }

    bool IsNumericObsColumn(const std::vector<QString>& column,
        size_t maxValuesToCheck = 5)
    {
        size_t checked = 0;

        for (const QString& value : column)
        {
            if (value.trimmed().isEmpty())
                continue;

            float parsed = 0.0f;

            if (!TryParseFloat(value, parsed))
                return false;

            ++checked;

            if (checked >= maxValuesToCheck)
                break;
        }

        return checked > 0;
    }

    static constexpr const char* COLOR_SUFFIX = "_color";

    static bool IsColorColumn(const QString& columnName)
    {
        return columnName.endsWith(COLOR_SUFFIX);
    }

    static QString BaseColumnNameForColorColumn(const QString& colorColumnName)
    {
        QString baseName = colorColumnName;

        if (baseName.endsWith(COLOR_SUFFIX))
            baseName.chop(QString(COLOR_SUFFIX).size());

        return baseName;
    }

    static bool TryParseColor(const QString& text, QColor& color)
    {
        const QString trimmed = text.trimmed();

        if (trimmed.isEmpty())
            return false;

        QColor parsed(trimmed);

        if (!parsed.isValid())
            return false;

        color = parsed;
        return true;
    }
}

class CollectMetadataStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "CollectMetadata";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        AnnotatedData combined;
        combined.obs.indexName = "cell_id";

        QHash<QString, size_t> rowByIndex;
        QHash<QString, size_t> columnByName;

        for (const AnnotatedData& table : ctx.normalizedTables)
            MergeObsTable(table, combined, rowByIndex, columnByName);

        ctx.metadataColorMaps.clear();
        ExtractColorMaps(combined.obs, ctx.metadataColorMaps);
        RemoveColorColumns(combined.obs);
        AssignDefaultBinaryMetadataColors(combined.obs, ctx.metadataColorMaps);

        combined.X.rowCount = combined.obs.index.size();
        combined.X.columnCount = 0;
        combined.X.values.clear();
        combined.X.imputed.clear();

        combined.var.indexName.clear();
        combined.var.index.clear();
        combined.var.columnNames.clear();
        combined.var.values.clear();

        SplitNumericObs(combined);

        ctx.metadata = std::move(combined);

        AnnotatedDataExporter::Export(ctx.metadata, "metadata_export.csv");

        qDebug() << "Collected metadata table with" << ctx.metadata.obs.index.size() << "rows and" << ctx.metadata.obs.columnNames.size() << "columns.";

        return true;
    }

private:
    static void MergeObsTable(const AnnotatedData& source, AnnotatedData& target, QHash<QString, size_t>& rowByIndex, QHash<QString, size_t>& columnByName)
    {
        ValidateAnnotationTable(source.obs);

        if (target.obs.indexName.isEmpty())
            target.obs.indexName = source.obs.indexName;
        else if (!source.obs.indexName.isEmpty() && target.obs.indexName != source.obs.indexName)
            qWarning() << "Merging obs tables with different index names:" << target.obs.indexName << "and" << source.obs.indexName;

        EnsureColumnsExist(source, target, columnByName);

        for (size_t sourceRow = 0; sourceRow < source.obs.index.size(); ++sourceRow)
        {
            const QString& indexValue = source.obs.index[sourceRow];

            if (indexValue.isEmpty())
            {
                qWarning() << "Skipping obs row with empty index.";
                continue;
            }

            const size_t targetRow = EnsureRowExists(indexValue, target, rowByIndex);

            MergeObsRow(source, sourceRow, target, targetRow, columnByName);
        }
    }

    static void ValidateAnnotationTable(const AnnotationTable& table)
    {
        if (table.values.size() != table.columnNames.size())
        {
            throw std::runtime_error(
                QString("Invalid AnnotationTable: values column count (%1) does not match columnNames count (%2).")
                .arg(table.values.size())
                .arg(table.columnNames.size())
                .toStdString()
            );
        }

        for (size_t col = 0; col < table.values.size(); ++col)
        {
            if (table.values[col].size() != table.index.size())
            {
                throw std::runtime_error(
                    QString("Invalid AnnotationTable: column %1 has %2 values, expected %3 rows.")
                    .arg(col)
                    .arg(table.values[col].size())
                    .arg(table.index.size())
                    .toStdString()
                );
            }
        }
    }

    static void EnsureColumnsExist(const AnnotatedData& source, AnnotatedData& target, QHash<QString, size_t>& columnByName)
    {
        for (const QString& columnName : source.obs.columnNames)
        {
            if (columnByName.contains(columnName))
                continue;

            const size_t newColumnIndex = target.obs.columnNames.size();

            target.obs.columnNames.push_back(columnName);
            columnByName.insert(columnName, newColumnIndex);

            // Column-major: one vector per column, with one value per existing row.
            target.obs.values.emplace_back(target.obs.index.size());
        }
    }

    static size_t EnsureRowExists(const QString& indexValue, AnnotatedData& target, QHash<QString, size_t>& rowByIndex)
    {
        if (rowByIndex.contains(indexValue))
            return rowByIndex.value(indexValue);

        const size_t newRowIndex = target.obs.index.size();

        target.obs.index.push_back(indexValue);

        // Column-major: append one cell to each metadata column.
        for (std::vector<QString>& column : target.obs.values)
            column.push_back(QString());

        rowByIndex.insert(indexValue, newRowIndex);

        return newRowIndex;
    }

    static void MergeObsRow(const AnnotatedData& source, size_t sourceRow, AnnotatedData& target, size_t targetRow, const QHash<QString, size_t>& columnByName)
    {
        for (size_t sourceCol = 0; sourceCol < source.obs.columnNames.size(); ++sourceCol)
        {
            const QString& columnName = source.obs.columnNames[sourceCol];

            if (sourceCol >= source.obs.values.size() || sourceRow >= source.obs.values[sourceCol].size())
            {
                qWarning() << "Skipping invalid metadata cell for column" << columnName << "row" << sourceRow;
                continue;
            }

            const QString& incomingValue = source.obs.values[sourceCol][sourceRow];

            const size_t targetCol = columnByName.value(columnName);
            QString& existingValue = target.obs.values[targetCol][targetRow];

            if (existingValue.isEmpty())
            {
                existingValue = incomingValue;
            }
            else if (!incomingValue.isEmpty() && existingValue != incomingValue)
            {
                qWarning() << "Conflicting metadata value for index" << target.obs.index[targetRow]
                    << "column" << columnName
                    << ". Keeping existing value:" << existingValue
                    << "ignoring incoming value:" << incomingValue;
            }
        }
    }

    static void SplitNumericObs(AnnotatedData& data)
    {
        constexpr const char* obsmKey = "numeric_obs";

        if (data.obs.values.size() != data.obs.columnNames.size())
        {
            throw std::runtime_error(
                QString("Cannot split numeric obs: obs.values has %1 columns, expected %2.")
                .arg(data.obs.values.size())
                .arg(data.obs.columnNames.size())
                .toStdString()
            );
        }

        std::vector<size_t> numericColumns;

        for (size_t col = 0; col < data.obs.values.size(); ++col)
        {
            const std::vector<QString>& column = data.obs.values[col];

            if (column.size() != data.obs.index.size())
            {
                throw std::runtime_error(
                    QString("Cannot split numeric obs: column %1 has %2 rows, expected %3.")
                    .arg(col)
                    .arg(column.size())
                    .arg(data.obs.index.size())
                    .toStdString()
                );
            }

            if (IsNumericObsColumn(column, 5))
                numericColumns.push_back(col);
        }

        if (numericColumns.empty())
            return;

        NamedNumericMatrix numericObs;
        for (size_t sourceCol : numericColumns)
        {
            numericObs.columnNames.push_back(data.obs.columnNames[sourceCol]);
        }
        numericObs.matrix.rowCount = data.obs.index.size();
        numericObs.matrix.columnCount = numericColumns.size();
        numericObs.matrix.values.reserve(numericObs.matrix.rowCount * numericObs.matrix.columnCount);
        numericObs.matrix.imputed.reserve(numericObs.matrix.rowCount * numericObs.matrix.columnCount);

        // Build obsm matrix row-major.
        for (size_t row = 0; row < numericObs.matrix.rowCount; ++row)
        {
            for (size_t outputCol = 0; outputCol < numericColumns.size(); ++outputCol)
            {
                const size_t sourceCol = numericColumns[outputCol];
                const QString& text = data.obs.values[sourceCol][row];

                float value = 0.0f;

                if (TryParseFloat(text, value))
                {
                    numericObs.matrix.values.push_back(value);
                    numericObs.matrix.imputed.push_back(0);
                }
                else
                {
                    numericObs.matrix.values.push_back(MISSING_VALUE);
                    numericObs.matrix.imputed.push_back(1);
                }
            }
        }

        data.obsm.insert(obsmKey, std::move(numericObs));

        // Remove numeric columns from obs, preserving column-major layout.
        std::vector<QString> newColumnNames;
        std::vector<std::vector<QString>> newValues;

        newColumnNames.reserve(data.obs.columnNames.size() - numericColumns.size());
        newValues.reserve(data.obs.values.size() - numericColumns.size());

        size_t numericCursor = 0;

        for (size_t col = 0; col < data.obs.columnNames.size(); ++col)
        {
            const bool isNumericColumn =
                numericCursor < numericColumns.size() &&
                numericColumns[numericCursor] == col;

            if (isNumericColumn)
            {
                ++numericCursor;
                continue;
            }

            newColumnNames.push_back(std::move(data.obs.columnNames[col]));
            newValues.push_back(std::move(data.obs.values[col]));
        }

        data.obs.columnNames = std::move(newColumnNames);
        data.obs.values = std::move(newValues);
    }

    static QHash<QString, size_t> BuildColumnIndexByName(const AnnotationTable& obs)
    {
        QHash<QString, size_t> columnByName;

        for (size_t col = 0; col < obs.columnNames.size(); ++col)
            columnByName.insert(obs.columnNames[col], col);

        return columnByName;
    }

    static void ExtractColorMaps(const AnnotationTable& obs, QHash<QString, QHash<QString, QColor>>& outputColorMaps)
    {
        const QHash<QString, size_t> columnByName = BuildColumnIndexByName(obs);

        for (size_t colorCol = 0; colorCol < obs.columnNames.size(); ++colorCol)
        {
            const QString& colorColumnName = obs.columnNames[colorCol];

            if (!IsColorColumn(colorColumnName))
                continue;

            const QString baseColumnName =
                BaseColumnNameForColorColumn(colorColumnName);

            if (!columnByName.contains(baseColumnName))
            {
                qWarning() << "Color column" << colorColumnName << "has no matching metadata column" << baseColumnName << "; ignoring it.";
                continue;
            }

            const size_t labelCol = columnByName.value(baseColumnName);

            const std::vector<QString>& labels = obs.values[labelCol];
            const std::vector<QString>& colors = obs.values[colorCol];

            QHash<QString, QColor>& colorMap = outputColorMaps[baseColumnName];

            const size_t rowCount =
                std::min(labels.size(), colors.size());

            for (size_t row = 0; row < rowCount; ++row)
            {
                const QString label = labels[row].trimmed();
                const QString colorText = colors[row].trimmed();

                if (label.isEmpty() || colorText.isEmpty())
                    continue;

                QColor color;

                if (!TryParseColor(colorText, color))
                {
                    qWarning() << "Invalid color" << colorText << "for metadata column" << baseColumnName << "label" << label;
                    continue;
                }

                if (colorMap.contains(label) && colorMap.value(label) != color)
                {
                    qWarning() << "Conflicting colors for metadata column" << baseColumnName << "label" << label << ". Keeping existing color" << colorMap.value(label).name() << "and ignoring" << color.name();
                    continue;
                }

                colorMap.insert(label, color);
            }

            qDebug() << "Extracted" << colorMap.size() << "colors for metadata column" << baseColumnName;
        }
    }

    enum class BooleanLikeValue
    {
        NotBooleanLike,
        FalseLike,
        TrueLike
    };

    static QString NormalizeBooleanLikeText(const QString& text)
    {
        QString normalized = text.trimmed().toLower();
        normalized.remove(' ');
        normalized.remove('_');
        normalized.remove('-');
        return normalized;
    }

    static BooleanLikeValue ClassifyBooleanLikeValue(const QString& text)
    {
        const QString normalized = NormalizeBooleanLikeText(text);

        if (normalized.isEmpty())
            return BooleanLikeValue::NotBooleanLike;

        static const QSet<QString> trueLikeValues = {
            "true",
            "t",
            "1",
            "yes",
            "y",
            "on",
            "present",
            "positive",
            "pos",
            "enabled",
            "enable",
            "included",
            "include",
            "valid",
            "pass",
            "passed"
        };

        static const QSet<QString> falseLikeValues = {
            "false",
            "f",
            "0",
            "no",
            "n",
            "off",
            "absent",
            "negative",
            "neg",
            "disabled",
            "disable",
            "excluded",
            "exclude",
            "invalid",
            "fail",
            "failed"
        };

        if (trueLikeValues.contains(normalized))
            return BooleanLikeValue::TrueLike;

        if (falseLikeValues.contains(normalized))
            return BooleanLikeValue::FalseLike;

        return BooleanLikeValue::NotBooleanLike;
    }

    static bool TryAssignBooleanSemanticColors(
        const QString& columnName,
        const QStringList& labels,
        QHash<QString, QColor>& colorMap)
    {
        if (labels.size() != 2)
            return false;

        const BooleanLikeValue first = ClassifyBooleanLikeValue(labels[0]);
        const BooleanLikeValue second = ClassifyBooleanLikeValue(labels[1]);

        const bool isBooleanPair =
            (first == BooleanLikeValue::TrueLike && second == BooleanLikeValue::FalseLike) ||
            (first == BooleanLikeValue::FalseLike && second == BooleanLikeValue::TrueLike);

        if (!isBooleanPair)
            return false;

        constexpr const char* trueColor = "#1984A3";
        constexpr const char* falseColor = "#F5C767";

        for (const QString& label : labels)
        {
            const BooleanLikeValue classification = ClassifyBooleanLikeValue(label);

            if (classification == BooleanLikeValue::TrueLike)
                colorMap.insert(label, QColor(trueColor));
            else if (classification == BooleanLikeValue::FalseLike)
                colorMap.insert(label, QColor(falseColor));
        }

        qDebug() << "Assigned boolean-semantic metadata colors for column" << columnName
            << ": true-like ->" << trueColor
            << ", false-like ->" << falseColor;

        return true;
    }

    static void AssignDefaultBinaryMetadataColors(
        const AnnotationTable& obs,
        QHash<QString, QHash<QString, QColor>>& colorMaps)
    {
        constexpr const char* defaultFirstColor = "#1984A3";
        constexpr const char* defaultSecondColor = "#F5C767";

        const std::array<QColor, 2> binaryColors = {
            QColor(defaultFirstColor),
            QColor(defaultSecondColor)
        };

        for (size_t col = 0; col < obs.columnNames.size(); ++col)
        {
            const QString& columnName = obs.columnNames[col];

            if (IsColorColumn(columnName))
                continue;

            if (col >= obs.values.size())
                continue;

            const std::vector<QString>& values = obs.values[col];

            QStringList labels;
            QSet<QString> seen;

            for (const QString& value : values)
            {
                const QString label = value.trimmed();

                if (label.isEmpty() || seen.contains(label))
                    continue;

                labels.push_back(label);
                seen.insert(label);

                if (labels.size() > 2)
                    break;
            }

            if (labels.size() != 2)
                continue;

            QHash<QString, QColor>& colorMap = colorMaps[columnName];

            if (TryAssignBooleanSemanticColors(columnName, labels, colorMap))
                continue;

            for (int i = 0; i < labels.size(); ++i)
            {
                const QString& label = labels[i];

                if (colorMap.contains(label))
                    continue;

                colorMap.insert(label, binaryColors[static_cast<size_t>(i)]);
            }

            qDebug() << "Assigned default binary metadata colors for column" << columnName
                << ":" << labels[0] << binaryColors[0].name()
                << "," << labels[1] << binaryColors[1].name();
        }
    }

    static void RemoveColorColumns(AnnotationTable& obs)
    {
        std::vector<QString> newColumnNames;
        std::vector<std::vector<QString>> newValues;

        newColumnNames.reserve(obs.columnNames.size());
        newValues.reserve(obs.values.size());

        for (size_t col = 0; col < obs.columnNames.size(); ++col)
        {
            if (IsColorColumn(obs.columnNames[col]))
                continue;

            newColumnNames.push_back(std::move(obs.columnNames[col]));
            newValues.push_back(std::move(obs.values[col]));
        }

        obs.columnNames = std::move(newColumnNames);
        obs.values = std::move(newValues);
    }
};

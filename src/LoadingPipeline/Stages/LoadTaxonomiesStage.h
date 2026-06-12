#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"
#include "AnnotatedData.h"

#include <QColor>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QSet>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QDebug>

#include <vector>

class LoadTaxonomiesStage final : public PipelineStage
{
public:
    QString Name() const override { return "LoadTaxonomies"; }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage:" << Name();

        if (ctx.config.taxonomies.isEmpty())
        {
            ctx.result.Info(Name(), "No taxonomy files configured; skipping taxonomy loading.");
            return true;
        }

        if (ctx.metadata.obs.index.empty())
        {
            ctx.result.Warning(Name(), "Metadata table is empty; taxonomy files cannot be applied.");
            return true;
        }

        for (auto it = ctx.config.taxonomies.begin(); it != ctx.config.taxonomies.end(); ++it)
        {
            const QString taxonomyKey = it.key();
            const config::TaxonomySource& source = it.value();

            if (!LoadAndApplyTaxonomy(ctx, taxonomyKey, source))
                return false;
        }

        return true;
    }

private:
    struct TaxonomyLevelColumns
    {
        QString levelName;
        int labelColumn = -1;
        int colorColumn = -1;
    };

    struct TaxonomyData
    {
        QString leafLevel;
        QStringList levels;

        // leaf label -> level name -> label at that level
        QHash<QString, QHash<QString, QString>> labelsByLeaf;

        // level name -> label -> color
        QHash<QString, QHash<QString, QColor>> colorMaps;
    };

    static bool LoadAndApplyTaxonomy(PipelineContext& ctx, const QString& taxonomyKey, const config::TaxonomySource& source)
    {
        if (source.path.trimmed().isEmpty())
        {
            ctx.result.Warning("LoadTaxonomies", QString("Taxonomy '%1' has no path; skipping it.").arg(taxonomyKey));
            return true;
        }

        TaxonomyData taxonomy;
        if (!ReadTaxonomyCsv(ctx, taxonomyKey, source, taxonomy))
            return false;

        ApplyTaxonomyToMetadata(ctx, taxonomyKey, source, taxonomy);
        MergeTaxonomyColorMaps(ctx, taxonomyKey, source, taxonomy);
        PropagateAliasColorsFromMatchedMetadata(ctx, taxonomyKey, source, taxonomy);

        ctx.result.Info(
            "LoadTaxonomies",
            QString("Loaded taxonomy '%1'.").arg(DisplayName(taxonomyKey, source)),
            QString("levels=%1, entries=%2").arg(taxonomy.levels.join(", ")).arg(taxonomy.labelsByLeaf.size()));

        return true;
    }

    static bool ReadTaxonomyCsv(PipelineContext& ctx, const QString& taxonomyKey, const config::TaxonomySource& source, TaxonomyData& taxonomy)
    {
        QFile file(source.path);

        if (!file.exists())
        {
            ctx.result.Error(
                "LoadTaxonomies",
                QString("Taxonomy file for '%1' does not exist.").arg(DisplayName(taxonomyKey, source)),
                source.path);
            return false;
        }

        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            ctx.result.Error(
                "LoadTaxonomies",
                QString("Could not open taxonomy file for '%1'.").arg(DisplayName(taxonomyKey, source)),
                source.path);
            return false;
        }

        QTextStream stream(&file);

        if (stream.atEnd())
        {
            ctx.result.Error(
                "LoadTaxonomies",
                QString("Taxonomy file for '%1' is empty.").arg(DisplayName(taxonomyKey, source)),
                source.path);
            return false;
        }

        const QStringList headers = SplitCsvLine(stream.readLine());
        const QVector<TaxonomyLevelColumns> levelColumns = InferLevelColumns(ctx, taxonomyKey, source, headers);

        if (levelColumns.isEmpty())
        {
            ctx.result.Error(
                "LoadTaxonomies",
                QString("Taxonomy file for '%1' has no non-color taxonomy columns.").arg(DisplayName(taxonomyKey, source)),
                source.path);
            return false;
        }

        taxonomy.leafLevel = levelColumns.front().levelName;

        for (const TaxonomyLevelColumns& level : levelColumns)
            taxonomy.levels.push_back(level.levelName);

        int rowNumber = 1;

        while (!stream.atEnd())
        {
            ++rowNumber;

            const QString line = stream.readLine();
            if (line.trimmed().isEmpty())
                continue;

            const QStringList fields = SplitCsvLine(line);
            const QString leafLabel = Field(fields, levelColumns.front().labelColumn).trimmed();

            if (leafLabel.isEmpty())
            {
                ctx.result.Warning(
                    "LoadTaxonomies",
                    QString("Skipping taxonomy row with empty leaf label in '%1'.").arg(DisplayName(taxonomyKey, source)),
                    QString("row=%1").arg(rowNumber));
                continue;
            }

            QHash<QString, QString>& labelsByLevel = taxonomy.labelsByLeaf[leafLabel];

            for (const TaxonomyLevelColumns& level : levelColumns)
            {
                const QString label = Field(fields, level.labelColumn).trimmed();

                if (label.isEmpty())
                    continue;

                const QString existingLabel = labelsByLevel.value(level.levelName);
                if (!existingLabel.isEmpty() && existingLabel != label)
                {
                    ctx.result.Warning(
                        "LoadTaxonomies",
                        QString("Conflicting taxonomy value for leaf '%1' at level '%2'. Keeping the first value.")
                        .arg(leafLabel, level.levelName),
                        QString("existing=%1, incoming=%2, row=%3").arg(existingLabel, label).arg(rowNumber));
                }
                else
                {
                    labelsByLevel[level.levelName] = label;
                }

                if (level.colorColumn >= 0)
                {
                    const QString colorText = Field(fields, level.colorColumn).trimmed();
                    if (!colorText.isEmpty())
                        AddColor(ctx, taxonomyKey, source, taxonomy, level.levelName, label, colorText, rowNumber);
                }
            }
        }

        return true;
    }

    static QVector<TaxonomyLevelColumns> InferLevelColumns(
        PipelineContext& ctx,
        const QString& taxonomyKey,
        const config::TaxonomySource& source,
        const QStringList& headers)
    {
        QVector<TaxonomyLevelColumns> levels;
        QHash<QString, int> headerIndexByLowerName;

        for (int i = 0; i < headers.size(); ++i)
            headerIndexByLowerName[headers[i].trimmed().toLower()] = i;

        for (int i = 0; i < headers.size(); ++i)
        {
            const QString header = headers[i].trimmed();
            if (header.isEmpty() || IsColorColumn(header))
                continue;

            TaxonomyLevelColumns level;
            level.levelName = header;
            level.labelColumn = i;
            level.colorColumn = FindColorColumnForLevel(header, headerIndexByLowerName);

            levels.push_back(level);
        }

        QSet<QString> matchedColorColumns;
        for (const TaxonomyLevelColumns& level : levels)
        {
            if (level.colorColumn >= 0)
                matchedColorColumns.insert(headers[level.colorColumn].trimmed().toLower());
        }

        for (const QString& header : headers)
        {
            const QString trimmed = header.trimmed();
            if (!IsColorColumn(trimmed))
                continue;

            if (!matchedColorColumns.contains(trimmed.toLower()))
            {
                ctx.result.Warning(
                    "LoadTaxonomies",
                    QString("Taxonomy '%1' has color column '%2' but no matching taxonomy label column.")
                    .arg(DisplayName(taxonomyKey, source), trimmed));
            }
        }

        return levels;
    }

    static int FindColorColumnForLevel(const QString& levelName, const QHash<QString, int>& headerIndexByLowerName)
    {
        const QString normalizedLevel = levelName.trimmed().toLower();

        const QStringList candidates = {
            QString("color_hex_%1").arg(normalizedLevel),
            QString("%1_color").arg(normalizedLevel),
            QString("color_%1").arg(normalizedLevel)
        };

        for (const QString& candidate : candidates)
        {
            const auto it = headerIndexByLowerName.find(candidate);
            if (it != headerIndexByLowerName.end())
                return it.value();
        }

        return -1;
    }

    static void AddColor(
        PipelineContext& ctx,
        const QString& taxonomyKey,
        const config::TaxonomySource& source,
        TaxonomyData& taxonomy,
        const QString& levelName,
        const QString& label,
        const QString& colorText,
        int rowNumber)
    {
        const QColor color(colorText);

        if (!color.isValid())
        {
            ctx.result.Warning(
                "LoadTaxonomies",
                QString("Invalid color '%1' in taxonomy '%2'.").arg(colorText, DisplayName(taxonomyKey, source)),
                QString("level=%1, label=%2, row=%3").arg(levelName, label).arg(rowNumber));
            return;
        }

        QHash<QString, QColor>& colorMap = taxonomy.colorMaps[levelName];

        if (colorMap.contains(label) && colorMap[label] != color)
        {
            ctx.result.Warning(
                "LoadTaxonomies",
                QString("Conflicting colors for taxonomy label '%1' at level '%2'. Keeping the first color.")
                .arg(label, levelName),
                QString("existing=%1, incoming=%2, row=%3")
                .arg(colorMap[label].name(), color.name()).arg(rowNumber));
            return;
        }

        colorMap[label] = color;
    }

    static void ApplyTaxonomyToMetadata(
        PipelineContext& ctx,
        const QString& taxonomyKey,
        const config::TaxonomySource& source,
        const TaxonomyData& taxonomy)
    {
        AnnotationTable& obs = ctx.metadata.obs;

        if (taxonomy.levels.isEmpty())
            return;

        const QString leafLevel = taxonomy.levels.front();
        const int leafColumn = FindColumnForTaxonomyLevel(obs, leafLevel);

        if (leafColumn < 0)
        {
            ctx.result.Warning(
                "LoadTaxonomies",
                QString("Could not apply taxonomy '%1': metadata has no leaf column matching '%2', '%2_name', or '%2_label'.")
                .arg(DisplayName(taxonomyKey, source), leafLevel));
            return;
        }

        const QString leafColumnName = obs.columnNames[leafColumn];
        const QString preferredAliasSuffix = AliasSuffixForLevelColumn(leafColumnName, leafLevel);

        QHash<QString, int> metadataColumnByLevel;
        for (const QString& level : taxonomy.levels)
        {
            int column = FindColumnForTaxonomyLevel(obs, level, preferredAliasSuffix);
            if (column < 0)
                column = AddColumn(obs, ColumnNameForLevelAlias(level, preferredAliasSuffix));

            metadataColumnByLevel[level] = column;
        }

        QSet<QString> warnedMissingLabels;

        for (size_t row = 0; row < obs.index.size(); ++row)
        {
            const QString leafLabel = obs.values[leafColumn][row].trimmed();

            if (leafLabel.isEmpty())
                continue;

            const auto leafIt = taxonomy.labelsByLeaf.find(leafLabel);
            if (leafIt == taxonomy.labelsByLeaf.end())
            {
                if (!warnedMissingLabels.contains(leafLabel))
                {
                    warnedMissingLabels.insert(leafLabel);
                    ctx.result.Warning(
                        "LoadTaxonomies",
                        QString("Taxonomy '%1' has no entry for metadata label '%2'.")
                        .arg(DisplayName(taxonomyKey, source), leafLabel));
                }
                continue;
            }

            const QHash<QString, QString>& labelsByLevel = leafIt.value();

            for (auto labelIt = labelsByLevel.begin(); labelIt != labelsByLevel.end(); ++labelIt)
            {
                const QString level = labelIt.key();
                const QString incomingValue = labelIt.value();
                const int metadataColumn = metadataColumnByLevel.value(level, -1);

                if (metadataColumn < 0 || incomingValue.isEmpty())
                    continue;

                QString& metadataValue = obs.values[metadataColumn][row];

                if (metadataValue.trimmed().isEmpty())
                {
                    metadataValue = incomingValue;
                }
                else if (metadataValue != incomingValue)
                {
                    ctx.result.Warning(
                        "LoadTaxonomies",
                        QString("Taxonomy '%1' conflicts with existing metadata for cell '%2'. Keeping existing metadata value.")
                        .arg(DisplayName(taxonomyKey, source), obs.index[row]),
                        QString("column=%1, existing=%2, taxonomy=%3")
                        .arg(level, metadataValue, incomingValue));
                }
            }
        }
    }

    static void MergeTaxonomyColorMaps(
        PipelineContext& ctx,
        const QString& taxonomyKey,
        const config::TaxonomySource& source,
        const TaxonomyData& taxonomy)
    {
        for (auto levelIt = taxonomy.colorMaps.begin(); levelIt != taxonomy.colorMaps.end(); ++levelIt)
        {
            const QString levelName = levelIt.key();
            const QHash<QString, QColor>& taxonomyColors = levelIt.value();

            for (const QString& metadataColumnName : TaxonomyLevelAliases(levelName))
            {
                QHash<QString, QColor>& metadataColors = ctx.metadataColorMaps[metadataColumnName];

                for (auto colorIt = taxonomyColors.begin(); colorIt != taxonomyColors.end(); ++colorIt)
                {
                    const QString label = colorIt.key();
                    const QColor color = colorIt.value();

                    if (metadataColors.contains(label) && metadataColors[label] != color)
                    {
                        ctx.result.Warning(
                            "LoadTaxonomies",
                            QString("Taxonomy '%1' overrides an existing metadata color for '%2' label '%3'.")
                            .arg(DisplayName(taxonomyKey, source), metadataColumnName, label),
                            QString("metadata=%1, taxonomy=%2").arg(metadataColors[label].name(), color.name()));
                    }

                    // Taxonomy colors are authoritative. They intentionally overwrite
                    // colors extracted earlier from metadata *_color columns.
                    metadataColors[label] = color;
                }
            }
        }
    }

    static void PropagateAliasColorsFromMatchedMetadata(
        PipelineContext& ctx,
        const QString& taxonomyKey,
        const config::TaxonomySource& source,
        const TaxonomyData& taxonomy)
    {
        AnnotationTable& obs = ctx.metadata.obs;

        for (const QString& levelName : taxonomy.levels)
        {
            const auto taxonomyColorsIt = taxonomy.colorMaps.find(levelName);
            if (taxonomyColorsIt == taxonomy.colorMaps.end())
                continue;

            const QHash<QString, QColor>& taxonomyColors = taxonomyColorsIt.value();
            const QStringList aliases = TaxonomyLevelAliases(levelName);

            QVector<int> aliasColumns;
            QVector<QString> aliasColumnNames;

            for (const QString& alias : aliases)
            {
                const int column = FindColumn(obs, alias);
                if (column >= 0)
                {
                    aliasColumns.push_back(column);
                    aliasColumnNames.push_back(alias);
                }
            }

            if (aliasColumns.size() < 2)
                continue;

            for (size_t row = 0; row < obs.index.size(); ++row)
            {
                QColor matchedColor;
                QString matchedColumnName;
                QString matchedValue;

                for (int i = 0; i < aliasColumns.size(); ++i)
                {
                    const QString value = obs.values[aliasColumns[i]][row].trimmed();
                    if (value.isEmpty())
                        continue;

                    const auto colorIt = taxonomyColors.find(value);
                    if (colorIt == taxonomyColors.end())
                        continue;

                    matchedColor = colorIt.value();
                    matchedColumnName = aliasColumnNames[i];
                    matchedValue = value;
                    break;
                }

                if (!matchedColor.isValid())
                    continue;

                for (int i = 0; i < aliasColumns.size(); ++i)
                {
                    const QString aliasValue = obs.values[aliasColumns[i]][row].trimmed();
                    if (aliasValue.isEmpty())
                        continue;

                    QHash<QString, QColor>& metadataColors = ctx.metadataColorMaps[aliasColumnNames[i]];

                    if (metadataColors.contains(aliasValue) && metadataColors[aliasValue] != matchedColor)
                    {
                        ctx.result.Warning(
                            "LoadTaxonomies",
                            QString("Taxonomy '%1' colors alias metadata value '%2' in column '%3' from matching %4='%5'.")
                            .arg(DisplayName(taxonomyKey, source), aliasValue, aliasColumnNames[i], matchedColumnName, matchedValue),
                            QString("metadata=%1, taxonomy=%2").arg(metadataColors[aliasValue].name(), matchedColor.name()));
                    }

                    // If one alias column matches the taxonomy label, all other alias
                    // values on the same row represent the same category and should
                    // receive the same authoritative taxonomy color.
                    metadataColors[aliasValue] = matchedColor;
                }
            }
        }
    }

    static QStringList SplitCsvLine(const QString& line)
    {
        QStringList fields;
        QString field;
        bool inQuotes = false;

        for (int i = 0; i < line.size(); ++i)
        {
            const QChar c = line[i];

            if (c == '"')
            {
                if (inQuotes && i + 1 < line.size() && line[i + 1] == '"')
                {
                    field.append('"');
                    ++i;
                }
                else
                {
                    inQuotes = !inQuotes;
                }
            }
            else if (c == ',' && !inQuotes)
            {
                fields.push_back(field);
                field.clear();
            }
            else
            {
                field.append(c);
            }
        }

        fields.push_back(field);
        return fields;
    }

    static QString Field(const QStringList& fields, int index)
    {
        if (index < 0 || index >= fields.size())
            return {};

        return fields[index];
    }

    static bool IsColorColumn(const QString& columnName)
    {
        const QString lower = columnName.trimmed().toLower();
        return lower.startsWith("color_hex_") || lower.endsWith("_color") || lower.startsWith("color_");
    }

    static int FindColumn(const AnnotationTable& table, const QString& columnName)
    {
        for (int i = 0; i < static_cast<int>(table.columnNames.size()); ++i)
        {
            if (table.columnNames[i] == columnName)
                return i;
        }

        return -1;
    }

    static QStringList TaxonomyLevelAliases(const QString& levelName)
    {
        return {
            levelName,
            QString("%1_name").arg(levelName),
            QString("%1_label").arg(levelName)
        };
    }

    static QString AliasSuffixForLevelColumn(const QString& columnName, const QString& levelName)
    {
        if (columnName == QString("%1_name").arg(levelName))
            return "_name";

        if (columnName == QString("%1_label").arg(levelName))
            return "_label";

        return {};
    }

    static QString ColumnNameForLevelAlias(const QString& levelName, const QString& aliasSuffix)
    {
        return aliasSuffix.isEmpty() ? levelName : QString("%1%2").arg(levelName, aliasSuffix);
    }

    static int FindColumnForTaxonomyLevel(
        const AnnotationTable& table,
        const QString& levelName,
        const QString& preferredAliasSuffix = {})
    {
        if (!preferredAliasSuffix.isEmpty())
        {
            const int preferredColumn = FindColumn(table, ColumnNameForLevelAlias(levelName, preferredAliasSuffix));
            if (preferredColumn >= 0)
                return preferredColumn;
        }

        for (const QString& alias : TaxonomyLevelAliases(levelName))
        {
            const int column = FindColumn(table, alias);
            if (column >= 0)
                return column;
        }

        return -1;
    }

    static int AddColumn(AnnotationTable& table, const QString& columnName)
    {
        table.columnNames.push_back(columnName);

        std::vector<QString> values;
        values.resize(table.index.size());
        table.values.push_back(std::move(values));

        return static_cast<int>(table.columnNames.size()) - 1;
    }

    static QString DisplayName(const QString& taxonomyKey, const config::TaxonomySource& source)
    {
        return source.displayName.trimmed().isEmpty() ? taxonomyKey : source.displayName;
    }
};

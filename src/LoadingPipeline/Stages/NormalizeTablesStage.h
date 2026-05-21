#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "AnnotatedData.h"

#include <QDir>

namespace
{
    void TrimWhitespace(AnnotationTable& table)
    {
        for (QString& header : table.columnNames)
            header = header.trimmed();

        for (QString& index : table.index)
            index = index.trimmed();

        for (auto& row : table.values)
        {
            for (QString& value : row)
                value = value.trimmed();
        }
    }

    size_t RemoveEmptyIndexRows(AnnotatedData& data)
    {
        std::vector<size_t> rowsToDelete;

        for (size_t row = 0; row < data.obs.index.size(); row++)
        {
            if (data.obs.index[row].isEmpty())
            {
                rowsToDelete.push_back(row);
            }

            data.RemoveRows(rowsToDelete);

            return rowsToDelete.size();
        }
    }

    std::vector<size_t> FindDuplicateIndexRows(const AnnotatedData& data)
    {
        std::vector<size_t> duplicateRows;
        QSet<QString> seen;

        for (size_t row = 0; row < data.obs.index.size(); row++)
        {
            const QString& value = data.obs.index[row];

            if (seen.contains(value))
                duplicateRows.push_back(row);
            else
                seen.insert(value);
        }

        return duplicateRows;
    }

    struct NormalizationStats
    {
        int emptyIndexRowsRemoved = 0;
        int duplicateRowsRemoved = 0;
    };

    bool NormalizeIndexedTable(AnnotatedData& data, const QString& tableName, const QString& indexColumn, PipelineResult& result, const QString& stageName, NormalizationStats& stats)
    {
        TrimWhitespace(data.obs);
        TrimWhitespace(data.var);

        const int emptyRemoved = RemoveEmptyIndexRows(data);

        if (emptyRemoved > 0)
        {
            result.Warning(stageName, "Removed rows with empty index values.", QString("%1: removed %2 rows with empty '%3'")
                .arg(tableName)
                .arg(emptyRemoved)
                .arg(indexColumn)
            );

            stats.emptyIndexRowsRemoved += emptyRemoved;
        }

        const std::vector<size_t> duplicates = FindDuplicateIndexRows(data);

        if (!duplicates.empty())
        {
            data.RemoveRows(duplicates);

            const int duplicateCount = static_cast<int>(duplicates.size());

            result.Warning(stageName, "Removed duplicate index rows; kept first occurrence.", QString("%1: removed %2 duplicate rows based on '%3'")
                .arg(tableName)
                .arg(duplicateCount)
                .arg(indexColumn)
            );

            stats.duplicateRowsRemoved += duplicateCount;
        }

        return true;
    }
}


class NormalizeTablesStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "NormalizeTables";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        ctx.normalizedTables.clear();

        int normalizedCount = 0;
        int skippedCount = 0;

        NormalizationStats stats;

        NormalizeOptionalTable(ctx, config::keys::sources::Rna , ctx.config.rna, normalizedCount, skippedCount, stats);
        NormalizeOptionalTable(ctx, config::keys::sources::Ephys, ctx.config.ephys, normalizedCount, skippedCount, stats);
        NormalizeOptionalTable(ctx, config::keys::sources::Morphology, ctx.config.morphology, normalizedCount, skippedCount, stats);
        NormalizeOptionalTable(ctx, config::keys::sources::Metadata, ctx.config.metadata, normalizedCount, skippedCount, stats);
        NormalizeOptionalTable(ctx, config::keys::embeddings::RnaUmap, ctx.config.rnaUmap, normalizedCount, skippedCount, stats);
        NormalizeOptionalTable(ctx, config::keys::embeddings::EphysUmap, ctx.config.ephysUmap, normalizedCount, skippedCount, stats);
        NormalizeOptionalTable(ctx, config::keys::embeddings::MorphoUmap, ctx.config.morphoUmap, normalizedCount, skippedCount, stats);

        if (normalizedCount == 0)
        {
            ctx.result.Warning(Name(), "No tables were normalized.", "No raw tables were available from the previous stage.");
        }

        ctx.result.Info(Name(), "Table normalization complete.", QString("normalized=%1, skipped=%2, empty_index_rows_removed=%3, duplicate_rows_removed=%4")
            .arg(normalizedCount)
            .arg(skippedCount)
            .arg(stats.emptyIndexRowsRemoved)
            .arg(stats.duplicateRowsRemoved)
        );

        return true;
    }

    void NormalizeOptionalTable(PipelineContext& ctx, const QString& sourceName, const std::optional<config::TableSource>& source, int& normalizedCount, int& skippedCount, NormalizationStats& stats) const
    {
        if (!source)
        {
            ++skippedCount;
            return;
        }

        if (!ctx.rawTables.contains(sourceName))
        {
            ctx.result.Info(Name(), "Raw table is not loaded; skipping normalization.", sourceName);

            ++skippedCount;
            return;
        }

        AnnotatedData data = ctx.rawTables[sourceName];

        const bool ok = NormalizeIndexedTable(data, sourceName, source->index, ctx.result, Name(), stats);

        if (!ok)
        {
            ++skippedCount;
            return;
        }

        ctx.normalizedTables.insert(sourceName, std::move(data));
        ++normalizedCount;
    }

};

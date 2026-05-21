#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "CSVLoader.h"

#include <QDir>

class LoadTablesStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "LoadTables";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        int loadedCount = 0;
        int skippedCount = 0;
        int failedCount = 0;
        qDebug() << Name() << "rna";
        LoadOptionalTable(ctx, config::keys::sources::Rna, ctx.config.rna, loadedCount, skippedCount, failedCount);
        qDebug() << Name() << "ephys";
        LoadOptionalTable(ctx, config::keys::sources::Ephys, ctx.config.ephys, loadedCount, skippedCount, failedCount);
        qDebug() << Name() << "morphology";
        LoadOptionalTable(ctx, config::keys::sources::Morphology , ctx.config.morphology, loadedCount, skippedCount, failedCount);
        qDebug() << Name() << "metadata";
        LoadOptionalTable(ctx, config::keys::sources::Metadata, ctx.config.metadata, loadedCount, skippedCount, failedCount);
        qDebug() << Name() << "rna_umap";
        LoadOptionalTable(ctx, config::keys::embeddings::RnaUmap, ctx.config.rnaUmap, loadedCount, skippedCount, failedCount);
        qDebug() << Name() << "ephys_umap";
        LoadOptionalTable(ctx, config::keys::embeddings::EphysUmap, ctx.config.ephysUmap, loadedCount, skippedCount, failedCount);
        qDebug() << Name() << "morpho_umap";
        LoadOptionalTable(ctx, config::keys::embeddings::MorphoUmap, ctx.config.morphoUmap, loadedCount, skippedCount, failedCount);

        if (loadedCount == 0)
        {
            ctx.result.Warning(Name(), "No source tables were loaded.", "The config may contain only assets/embeddings, or configured source files may be missing.");
        }

        ctx.result.Info(Name(), "Table loading complete.", QString("loaded=%1, skipped=%2, failed=%3")
            .arg(loadedCount)
            .arg(skippedCount)
            .arg(failedCount)
        );

        return true;
    }

private:
    void LoadOptionalTable(PipelineContext& ctx, const QString& sourceName, const std::optional<config::TableSource>& source, int& loadedCount, int& skippedCount, int& failedCount) const
    {
        if (!source)
        {
            ctx.result.Info(Name(), "Optional table source is not configured.", sourceName);
            ++skippedCount;
            return;
        }

        if (source->path.isEmpty())
        {
            ctx.result.Warning(Name(), "Configured table source has an empty path.", sourceName);
            ++failedCount;
            return;
        }

        const QFileInfo fileInfo(source->path);

        if (!fileInfo.exists())
        {
            ctx.result.Warning(Name(), "Configured table file does not exist; skipping table.", QString("%1 = %2").arg(sourceName, source->path));
            ++skippedCount;
            return;
        }

        if (!fileInfo.isFile())
        {
            ctx.result.Warning(Name(), "Configured table path is not a file; skipping table.", QString("%1 = %2").arg(sourceName, source->path));
            ++skippedCount;
            return;
        }

        if (!fileInfo.isReadable())
        {
            ctx.result.Warning(Name(), "Configured table file is not readable; skipping table.", QString("%1 = %2").arg(sourceName, source->path));
            ++skippedCount;
            return;
        }

        AnnotatedData data;

        try
        {
            CsvLoader csvLoader;
            csvLoader.Load(source.value(), data);
        }
        catch (const std::exception& e)
        {
            ctx.result.Warning(Name(), "Failed to read table file.", QString("%1 = %2; error=%3")
                .arg(sourceName, source->path, e.what())
            );
            ++failedCount;
            return;
        }
        catch (...)
        {
            ctx.result.Warning(Name(), "Failed to read table file due to an unknown error.", QString("%1 = %2").arg(sourceName, source->path));
            ++failedCount;
            return;
        }

        //if (data.X.rowCount == 0 || data.X.columnCount == 0)
        //{
        //    ctx.result.Warning(Name(), "Table is empty after loading; skipping table.", QString("%1 = %2").arg(sourceName, source->path));
        //    ++failedCount;
        //    return;
        //}
        
        if (data.obs.indexName != source->index)
        {
            ctx.result.Warning(Name(), "Loaded table does not contain configured index column.", QString("%1: missing column '%2' in %3")
                .arg(sourceName, source->index, source->path)
            );
            ++failedCount;
            return;
        }

        ctx.rawTables.insert(sourceName, std::move(data));

        ctx.result.Info(Name(), "Loaded table.", QString("%1 from %2 using index column '%3'")
            .arg(sourceName, source->path, source->index)
        );

        ++loadedCount;
    }
};

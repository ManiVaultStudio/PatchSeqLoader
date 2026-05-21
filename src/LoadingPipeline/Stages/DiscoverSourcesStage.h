#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "Config/Config.h"

#include <QDir>
#include <QFileInfo>

class DiscoverSourcesStage final : public PipelineStage
{
public:
    QString Name() const override { return "DiscoverSources"; }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        int discoveredFileCount = 0;
        int missingFileCount = 0;
        int discoveredDirectoryCount = 0;
        int missingDirectoryCount = 0;

        CheckTableSource(ctx, "rna", ctx.config.rna, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, "ephys", ctx.config.ephys, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, "morphology", ctx.config.morphology, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, "metadata", ctx.config.metadata, discoveredFileCount, missingFileCount);

        CheckTableSource(ctx, "rna_umap", ctx.config.rnaUmap, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, "ephys_umap", ctx.config.ephysUmap, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, "morpho_umap", ctx.config.morphoUmap, discoveredFileCount, missingFileCount);

        for (auto it = ctx.config.extraEmbeddings.begin(); it != ctx.config.extraEmbeddings.end(); ++it)
        {
            CheckTableSource(ctx, it.key(), it.value(), discoveredFileCount, missingFileCount);
        }

        CheckDirectory(ctx, "assets.directories.morphology_reconstruction", ctx.config.assetDirectories.morphologyReconstruction, discoveredDirectoryCount, missingDirectoryCount);
        CheckDirectory(ctx, "assets.directories.ephys_traces", ctx.config.assetDirectories.ephysTraces, discoveredDirectoryCount, missingDirectoryCount);

        if (discoveredFileCount == 0 && discoveredDirectoryCount == 0)
            ctx.result.Warning(Name(), "No configured source files or asset directories were found.", "The config loaded successfully, but there may be nothing to load.");

        ctx.result.Info(Name(), "Source discovery complete.", QString("files found=%1, files missing=%2, directories found=%3, directories missing=%4")
            .arg(discoveredFileCount)
            .arg(missingFileCount)
            .arg(discoveredDirectoryCount)
            .arg(missingDirectoryCount)
        );

        return true;
    }

private:
    void CheckTableSource(PipelineContext& ctx, const QString& label, const std::optional<config::TableSource>& source, int& foundCount, int& missingCount) const
    {
        if (!source)
        {
            ctx.result.Info(Name(), "Optional table source is not configured.", label);
            return;
        }

        CheckFile(ctx, QString("sources.%1.path").arg(label), source->path, foundCount, missingCount);
    }

    void CheckFile(PipelineContext& ctx, const QString& configPath, const QString& filePath, int& foundCount, int& missingCount) const
    {
        if (filePath.isEmpty())
        {
            ctx.result.Warning(Name(), "Configured file path is empty.", configPath);
            ++missingCount;
            return;
        }

        const QFileInfo info(filePath);

        if (!info.exists())
        {
            ctx.result.Warning(Name(), "Configured file does not exist.", QString("%1 = %2").arg(configPath, filePath));
            ++missingCount;
            return;
        }

        if (!info.isFile())
        {
            ctx.result.Warning(Name(), "Configured path exists but is not a file.", QString("%1 = %2").arg(configPath, filePath));
            ++missingCount;
            return;
        }

        if (!info.isReadable())
        {
            ctx.result.Warning(Name(), "Configured file exists but is not readable.", QString("%1 = %2").arg(configPath, filePath));
            ++missingCount;
            return;
        }

        ctx.result.Info(Name(), "Discovered file.", QString("%1 = %2").arg(configPath, filePath));

        ++foundCount;
    }

    void CheckDirectory(PipelineContext& ctx, const QString& configPath, const QString& directoryPath, int& foundCount, int& missingCount) const
    {
        if (directoryPath.isEmpty())
        {
            ctx.result.Info(Name(), "Optional asset directory is not configured.", configPath);
            return;
        }

        const QFileInfo info(directoryPath);

        if (!info.exists())
        {
            ctx.result.Warning(Name(), "Configured directory does not exist.", QString("%1 = %2").arg(configPath, directoryPath));
            ++missingCount;
            return;
        }

        if (!info.isDir())
        {
            ctx.result.Warning(Name(), "Configured path exists but is not a directory.", QString("%1 = %2").arg(configPath, directoryPath));
            ++missingCount;
            return;
        }

        if (!info.isReadable())
        {
            ctx.result.Warning(Name(), "Configured directory exists but is not readable.", QString("%1 = %2").arg(configPath, directoryPath));
            ++missingCount;
            return;
        }

        ctx.result.Info(Name(), "Discovered directory.", QString("%1 = %2").arg(configPath, directoryPath));

        ++foundCount;
    }
};

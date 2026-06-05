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

        CheckTableSource(ctx, config::keys::sources::Rna, ctx.config.rna, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, config::keys::sources::Ephys, ctx.config.ephys, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, config::keys::sources::Morphology, ctx.config.morphology, discoveredFileCount, missingFileCount);
        CheckTableSource(ctx, config::keys::sources::Metadata, ctx.config.metadata, discoveredFileCount, missingFileCount);

        CheckTableSource(ctx, config::keys::embeddings::RnaUmap, ctx.config.rnaUmap, discoveredFileCount, missingFileCount, config::keys::Embeddings);
        CheckTableSource(ctx, config::keys::embeddings::EphysUmap, ctx.config.ephysUmap, discoveredFileCount, missingFileCount, config::keys::Embeddings);
        CheckTableSource(ctx, config::keys::embeddings::MorphoUmap, ctx.config.morphoUmap, discoveredFileCount, missingFileCount, config::keys::Embeddings);

        for (auto it = ctx.config.extraEmbeddings.begin(); it != ctx.config.extraEmbeddings.end(); ++it)
            CheckTableSource(ctx, it.key(), it.value(), discoveredFileCount, missingFileCount, config::keys::Embeddings);

        CheckEphysTracesSource(ctx, ctx.config.ephysTraces, discoveredFileCount, missingFileCount, discoveredDirectoryCount, missingDirectoryCount);
        CheckMorphologyReconstructionsSource(ctx, ctx.config.morphologyReconstructions, discoveredDirectoryCount, missingDirectoryCount);

        if (discoveredFileCount == 0 && discoveredDirectoryCount == 0)
            ctx.result.Warning(Name(), "No configured source files or directories were found.", "The config loaded successfully, but there may be nothing to load.");

        ctx.result.Info(Name(), "Source discovery complete.", QString("files found=%1, files missing=%2, directories found=%3, directories missing=%4")
            .arg(discoveredFileCount)
            .arg(missingFileCount)
            .arg(discoveredDirectoryCount)
            .arg(missingDirectoryCount)
        );

        return true;
    }

private:
    void CheckTableSource(PipelineContext& ctx, const QString& label, const std::optional<config::TableSource>& source, int& foundCount, int& missingCount, const QString& parentPath = config::keys::Sources) const
    {
        if (!source)
        {
            ctx.result.Info(Name(), "Optional table source is not configured.", QString("%1.%2").arg(parentPath, label));
            return;
        }

        CheckFile(ctx, QString("%1.%2.%3").arg(parentPath, label, config::keys::Path), source->path, foundCount, missingCount);
    }

    void CheckTableSource(PipelineContext& ctx, const QString& label, const config::TableSource& source, int& foundCount, int& missingCount, const QString& parentPath = config::keys::Sources) const
    {
        CheckFile(ctx, QString("%1.%2.%3").arg(parentPath, label, config::keys::Path), source.path, foundCount, missingCount);
    }

    void CheckEphysTracesSource(PipelineContext& ctx, const std::optional<config::EphysTracesSource>& source, int& foundFileCount, int& missingFileCount, int& foundDirectoryCount, int& missingDirectoryCount) const
    {
        const QString sourcePath = QString("%1.%2").arg(config::keys::Sources, config::keys::sources::EphysTraces);

        if (!source)
        {
            ctx.result.Info(Name(), "Optional directory-backed source is not configured.", sourcePath);
            return;
        }

        CheckDirectory(ctx, QString("%1.%2").arg(sourcePath, config::keys::Directory), source->directory, foundDirectoryCount, missingDirectoryCount);

        if (!source->failedSweepsPath.isEmpty())
            CheckFile(ctx, QString("%1.%2").arg(sourcePath, config::keys::FailedSweepsPath), source->failedSweepsPath, foundFileCount, missingFileCount);
    }

    void CheckMorphologyReconstructionsSource(PipelineContext& ctx, const std::optional<config::MorphologyReconstructionsSource>& source, int& foundDirectoryCount, int& missingDirectoryCount) const
    {
        const QString sourcePath = QString("%1.%2").arg(config::keys::Sources, config::keys::sources::MorphologyReconstructions);

        if (!source)
        {
            ctx.result.Info(Name(), "Optional directory-backed source is not configured.", sourcePath);
            return;
        }

        CheckDirectory(ctx, QString("%1.%2").arg(sourcePath, config::keys::Directory), source->directory, foundDirectoryCount, missingDirectoryCount);
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
            ctx.result.Warning(Name(), "Configured directory path is empty.", configPath);
            ++missingCount;
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

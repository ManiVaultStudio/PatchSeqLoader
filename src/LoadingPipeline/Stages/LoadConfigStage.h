#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "Config/Config.h"

#include <QString>
#include <QFileInfo>

class LoadConfigStage final : public PipelineStage
{
public:
    QString Name() const override { return "LoadConfig"; }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        if (ctx.configPath.isEmpty())
        {
            ctx.result.Error(Name(), "No config file path was provided.");
            return false;
        }

        if (!ctx.config.Load(ctx.configPath))
        {
            ctx.result.Error(Name(), "Failed to load config file.", ctx.configPath);
            return false;
        }

        ctx.datasetRoot = QFileInfo(ctx.configPath).absolutePath();

        ctx.result.Info(Name(), "Loaded config file.", QString("dataset_id=%1").arg(ctx.config.datasetId));

        return true;
    }
};

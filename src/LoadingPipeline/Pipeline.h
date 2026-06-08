#pragma once

#include "PipelineStage.h"
#include "PipelineContext.h"

#include <memory>
#include <vector>

class Pipeline
{
public:
    void Add(std::unique_ptr<PipelineStage> stage)
    {
        _stages.push_back(std::move(stage));
    }

    QStringList GetStageNames()
    {
        QStringList stageNames;
        for (size_t i = 0; i < _stages.size(); i++)
            stageNames.push_back(_stages[i]->Name());
        return stageNames;
    }

    PipelineResult Run(PipelineContext& ctx)
    {
        if (!ctx.task)
        {
            ctx.result.Error("Pipeline", "Pipeline task not initialized, aborting...");
            return ctx.result;
        }

        ctx.task->setSubtasks(GetStageNames());

        for (size_t i = 0; i < _stages.size(); ++i)
        {
            auto& stage = _stages[i];

            ctx.task->setSubtaskStarted(stage->Name(), "Running " + stage->Name()); QApplication::processEvents();

            const bool success = stage->Run(ctx);

            if (!success) {
                ctx.result.Error(stage->Name(), "Pipeline stopped at stage: " + stage->Name());
                return ctx.result;
            }

            ctx.task->setSubtaskFinished(stage->Name(), "Finished " + stage->Name()); QApplication::processEvents();
        }

        ctx.task->setFinished();
        ctx.task->setDescription("Patch-seq load complete");
        QApplication::processEvents();

        return ctx.result;
    }

private:
    std::vector<std::unique_ptr<PipelineStage>> _stages;
};

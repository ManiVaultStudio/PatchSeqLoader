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

    PipelineResult Run(PipelineContext& ctx)
    {
        const int total = static_cast<int>(_stages.size());

        for (int i = 0; i < total; ++i) {
            auto& stage = _stages[i];

            if (ctx.task) {
                ctx.task->setProgress(static_cast<float>(i) / total);
                ctx.task->setDescription("Running " + stage->Name());
            }

            const bool success = stage->Run(ctx);

            if (!success) {
                ctx.result.Error(
                    stage->Name(),
                    "Pipeline stopped at stage: " + stage->Name()
                );
                return ctx.result;
            }
        }

        if (ctx.task) {
            ctx.task->setProgress(1.0f);
            ctx.task->setDescription("Patch-seq load complete");
        }

        return ctx.result;
    }

private:
    std::vector<std::unique_ptr<PipelineStage>> _stages;
};

#pragma once

#include "PipelineContext.h"

class PipelineStage
{
public:
    virtual ~PipelineStage() = default;

    virtual QString Name() const = 0;

    virtual bool Run(PipelineContext& ctx) = 0;
};

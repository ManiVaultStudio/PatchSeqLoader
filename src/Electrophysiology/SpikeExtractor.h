#pragma once

#include "EphysData/TimeSeries.h"

class ActionPotential;

class SpikeExtractor
{
public:
    ActionPotential* DetectActionPotential(const TimeSeries& stim, const TimeSeries& acq);
private:

};

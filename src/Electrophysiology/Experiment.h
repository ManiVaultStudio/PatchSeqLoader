#pragma once

#include "Recording.h"

#include <vector>

class Experiment
{
public:
    const std::vector<Recording>& getAcquisitions() { return _acquisitions; }

    void addAcquisition(Recording& recording);
    void addStimulus(Recording& recording);

private:
    std::vector<Recording> _acquisitions;
    std::vector<Recording> _stimuli;
};

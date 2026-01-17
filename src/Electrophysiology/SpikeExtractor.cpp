#include "SpikeExtractor.h"

#include <EphysData/ActionPotential.h>

#include <fstream>

ActionPotential* SpikeExtractor::DetectActionPotential(const TimeSeries& stim, const TimeSeries& acq)
{
    if (stim.ySeries.empty())
        ;// throw exception;

    // Find first index of stimulus
    int stimIndex = 0;

    float prevValue = stim.ySeries[0];
    for (int i = 10000; i < stim.ySeries.size(); i++)
    {
        if (stim.ySeries[i] > prevValue)
        {
            stimIndex = i;
            break;
        }
    }

    // Find first peak
    int peakIndex = 0;

    {
        float maxValue = 0;
        for (int i = 0; i < acq.ySeries.size(); i++)
        {
            if (acq.ySeries[i] > maxValue)
            {
                maxValue = acq.ySeries[i];
                peakIndex = i;
            }
        }
    }

    // Isolate action potential
    float dydx = 0;
    bool up = true;
    prevValue = acq.ySeries[stimIndex];

    std::vector<float> actionPotential;
    for (int i = peakIndex - 70; i < acq.ySeries.size(); i++)
    {
        float y = acq.ySeries[i];
        dydx = (y - prevValue);
        if (dydx >= 0)
        {
            if (up == false)
            {
                break;
            }
            up = true;
        }
        if (dydx < 0 && i > peakIndex)
            up = false;

        actionPotential.push_back(y);

        prevValue = y;
    }

    std::vector<float> timeSeries(actionPotential.size());
    std::iota(timeSeries.begin(), timeSeries.end(), 0);
    for (int i = 0; i < timeSeries.size(); i++)
    {
        timeSeries[i] *= 0.02f;
    }
    //std::ofstream file("output.csv");
    //for (float v : actionPotential)
    //{
    //    file << v << "\n";
    //}

    //file.close();

    return new ActionPotential(timeSeries, actionPotential, peakIndex - stimIndex);
}

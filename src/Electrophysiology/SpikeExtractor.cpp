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

ActionPotential* SpikeExtractor::ExtractActionPotential(
    const TimeSeries& acq,
    int thresholdIndex)
{
    if (acq.xSeries.empty() ||
        acq.ySeries.empty() ||
        acq.xSeries.size() != acq.ySeries.size())
    {
        return nullptr;
    }

    if (thresholdIndex < 0 ||
        thresholdIndex >= static_cast<int>(acq.ySeries.size()))
    {
        return nullptr;
    }

    constexpr float PRE_THRESHOLD_SECONDS = 0.001f;
    constexpr float POST_THRESHOLD_SECONDS = 0.003f;

    const float thresholdTime =
        acq.xSeries[thresholdIndex];

    const float startTime =
        thresholdTime - PRE_THRESHOLD_SECONDS;

    const float endTime =
        thresholdTime + POST_THRESHOLD_SECONDS;

    auto startIt = std::lower_bound(
        acq.xSeries.begin(),
        acq.xSeries.end(),
        startTime);

    auto endIt = std::upper_bound(
        acq.xSeries.begin(),
        acq.xSeries.end(),
        endTime);

    const size_t startIndex =
        std::distance(acq.xSeries.begin(), startIt);

    const size_t endIndex =
        std::distance(acq.xSeries.begin(), endIt);

    if (startIndex >= endIndex)
        return nullptr;

    std::vector<float> apTime;
    std::vector<float> apVoltage;

    apTime.reserve(endIndex - startIndex);
    apVoltage.reserve(endIndex - startIndex);

    const float t0 = acq.xSeries[startIndex];

    for (size_t i = startIndex; i < endIndex; ++i)
    {
        // Assuming xSeries is seconds; output AP time in ms.
        apTime.push_back(
            (acq.xSeries[i] - t0) * 1000.0f);

        apVoltage.push_back(acq.ySeries[i]);
    }

    // Find the AP peak inside the extracted region.
    const auto peakIt = std::max_element(
        acq.ySeries.begin() + startIndex,
        acq.ySeries.begin() + endIndex);

    const int peakIndex =
        static_cast<int>(
            std::distance(acq.ySeries.begin(), peakIt));

    return new ActionPotential(
        apTime,
        apVoltage,
        peakIndex - thresholdIndex);
}

#include "SpikeDetector.h"

#include "EphysData/TimeSeries.h"

std::vector<int> DetectSpikes(const TimeSeries& timeSeries)
{
    int windowSize = 10;
    float threshold = 3.0f;

    const std::vector<float>& y = timeSeries.ySeries;
    std::vector<int> spikes;

    if (y.size() < windowSize * 2)
        return spikes;

    for (size_t i = windowSize; i < y.size() - windowSize; ++i)
    {
        float mean = 0.0f;
        float var = 0.0f;

        // Compute mean
        for (int j = -windowSize; j <= windowSize; ++j)
            mean += y[i + j];

        mean /= (2 * windowSize + 1);

        // Compute variance
        for (int j = -windowSize; j <= windowSize; ++j)
        {
            float diff = y[i + j] - mean;
            var += diff * diff;
        }

        float stddev = std::sqrt(var / (2 * windowSize + 1));

        if ((y[i] - mean) > threshold * stddev)
            spikes.push_back(i);
    }

    return spikes;
}

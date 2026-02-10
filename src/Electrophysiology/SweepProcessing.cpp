#include "SweepProcessing.h"

#include <algorithm>

namespace
{
    constexpr float EPSILON = 1e-6f;

    void FindNonzeroRanges(const std::vector<float>& v, std::vector<Envelope>& envelopes)
    {
        bool inRange = false;
        size_t startIndex = 0;

        for (size_t i = 0; i < v.size(); i++)
        {
            if (std::isnan(v[i]))
                break;

            if (abs(v[i]) > EPSILON)
            {
                if (!inRange)
                {
                    startIndex = i;
                    inRange = true;
                }
            }
            else
            {
                if (inRange)
                {
                    envelopes.emplace_back(startIndex, i, 0);
                    inRange = false;
                }
            }
        }

        // Handle trailing non-zero range
        if (inRange)
            envelopes.emplace_back(startIndex, v.size(), 0);
    }

    void ComputeEnvelopeAreas(std::vector<Envelope>& envelopes, const std::vector<float>& v)
    {
        for (Envelope& env : envelopes)
        {
            env.area = 0;
            for (int i = env.startIndex; i < env.endIndex; i++)
            {
                env.area += abs(v[i]);
            }
        }
    }
}

std::vector<Envelope> ComputeStimulusEnvelopes(Stimulus stimulus)
{
    TimeSeries& ts = stimulus.GetRecording().GetData();

    std::vector<Envelope> envelopes;
    FindNonzeroRanges(ts.ySeries, envelopes);

    ComputeEnvelopeAreas(envelopes, ts.ySeries);

    if (envelopes.size() <= 1)
        return envelopes;

    //float avgPeak = envelopes[0].area / (envelopes[0].endIndex - envelopes[0].startIndex);

    if (envelopes[0].area < envelopes[1].area / 4)
        envelopes.erase(envelopes.begin()); // Remove test spike (which is around 50mV for ~2000 steps)
    else
        qDebug() << "Interesting Envelope area: " << envelopes[0].area << envelopes[1].area;

    return envelopes;
}

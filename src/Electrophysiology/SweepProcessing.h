#pragma once

#include <EphysData/Recording.h>

class Envelope
{
public:
    int startIndex;
    int endIndex;
    float area;
    float peak;

    bool operator<(const Envelope& other) const
    {
        return area > other.area;
    }
};

std::vector<Envelope> ComputeStimulusEnvelopes(Stimulus stimulus);

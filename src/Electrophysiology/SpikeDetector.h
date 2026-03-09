#pragma once

#include <vector>

class TimeSeries;

std::vector<int> DetectSpikes(const TimeSeries& timeSeries);

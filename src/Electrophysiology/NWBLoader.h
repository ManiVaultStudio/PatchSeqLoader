#pragma once

#include <QString>

class Experiment;

class NWBLoader
{
public:
    void LoadNWB(QString fileName, Experiment& experiment);
private:

};

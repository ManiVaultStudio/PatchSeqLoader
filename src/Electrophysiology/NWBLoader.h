#pragma once

#include <QString>
#include <QHash>

class Experiment;

class LoadInfo
{
public:
    QString failedSweepPath;
    // List of loaded vs. not loaded stim descriptions
    QHash<QString, int> ignoredStimsets;
    QHash<QString, int> loadedStimsets;
};

class NWBLoader
{
public:
    void LoadNWB(QString fileName, Experiment& experiment, LoadInfo& info);
private:

};

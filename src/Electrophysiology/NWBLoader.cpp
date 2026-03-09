#include "NWBLoader.h"

#include "StimulusCodeMap.h"

#include "EphysData/Experiment.h"
#include "EphysData/ActionPotential.h"

#include <QDebug>
#include <iostream>
#include <string>
#include <fstream>
#include <regex>

#include "LEAD/NWBFile.h"
#include "Electrophysiology/SpikeExtractor.h"
#include "Electrophysiology/SweepProcessing.h"
#include "Electrophysiology/FailedSweepDetector.h"
#include "Electrophysiology/SpikeDetector.h"

///
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
///

#include <windows.h>

using namespace H5;
int failIndex = 0;
namespace
{
    void exportToCSV(const std::vector<float>& x, const std::vector<float>& y, const std::string& filename) {
        if (x.size() != y.size()) {
            std::cerr << "Error: x and y vectors must be the same size." << std::endl;
            return;
        }

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: could not open file " << filename << " for writing." << std::endl;
            return;
        }

        file << "x,y\n"; // CSV header

        for (size_t i = 0; i < x.size(); ++i) {
            file << x[i] << "," << y[i] << "\n";
        }

        file.close();
        std::cout << "Data exported to " << filename << std::endl;
    }

    int extractSweepNumber(const QString& input)
    {
        QRegularExpression re("data_(\\d{5})_");
        QRegularExpressionMatch match = re.match(input);
        if (match.hasMatch()) {
            QString numberStr = match.captured(1);
            return numberStr.toInt();
        }
        return -1;  // or throw, depending on your use case
    }
}

class RecordingPair
{
public:
    LEAD::Group acquisition;
    LEAD::Group stimulus;
};

namespace
{
    QString ExtractFileId(const QString& path)
    {
        QFileInfo fi(path);
        return fi.completeBaseName();   // filename without extension
    }

    void ExtractRecordings(const Groups& groups, QHash<QString, RecordingPair>& recordingPairs)
    {
        // Find acquisitions and stimuli
        for (int i = 0; i < groups.size(); i++)
        {
            QString groupPath = QString::fromStdString(groups[i].GetName());

            // Extract the group name by taking the section beyond the last forward slash
            QString groupName = groupPath.section('/', -1);

            // Extract the number from e.g. data_00011_AD0
            QRegularExpression re("^.*_(\\d+)_.*$");
            QRegularExpressionMatch match = re.match(groupName);
            if (match.hasMatch()) {
                QString number = match.captured(1);

                if (!recordingPairs.contains(number))
                {
                    recordingPairs[number] = RecordingPair();
                }

                if (groupPath.startsWith("acquisition/"))
                    recordingPairs[number].acquisition = groups[i];

                if (groupPath.startsWith("stimulus/presentation/"))
                    recordingPairs[number].stimulus = groups[i];
            }
        }
    }

    LEAD::Dataset FindDataset(const Datasets& datasets, std::string datasetName)
    {
        for (int i = 0; i < datasets.size(); i++)
        {
            if (datasets[i].GetName() == datasetName)
                return datasets[i];
        }
    }

    bool FindStimulusDescription(NWBFile& nwbFile, LEAD::Group& group, QString& stimDescription)
    {
        // Check whether the acquisition group has a stimulus description dataset
        std::string stimulusDescriptionDatasetName = group.GetName() + "/stimulus_description";
        if (nwbFile.DatasetExists(stimulusDescriptionDatasetName))
        {
            // Attempt to load the stimulus description dataset
            std::vector<std::string> stimDescriptions;
            nwbFile.OpenStringDataset(stimulusDescriptionDatasetName, stimDescriptions);
            if (stimDescriptions.empty())
            {
                qWarning() << "Found stimulus description dataset, but couldn't load it";
                return false;
            }
            else
            {
                stimDescription = QString::fromStdString(stimDescriptions[0]);
                return true;
            }
        }
        else
        {
            // No stimulus description dataset, check if theres a stimulus description attribute
            group.LoadAllAttributes(nwbFile.GetFileId());
            for (const auto& attr : group.GetAttributes())
            {
                if (attr.GetName() == "stimulus_description")
                {
                    stimDescription = QString::fromStdString(attr.GetValue());
                    return true;
                }
            }
        }

        qWarning() << "No stimulus description found.";
        return false;
    }

    int ExtractStimScale(std::string input)
    {
        std::regex re(R"(Stim Scale Factor:\s*([0-9]+)\.)");
        std::smatch match;

        if (std::regex_search(input, match, re)) {
            int stimScaleFactor = std::stoi(match[1]);
            return stimScaleFactor;
        }
        return -1;
    }

    void ReadTimeseries(NWBFile& file, std::string groupName, Recording& recording)
    {
        //std::cout << "TIMESERIES " << groupName << std::endl;
        std::vector<hsize_t> dims;
        file.OpenFloatDataset(groupName + "/data", recording.GetData().ySeries, dims);

        recording.GetData().xSeries.resize(recording.GetData().ySeries.size());
        std::iota(recording.GetData().xSeries.begin(), recording.GetData().xSeries.end(), 0);

        // Read sampling rate for xSeries
        std::string rateDatasetName = groupName + "/starting_time";
        LEAD::Dataset rateDataset = FindDataset(file.GetDatasets(), rateDatasetName);
        rateDataset.LoadAllAttributes(file.GetFileId());
        float rate = 1;

        for (const auto& attr : rateDataset.GetAttributes())
        {
            if (attr.GetName() == "rate")
                rate = QString::fromStdString(attr.GetValue()).toFloat();
        }

        float timeStep = (1 / rate);
        std::transform(recording.GetData().xSeries.begin(), recording.GetData().xSeries.end(), recording.GetData().xSeries.begin(), [timeStep](auto& c) { return c * timeStep; });
        //recording.GetData().downsample();

        /////
        //std::cout << file.GetFileName() << std::endl;
        //if (file.GetFileName().find("H19.03.302.11.14.02.05") != std::string::npos)
        //{
        //    QString fileName = QString::fromStdString(groupName);
        //    fileName = fileName.replace("/", "_");

        //    exportToCSV(recording.GetData().xSeries, recording.GetData().ySeries, file.GetFileName() + "-" + fileName.toStdString() + ".csv");
        //}
        /////
    }

    bool DetectFailedAcquisition(const Stimulus& stimulus, const Recording& acquisition)
    {
        const TimeSeries& stimSeries = stimulus.GetRecording().GetData();

        int lastSignal = -1;
        int failCount = 0;
        for (int i = 0; i < stimSeries.ySeries.size(); i++)
        {
            if (abs(stimSeries.ySeries[i]) > 0.001f)
                lastSignal = i;
            //if (stimSeries.ySeries[i] > 0.001f && acquisition.GetData().ySeries[i] == 0)
            //    failCount++;
            bool acqFail = abs(acquisition.GetData().ySeries[i]) < 0.001f || std::isnan(acquisition.GetData().ySeries[i]);

            if (acqFail && lastSignal != -1 && i - lastSignal < 2000)
                failCount++;
            //if (abs(stimSeries.ySeries[i]) > 0.001f && abs(acquisition.GetData().ySeries[i]) < 0.001f)
            //    failCount++;
        }
        //qDebug() << "Fail count: " << failCount;
        return failCount > 50;
        //if (failCount > 30)
        //{
        //    return true;
        //}

        //float refValue = recording.GetData().ySeries[recording.GetData().ySeries.size() / 2];
        //for (int i = recording.GetData().ySeries.size() / 2; i < recording.GetData().ySeries.size(); i++)
        //{
        //    if (std::abs(recording.GetData().ySeries[i] - refValue) > 0.001f)
        //        return false;
        //}
        //return true;
    }
}

void NWBLoader::LoadNWB(QString filePath, Experiment& experiment, LoadInfo& info)
{
    NWBFile nwbFile;
    nwbFile.Load(filePath.toStdString());
    Groups groups = nwbFile.GetGroups();

    nwbFile.Open(filePath.toStdString());
    bool written = false;

    float totalSize = 0;
    qDebug() << "Filepath: " << filePath;

    QString fileName = ExtractFileId(filePath);

    QHash<QString, QVector<int>> failedSweepDict = LoadFailedSweeps(":met_loader/failed_sweeps.json");

    QVector<int> failedSweeps;
    if (failedSweepDict.contains(fileName))
        failedSweeps = failedSweepDict[fileName];

    QHash<QString, RecordingPair> recordingPairs;
    ExtractRecordings(groups, recordingPairs);

    for (RecordingPair& recordingPair : recordingPairs)
    {
        // Find stimulus description
        QString stimDescription;
        bool stimDescriptionFound = FindStimulusDescription(nwbFile, recordingPair.acquisition, stimDescription);

        if (!stimDescriptionFound)
            continue;

        // There is a stimulus description, chop it, and determine if we should load the associated recordings
        stimDescription.chop(5); // Trim _DA_0

        if (!USEFUL_STIM_CODES.contains(stimDescription))
        {
            //qWarning() << "Not loading recordings because stimulus description was: " << stimDescription; // TEMP
            if (!info.ignoredStimsets.contains(stimDescription))
                info.ignoredStimsets[stimDescription] = 1;
            else
                info.ignoredStimsets[stimDescription]++;
            continue;
        }
        else
        {
            if (!info.loadedStimsets.contains(stimDescription))
                info.loadedStimsets[stimDescription] = 1;
            else
                info.loadedStimsets[stimDescription]++;
        }

        //if (STIMULUS_CODE_NAME_MAP.contains(stimDescription))
        //    stimDescription = STIMULUS_CODE_NAME_MAP[stimDescription];

        // Extract sweep number
        int acqSweepNumber = extractSweepNumber(QString::fromStdString(recordingPair.acquisition.GetName()));
        int stimSweepNumber = extractSweepNumber(QString::fromStdString(recordingPair.stimulus.GetName()));

        if (stimSweepNumber != acqSweepNumber)
        {
            qCritical() << "[ERROR] Stimulus has sweep number: " << stimSweepNumber << " but acquisition: " << acqSweepNumber;
            continue;
        }

        // Check if sweep failed QC
        bool failedSweep = false;
        for (const int& sweepNum : failedSweeps)
            if (sweepNum == acqSweepNumber)
                failedSweep = true;

        if (failedSweep)
        {
            qDebug() << "[" << fileName << "]" << "Discarding failed sweep:" << acqSweepNumber;
            continue;
        }

        // Recordings should be loaded, so store stimulus description in both recordings
        Sweep sweep;

        sweep.SetSweepNumber(stimSweepNumber);
        sweep.stimulus.SetStimulusDescription(stimDescription);

        // Load associated timeseries
        // ACQUISITION
        {
            // Only load attributes for the acquisition if they have not been previously loaded while checking for a stim description FIXME
            if (recordingPair.acquisition.GetAttributes().empty())
                recordingPair.acquisition.LoadAllAttributes(nwbFile.GetFileId());

            // Load all attributes
            for (int j = 0; j < recordingPair.acquisition.GetAttributes().size(); j++)
            {
                const LEAD::Attribute& attribute = recordingPair.acquisition.GetAttributes()[j];

                sweep.acquisition.AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
            }

            ReadTimeseries(nwbFile, recordingPair.acquisition.GetName(), sweep.acquisition);

            totalSize += ((sweep.acquisition.GetData().xSeries.size() + sweep.acquisition.GetData().ySeries.size()) * sizeof(float)) / 1000000.0f;

            for (auto it = sweep.acquisition.GetAttributes().constBegin(); it != sweep.acquisition.GetAttributes().constEnd(); ++it)
            {
                totalSize += it.value().size() / 1000000.0f;
                //qDebug() << "Attribute size: " << it.key() << " " << it.value().size();
            }
        }

        // STIMULUS
        {
            recordingPair.stimulus.LoadAllAttributes(nwbFile.GetFileId());

            // Load all attributes
            for (int j = 0; j < recordingPair.stimulus.GetAttributes().size(); j++)
            {
                const LEAD::Attribute& attribute = recordingPair.stimulus.GetAttributes()[j];

                sweep.stimulus.GetRecording().AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
            }

            ReadTimeseries(nwbFile, recordingPair.stimulus.GetName(), sweep.stimulus.GetRecording());
            
            // Extract action potential before downsampling
            if (stimDescription.contains("Rheo"))
            {
                SpikeExtractor extractor;
                ActionPotential* actionPotential = extractor.DetectActionPotential(sweep.stimulus.GetRecording().GetData(), sweep.acquisition.GetData());
                experiment.setActionPotential(actionPotential);

                for (const auto& attribute : recordingPair.stimulus.GetAttributes())
                {
                    if (attribute.GetName().find("comment") != std::string::npos)
                    {
                        int stimScale = ExtractStimScale(attribute.GetValue());

                        if (true)
                        //if (stimScale == 100)
                        {

                        }
                    }
                }
            }
        }

        // Detect failed acquisitions
        //if (stimDescription.contains("Ramp"))
        //    qDebug() << "Ramp";
        //if (fileName.contains("H23.06.351.11.56.01.06"))
        //{
        //    qDebug() << "Test" << stimSweepNumber;
        //}
        //if (failIndex >= 2302)
        //{
        //    qDebug() << "Test" << stimSweepNumber;
        //}

        //bool failed = DetectFailedAcquisition(sweep.stimulus, sweep.acquisition);
        //if (failed)
        //{
        //    qDebug() << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Discarded failed acq " << failIndex << stimSweepNumber << stimDescription;

        //    if (stimDescription.contains("ramp", Qt::CaseInsensitive))
        //    {
        //        exportToCSV(sweep.acquisition.GetData().xSeries, sweep.acquisition.GetData().ySeries, std::to_string(failIndex) + "_acq.csv");
        //        exportToCSV(sweep.stimulus.GetRecording().GetData().xSeries, sweep.stimulus.GetRecording().GetData().ySeries, std::to_string(failIndex) + "_stim.csv");
        //        failIndex += 1;
        //    }

        //    continue;
        //}

        // Downsample the recording
        sweep.acquisition.GetData().downsample();
        sweep.stimulus.GetRecording().GetData().downsample();

        std::vector<Envelope> stimEnvelopes = ComputeStimulusEnvelopes(sweep.stimulus);
        if (stimEnvelopes.empty())
        {
            qDebug() << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Empty Envelopes";
            continue;
        }

        std::pair<int, int> stimRange = { stimEnvelopes[0].startIndex, stimEnvelopes[stimEnvelopes.size() - 1].endIndex };
        //std::pair<int, int> stimRange = sweep.stimulus.GetRecording().GetData().FindStimulusRange();
        if (stimRange.first != -1)
        {
            sweep.stimulus.GetRecording().GetData().trim(stimRange.first, stimRange.second);
            sweep.acquisition.GetData().trim(stimRange.first, stimRange.second);
        }
        else
            continue;

        sweep.stimulus.GetRecording().GetData().computeExtents();
        sweep.stimulus.CalculateStimulusAmplitude();
        sweep.stimulus.DetectStimulusType();
        sweep.acquisition.GetData().computeExtents();

        // Detect spikes
        std::vector<int> spikeIndices = DetectSpikes(sweep.acquisition.GetData());
        sweep.acquisition.AddAttribute("NumSpikes", QString::number(static_cast<qulonglong>(spikeIndices.size())));

        experiment.AddSweep(std::move(sweep));
    }

    //for (int i = 0; i < groups.size(); i++)
    //{
    //    QString groupName = QString::fromStdString(groups[i].GetName());

    //    Recording recording;

    //    // Load all acquisitions
    //    if (groupName.startsWith("acquisition/") && groupName.contains("data"))
    //    {
    //        // Check stimulus description whether we should load the dataset
    //        // Determine if the acquisition has a stimulus description, and if so load it
    //        std::string stimulusDescriptionDatasetName = (groupName + "/stimulus_description").toStdString();
    //        if (nwbFile.DatasetExists(stimulusDescriptionDatasetName))
    //        {
    //            std::vector<std::string> stimDescriptions;
    //            nwbFile.OpenStringDataset(stimulusDescriptionDatasetName, stimDescriptions);
    //            if (!stimDescriptions.empty())
    //            {
    //                QString stimDescription = QString::fromStdString(stimDescriptions[0]);
    //                stimDescription.chop(5); // Trim _DA_0

    //                if (!USEFUL_STIM_CODES.contains(stimDescription))
    //                {
    //                    loadedDatasets.push_back(false);
    //                    continue;
    //                }
    //                
    //                loadedDatasets.push_back(true);

    //                if (STIMULUS_CODE_NAME_MAP.contains(stimDescription))
    //                    stimDescription = STIMULUS_CODE_NAME_MAP[stimDescription];

    //                recording.SetStimulusDescription(stimDescription);
    //            }
    //        }

    //        //std::cout << i << ": " << groups[i].GetName() << std::endl;
    //        groups[i].LoadAllAttributes(nwbFile.GetFileId());
    //        
    //        // Load all attributes
    //        for (int j = 0; j < groups[i].GetAttributes().size(); j++)
    //        {
    //            const LEAD::Attribute& attribute = groups[i].GetAttributes()[j];

    //            recording.AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
    //        }

    //        ReadTimeseries(nwbFile, groupName.toStdString(), recording);

    //        totalSize += ((recording.GetData().xSeries.size() + recording.GetData().ySeries.size()) * sizeof(float)) / 1000000.0f;

    //        for (auto it = recording.GetAttributes().constBegin(); it != recording.GetAttributes().constEnd(); ++it)
    //        {
    //            totalSize += it.value().size() / 1000000.0f;
    //            //qDebug() << "Attribute size: " << it.key() << " " << it.value().size();
    //        }

    //        experiment.addAcquisition(std::move(recording));
    //    }
    //    // Load all stimuli
    //    if (groupName.startsWith("stimulus/presentation/") && groupName.contains("data"))
    //    {
    //        //std::cout << i << ": " << groups[i].GetName() << std::endl;
    //        groups[i].LoadAllAttributes(nwbFile.GetFileId());

    //        // Load all attributes
    //        for (int j = 0; j < groups[i].GetAttributes().size(); j++)
    //        {
    //            const LEAD::Attribute& attribute = groups[i].GetAttributes()[j];

    //            recording.AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
    //        }

    //        ReadTimeseries(nwbFile, groupName.toStdString(), recording);

    //        totalSize += ((recording.GetData().xSeries.size() + recording.GetData().ySeries.size()) * sizeof(float)) / 1000000.0f;

    //        for (auto it = recording.GetAttributes().constBegin(); it != recording.GetAttributes().constEnd(); ++it)
    //        {
    //            totalSize += it.value().size() / 1000000.0f;
    //            //qDebug() << "Attribute2 size: " << it.key() << " " << it.value().size();
    //        }

    //        experiment.addStimulus(std::move(recording));
    //    }
    //}
    nwbFile.Close();
    qDebug() << ">>>>>>>> NUM SWEEPS: " << experiment.GetSweeps().size();
    std::cout << "Size: " << totalSize << "MB" << std::endl;

    //if (experiment.getAcquisitions().size() > 0)
    //{
    //    const QHash<QString, QString>& attrs = experiment.getAcquisitions()[0].GetAttributes();
    //    for (auto it = attrs.constBegin(); it != attrs.constEnd(); ++it)
    //    {
    //        qDebug() << "Attr: " << it.key() << "Value: " << it.value();
    //    }
    //}

    //dataset.iterateAttrs((H5::attr_operator_t) attr_op);
}

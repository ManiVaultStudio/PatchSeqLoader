#define NOMINMAX

#include "NWBLoader.h"

#include "StimulusCodeMap.h"

#include "EphysData/Experiment.h"
#include "EphysData/ActionPotential.h"

#include <QDebug>
#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <limits>

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

        for (size_t i = 0; i < x.size(); i += 10) {
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

        float timeStep = (1.0f / rate);
        std::transform(recording.GetData().xSeries.begin(), recording.GetData().xSeries.end(), recording.GetData().xSeries.begin(), [timeStep](auto& c) { return c * timeStep; });
        //recording.GetData().downsample();

        /////
        //std::cout << file.GetFileName() << std::endl;
        //if (file.GetFileName().find("H19.03.302.11.14.02.05") != std::string::npos)
        //{
            //QString fileName = QString::fromStdString(groupName);
            //fileName = fileName.replace("/", "_");

            //exportToCSV(recording.GetData().xSeries, recording.GetData().ySeries, file.GetFileName() + "-" + fileName.toStdString() + ".csv");
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
        return failCount > 1800;
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

    // Open and load file hierarchy as groups and datasets, then close it again
    nwbFile.Load(filePath.toStdString());
    Groups groups = nwbFile.GetGroups();

    // Reopen file
    nwbFile.Open(filePath.toStdString());

    float totalSize = 0;
    qDebug() << "Filepath: " << filePath;

    // Get name of the file (e.g. QN24.26.017.15.06A.06)
    QString fileName = ExtractFileId(filePath);
    experiment.SetName(fileName.toStdString());

    // Load list of sweeps marked as failed (FIXME: Probably should be optional)
    QHash<QString, QVector<int>> failedSweepDict = LoadFailedSweeps(info.failedSweepPath);

    // Get indices of failed sweeps for this file
    QVector<int> failedSweeps;
    if (failedSweepDict.contains(fileName))
        failedSweeps = failedSweepDict[fileName];

    // Extract and store stimulus and acquisition as a recording pair
    QHash<QString, RecordingPair> recordingPairs;
    ExtractRecordings(groups, recordingPairs);

    // For every recording pair, load and process all the data
    for (RecordingPair& recordingPair : recordingPairs)
    {
        // Find stimulus description (e.g. X4PS_SupraThresh_DA_1)
        QString stimDescription;
        bool stimDescriptionFound = FindStimulusDescription(nwbFile, recordingPair.acquisition, stimDescription);

        if (!stimDescriptionFound)
            continue;

        // There is a stimulus description, chop it, and determine if we should load the associated recordings
        stimDescription.chop(5); // Trim e.g. _DA_0

        // Check if stim description is part of ignored stimsets or not
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

        // Extract sweep number (e.g. data_00083 -> 83)
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
        sweep.stimulus.DetectStimulusType();

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

                sweep.acquisition.GetRecording().AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
            }

            ReadTimeseries(nwbFile, recordingPair.acquisition.GetName(), sweep.acquisition.GetRecording());

            totalSize += ((sweep.acquisition.GetRecording().GetData().xSeries.size() + sweep.acquisition.GetRecording().GetData().ySeries.size()) * sizeof(float)) / 1000000.0f;

            for (auto it = sweep.acquisition.GetRecording().GetAttributes().constBegin(); it != sweep.acquisition.GetRecording().GetAttributes().constEnd(); ++it)
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
            
            //// Extract action potential before downsampling
            //if (stimDescription.contains("Rheo"))
            //{
            //    SpikeExtractor extractor;
            //    ActionPotential* actionPotential = extractor.DetectActionPotential(sweep.stimulus.GetRecording().GetData(), sweep.acquisition.GetRecording().GetData());
            //    experiment.setActionPotential(actionPotential);

            //    for (const auto& attribute : recordingPair.stimulus.GetAttributes())
            //    {
            //        if (attribute.GetName().find("comment") != std::string::npos)
            //        {
            //            int stimScale = ExtractStimScale(attribute.GetValue());

            //            if (true)
            //            //if (stimScale == 100)
            //            {

            //            }
            //        }
            //    }
            //}
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

        //// Downsample the recording
        //sweep.acquisition.GetRecording().GetData().Downsample();
        //sweep.stimulus.GetRecording().GetData().Downsample();

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
            constexpr float SWEEP_PADDING_SECONDS = 0.1f;

            sweep.stimulus.GetRecording().GetData().Trim(stimRange.first, stimRange.second, SWEEP_PADDING_SECONDS);
            sweep.acquisition.GetRecording().GetData().Trim(stimRange.first, stimRange.second, SWEEP_PADDING_SECONDS);
        }
        else
            continue;

        sweep.stimulus.CalculateStimulusAmplitude();

        sweep.DetectSpikes();
        if (experiment.GetName() == "QM24.26.011.14v.02A.04")
        {
            std::cout << "beep" << std::endl;
            if (sweep.GetSweepProperties().spikeIndices.size() == 1)
                sweep.DetectSpikes();
        }
        sweep.stimulus.GetRecording().GetData().ComputeExtents();
        sweep.acquisition.GetRecording().GetData().ComputeExtents();

        experiment.AddSweep(std::move(sweep));
    }

    // Find Rheobase
    double minRheobaseAmplitude = std::numeric_limits<double>::max();
    int rheobaseIndex = -1;
    for (int i = 0; i < experiment.GetSweeps().size(); i++)
    {
        const Sweep& sweep = experiment.GetSweeps()[i];

        // Only regard long square stimuli for computing rheobase
        if (sweep.stimulus.GetStimulusType() != StimulusType::LongSquare)
            continue;

        // Find the sweep with the minimum stimulus amplitude required to produce a spike
        const int spikeCount = sweep.GetSweepProperties().GetSpikeCount();
        const double stimulusAmplitude = sweep.stimulus.GetStimulusAmplitude();

        if (spikeCount > 0 && stimulusAmplitude < minRheobaseAmplitude)
        {
            rheobaseIndex = i;
            minRheobaseAmplitude = stimulusAmplitude;
        }
    }

    if (rheobaseIndex != -1)
    {
        const Sweep& rheobase = experiment.GetSweeps()[rheobaseIndex];
        SpikeExtractor extractor;

            // Temp code
            if (experiment.GetName() == "QM24.26.011.14v.02A.04")
            {
            const auto& acqData = rheobase.acquisition.GetRecording().GetData();
            const auto& spikeIndices = rheobase.GetSweepProperties().spikeIndices;

            // Print spike indices
            std::cout << "Spike indices: ";
            for (int index : spikeIndices)
            {
                std::cout << index << " ";
            }
            std::cout << std::endl;

            // Save acquisition y-series to disk
            {
                std::ofstream file(
                    experiment.GetName() + "_rheobase_acq_y.csv"
                );

                if (file.is_open())
                {
                    file << "index,y\n";

                    for (size_t i = 0; i < acqData.ySeries.size(); ++i)
                    {
                        file << acqData.xSeries[i] << "," << acqData.ySeries[i] << "\n";
                    }

                    file.close();

                    std::cout
                        << "Saved rheobase acquisition to "
                        << experiment.GetName()
                        << "_rheobase_acq_y.csv"
                        << std::endl;
                }
                else
                {
                    std::cerr << "Failed to open rheobase acquisition output file."
                        << std::endl;
                }
            }
        }
        //

        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << experiment.GetName();
        
        ActionPotential* actionPotential = extractor.ExtractActionPotential(rheobase.acquisition.GetRecording().GetData(), rheobase.GetSweepProperties().spikeIndices[0]);
        experiment.setActionPotential(actionPotential);
    }

    for (Sweep& sweep : experiment.GetSweeps())
    {
        sweep.stimulus.GetRecording().GetData().Downsample();
        sweep.acquisition.GetRecording().GetData().Downsample();
    }

    nwbFile.Close();
    qDebug() << ">>>>>>>> NUM SWEEPS: " << experiment.GetSweeps().size();
    std::cout << "Size: " << totalSize << "MB" << std::endl;
}

#include "NWBLoader.h"

#include "StimulusCodeMap.h"

#include "EphysData/Experiment.h"

#include <QDebug>
#include <iostream>
#include <string>
#include <fstream>

#include "LEAD/NWBFile.h"

///
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
///

#include <windows.h>

using namespace H5;

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

    void ReadTimeseries(NWBFile& file, std::string groupName, Recording& recording)
    {
        //std::cout << "TIMESERIES " << groupName << std::endl;
        std::vector<hsize_t> dims;
        file.OpenFloatDataset(groupName + "/data", recording.GetData().ySeries, dims);

        recording.GetData().xSeries.resize(recording.GetData().ySeries.size());
        std::iota(recording.GetData().xSeries.begin(), recording.GetData().xSeries.end(), 0);

        //QString fileName = QString::fromStdString(groupName);
        //fileName = fileName.replace("/", "_");

        //exportToCSV(recording.data.xSeries, recording.data.ySeries, file.getFileName() + "-" + fileName.toStdString() + ".csv");

        std::transform(recording.GetData().xSeries.begin(), recording.GetData().xSeries.end(), recording.GetData().xSeries.begin(), [](auto& c) { return c / 1000.0f; });
        recording.GetData().downsample();
        recording.GetData().trim();
        recording.GetData().computeExtents();
    }
}

void NWBLoader::LoadNWB(QString fileName, Experiment& experiment, LoadInfo& info)
{
    NWBFile nwbFile;
    nwbFile.Load(fileName.toStdString());
    Groups groups = nwbFile.GetGroups();

    nwbFile.Open(fileName.toStdString());
    bool written = false;

    float totalSize = 0;

    QHash<QString, RecordingPair> recordingPairs;
    ExtractRecordings(groups, recordingPairs);

    for (RecordingPair& recordingPair : recordingPairs) {
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

        // Recordings should be loaded, so store stimulus description in both recordings
        Recording acquisition;
        Recording stimulus;

        acquisition.SetSweepNumber(acqSweepNumber);
        stimulus.SetSweepNumber(stimSweepNumber);

        acquisition.SetStimulusDescription(stimDescription);
        stimulus.SetStimulusDescription(stimDescription);

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

                acquisition.AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
            }

            ReadTimeseries(nwbFile, recordingPair.acquisition.GetName(), acquisition);

            totalSize += ((acquisition.GetData().xSeries.size() + acquisition.GetData().ySeries.size()) * sizeof(float)) / 1000000.0f;

            for (auto it = acquisition.GetAttributes().constBegin(); it != acquisition.GetAttributes().constEnd(); ++it)
            {
                totalSize += it.value().size() / 1000000.0f;
                //qDebug() << "Attribute size: " << it.key() << " " << it.value().size();
            }

            experiment.addAcquisition(std::move(acquisition));
        }

        // STIMULUS
        {
            recordingPair.stimulus.LoadAllAttributes(nwbFile.GetFileId());

            // Load all attributes
            for (int j = 0; j < recordingPair.stimulus.GetAttributes().size(); j++)
            {
                const LEAD::Attribute& attribute = recordingPair.stimulus.GetAttributes()[j];

                stimulus.AddAttribute(QString::fromStdString(attribute.GetName()), QString::fromStdString(attribute.GetValue()));
            }

            ReadTimeseries(nwbFile, recordingPair.stimulus.GetName(), stimulus);

            totalSize += ((stimulus.GetData().xSeries.size() + stimulus.GetData().ySeries.size()) * sizeof(float)) / 1000000.0f;

            for (auto it = stimulus.GetAttributes().constBegin(); it != stimulus.GetAttributes().constEnd(); ++it)
            {
                totalSize += it.value().size() / 1000000.0f;
                //qDebug() << "Attribute size: " << it.key() << " " << it.value().size();
            }

            experiment.addStimulus(std::move(stimulus));
        }
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
    qDebug() << ">>>>>>>> NUM ACQUISITIONS: " << experiment.getAcquisitions().size();
    qDebug() << ">>>>>>>> NUM STIMULI: " << experiment.getStimuli().size();
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

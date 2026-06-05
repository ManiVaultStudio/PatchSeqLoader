#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "Electrophysiology/NWBLoader.h"
#include "Morphology/SWCLoader.h"

#include <EphysData/EphysData.h>
#include <CellMorphologyData/CellMorphologyData.h>

#include <util/Timer.h>

class LoadAssetsStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "CollectMetadata";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();



        return true;
    }

private:
    void LoadEphysTraces(PipelineContext& ctx, QDir dir)
    {
        // Create ephys dataset
        ctx.ephysTraces = mv::data().createDataset<EphysExperiments>("Electrophysiology Data", "Ephys Traces", mv::Dataset<DatasetImpl>(), "", false);
        ctx.ephysTraces->setProperty("PatchSeqType", "EphysTraces");

        // Find all .nwb files in given directory
        QDir ephysTracesDir(dir);

        QStringList nwbFiles = ephysTracesDir.entryList(QStringList() << "*.nwb" << "*.NWB", QDir::Files);

        // Map metadata cell names to cell_ids
        std::vector<QString> metaCellSpecimenNames = _metadataDf[CELL_NAME_TAG];
        std::vector<QString> metaCellIds = _metadataDf[CELL_ID_TAG];

        std::unordered_map<QString, QString> specimenNameToCellId;
        for (size_t i = 0; i < metaCellSpecimenNames.size(); ++i)
            specimenNameToCellId[metaCellSpecimenNames[i]] = metaCellIds[i];

        qDebug() << "Found" << nwbFiles.size() << "NWB files, attempting to load them..";

        // Load NWB files and add them to dataset
        LoadInfo loadInfo;
        loadInfo.failedSweepPath = ctx.config.ephysTraces.failedSweepPath;
        NWBLoader loader;
        for (int i = 0; i < nwbFiles.size(); i++)
        {
            Experiment experiment;

            QString fileName = nwbFiles[i];

            qDebug() << "Loading NWB file: " << ephysTracesDir.filePath(fileName);

            loader.LoadNWB(ephysTracesDir.filePath(fileName), experiment, loadInfo);

            QString specimenName = fileName;
            specimenName.chop(4); // Cut off the .nwb part
            qDebug() << "Specimen name: " << specimenName;
            if (specimenNameToCellId.find(specimenName) == specimenNameToCellId.end())
                continue;

            _ephysTraceCellIds.push_back(specimenNameToCellId[specimenName]);
            qDebug() << "Cell ID: " << _ephysTraceCellIds[_ephysTraceCellIds.size() - 1];
            ctx.ephysTraces->addExperiment(std::move(experiment));
        }

        qDebug() << "Ignored stimsets: " << loadInfo.ignoredStimsets;
        qDebug() << "Loaded stimsets: " << loadInfo.loadedStimsets;

        events().notifyDatasetAdded(ctx.ephysTraces);
        events().notifyDatasetDataChanged(ctx.ephysTraces);
    }

    void LoadMorphologyCells(PipelineContext& ctx, QDir dir)
    {
        Timer timer("SWC Morphology Loading");

        // Load morphology cells
        ctx.cellMorphologies = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", "Cell Morphologies", mv::Dataset<DatasetImpl>(), "", false);
        ctx.cellMorphologies->setProperty("PatchSeqType", "Morphologies");
#if defined(DALLEYLEE) || defined(WALEBOER)
        ctx.cellMorphologies->setProperty("isCortical", true);
#endif

        QDir morphologyDir(dir);

        QStringList swcFiles = morphologyDir.entryList(QStringList() << "*.swc" << "*.SWC", QDir::Files);

        QStringList cellIds;
        std::vector<CellMorphology> cellMorphologies(swcFiles.size());

        SWCLoader loader;
        for (int i = 0; i < cellMorphologies.size(); i++)
        {
            QString swcFile = swcFiles[i];
            CellMorphology& cellMorphology = cellMorphologies[i];

            loader.LoadSWC(morphologyDir.filePath(swcFile), cellMorphology);

            cellMorphology.findCentroid();
            cellMorphology.findExtents();
            cellMorphology.cellTypeColor.set(0.11f, 0.79f, 0);

            cellIds.append(QFileInfo(swcFile).baseName());
        }

        ctx.cellMorphologies->setCellIdentifiers(cellIds);
        ctx.cellMorphologies->setData(cellMorphologies);

        events().notifyDatasetAdded(ctx.cellMorphologies);
        events().notifyDatasetDataChanged(ctx.cellMorphologies);
        events().notifyDatasetDataDimensionsChanged(ctx.cellMorphologies);
    }
};

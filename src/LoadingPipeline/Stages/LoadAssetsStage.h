#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "Electrophysiology/NWBLoader.h"
#include "Morphology/SWCLoader.h"

#include <EphysData/EphysData.h>
#include <CellMorphologyData/CellMorphologyData.h>

#include <util/Timer.h>

#include <QDir>
#include <QFileInfo>
#include <QHash>
#include <QSet>

class LoadAssetsStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "LoadAssets";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        int loadedCount = 0;
        int skippedCount = 0;
        int failedCount = 0;

        if (ctx.metadata.obs.index.empty())
        {
            ctx.result.Warning(Name(), "Metadata table is empty; asset sources cannot be linked to cells.", "Run CollectMetadataStage before LoadAssetsStage.");
        }

        if (ctx.config.ephysTraces)
        {
            if (LoadEphysTraces(ctx, *ctx.config.ephysTraces))
                ++loadedCount;
            else
                ++failedCount;
        }
        else
        {
            ctx.result.Info(Name(), "Ephys traces source is not configured.", "sources.ephys_traces");
            ++skippedCount;
        }

        if (ctx.config.morphologyReconstructions)
        {
            if (LoadMorphologyReconstructions(ctx, *ctx.config.morphologyReconstructions))
                ++loadedCount;
            else
                ++failedCount;
        }
        else
        {
            ctx.result.Info(Name(), "Morphology reconstructions source is not configured.", "sources.morphology_reconstructions");
            ++skippedCount;
        }

        ctx.result.Info(Name(), "Asset loading complete.", QString("loaded=%1, skipped=%2, failed=%3")
            .arg(loadedCount)
            .arg(skippedCount)
            .arg(failedCount));

        return failedCount == 0;
    }

private:
    QHash<QString, QString> BuildFilenameValueToCellIdMap(const AnnotationTable& obs, const QString& filenameMetadataColumn) const
    {
        QHash<QString, QString> result;

        if (filenameMetadataColumn.isEmpty())
            return result;

        if (filenameMetadataColumn == obs.indexName)
        {
            for (const QString& cellId : obs.index)
                result.insert(cellId.trimmed(), cellId.trimmed());

            return result;
        }

        const int columnIndex = FindColumn(obs, filenameMetadataColumn);

        if (columnIndex < 0)
            return result;

        const std::vector<QString>& filenameValues = obs.values[static_cast<size_t>(columnIndex)];

        const size_t rowCount = std::min(obs.index.size(), filenameValues.size());

        for (size_t row = 0; row < rowCount; ++row)
        {
            const QString filenameValue = filenameValues[row].trimmed();
            const QString cellId = obs.index[row].trimmed();

            if (filenameValue.isEmpty() || cellId.isEmpty())
                continue;

            result.insert(filenameValue, cellId);
        }

        return result;
    }

    int FindColumn(const AnnotationTable& table, const QString& columnName) const
    {
        for (size_t col = 0; col < table.columnNames.size(); ++col)
        {
            if (table.columnNames[col] == columnName)
                return static_cast<int>(col);
        }

        return -1;
    }

    QString FileStem(const QString& fileName) const
    {
        return QFileInfo(fileName).completeBaseName();
    }

    bool ValidateDirectory(PipelineContext& ctx, const QString& configPath, const QString& directoryPath) const
    {
        if (directoryPath.isEmpty())
        {
            ctx.result.Warning(Name(), "Configured asset directory path is empty.", configPath);
            return false;
        }

        const QFileInfo info(directoryPath);

        if (!info.exists())
        {
            ctx.result.Warning(Name(), "Configured asset directory does not exist.", QString("%1 = %2").arg(configPath, directoryPath));
            return false;
        }

        if (!info.isDir())
        {
            ctx.result.Warning(Name(), "Configured asset path exists but is not a directory.", QString("%1 = %2").arg(configPath, directoryPath));
            return false;
        }

        if (!info.isReadable())
        {
            ctx.result.Warning(Name(), "Configured asset directory is not readable.", QString("%1 = %2").arg(configPath, directoryPath));
            return false;
        }

        return true;
    }

    bool LoadEphysTraces(PipelineContext& ctx, const config::EphysTracesSource& source) const
    {
        if (!ValidateDirectory(ctx, "sources.ephys_traces.directory", source.directory))
            return false;

        const QHash<QString, QString> filenameValueToCellId = BuildFilenameValueToCellIdMap(ctx.metadata.obs, source.filenameMetadataColumn);

        if (filenameValueToCellId.isEmpty())
        {
            ctx.result.Warning(Name(), "Could not build filename-to-cell_id map for ephys traces.", QString("filename_metadata_column='%1'").arg(source.filenameMetadataColumn));
            return false;
        }

        QDir ephysTracesDir(source.directory);
        const QStringList nwbFiles = ephysTracesDir.entryList(QStringList() << "*.nwb" << "*.NWB", QDir::Files);

        if (nwbFiles.isEmpty())
        {
            ctx.result.Warning(Name(), "No NWB files found in ephys traces directory.", source.directory);
            return false;
        }

        const QString displayName = source.displayName.isEmpty() ? QString("Ephys Traces") : source.displayName;

        ctx.ephysTraces = mv::data().createDataset<EphysExperiments>("Electrophysiology Data", displayName, mv::Dataset<mv::DatasetImpl>(), "", false);

        ctx.ephysTraces->setProperty("PatchSeqType", "EphysTraces");

        LoadInfo loadInfo;
        loadInfo.failedSweepPath = source.failedSweepsPath;

        NWBLoader loader;

        int loaded = 0;
        int skipped = 0;
        int failed = 0;

        qDebug() << "Found" << nwbFiles.size() << "NWB files, attempting to load them from" << source.directory;

        ctx.ephysTraceCellIds.clear();
        for (const QString& fileName : nwbFiles)
        {
            const QString fileStem = FileStem(fileName);

            if (!filenameValueToCellId.contains(fileStem))
            {
                ++skipped;
                ctx.result.Info(Name(), "Skipping NWB file because it does not match metadata.", fileName);
                continue;
            }

            const QString filePath = ephysTracesDir.filePath(fileName);

            try
            {
                Experiment experiment;

                qDebug() << "Loading NWB file:" << filePath;

                const QString cellId = filenameValueToCellId.value(fileStem);

                loader.LoadNWB(filePath, experiment, loadInfo);

                ctx.ephysTraces->addExperiment(std::move(experiment));
                ctx.ephysTraceCellIds.push_back(cellId);

                ++loaded;
            }
            catch (const std::exception& e)
            {
                ++failed;
                ctx.result.Warning(Name(), "Failed to load NWB file.", QString("%1; error=%2").arg(filePath, e.what()));
            }
            catch (...)
            {
                ++failed;
                ctx.result.Warning(Name(), "Failed to load NWB file due to an unknown error.", filePath);
            }
        }

        if (loaded == 0)
        {
            ctx.result.Warning(Name(), "No ephys traces were loaded.", QString("files=%1, skipped=%2, failed=%3").arg(nwbFiles.size()).arg(skipped).arg(failed));
            return false;
        }

        qDebug() << "Ignored stimsets:" << loadInfo.ignoredStimsets;
        qDebug() << "Loaded stimsets:" << loadInfo.loadedStimsets;

        mv::events().notifyDatasetAdded(ctx.ephysTraces);
        mv::events().notifyDatasetDataChanged(ctx.ephysTraces);
        mv::events().notifyDatasetDataDimensionsChanged(ctx.ephysTraces);

        ctx.result.Info(Name(), "Loaded ephys traces.", QString("loaded=%1, skipped=%2, failed=%3").arg(loaded).arg(skipped).arg(failed));

        return true;
    }

    bool LoadMorphologyReconstructions(PipelineContext& ctx, const config::MorphologyReconstructionsSource& source) const
    {
        if (!ValidateDirectory(ctx, "sources.morphology_reconstructions.directory", source.directory))
            return false;

        const QHash<QString, QString> filenameValueToCellId = BuildFilenameValueToCellIdMap(ctx.metadata.obs, source.filenameMetadataColumn);

        if (filenameValueToCellId.isEmpty())
        {
            ctx.result.Warning(Name(), "Could not build filename-to-cell_id map for morphology reconstructions.", QString("filename_metadata_column='%1'").arg(source.filenameMetadataColumn));
            return false;
        }

        QDir morphologyDir(source.directory);
        const QStringList swcFiles = morphologyDir.entryList(QStringList() << "*.swc" << "*.SWC", QDir::Files);

        if (swcFiles.isEmpty())
        {
            ctx.result.Warning(Name(), "No SWC files found in morphology reconstruction directory.", source.directory);
            return false;
        }

        Timer timer("SWC Morphology Loading");

        const QString displayName = source.displayName.isEmpty() ? QString("Cell Morphologies") : source.displayName;

        ctx.cellMorphologies = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", displayName, mv::Dataset<mv::DatasetImpl>(), "", false);

        ctx.cellMorphologies->setProperty("PatchSeqType", "Morphologies");

#if defined(DALLEYLEE) || defined(WALEBOER)
        ctx.cellMorphologies->setProperty("isCortical", true);
#endif

        QStringList cellIds;
        std::vector<CellMorphology> loadedMorphologies;

        cellIds.reserve(swcFiles.size());
        loadedMorphologies.reserve(static_cast<size_t>(swcFiles.size()));

        SWCLoader loader;

        int loaded = 0;
        int skipped = 0;
        int failed = 0;

        qDebug() << "Found" << swcFiles.size() << "SWC files, attempting to load them from" << source.directory;
        
        ctx.morphologyCellIds.clear();
        for (const QString& fileName : swcFiles)
        {
            const QString fileStem = FileStem(fileName);

            if (!filenameValueToCellId.contains(fileStem))
            {
                ++skipped;
                ctx.result.Info(Name(), "Skipping SWC file because it does not match metadata.", fileName);
                continue;
            }

            const QString filePath = morphologyDir.filePath(fileName);

            try
            {
                CellMorphology cellMorphology;

                qDebug() << "Loading SWC file:" << filePath;

                loader.LoadSWC(filePath, cellMorphology);

                cellMorphology.findCentroid();
                cellMorphology.findExtents();
                cellMorphology.cellTypeColor.set(0.11f, 0.79f, 0.0f);

                cellIds.append(filenameValueToCellId.value(fileStem));
                loadedMorphologies.push_back(std::move(cellMorphology));

                ++loaded;
            }
            catch (const std::exception& e)
            {
                ++failed;
                ctx.result.Warning(Name(), "Failed to load SWC file.", QString("%1; error=%2").arg(filePath, e.what()));
            }
            catch (...)
            {
                ++failed;
                ctx.result.Warning(Name(), "Failed to load SWC file due to an unknown error.", filePath);
            }
        }

        if (loaded == 0)
        {
            ctx.result.Warning(Name(), "No morphology reconstructions were loaded.", QString("files=%1, skipped=%2, failed=%3").arg(swcFiles.size()).arg(skipped).arg(failed));
            return false;
        }

        ctx.cellMorphologies->setCellIdentifiers(cellIds);
        ctx.cellMorphologies->setData(loadedMorphologies);

        std::vector<QString> morphologyCellIds;
        for (int i = 0; i < cellIds.size(); i++)
            morphologyCellIds.push_back(cellIds[i]);
        ctx.morphologyCellIds = morphologyCellIds;

        mv::events().notifyDatasetAdded(ctx.cellMorphologies);
        mv::events().notifyDatasetDataChanged(ctx.cellMorphologies);
        mv::events().notifyDatasetDataDimensionsChanged(ctx.cellMorphologies);

        ctx.result.Info(Name(), "Loaded morphology reconstructions.", QString("loaded=%1, skipped=%2, failed=%3").arg(loaded).arg(skipped).arg(failed));

        return true;
    }

//    void LoadEphysTraces(PipelineContext& ctx, QDir dir)
//    {
//        // Create ephys dataset
//        ctx.ephysTraces = mv::data().createDataset<EphysExperiments>("Electrophysiology Data", "Ephys Traces", mv::Dataset<DatasetImpl>(), "", false);
//        ctx.ephysTraces->setProperty("PatchSeqType", "EphysTraces");
//
//        // Find all .nwb files in given directory
//        QDir ephysTracesDir(dir);
//
//        QStringList nwbFiles = ephysTracesDir.entryList(QStringList() << "*.nwb" << "*.NWB", QDir::Files);
//
//        // Map metadata cell names to cell_ids
//        std::vector<QString> metaCellSpecimenNames = _metadataDf[CELL_NAME_TAG];
//        std::vector<QString> metaCellIds = _metadataDf[CELL_ID_TAG];
//
//        std::unordered_map<QString, QString> specimenNameToCellId;
//        for (size_t i = 0; i < metaCellSpecimenNames.size(); ++i)
//            specimenNameToCellId[metaCellSpecimenNames[i]] = metaCellIds[i];
//
//        qDebug() << "Found" << nwbFiles.size() << "NWB files, attempting to load them..";
//
//        // Load NWB files and add them to dataset
//        LoadInfo loadInfo;
//        loadInfo.failedSweepPath = ctx.config.ephysTraces.failedSweepPath;
//        NWBLoader loader;
//        for (int i = 0; i < nwbFiles.size(); i++)
//        {
//            Experiment experiment;
//
//            QString fileName = nwbFiles[i];
//
//            qDebug() << "Loading NWB file: " << ephysTracesDir.filePath(fileName);
//
//            loader.LoadNWB(ephysTracesDir.filePath(fileName), experiment, loadInfo);
//
//            QString specimenName = fileName;
//            specimenName.chop(4); // Cut off the .nwb part
//            qDebug() << "Specimen name: " << specimenName;
//            if (specimenNameToCellId.find(specimenName) == specimenNameToCellId.end())
//                continue;
//
//            _ephysTraceCellIds.push_back(specimenNameToCellId[specimenName]);
//            qDebug() << "Cell ID: " << _ephysTraceCellIds[_ephysTraceCellIds.size() - 1];
//            ctx.ephysTraces->addExperiment(std::move(experiment));
//        }
//
//        qDebug() << "Ignored stimsets: " << loadInfo.ignoredStimsets;
//        qDebug() << "Loaded stimsets: " << loadInfo.loadedStimsets;
//
//        events().notifyDatasetAdded(ctx.ephysTraces);
//        events().notifyDatasetDataChanged(ctx.ephysTraces);
//    }
//
//    void LoadMorphologyCells(PipelineContext& ctx, QDir dir)
//    {
//        Timer timer("SWC Morphology Loading");
//
//        // Load morphology cells
//        ctx.cellMorphologies = mv::data().createDataset<CellMorphologies>("Cell Morphology Data", "Cell Morphologies", mv::Dataset<DatasetImpl>(), "", false);
//        ctx.cellMorphologies->setProperty("PatchSeqType", "Morphologies");
//#if defined(DALLEYLEE) || defined(WALEBOER)
//        ctx.cellMorphologies->setProperty("isCortical", true);
//#endif
//
//        QDir morphologyDir(dir);
//
//        QStringList swcFiles = morphologyDir.entryList(QStringList() << "*.swc" << "*.SWC", QDir::Files);
//
//        QStringList cellIds;
//        std::vector<CellMorphology> cellMorphologies(swcFiles.size());
//
//        SWCLoader loader;
//        for (int i = 0; i < cellMorphologies.size(); i++)
//        {
//            QString swcFile = swcFiles[i];
//            CellMorphology& cellMorphology = cellMorphologies[i];
//
//            loader.LoadSWC(morphologyDir.filePath(swcFile), cellMorphology);
//
//            cellMorphology.findCentroid();
//            cellMorphology.findExtents();
//            cellMorphology.cellTypeColor.set(0.11f, 0.79f, 0);
//
//            cellIds.append(QFileInfo(swcFile).baseName());
//        }
//
//        ctx.cellMorphologies->setCellIdentifiers(cellIds);
//        ctx.cellMorphologies->setData(cellMorphologies);
//
//        events().notifyDatasetAdded(ctx.cellMorphologies);
//        events().notifyDatasetDataChanged(ctx.cellMorphologies);
//        events().notifyDatasetDataDimensionsChanged(ctx.cellMorphologies);
//    }
};

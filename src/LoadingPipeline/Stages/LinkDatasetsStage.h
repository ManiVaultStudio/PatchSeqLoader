#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "Config/ConfigSchema.h"

#include <PointData/PointData.h>
#include <SelectionGroup.h>

#include <QString>
#include <optional>

class LinkDatasetsStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "LinkDatasets";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        int linkedCount = 0;

        AddDatasetToSelectionGroupIfAvailable(ctx, config::keys::sources::Rna, linkedCount);
        AddDatasetToSelectionGroupIfAvailable(ctx, config::keys::sources::Ephys, linkedCount);
        AddDatasetToSelectionGroupIfAvailable(ctx, config::keys::sources::Morphology, linkedCount);

        AddEmbeddingToSelectionGroupIfAvailable(ctx, config::keys::embeddings::RnaUmap, linkedCount);
        AddEmbeddingToSelectionGroupIfAvailable(ctx, config::keys::embeddings::EphysUmap, linkedCount);
        AddEmbeddingToSelectionGroupIfAvailable(ctx, config::keys::embeddings::MorphoUmap, linkedCount);

        AddMetadataToSelectionGroupIfAvailable(ctx, config::keys::sources::Metadata, linkedCount);

        // Link ephys and morpho reconstructions
        AddEphysTracesToSelectionGroupIfAvailable(ctx, linkedCount);
        AddMorphologiesToSelectionGroupIfAvailable(ctx, linkedCount);

        events().addSelectionGroup(ctx.selectionGroup);

        ctx.result.Info(Name(), "Dataset linking complete.", QString("linked=%1").arg(linkedCount));

        return true;
    }

private:
    void AddDatasetToSelectionGroupIfAvailable(PipelineContext& ctx, const QString& sourceName, int& linkedCount)
    {
        if (!ctx.featureDatasets.contains(sourceName))
            return;

        if (!ctx.normalizedTables.contains(sourceName))
        {
            ctx.result.Warning(Name(), "Dataset is available for linking, but no source table", QString("%1").arg(sourceName));
            return;
        }

        const AnnotatedData& data = ctx.normalizedTables[sourceName];
        const mv::Dataset<Points>& dataset = ctx.featureDatasets[sourceName];

        if (data.obs.index.size() != dataset->getNumPoints())
        {
            ctx.result.Warning(Name(), "Dataset index size is different than number of points.", QString("%1").arg(sourceName));
            return;
        }

        //BiMap biMap;
        //std::vector<uint32_t> indices(dataset->getNumPoints());
        //std::iota(indices.begin(), indices.end(), 0);
        //biMap.addKeyValuePairs(data.obs.index, indices);

        ctx.selectionGroup.addDataset(dataset, data.obs.index);

        ++linkedCount;
    }

    void AddEmbeddingToSelectionGroupIfAvailable(PipelineContext& ctx, const QString& sourceName, int& linkedCount)
    {
        if (!ctx.embeddingDatasets.contains(sourceName))
            return;

        if (!ctx.normalizedTables.contains(sourceName))
        {
            ctx.result.Warning(Name(), "Dataset is available for linking, but no source table", QString("%1").arg(sourceName));
            return;
        }

        const AnnotatedData& data = ctx.normalizedTables[sourceName];
        const mv::Dataset<Points>& dataset = ctx.embeddingDatasets[sourceName];

        if (data.obs.index.size() != dataset->getNumPoints())
        {
            ctx.result.Warning(Name(), "Dataset index size is different than number of points.", QString("%1").arg(sourceName));
            return;
        }

        ctx.selectionGroup.addDataset(dataset, data.obs.index);

        ++linkedCount;
    }

    void AddMetadataToSelectionGroupIfAvailable(PipelineContext& ctx, const QString& sourceName, int& linkedCount)
    {
        if (!ctx.textDatasets.contains(sourceName))
        {
            ctx.result.Warning(Name(), "Context text datasets does not contain", QString("%1").arg(sourceName));
            return;
        }

        //if (!ctx.normalizedTables.contains(sourceName))
        //{
        //    ctx.result.Warning(Name(), "Dataset is available for linking, but no source table", QString("%1").arg(sourceName));
        //    return;
        //}

        const AnnotatedData& data = ctx.metadata;
        const mv::Dataset<Text>& dataset = ctx.textDatasets[sourceName];

        ctx.selectionGroup.addDataset(dataset, data.obs.index);

        ++linkedCount;
    }

    void AddEphysTracesToSelectionGroupIfAvailable(PipelineContext& ctx, int& linkedCount)
    {
        if (!ctx.ephysTraces.isValid())
            return;

        if (ctx.ephysTraceCellIds.empty())
        {
            ctx.result.Warning(Name(), "Ephys traces dataset is available for linking, but no cell IDs were stored.");
            return;
        }

        ctx.selectionGroup.addDataset(ctx.ephysTraces, ctx.ephysTraceCellIds);

        ++linkedCount;
    }

    void AddMorphologiesToSelectionGroupIfAvailable(PipelineContext& ctx, int& linkedCount)
    {
        if (!ctx.cellMorphologies.isValid())
            return;

        if (ctx.morphologyCellIds.empty())
        {
            ctx.result.Warning(Name(), "Morphology dataset is available for linking, but no cell IDs were stored.");
            return;
        }

        ctx.selectionGroup.addDataset(ctx.cellMorphologies, ctx.morphologyCellIds);

        ++linkedCount;
    }
};

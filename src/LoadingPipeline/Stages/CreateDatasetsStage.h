#pragma once

#include "LoadingPipeline/PipelineStage.h"
#include "LoadingPipeline/PipelineContext.h"

#include "Config/ConfigSchema.h"

#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>

#include <QString>
#include <optional>

namespace
{
    std::map<QString, std::vector<unsigned int>> MakeClustersFromList(std::vector<QString> list)
    {
        std::map<QString, std::vector<unsigned int>> clusterData;

        for (int i = 0; i < list.size(); i++)
        {
            const QString& str = list[i];
            clusterData[str].push_back(i);
        }
        return clusterData;
    }
}

class CreateDatasetsStage final : public PipelineStage
{
public:
    QString Name() const override
    {
        return "CreateDatasets";
    }

    bool Run(PipelineContext& ctx) override
    {
        qDebug() << "Running pipeline stage: " << Name();

        int createdCount = 0;

        CreateTableDatasetIfAvailable(ctx, config::keys::sources::Rna, ctx.config.rna, createdCount);
        CreateTableDatasetIfAvailable(ctx, config::keys::sources::Ephys, ctx.config.ephys, createdCount);
        CreateTableDatasetIfAvailable(ctx, config::keys::sources::Morphology, ctx.config.morphology, createdCount);

        CreateEmbeddingDatasetIfAvailable(ctx, config::keys::embeddings::RnaUmap, ctx.config.rnaUmap, ctx.featureDatasets[config::keys::sources::Rna], createdCount);
        CreateEmbeddingDatasetIfAvailable(ctx, config::keys::embeddings::EphysUmap, ctx.config.ephysUmap, ctx.featureDatasets[config::keys::sources::Ephys], createdCount);
        CreateEmbeddingDatasetIfAvailable(ctx, config::keys::embeddings::MorphoUmap, ctx.config.morphoUmap, ctx.featureDatasets[config::keys::sources::Morphology], createdCount);

        CreateMetadataDataset(ctx, config::keys::sources::Metadata, ctx.metadata, ctx.config.metadata.value());

        //for (auto it = ctx.config.extraEmbeddings.begin(); it != ctx.config.extraEmbeddings.end(); ++it)
        //    CreateExtraEmbeddingDataset(ctx, it.key(), it.value(), createdCount);

        ctx.result.Info(Name(), "Dataset creation complete.", QString("created=%1").arg(createdCount));

        return true;
    }

private:
    void CreateTableDatasetIfAvailable(PipelineContext& ctx, const QString& sourceName, const std::optional<config::TableSource>& source, int& createdCount) const
    {
        if (!source)
            return;

        if (!ctx.normalizedTables.contains(sourceName))
            return;

        AnnotatedData& data = ctx.normalizedTables[sourceName];

        if (sourceName == config::keys::sources::Metadata)
            CreateMetadataDataset(ctx, sourceName, data, source.value());
        else
            CreateFeatureDataset(ctx, sourceName, data, source.value());

        ++createdCount;
    }

    void CreateEmbeddingDatasetIfAvailable(PipelineContext& ctx, const QString& embeddingName, const std::optional<config::TableSource>& embedding, const mv::Dataset<Points>& parent, int& createdCount) const
    {
        if (!embedding)
            return;

        if (!ctx.normalizedTables.contains(embeddingName))
            return;

        AnnotatedData& data = ctx.normalizedTables[embeddingName];

        CreateEmbeddingDataset(ctx, embeddingName, data, embedding.value(), parent);

        ++createdCount;
    }

    //void CreateExtraEmbeddingDataset(PipelineContext& ctx, const QString& embeddingName, const config::Embedding& embedding, int& createdCount) const
    //{
    //    if (!ctx.embeddingTables.contains(embeddingName))
    //        return;

    //    const AnnotatedData& table = ctx.embeddingTables[embeddingName];

    //    createEmbeddingDataset(ctx, embeddingName, table, embedding);
    //    ++createdCount;
    //}

    void CreateMetadataDataset(PipelineContext& ctx, const QString& sourceName, AnnotatedData& data, const config::TableSource& config) const
    {
        auto textDataset = mv::data().createDataset<Text>("Text", config.displayName, mv::Dataset<mv::DatasetImpl>(), "", false);
        textDataset->setProperty("PatchSeqType", sourceName);
        textDataset->addColumn(data.obs.indexName, data.obs.index);
        for (int i = 0; i < data.obs.columnNames.size(); i++)
            textDataset->addColumn(data.obs.columnNames[i], data.obs.values[i]);

        mv::events().notifyDatasetAdded(textDataset);
        mv::events().notifyDatasetDataChanged(textDataset);
        mv::events().notifyDatasetDataDimensionsChanged(textDataset);

        ctx.textDatasets.insert(sourceName, textDataset);

        for (int i = 0; i < data.obs.values.size(); i++)
        {
            mv::Dataset<Clusters> clusterData = mv::data().createDataset<Clusters>("Cluster", data.obs.columnNames[i], textDataset);

            const std::vector<QString>& clusterAsList = data.obs.values[i];
            std::map<QString, std::vector<unsigned int>> clusterMap = MakeClustersFromList(clusterAsList);

            const auto colorMapIt = ctx.metadataColorMaps.find(data.obs.columnNames[i]);

            const bool hasColorMap = colorMapIt != ctx.metadataColorMaps.end();

            for (auto& kv : clusterMap)
            {
                Cluster cluster;

                cluster.setName(kv.first);
                cluster.setIndices(kv.second);

                if (hasColorMap)
                {
                    const QHash<QString, QColor>& colorMap = colorMapIt.value();
                    const auto colorIt = colorMap.find(kv.first);

                    if (colorIt != colorMap.end())
                        cluster.setColor(colorIt.value());
                }

                clusterData->addCluster(cluster);
            }

            if (!hasColorMap)
                Cluster::colorizeClusters(clusterData->getClusters());

            mv::events().notifyDatasetDataChanged(clusterData);
            mv::events().notifyDatasetDataDimensionsChanged(clusterData);
        }
    }

    void CreateFeatureDataset(PipelineContext& ctx, const QString& sourceName, const AnnotatedData& data, const config::TableSource& config) const
    {
        auto featureDataset = mv::data().createDataset<Points>("Points", config.displayName, mv::Dataset<mv::DatasetImpl>(), "", false);
        featureDataset->setProperty("PatchSeqType", sourceName);
        featureDataset->setData(std::move(data.X.values), data.X.columnCount);
        featureDataset->setDimensionNames(data.var.index);

        mv::events().notifyDatasetAdded(featureDataset);
        mv::events().notifyDatasetDataChanged(featureDataset);
        mv::events().notifyDatasetDataDimensionsChanged(featureDataset);

        ctx.featureDatasets.insert(sourceName, featureDataset);
    }

    void CreateEmbeddingDataset(PipelineContext& ctx, const QString& embeddingName, const AnnotatedData& data, const config::TableSource& embedding, const mv::Dataset<Points>& parent) const
    {
        // TODO: create 2D Points/embedding dataset here.
        auto embeddingDataset = mv::data().createDataset<Points>("Points", embedding.displayName, parent, "", false);
        embeddingDataset->setProperty("PatchSeqType", embeddingName);
        embeddingDataset->setData(data.X.values, data.X.columnCount);
        embeddingDataset->setDimensionNames(data.var.index);

        mv::events().notifyDatasetAdded(embeddingDataset);
        mv::events().notifyDatasetDataChanged(embeddingDataset);
        mv::events().notifyDatasetDataDimensionsChanged(embeddingDataset);

        ctx.embeddingDatasets.insert(embeddingName, embeddingDataset);
    }
};

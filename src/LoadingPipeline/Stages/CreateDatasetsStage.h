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
    QStringList OrderedClusterLabels(const std::map<QString, std::vector<unsigned int>>& clusterMap, const QStringList& priority)
    {
        QStringList ordered;
        QSet<QString> added;

        // Add configured cluster labels first, if they actually exist.
        for (const QString& label : priority)
        {
            if (added.contains(label))
                continue;

            if (clusterMap.find(label) != clusterMap.end())
            {
                ordered.push_back(label);
                added.insert(label);
            }
        }

        // Add all remaining labels in the map's default order.
        for (const auto& kv : clusterMap)
        {
            if (!added.contains(kv.first))
            {
                ordered.push_back(kv.first);
                added.insert(kv.first);
            }
        }

        return ordered;
    }

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

    QStringList OrderedMetadataColumns(const AnnotationTable& obs, const QStringList& priority)
    {
        QStringList ordered;
        QSet<QString> added;

        auto findColumn = [&](const QString& name) -> int
            {
                for (int i = 0; i < static_cast<int>(obs.columnNames.size()); ++i)
                {
                    if (obs.columnNames[i] == name)
                        return i;
                }

                return -1;
            };

        // First add priority columns that actually exist.
        for (const QString& columnName : priority)
        {
            if (added.contains(columnName))
                continue;

            if (findColumn(columnName) >= 0)
            {
                ordered.push_back(columnName);
                added.insert(columnName);
            }
        }

        // Then append every remaining metadata column in its existing order.
        for (const QString& columnName : obs.columnNames)
        {
            if (!added.contains(columnName))
            {
                ordered.push_back(columnName);
                added.insert(columnName);
            }
        }

        return ordered;
    }

    int FindObsColumn(const AnnotationTable& obs, const QString& columnName)
    {
        for (int i = 0; i < static_cast<int>(obs.columnNames.size()); ++i)
        {
            if (obs.columnNames[i] == columnName)
                return i;
        }

        return -1;
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
        auto textDataset = mv::data().createDataset<Text>(
            "Text",
            config.displayName,
            mv::Dataset<mv::DatasetImpl>(),
            "",
            false);

        textDataset->setProperty("PatchSeqType", sourceName);

        textDataset->addColumn(data.obs.indexName, data.obs.index);

        // Keep the text table itself in the actual metadata column order.
        for (int i = 0; i < static_cast<int>(data.obs.columnNames.size()); ++i)
            textDataset->addColumn(data.obs.columnNames[i], data.obs.values[i]);

        mv::events().notifyDatasetAdded(textDataset);
        mv::events().notifyDatasetDataChanged(textDataset);
        mv::events().notifyDatasetDataDimensionsChanged(textDataset);

        ctx.textDatasets.insert(sourceName, textDataset);

        const QStringList orderedColumnNames =
            OrderedMetadataColumns(data.obs, config.obsPriority);

        for (const QString& columnName : orderedColumnNames)
        {
            const int columnIndex = FindObsColumn(data.obs, columnName);

            if (columnIndex < 0)
                continue;

            mv::Dataset<Clusters> clusterData =
                mv::data().createDataset<Clusters>("Cluster", columnName, textDataset);

            const std::vector<QString>& clusterAsList = data.obs.values[columnIndex];

            std::map<QString, std::vector<unsigned int>> clusterMap =
                MakeClustersFromList(clusterAsList);

            const auto colorMapIt = ctx.metadataColorMaps.find(columnName);
            const bool hasColorMap = colorMapIt != ctx.metadataColorMaps.end();

            const QStringList clusterPriority = config.clusterPriority.value(columnName);

            const QStringList orderedClusterLabels = OrderedClusterLabels(clusterMap, clusterPriority);

            for (const QString& label : orderedClusterLabels)
            {
                const auto clusterIt = clusterMap.find(label);

                if (clusterIt == clusterMap.end())
                    continue;

                Cluster cluster;

                cluster.setName(label);
                cluster.setIndices(clusterIt->second);

                if (hasColorMap)
                {
                    const QHash<QString, QColor>& colorMap = colorMapIt.value();
                    const auto colorIt = colorMap.find(label);

                    if (colorIt != colorMap.end())
                        cluster.setColor(colorIt.value());
                }

                clusterData->addCluster(cluster);
            }

            if (!hasColorMap)
                Cluster::colorizeClusters(clusterData->getClusters());

            // Set empty label clusters to grey
            for (Cluster& cluster : clusterData->getClusters())
            {
                if (cluster.getName().isEmpty())
                {
                    cluster.setColor(QColor(220, 220, 220));
                    cluster.setName("Unassigned");
                }
            }

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

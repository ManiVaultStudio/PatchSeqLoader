#pragma once

#include "Config/ConfigSchema.h"

#include <QString>
#include <QStringList>
#include <QMap>
#include <optional>

namespace config
{
    enum class ObsColumnMode
    {
        Explicit,
        All
    };

    struct TableSource
    {
        QString path;
        QString displayName;
        QString index;

        ObsColumnMode obsColumnMode = ObsColumnMode::Explicit;
        QStringList obsColumns; // Optional

        // Preferred order for creating metadata-derived datasets.
        // Columns not listed here are appended afterward in their existing order
        QStringList obsPriority;

        // Metadata column -> preferred cluster/category order.
        // Unlisted cluster labels are appended afterward.
        QMap<QString, QStringList> clusterPriority;
    };

    struct EphysTracesSource
    {
        QString directory;
        QString failedSweepsPath;
        QString displayName;
        QString filenameMetadataColumn;
    };

    struct MorphologyReconstructionsSource
    {
        QString directory;
        QString displayName;
        QString filenameMetadataColumn;
    };

    struct TaxonomySource
    {
        QString path;
        QString displayName;
    };

    struct File
    {
    public:
        bool Load(QString filePath);

    public:
        /** String identifying which file schema we are following */
        QString format;
        /** Version of the given schema */
        QString version;
        /** Name of the dataset this config is for */
        QString datasetId;
        /** Determines whether the data is from the brain cortex, which influences our drawing of layers */
        bool isCortical = false;

        /** Table source for transcriptomic gene-expression data */
        std::optional<config::TableSource> rna;
        /** Table source for electrophysiological feature data */
        std::optional<config::TableSource> ephys;
        /** Table source for morphological feature data */
        std::optional<config::TableSource> morphology;
        /** Table source for associated metadata */
        std::optional<config::TableSource> metadata;

        /** Source properties for electrophysiological experiments */
        std::optional<config::EphysTracesSource> ephysTraces;
        /** Source properties for morphological reconstructions */
        std::optional<config::MorphologyReconstructionsSource> morphologyReconstructions;

        /** Table source for transcriptomic gene-expression UMAP */
        std::optional<config::TableSource> rnaUmap;
        /** Table source for electrophysiological feature UMAP */
        std::optional<config::TableSource> ephysUmap;
        /** Table source for morphological feature UMAP */
        std::optional<config::TableSource> morphoUmap;

        /** Room for additional UMAPs or other embeddings */
        QMap<QString, config::TableSource> extraEmbeddings;

        /** Source properties of file describing global taxonomy and coloring */
        QMap<QString, TaxonomySource> taxonomies;
    };

} // namespace config

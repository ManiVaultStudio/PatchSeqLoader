#pragma once

#include <QString>

namespace config::keys
{
    // Top-level
    inline const QString Format = "format";
    inline const QString Version = "version";
    inline const QString DatasetId = "dataset_id";

    inline const QString Sources = "sources";
    inline const QString Assets = "assets";
    inline const QString Directories = "directories";
    inline const QString Embeddings = "embeddings";

    // Common object fields
    inline const QString Path = "path";
    inline const QString Index = "index";
    inline const QString ObsColumns = "obs_columns";
    inline const QString DisplayName = "display_name";
    inline const QString XColumn = "x_column";
    inline const QString YColumn = "y_column";

    // Source keys
    namespace sources
    {
        inline const QString Rna = "rna";
        inline const QString Ephys = "ephys";
        inline const QString Morphology = "morphology";
        inline const QString Metadata = "metadata";
    }

    // Asset directory keys
    namespace assets
    {
        inline const QString MorphologyReconstruction = "morphology_reconstruction";
        inline const QString EphysTraces = "ephys_traces";
    }

    // Embedding keys
    namespace embeddings
    {
        inline const QString RnaUmap = "rna_umap";
        inline const QString EphysUmap = "ephys_umap";
        inline const QString MorphoUmap = "morpho_umap";
    }

    // Format values
    namespace values
    {
        inline const QString CytosplorePatchSeqConfig = "cytosplore-patchseq-config";
    }
}

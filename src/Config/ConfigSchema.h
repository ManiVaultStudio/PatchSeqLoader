#pragma once

#include <QString>

namespace config::keys
{
    // Top-level
    inline const QString Format = "format";
    inline const QString Version = "version";
    inline const QString DatasetId = "dataset_id";

    inline const QString Sources = "sources";
    inline const QString Embeddings = "embeddings";

    // Common object fields
    inline const QString Path = "path";
    inline const QString Directory = "directory";
    inline const QString Index = "index";
    inline const QString ObsColumns = "obs_columns";
    inline const QString DisplayName = "display_name";

    // Directory-backed source fields
    inline const QString FailedSweepsPath = "failed_sweeps_path";
    inline const QString FilenameMetadataColumn = "filename_metadata_column";

    // Source keys
    namespace sources
    {
        inline const QString Rna = "rna";
        inline const QString Ephys = "ephys";
        inline const QString Morphology = "morphology";
        inline const QString Metadata = "metadata";
        inline const QString EphysTraces = "ephys_traces";
        inline const QString MorphologyReconstructions = "morphology_reconstructions";
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

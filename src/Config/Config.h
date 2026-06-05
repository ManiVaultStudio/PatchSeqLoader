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
        QStringList obsColumns;
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

    struct File
    {
    public:
        bool Load(QString filePath);

    public:
        QString format;
        QString version;
        QString datasetId;

        std::optional<config::TableSource> rna;
        std::optional<config::TableSource> ephys;
        std::optional<config::TableSource> morphology;
        std::optional<config::TableSource> metadata;

        std::optional<config::EphysTracesSource> ephysTraces;
        std::optional<config::MorphologyReconstructionsSource> morphologyReconstructions;

        std::optional<config::TableSource> rnaUmap;
        std::optional<config::TableSource> ephysUmap;
        std::optional<config::TableSource> morphoUmap;

        QMap<QString, config::TableSource> extraEmbeddings;
    };

} // namespace config

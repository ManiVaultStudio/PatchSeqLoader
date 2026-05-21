#pragma once

#include "Config/ConfigSchema.h"

#include <QString>
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

    //struct Embedding
    //{
    //    QString path;
    //    QString displayName;
    //    QString index;
    //    QString xColumn;
    //    QString yColumn;
    //};

    struct AssetDirectories
    {
        QString morphologyReconstruction;
        QString ephysTraces;
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

        config::AssetDirectories assetDirectories;

        std::optional<config::TableSource> rnaUmap;
        std::optional<config::TableSource> ephysUmap;
        std::optional<config::TableSource> morphoUmap;

        QMap<QString, config::TableSource> extraEmbeddings;
    };

} // namespace config

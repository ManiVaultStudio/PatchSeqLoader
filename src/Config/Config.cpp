#include "Config.h"

#include "ConfigSchema.h"

#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QJsonValue>
#include <QJsonArray>
#include <QSet>
#include <QDebug>

namespace
{
    QString ReadString(const QJsonObject& object, const QString& key, const QString& fallback = {})
    {
        const QJsonValue value = object.value(key);

        if (!value.isString())
            return fallback;

        return value.toString();
    }

    bool HasObject(const QJsonObject& object, const QString& key)
    {
        return object.contains(key) && object.value(key).isObject();
    }

    std::optional<QStringList> ReadStringArray(const QJsonObject& object, const QString& key)
    {
        const QJsonValue value = object.value(key);

        if (value.isUndefined())
            return std::nullopt;

        if (!value.isArray())
        {
            qWarning() << "Config field" << key << "must be an array of strings.";
            return std::nullopt;
        }

        QStringList strings;

        const QJsonArray array = value.toArray();

        for (const QJsonValue& item : array)
        {
            if (!item.isString())
            {
                qWarning() << "Config field" << key << "contains a non-string value.";
                return std::nullopt;
            }

            strings.push_back(item.toString());
        }

        return strings;
    }

    bool ReadObsColumns(const QJsonObject& object, const QString& sourceKey, config::TableSource& source)
    {
        const std::optional<QStringList> obsColumns = ReadStringArray(object, config::keys::ObsColumns);

        if (!obsColumns)
        {
            qWarning()
                << "Config source" << sourceKey
                << "has missing or invalid field:"
                << QString("%1.%2.%3")
                .arg(config::keys::Sources,
                    sourceKey,
                    config::keys::ObsColumns);

            return false;
        }

        const bool containsWildcard = obsColumns->contains("*");

        if (containsWildcard && obsColumns->size() > 1)
        {
            qWarning()
                << "Config source" << sourceKey
                << "has invalid obs_columns. The wildcard \"*\" cannot be mixed with explicit columns.";

            return false;
        }

        if (containsWildcard)
        {
            source.obsColumnMode = config::ObsColumnMode::All;
            source.obsColumns.clear();
        }
        else
        {
            source.obsColumnMode = config::ObsColumnMode::Explicit;
            source.obsColumns = *obsColumns;
        }

        return true;
    }

    std::optional<config::TableSource> ReadTableSource(const QJsonObject& sources, const QString& key)
    {
        if (!sources.contains(key))
        {
            qInfo() << "Config source" << key << "is not present; skipping it.";
            return std::nullopt;
        }

        if (!sources.value(key).isObject())
        {
            qWarning() << "Config source" << key << "exists but is not an object. Expected:" << QString("\"%1\": { \"%2\": \"...\", \"%3\": \"...\", \"%4\": \"...\", \"%5\": [...] }")
                .arg(key,
                    config::keys::Path,
                    config::keys::DisplayName,
                    config::keys::Index,
                    config::keys::ObsColumns);

            return std::nullopt;
        }

        const QJsonObject object = sources.value(key).toObject();

        config::TableSource source;
        source.path = ReadString(object, config::keys::Path);
        source.displayName = ReadString(object, config::keys::DisplayName);
        source.index = ReadString(object, config::keys::Index);

        if (source.path.isEmpty())
        {
            qWarning() << "Config source" << key << "has no path. Expected field:" << QString("%1.%2.%3")
                .arg(config::keys::Sources, key, config::keys::Path);

            return std::nullopt;
        }

        if (source.displayName.isEmpty())
        {
            qWarning() << "Config source" << key << "has no display name. Expected field:" << QString("%1.%2.%3")
                .arg(config::keys::Sources, key, config::keys::DisplayName);

            return std::nullopt;
        }

        if (source.index.isEmpty())
        {
            qWarning() << "Config source" << key << "has no index column. Expected field:" << QString("%1.%2.%3")
                .arg(config::keys::Sources, key, config::keys::Index);

            return std::nullopt;
        }

        if (!ReadObsColumns(object, key, source))
            return std::nullopt;

        if (source.obsColumnMode == config::ObsColumnMode::All)
            qInfo() << "Loaded config source" << key << "from" << source.path << "using index column" << source.index << "with all non-index columns treated as obs columns.";
        else
            qInfo() << "Loaded config source" << key << "from" << source.path << "using index column" << source.index << "with" << source.obsColumns.size() << "explicit obs columns.";

        return source;
    }

    //std::optional<config::TableSource> ReadEmbedding(const QJsonObject& embeddings, const QString& key)
    //{
    //    if (!HasObject(embeddings, key))
    //        return std::nullopt;

    //    const QJsonObject object = embeddings.value(key).toObject();

    //    config::TableSource embedding;
    //    embedding.path = ReadString(object, config::keys::Path);
    //    embedding.displayName = ReadString(object, config::keys::DisplayName);
    //    embedding.index = ReadString(object, config::keys::Index);

    //    if (embedding.path.isEmpty())
    //    {
    //        qWarning() << "Embedding" << key << "has no path; ignoring it. Expected field:" << QString("%1.%2.%3").arg(config::keys::Embeddings, key, config::keys::Path);
    //        return std::nullopt;
    //    }

    //    if (embedding.index.isEmpty())
    //    {
    //        qWarning() << "Embedding" << key << "has no index column; ignoring it. Expected field:" << QString("%1.%2.%3").arg(config::keys::Embeddings, key, config::keys::Index);
    //        return std::nullopt;
    //    }

    //    if (embedding.displayName.isEmpty())
    //        embedding.displayName = key;

    //    if (!ReadObsColumns(object, key, embedding))
    //        return std::nullopt;

    //    if (embedding.obsColumnMode == config::ObsColumnMode::All)
    //        qInfo() << "Loaded config source" << key << "from" << embedding.path << "using index column" << embedding.index << "with all non-index columns treated as obs columns.";
    //    else
    //        qInfo() << "Loaded config source" << key << "from" << embedding.path << "using index column" << embedding.index << "with" << source.obsColumns.size() << "explicit obs columns.";


    //    return embedding;
    //}

    config::AssetDirectories ReadAssetDirectories(const QJsonObject& root)
    {
        config::AssetDirectories directories;

        if (!HasObject(root, config::keys::Assets))
            return directories;

        const QJsonObject assets = root.value(config::keys::Assets).toObject();

        if (!HasObject(assets, config::keys::Directories))
            return directories;

        const QJsonObject dirs = assets.value(config::keys::Directories).toObject();

        directories.morphologyReconstruction = ReadString(dirs, config::keys::assets::MorphologyReconstruction);
        directories.ephysTraces = ReadString(dirs, config::keys::assets::EphysTraces);

        return directories;
    }

    void ValidateKnownFormat(const QString& format)
    {
        if (format != config::keys::values::CytosplorePatchSeqConfig)
        {
            qWarning()
                << "Unexpected config format:"
                << format
                << "Expected:" << config::keys::values::CytosplorePatchSeqConfig;
        }
    }

    void WarnUnknownSourceKeys(const QJsonObject& sources)
    {
        const QSet<QString> knownKeys = {
            config::keys::sources::Rna,
            config::keys::sources::Ephys,
            config::keys::sources::Morphology,
            config::keys::sources::Metadata
        };

        for (auto it = sources.begin(); it != sources.end(); ++it)
        {
            if (!knownKeys.contains(it.key()))
            {
                qWarning()
                    << "Unknown source key in config:"
                    << it.key()
                    << "Known source keys are:"
                    << QStringList(knownKeys.begin(), knownKeys.end()).join(", ");
            }
        }
    }

    QString ResolveRelativePath(const QString& configFilePath, const QString& path)
    {
        if (path.isEmpty())
            return path;

        QFileInfo info(path);

        if (info.isAbsolute())
            return path;

        const QFileInfo configInfo(configFilePath);
        const QDir configDir = configInfo.absoluteDir();

        return configDir.filePath(path);
    }

    void ResolveTableSourcePath(const QString& configFilePath, std::optional<config::TableSource>& source)
    {
        if (!source)
            return;

        source->path = ResolveRelativePath(configFilePath, source->path);
    }
}

namespace config
{
    bool File::Load(QString filePath)
    {
        // Reset existing state so reusing the same config::File object is safe.
        format.clear();
        version.clear();
        datasetId.clear();

        rna.reset();
        ephys.reset();
        morphology.reset();
        metadata.reset();

        assetDirectories = config::AssetDirectories{};

        rnaUmap.reset();
        ephysUmap.reset();
        morphoUmap.reset();

        extraEmbeddings.clear();

        QFile file(filePath);

        if (!file.exists())
        {
            qWarning() << "Config file does not exist:" << filePath;
            return false;
        }

        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            qWarning() << "Could not open config file:" << filePath;
            return false;
        }

        const QByteArray bytes = file.readAll();

        QJsonParseError parseError;
        const QJsonDocument document = QJsonDocument::fromJson(bytes, &parseError);

        if (parseError.error != QJsonParseError::NoError)
        {
            qWarning() << "Could not parse config JSON:" << parseError.errorString() << "at offset" << parseError.offset;
            return false;
        }

        if (!document.isObject())
        {
            qWarning() << "Config root must be a JSON object.";
            return false;
        }

        const QJsonObject root = document.object();

        format = ReadString(root, config::keys::Format);
        version = ReadString(root, config::keys::Version);
        datasetId = ReadString(root, config::keys::DatasetId);

        ValidateKnownFormat(format);

        if (format != config::keys::values::CytosplorePatchSeqConfig)
        {
            qWarning() << "Unsupported config format:" << format << "Expected:" << config::keys::values::CytosplorePatchSeqConfig;
            return false;
        }

        if (version.isEmpty())
        {
            qWarning() << "Config is missing version.";
            return false;
        }

        if (datasetId.isEmpty())
        {
            qWarning() << "Config is missing dataset_id.";
            return false;
        }

        if (HasObject(root, config::keys::Sources))
        {
            const QJsonObject sources = root.value(config::keys::Sources).toObject();

            WarnUnknownSourceKeys(sources);

            rna = ReadTableSource(sources, config::keys::sources::Rna);
            ephys = ReadTableSource(sources, config::keys::sources::Ephys);
            morphology = ReadTableSource(sources, config::keys::sources::Morphology);
            metadata = ReadTableSource(sources, config::keys::sources::Metadata);
        }
        else
        {
            qWarning() << "Config has no sources object.";
            return false;
        }

        if (!metadata)
        {
            qWarning() << "Config is missing required source: sources.metadata.";
            return false;
        }

        assetDirectories = ReadAssetDirectories(root);

        if (HasObject(root, config::keys::Embeddings))
        {
            const QJsonObject embeddings = root.value(config::keys::Embeddings).toObject();

            rnaUmap = ReadTableSource(embeddings, config::keys::embeddings::RnaUmap);
            ephysUmap = ReadTableSource(embeddings, config::keys::embeddings::EphysUmap);
            morphoUmap = ReadTableSource(embeddings, config::keys::embeddings::MorphoUmap);

            //rnaUmap = ReadEmbedding(embeddings, config::keys::embeddings::RnaUmap);
            //ephysUmap = ReadEmbedding(embeddings, config::keys::embeddings::EphysUmap);
            //morphoUmap = ReadEmbedding(embeddings, config::keys::embeddings::MorphoUmap);

            const QSet<QString> reservedEmbeddingNames = {
                config::keys::embeddings::RnaUmap,
                config::keys::embeddings::EphysUmap,
                config::keys::embeddings::MorphoUmap
            };

            for (auto it = embeddings.begin(); it != embeddings.end(); ++it)
            {
                const QString key = it.key();

                if (reservedEmbeddingNames.contains(key))
                    continue;

                if (!it.value().isObject())
                {
                    qWarning() << "Extra embedding" << key << "is not an object; ignoring it.";
                    continue;
                }

                QJsonObject temp;
                temp.insert(key, it.value());

                const std::optional<config::TableSource> embedding = ReadTableSource(temp, key);

                if (embedding)
                    extraEmbeddings.insert(key, *embedding);
            }
        }
        
        // Resolve relative paths against the config file location.
        ResolveTableSourcePath(filePath, rna);
        ResolveTableSourcePath(filePath, ephys);
        ResolveTableSourcePath(filePath, morphology);
        ResolveTableSourcePath(filePath, metadata);

        ResolveTableSourcePath(filePath, rnaUmap);
        ResolveTableSourcePath(filePath, ephysUmap);
        ResolveTableSourcePath(filePath, morphoUmap);

        for (auto it = extraEmbeddings.begin(); it != extraEmbeddings.end(); ++it)
            it->path = ResolveRelativePath(filePath, it->path);

        assetDirectories.morphologyReconstruction =
            ResolveRelativePath(filePath, assetDirectories.morphologyReconstruction);

        assetDirectories.ephysTraces =
            ResolveRelativePath(filePath, assetDirectories.ephysTraces);

        return true;
    }

} // namespace config

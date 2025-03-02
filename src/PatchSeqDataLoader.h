#pragma once

#include <LoaderPlugin.h>

#include "DataFrame.h"
#include "PatchSeqFilePaths.h"

#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>
#include <TextData/TextData.h>
#include <CellMorphologyData/CellMorphologyData.h>
#include <EphysData/EphysData.h>

#include <Task.h>

using namespace mv::plugin;

// =============================================================================
// Loading input box
// =============================================================================

class PatchSeqDataLoader;

// =============================================================================
// View
// =============================================================================

class PatchSeqDataLoader : public LoaderPlugin
{
    Q_OBJECT
public:
    PatchSeqDataLoader(const PluginFactory* factory) :
        LoaderPlugin(factory),
        _task(this, "Loading patch-seq data", mv::Task::Status::Idle, true)
    { }

    ~PatchSeqDataLoader(void) override;

    void init() override;

    void addTaxonomyClustersForDf(DataFrame& df, DataFrame& metadata, DataFrame& taxonomyDf, QString name, mv::Dataset<mv::DatasetImpl> parent);
    void createClusterData(std::vector<QString> stringList, QString dataName, mv::Dataset<mv::DatasetImpl> parent);

    void loadData() Q_DECL_OVERRIDE;

private:
    void loadGeneExpressionData(QString filePath, const DataFrame& metadata);
    void loadEphysData(QString filePath, const DataFrame& metadata);
    void loadMorphologyData(QString filePath, const DataFrame& metadata);
    void loadMorphologyCells(QDir dir);
private:
    DataFrame _taxonomyDf;

    // Metadata
    DataFrame _metadata;

    // Gene expressions
    DataFrame _transcriptomicsDf;
    Dataset<Points> _geneExpressionData;

    // Electrophysiology
    DataFrame _ephysDf;
    Dataset<Points> _ephysData;

    // Ephys traces
    Dataset<EphysExperiments> _ephysTraces;

    // Morphology
    DataFrame _morphologyDf;
    Dataset<Points> _morphoData;
    DataFrame _morphoMetadata;

    // Cell morphology
    Dataset<CellMorphologies> _cellMorphoData;

    mv::ModalTask _task;
};


// =============================================================================
// Factory
// =============================================================================

class PatchSeqDataLoaderFactory : public LoaderPluginFactory
{
    Q_INTERFACES(mv::plugin::LoaderPluginFactory mv::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "studio.manivault.PatchSeqDataLoader"
                      FILE  "PatchSeqDataLoader.json")

public:
    PatchSeqDataLoaderFactory(void) {}
    ~PatchSeqDataLoaderFactory(void) override {}

    /**
     * Get plugin icon
     * @param color Icon color for flat (font) icons
     * @return Icon
     */
    QIcon getIcon(const QColor& color = Qt::black) const override;

    LoaderPlugin* produce() override;

    mv::DataTypes supportedDataTypes() const override;
};

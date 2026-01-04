#pragma once

#include <LoaderPlugin.h>

#include "DataFrame.h"
#include "PatchSeqFilePaths.h"
#include "ColorTaxonomy.h"

#include <PointData/PointData.h>
#include <ClusterData/ClusterData.h>
#include <TextData/TextData.h>
#include <CellMorphologyData/CellMorphologyData.h>
#include <EphysData/EphysData.h>

#include <Task.h>
#include <SelectionGroup.h>

#include <QString>
#include <QColor>

using namespace mv::plugin;

// =============================================================================
// Loading input box
// =============================================================================

class PatchSeqDataLoader;

namespace mv
{
    class BiMap;
}

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

    void addTaxonomyClustersForDf(DataFrame& df, DataFrame& metadata, TaxonomyLevel level, QString name, mv::Dataset<mv::DatasetImpl> parent, QString metaLabel);
    void createClusterData(std::vector<QString> stringList, QString dataName, mv::Dataset<mv::DatasetImpl> parent);

    void loadData() Q_DECL_OVERRIDE;

private:
    void loadGeneExpressionData(QString filePath, const DataFrame& metadata);
    void loadEphysData(QString filePath, const DataFrame& metadata);
    void loadMorphologyData(QString filePath, const DataFrame& metadata);
    void loadMorphologyCells(QDir dir);
    void loadEphysTraces(QDir dir);
    void loadUMap(QString filePath, mv::Dataset<Points> parent, QString datasetName, BiMap& bimap);

private:
    DataFrame _taxonomyDf;
    QHash<QString, QColor> _cellTypeColors;
    ColorTaxonomy _colorTaxonomy;

    KeyBasedSelectionGroup _selectionGroup;

    // Metadata
    DataFrame _metadataDf;
    Dataset<Text> _metadata;

    // Gene expressions
    DataFrame _transcriptomicsDf;
    Dataset<Points> _geneExpressionData;
    DataFrame _gexprMetadata;

    // Electrophysiology
    DataFrame _ephysDf;
    Dataset<Points> _ephysData;
    DataFrame _ephysMetadata;

    // Ephys traces
    Dataset<EphysExperiments> _ephysTraces;
    std::vector<QString> _ephysTraceCellIds;

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
    PatchSeqDataLoaderFactory(void)
    {
        setIconByName("database");
    }
    ~PatchSeqDataLoaderFactory(void) override {}

    LoaderPlugin* produce() override;

    mv::DataTypes supportedDataTypes() const override;
};

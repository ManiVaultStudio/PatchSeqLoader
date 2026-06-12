#include "PatchSeqDataLoader.h"

#include "InputDialog.h"

#include "LoadingPipeline/Stages/LoadConfigStage.h"
#include "LoadingPipeline/Stages/DiscoverSourcesStage.h"
#include "LoadingPipeline/Stages/LoadTablesStage.h"
#include "LoadingPipeline/Stages/NormalizeTablesStage.h"
#include "LoadingPipeline/Stages/CollectMetadataStage.h"
#include "LoadingPipeline/Stages/LoadTaxonomiesStage.h"
#include "LoadingPipeline/Stages/LoadAssetsStage.h"
#include "LoadingPipeline/Stages/CreateDatasetsStage.h"
#include "LoadingPipeline/Stages/LinkDatasetsStage.h"

#include <Set.h>

Q_PLUGIN_METADATA(IID "studio.manivault.PatchSeqDataLoader")

using namespace mv;
using namespace mv::gui;

namespace
{
    void ShowPipelineErrors(const PipelineResult& result)
    {
        for (const PipelineIssue& issue : result.issues)
        {
            const QString prefix =
                issue.severity == PipelineIssue::Severity::Error ? "ERROR" :
                issue.severity == PipelineIssue::Severity::Warning ? "WARNING" : "INFO";

            qWarning().noquote() << QString("[%1] %2: %3 %4")
                .arg(prefix)
                .arg(issue.stage)
                .arg(issue.message)
                .arg(issue.detail);
        }
    }
}

// =============================================================================
// View
// =============================================================================

PatchSeqDataLoader::~PatchSeqDataLoader(void)
{

}

void PatchSeqDataLoader::init()
{

}

void PatchSeqDataLoader::loadData()
{
    Q_INIT_RESOURCE(met_loader_resources);

    auto configPath = QString("D:/Dropbox/Julian/Patchseq/Projects/config_basal_ganglia.json");//askUserForConfigPath();
    if (configPath.isEmpty())
        return;

    Pipeline pipeline;
    pipeline.Add(std::make_unique<LoadConfigStage>());
    pipeline.Add(std::make_unique<DiscoverSourcesStage>());
    pipeline.Add(std::make_unique<LoadTablesStage>());
    pipeline.Add(std::make_unique<NormalizeTablesStage>());
    pipeline.Add(std::make_unique<CollectMetadataStage>());
    pipeline.Add(std::make_unique<LoadTaxonomiesStage>());
    pipeline.Add(std::make_unique<LoadAssetsStage>());
    pipeline.Add(std::make_unique<CreateDatasetsStage>());
    pipeline.Add(std::make_unique<LinkDatasetsStage>());

    PipelineContext ctx;
    ctx.configPath = configPath;
    ctx.task = &_task;
    ctx.task->setRunning();

    const auto result = pipeline.Run(ctx);
    qDebug() << "End of pipeline";

    //if (!result.Ok())
    ShowPipelineErrors(result);
}

// =============================================================================
// Factory
// =============================================================================

LoaderPlugin* PatchSeqDataLoaderFactory::produce()
{
    return new PatchSeqDataLoader(this);
}

DataTypes PatchSeqDataLoaderFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;
    supportedTypes.append(PointType);
    return supportedTypes;
}

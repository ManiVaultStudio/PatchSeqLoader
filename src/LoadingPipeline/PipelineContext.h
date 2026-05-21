// PipelineContext.hpp
#pragma once

#include "Config/Config.h"
#include "AnnotatedData.h"

//#include "PatchSeqConfig.hpp"
//#include "PatchSeqCellRecord.hpp"
//#include "PatchSeqLoadedData.hpp"

#include <PointData/PointData.h>
#include <TextData/TextData.h>
#include <SelectionGroup.h>

#include <QString>
#include <QMap>
#include <QVector>

struct PipelineIssue
{
    enum class Severity
    {
        Info,
        Warning,
        Error
    };

    Severity severity = Severity::Info;
    QString stage;
    QString message;
    QString detail;
};

struct PipelineResult
{
    QVector<PipelineIssue> issues;

    bool Ok() const
    {
        for (const auto& issue : issues) {
            if (issue.severity == PipelineIssue::Severity::Error)
                return false;
        }
        return true;
    }

    void Error(QString stage, QString message, QString detail = {})
    {
        issues.push_back({ PipelineIssue::Severity::Error, stage, message, detail });
    }

    void Warning(QString stage, QString message, QString detail = {})
    {
        issues.push_back({ PipelineIssue::Severity::Warning, stage, message, detail });
    }

    void Info(QString stage, QString message, QString detail = {})
    {
        issues.push_back({ PipelineIssue::Severity::Info, stage, message, detail });
    }
};

struct PipelineContext
{
    QString configPath;
    QString datasetRoot;

    config::File config;
    //PatchSeqConfig config;

    // Raw input tables, keyed by source name:
    // "metadata", "ephys", "morphology", "rna", etc.
    QMap<QString, AnnotatedData> rawTables;

    // Cleaned, normalized tables.
    QMap<QString, AnnotatedData> normalizedTables;

    QMap<QString, mv::Dataset<Points>> featureDatasets;
    QMap<QString, mv::Dataset<Text>> textDatasets;
    QMap<QString, mv::Dataset<Points>> embeddingDatasets;

    KeyBasedSelectionGroup selectionGroup;

    //// Optional loaded assets.
    //QVector<CellMorphology> morphologyCells;
    //QVector<EphysExperiment> ephysExperiments;

    //// Manivault/viewer datasets created near the end.
    //PatchSeqLoadedData loaded;

    PipelineResult result;

    // Optional task/progress bridge.
    mv::Task* task = nullptr;
};

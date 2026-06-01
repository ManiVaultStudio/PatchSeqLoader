#pragma once

#include "AnnotatedData.h"

#include <QString>

class AnnotatedDataExporter
{
public:
    static void Export(const AnnotatedData& data, const QString& filePath);

private:
    static void ValidateObs(const AnnotationTable& obs);
    static QString EscapeCsvField(const QString& value);
};

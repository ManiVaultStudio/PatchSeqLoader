#pragma once

#include <QString>

class QDir;

class PatchSeqFilePaths
{
public:
    void locateFilePaths(QDir dir);

    bool allFilesLocated();

    bool hasGeneExpressions() { return !gexprFilePath.isNull() && !gexprFilePath.isEmpty(); }
    bool hasEphysFeatures() { return !ephysFilePath.isNull() && !ephysFilePath.isEmpty(); }
    bool hasMorphologyFeatures() { return !morphoFilePath.isNull() && !morphoFilePath.isEmpty(); }

    bool hasEphysUMap() { return !ephysUMapFilePath.isNull() && !ephysUMapFilePath.isEmpty(); }
    bool hasMorphoUMap() { return !morphoUMapFilePath.isNull() && !morphoUMapFilePath.isEmpty(); }
    bool hasTxUMap() { return !txUMapFilePath.isNull() && !txUMapFilePath.isEmpty(); }
    bool hasMEUMap() { return !meUMapFilePath.isNull() && !meUMapFilePath.isEmpty(); }

    bool hasEphysTraces() { return !ephysTracesDir.isNull() && !ephysTracesDir.isEmpty(); }
    bool hasMorphologies() { return !morphologiesDir.isNull() && !morphologiesDir.isEmpty(); }

public:
    QString gexprFilePath;
    QString ephysFilePath;
    QString morphoFilePath;
    QString metadataFilePath;
    QString annotationFilePath;

    QString ephysUMapFilePath;
    QString morphoUMapFilePath;
    QString txUMapFilePath;
    QString meUMapFilePath;

    QString ephysTracesDir;
    QString morphologiesDir;
};

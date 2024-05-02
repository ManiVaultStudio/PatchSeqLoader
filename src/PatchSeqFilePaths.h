#pragma once

#include <QString>

class QDir;

class PatchSeqFilePaths
{
public:
    void locateFilePaths(QDir dir);

    bool allFilesLocated();

public:
    QString gexprFilePath;
    QString ephysFilePath;
    QString morphoFilePath;
    QString metadataFilePath;
    QString annotationFilePath;
};

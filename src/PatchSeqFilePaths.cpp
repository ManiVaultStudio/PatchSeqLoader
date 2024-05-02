#include "PatchSeqFilePaths.h"

#include <QDir>
#include <QStringList>
#include <QDebug>

void PatchSeqFilePaths::locateFilePaths(QDir dir)
{
    QStringList csvFiles = dir.entryList(QStringList() << "*.csv" << "*.CSV", QDir::Files);

    for (QString filePath : csvFiles)
    {
        if (filePath.contains("IDs"))
            gexprFilePath = dir.filePath(filePath);
        if (filePath.contains("ephys"))
            ephysFilePath = dir.filePath(filePath);
        if (filePath.contains("morpho"))
            morphoFilePath = dir.filePath(filePath);
        if (filePath.contains("metadata"))
            metadataFilePath = dir.filePath(filePath);
        if (filePath.contains("annotation"))
            annotationFilePath = dir.filePath(filePath);
    }

    qDebug() << "Located gene expression file: " << gexprFilePath;
    qDebug() << "Located electrophysiology file: " << ephysFilePath;
    qDebug() << "Located morphology data file: " << morphoFilePath;
    qDebug() << "Located metadata file: " << metadataFilePath;
    qDebug() << "Located annotation file: " << annotationFilePath;

    // Don't try to load a file if the dialog was cancelled or the file name is empty
    if (metadataFilePath.isNull() || metadataFilePath.isEmpty())
        return;
}

bool PatchSeqFilePaths::allFilesLocated()
{
    if (gexprFilePath.isNull() || gexprFilePath.isEmpty())
        return false;
    if (ephysFilePath.isNull() || ephysFilePath.isEmpty())
        return false;
    if (morphoFilePath.isNull() || morphoFilePath.isEmpty())
        return false;
    if (metadataFilePath.isNull() || metadataFilePath.isEmpty())
        return false;
    if (annotationFilePath.isNull() || annotationFilePath.isEmpty())
        return false;
    return true;
}

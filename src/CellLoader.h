#pragma once

#include <QString>

#include <string>

class CellMorphology;

void loadCellContentsFromFile(QString filePath, std::string& result);

void readCell(const std::string& contents, CellMorphology& cellMorphology);

#include "SWCLoader.h"

#include "CellMorphologyData/CellMorphology.h"

#include <graphics/Vector3f.h>

#include <fstream>
#include <sstream>
#include <iostream>

#include <vector>
#include <limits>

namespace
{
    void loadCellContentsFromFile(QString filePath, std::string& result)
    {
        std::ifstream file(filePath.toStdString());

        if (file) {
            std::ostringstream ss;
            ss << file.rdbuf(); // reading data
            result = ss.str();
        }
    }
}

void SWCLoader::LoadSWC(QString filePath, CellMorphology& cellMorphology)
{
    std::string fileContents;
    loadCellContentsFromFile(filePath, fileContents);

    std::istringstream fileStream(fileContents);
    std::string line;

    // #n type x y z radius parent
    while (getline(fileStream, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;
        std::string value;

        if (line.at(0) == '#')
            continue;

        int i = 0;
        while (getline(ss, value, ' ')) {
            //row.push_back(value);
            switch (i)
            {
            case 0: cellMorphology.ids.push_back(stoi(value)); break;
            case 1: cellMorphology.types.push_back(stoi(value)); break;
            case 2:
            {
                float x = stof(value);
                getline(ss, value, ' ');
                float y = stof(value);
                getline(ss, value, ' ');
                float z = stof(value);
                cellMorphology.positions.emplace_back(x, y, z);
                break;
            }
            case 3: cellMorphology.radii.push_back(stof(value)); break;
            case 4: cellMorphology.parents.push_back(stoi(value)); break;
            }
            i++;
        }
    }

    // Compute id to index map
    for (int i = 0; i < cellMorphology.ids.size(); i++)
    {
        int id = cellMorphology.ids[i];
        cellMorphology.idMap[id] = i;
    }

    // Find soma
    for (int i = 0; i < cellMorphology.ids.size(); i++)
    {
        mv::Vector3f position = cellMorphology.positions.at(i);

        if (cellMorphology.types.at(i) == 1) // Soma
        {
            cellMorphology.somaPosition = position;
            break;
        }
    }
}

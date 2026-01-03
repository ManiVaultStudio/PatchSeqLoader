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

    // Header format: #n type x y z radius parent

    // Read file in line-by-line and tokenize it
    while (getline(fileStream, line))
    {
        std::stringstream ss(line);

        std::string token;

        // Skip comments
        if (line.at(0) == '#')
            continue;

        // Tokenize based on space delimiter
        int i = 0;
        while (getline(ss, token, ' '))
        {
            switch (i)
            {
            case 0: cellMorphology.ids.push_back(stoi(token)); break;
            case 1: cellMorphology.types.push_back(stoi(token)); break;
            case 2:
            {
                float x = stof(token);
                getline(ss, token, ' ');
                float y = stof(token);
                getline(ss, token, ' ');
                float z = stof(token);
                cellMorphology.positions.emplace_back(x, y, z);
                break;
            }
            case 3: cellMorphology.radii.push_back(stof(token)); break;
            case 4: cellMorphology.parents.push_back(stoi(token)); break;
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

        if (cellMorphology.types.at(i) == (int) CellMorphology::Type::Soma)
        {
            cellMorphology.somaPosition = position;
            break;
        }
    }
}

#include "SWCLoader.h"

#include "CellMorphologyData/CellMorphology.h"

#include <graphics/Vector3f.h>

#include <fstream>
#include <sstream>
#include <iostream>

#include <vector>
#include <limits>

void loadCellContentsFromFile(QString filePath, std::string& result)
{
    std::ifstream file(filePath.toStdString());

    if (file) {
        std::ostringstream ss;
        ss << file.rdbuf(); // reading data
        result = ss.str();
    }
}

void readCell(const std::string& contents, CellMorphology& cellMorphology)
{
    std::istringstream fileStream(contents);
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

    // Print headers
    //for (const auto& header : headers) {
    //    std::cout << header << "\t";
    //}

    // Print data
    //for (const auto& position : _positions) {
    //    std::cout << position.str() << std::endl;
    //}

    //// Find centroid and extents
    //mv::Vector3f avgPos;
    //for (const auto& pos : cellMorphology.positions)
    //    avgPos += pos;
    //avgPos /= cellMorphology.positions.size();

    //// Center cell positions
    //for (auto& pos : cellMorphology.positions)
    //    pos -= avgPos;

    //// Find cell position ranges
    //mv::Vector3f minV(std::numeric_limits<float>::max());
    //mv::Vector3f maxV(-std::numeric_limits<float>::max());
    //for (const auto& pos : cellMorphology.positions)
    //{
    //    if (pos.x < minV.x) minV.x = pos.x;
    //    if (pos.y < minV.y) minV.y = pos.y;
    //    if (pos.z < minV.z) minV.z = pos.z;
    //    if (pos.x > maxV.x) maxV.x = pos.x;
    //    if (pos.y > maxV.y) maxV.y = pos.y;
    //    if (pos.z > maxV.z) maxV.z = pos.z;
    //}
    //mv::Vector3f range = (maxV - minV);
    //float maxRange = std::max(std::max(range.x, range.y), range.z);
    //// Rescale positions
    //for (auto& pos : cellMorphology.positions)
    //{
    //    pos /= maxRange;
    //}

    //std::cout << maxRange << std::endl;
    //std::cout << minV << std::endl;
    //std::cout << maxV << std::endl;
}

#include "ColorTaxonomy.h"

// Group,color_hex_group,Subclass,color_hex_subclass,Class,color_hex_class,Neighborhood,color_hex_neighborhood
#define GROUP_LABEL "Group"
#define SUBCLASS_LABEL "Subclass"
#define CLASS_LABEL "Class"
#define NEIGHBORHOOD_LABEL "Neighborhood"

#define GROUP_COLORS "color_hex_group"
#define SUBCLASS_COLORS "color_hex_subclass"
#define CLASS_COLORS "color_hex_class"
#define NEIGHBORHOOD_COLORS "color_hex_neighborhood"

void ColorTaxonomy::processFromDataFrame(const DataFrame& dataFrame)
{
    std::vector<QString> groupLabels = dataFrame[GROUP_LABEL];
    std::vector<QString> subclassLabels = dataFrame[SUBCLASS_LABEL];
    std::vector<QString> classLabels = dataFrame[CLASS_LABEL];
    std::vector<QString> neighborhoodLabels = dataFrame[NEIGHBORHOOD_LABEL];

    std::vector<QString> groupColors = dataFrame[GROUP_COLORS];
    std::vector<QString> subclassColors = dataFrame[SUBCLASS_COLORS];
    std::vector<QString> classColors = dataFrame[CLASS_COLORS];
    std::vector<QString> neighborhoodColors = dataFrame[NEIGHBORHOOD_COLORS];

    for (int i = 0; i < groupLabels.size(); i++)
        _groupColors[groupLabels[i]] = groupColors[i];
    
    for (int i = 0; i < subclassLabels.size(); i++)
        _subclassColors[subclassLabels[i]] = subclassColors[i];

    for (int i = 0; i < classLabels.size(); i++)
        _classColors[classLabels[i]] = classColors[i];

    for (int i = 0; i < neighborhoodLabels.size(); i++)
        _neighborhoodColors[neighborhoodLabels[i]] = neighborhoodColors[i];
}

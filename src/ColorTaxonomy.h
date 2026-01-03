#pragma once

#include "DataFrame.h"

#include <QHash>
#include <QString>
#include <QColor>

enum class TaxonomyLevel
{
    GROUP, SUBCLASS, CLASS, NEIGHBORHOOD
};

class ColorTaxonomy
{
public:
    const QColor& getGroupColor(const QString& label) { return _groupColors[label]; }
    const QColor& getSubclassColor(const QString& label) { return _subclassColors[label]; }
    const QColor& getClassColor(const QString& label) { return _classColors[label]; }
    const QColor& getNeighborhoodColor(const QString& label) { return _neighborhoodColors[label]; }

public:
    void processFromDataFrame(const DataFrame& dataFrame);

private:
    QHash<QString, QColor> _groupColors;
    QHash<QString, QColor> _subclassColors;
    QHash<QString, QColor> _classColors;
    QHash<QString, QColor> _neighborhoodColors;
};

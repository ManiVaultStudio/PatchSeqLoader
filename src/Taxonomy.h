#pragma once

struct Taxonomy
{
    QString key;
    QString displayName;

    // Non-color taxonomy columns, e.g. Group, Subclass, Class, Neighborhood.
    QStringList levels;

    // level -> label -> color
    QHash<QString, QHash<QString, QColor>> colorMaps;

    // leaf level -> leaf label -> level -> label
    //
    // In practice, use the first non-color column as the leaf/key column.
    // For your file: Group -> "Astrocyte" -> { Group, Subclass, Class, Neighborhood }
    QHash<QString, QHash<QString, QString>> labelsByLeaf;
};

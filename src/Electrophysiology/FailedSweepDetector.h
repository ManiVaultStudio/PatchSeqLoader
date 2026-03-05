#pragma once

#include <QHash>
#include <QString>
#include <QVector>

QHash<QString, QVector<int>> LoadFailedSweeps(const QString& jsonPath);

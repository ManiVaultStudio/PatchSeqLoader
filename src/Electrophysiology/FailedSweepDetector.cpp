#include <QFile>
#include <QHash>
#include <QString>
#include <QVector>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QDebug>

QHash<QString, QVector<int>> LoadFailedSweeps(const QString& jsonPath)
{
    QHash<QString, QVector<int>> result;

    QFile file(jsonPath);
    if (!file.open(QIODevice::ReadOnly))
    {
        qWarning() << "Failed to open file:" << jsonPath;
        return result;
    }

    QByteArray data = file.readAll();
    file.close();

    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(data, &parseError);

    if (parseError.error != QJsonParseError::NoError)
    {
        qWarning() << "JSON parse error:" << parseError.errorString();
        return result;
    }

    if (!doc.isObject())
    {
        qWarning() << "JSON root is not an object.";
        return result;
    }

    QJsonObject root = doc.object();

    for (auto it = root.begin(); it != root.end(); ++it)
    {
        QString fileName = it.key();
        QJsonObject entry = it.value().toObject();

        QVector<int> sweeps;

        if (entry.contains("failed_sweeps"))
        {
            QJsonArray arr = entry["failed_sweeps"].toArray();

            for (const QJsonValue& v : arr)
                sweeps.append(v.toInt());
        }

        result.insert(fileName, sweeps);
    }

    return result;
}

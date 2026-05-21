#include "CSVLoader.h"

#include "csv.h"

#include <QFileInfo>

#include <stdexcept>

// Magic number that represents a missing value, to be imputed
constexpr float MISSING_VALUE = 1234567.0f;
constexpr char DELIMITER = ',';

namespace
{
    namespace
    {
        bool IsMetadataColumn(const QString& token, const QStringList& metadataHeaders, bool allColumnsAreMetadata)
        {
            if (allColumnsAreMetadata)
                return true;

            return metadataHeaders.contains(token);
        }
    }

    // Split header-line of CSV into stringlist of headers, trims headers and handles quoted headers
    QStringList SplitHeaders(const QString& line)
    {
        QStringList headers;
        QString token;

        bool inQuotes = false;

        for (int i = 0; i < line.size(); ++i)
        {
            const QChar ch = line[i];

            if (ch == '"')
            {
                // Escaped quote: "" inside a quoted field means one literal "
                if (inQuotes && i + 1 < line.size() && line[i + 1] == '"')
                {
                    token += '"';
                    ++i;
                }
                else
                    inQuotes = !inQuotes;
            }
            else if (ch == ',' && !inQuotes)
            {
                headers.push_back(token.trimmed());
                token.clear();
            }
            else
                token += ch;
        }

        if (inQuotes)
            throw std::runtime_error("Malformed CSV line: quoted field was not closed.");

        headers.push_back(token.trimmed());

        return headers;
    }

    void LineRead(char* line, QString& indexValue, std::vector<QString>& metadataRow, std::vector<float>& dataRow, const CsvLoadContext& ctx)
    {
        int colIndex = 0;

        char* p = line;
        while (true)
        {
            char* p2 = strchr(p, DELIMITER);
            if (p2 != NULL)
                *p2 = '\0';

            if (ctx.columnKinds[colIndex] == CsvColumnKind::Data)
            {
                if (*p == '\0')
                    dataRow.push_back(MISSING_VALUE);
                else
                    dataRow.push_back(atof(p));
            }
            else if (ctx.columnKinds[colIndex] == CsvColumnKind::Metadata)
            {
                metadataRow.push_back(p);
                metadataRow[metadataRow.size() - 1].replace("\"", "");
            }
            else if (ctx.columnKinds[colIndex] == CsvColumnKind::Index)
            {
                indexValue = p;
            }

            colIndex++;

            if (p2 == NULL)
                break;
            p = p2 + 1;
        }
    }

    // Read the first line of the CSV file and split it into metadata headers and data headers
    void ReadHeader(QString filePath, const QStringList& metadataHeaders, const QString& index, CsvLoadContext& ctx, AnnotatedData& data)
    {
        QFile inputFile(filePath);

        // Check if file can open, if not throw an exception
        if (!inputFile.open(QIODevice::ReadOnly))
            throw std::runtime_error(QString("Failed to open file at location: %1").arg(filePath).toStdString());

        QTextStream in(&inputFile);

        if (in.atEnd())
            throw std::runtime_error(QString("CSV file is empty: %1").arg(filePath).toStdString());

        const QString headerLine = in.readLine();
        const QStringList tokens = SplitHeaders(headerLine);

        data.obs.columnNames.clear();

        for (const QString& token : tokens)
        {
            if (token == index)
            {
                data.obs.indexName = token;
                ctx.columnKinds.push_back(CsvColumnKind::Index);
            }
            else if (IsMetadataColumn(token, metadataHeaders, ctx.allColumnsAreMetadata))
            {
                data.obs.columnNames.push_back(token);
                ctx.columnKinds.push_back(CsvColumnKind::Metadata);
            }
            else
            {
                data.var.index.push_back(token);
                ctx.columnKinds.push_back(CsvColumnKind::Data);
            }
        }
    }

    void ReadBody(QString filePath, CsvLoadContext& ctx, AnnotatedData& data)
    {
        std::vector<QString> metadataRow;
        std::vector<float> dataRow;

        // Open file again
        io::LineReader fin(filePath.toStdString());

        // Define variables
        char* token;
        int lineCount = 0;

        // Skip header
        fin.next_line();

        data.obs.values.resize(data.obs.columnNames.size());
        // Process data line-by-line
        while (char* line = fin.next_line())
        {
            metadataRow.clear();
            dataRow.clear();

            QString indexValue;
            metadataRow.reserve(data.obs.columnNames.size());
            dataRow.reserve(data.var.index.size());

            LineRead(line, indexValue, metadataRow, dataRow, ctx);

            data.obs.index.push_back(indexValue);
            for (int i = 0; i < metadataRow.size(); i++)
                data.obs.values[i].push_back(metadataRow[i]); // Column-major storage
            data.X.values.insert(data.X.values.end(), dataRow.begin(), dataRow.end());

            lineCount++;
        }
        data.X.rowCount = lineCount;
        data.X.columnCount = data.var.index.size();
    }
}

void CsvLoader::Load(const QString& filePath, const QStringList& metadataHeaders, const QString& index, AnnotatedData& data)
{
    // Check if the file exists
    QFileInfo fileInfo(filePath);
    if (!fileInfo.exists())
        throw std::runtime_error(QString("File was not found at location: %1").arg(filePath).toStdString());

    CsvLoadContext ctx;
    ReadHeader(filePath, metadataHeaders, index, ctx, data);
    ReadBody(filePath, ctx, data);
}

void CsvLoader::Load(const config::TableSource& config, AnnotatedData& data)
{
    // Check if the file exists
    QFileInfo fileInfo(config.path);
    if (!fileInfo.exists())
        throw std::runtime_error(QString("File was not found at location: %1").arg(config.path).toStdString());

    CsvLoadContext ctx;
    ctx.allColumnsAreMetadata = config.obsColumnMode == config::ObsColumnMode::All;
    ReadHeader(config.path, config.obsColumns, config.index, ctx, data);
    ReadBody(config.path, ctx, data);
}

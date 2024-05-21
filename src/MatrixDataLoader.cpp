#include "MatrixDataLoader.h"

#include <LoaderPlugin.h>
#include <util/Timer.h>
#include <Task.h>

#include <QFile>
#include <QTextStream>
#include <QStringList>

#include "csv.h"

constexpr char DELIMITER = ',';

namespace
{
    void ReadHeader(QString fileName, DataFrame& df, MatrixData& matrix, int numMetaColumns)
    {
        QFile inputFile(fileName);

        // Check if file can open, if not throw an exception
        if (!inputFile.open(QIODevice::ReadOnly))
        {
            throw mv::plugin::DataLoadException(fileName, "Failed to open file at location.");
        }

        QTextStream in(&inputFile);

        if (!in.atEnd())
        {
            // Read a single line from the file (hopefully the header)
            QString line = in.readLine();

            // Split the header line up into tokens
            QStringList tokens = line.split(",");

            // Determine the number of metadata columns and matrix columns
            matrix.numCols = (tokens.size() - numMetaColumns);

            // Get rid of any extra double quotes
            for (int i = 0; i < tokens.size(); i++)
            {
                QString& token = tokens[i];
                token.replace("\"", "");

                // Set the metadata columns in the dataframe or matrix
                if (i < numMetaColumns)
                    df.addHeader(token);
                else
                    matrix.headers.push_back(token);
            }
        }

        inputFile.close();
    }

    void FastLineRead(char* line, std::vector<QString>& metadataRow, std::vector<float>& dataRow, int numMetaColumns)
    {
        char* token;
        int colIndex = 0;

        // Read metadata part
        for (int i = 0; i < numMetaColumns; i++)
        {
            token = strtok(i == 0 ? line : NULL, &DELIMITER);
            metadataRow[i] = token;
            metadataRow[i].replace("\"", "");
        }

        // Read the data part
        token = strtok(NULL, &DELIMITER);

        while (token != NULL)
        {
            dataRow[colIndex++] = atof(token);

            token = strtok(NULL, &DELIMITER);
        }
    }

    void MissingValueLineRead(char* line, std::vector<QString>& metadataRow, std::vector<float>& dataRow, int numMetaColumns)
    {
        int colIndex = 0;

        char* p = line;
        while (true)
        {
            char* p2 = strchr(p, DELIMITER);
            if (p2 != NULL)
                *p2 = '\0';

            if (colIndex < numMetaColumns)
            {
                metadataRow[colIndex] = p;
                metadataRow[colIndex].replace("\"", "");
            }
            else
            {
                if (*p == '\0')
                    dataRow[colIndex - numMetaColumns] = MISSING_VALUE;
                else
                    dataRow[colIndex - numMetaColumns] = atof(p);
            }
            colIndex++;

            if (p2 == NULL)
                break;
            p = p2 + 1;
        }
    }

    void ReadBody(QString fileName, DataFrame& df, MatrixData& matrix, int numMetaColumns, bool handleMissingValues)
    {
        std::vector<QString> metadataRow(numMetaColumns);
        std::vector<float> dataRow(matrix.numCols);

        // Open file again
        io::LineReader fin(fileName.toStdString());

        // Define variables
        char* token;
        int lineCount = 0;

        // Skip header
        fin.next_line();
        
        // Process data line-by-line
        while (char* line = fin.next_line())
        {
            if (handleMissingValues)
                MissingValueLineRead(line, metadataRow, dataRow, numMetaColumns);
            else
                FastLineRead(line, metadataRow, dataRow, numMetaColumns);

            df.getData().push_back(metadataRow);
            matrix.data.insert(matrix.data.end(), dataRow.begin(), dataRow.end());

            lineCount++;
        }
        matrix.numRows = lineCount;
    }
}

void MatrixDataLoader::LoadMatrixData(QString fileName, DataFrame& df, MatrixData& matrix, int numMetaCols)
{
    // Check if the file exists
    QFileInfo fileInfo(fileName);
    if (!fileInfo.exists())
    {
        throw mv::plugin::DataLoadException(fileName, "File was not found at location.");
    }

    // Measure time
    Timer timer("Data Load [" + fileName + "]");

    // File is open
    ReadHeader(fileName, df, matrix, numMetaCols);
    ReadBody(fileName, df, matrix, numMetaCols, _handleMissingValues);
}

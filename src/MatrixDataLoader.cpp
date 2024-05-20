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

    void ReadBody(QString fileName, DataFrame& df, MatrixData& matrix, int numMetaColumns)
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
            int colIndex = 0;

            // Read metadata part
            for (int i = 0; i < numMetaColumns; i++)
            {
                token = strtok(i == 0 ? line : NULL, &DELIMITER);
                metadataRow[i] = token;
                metadataRow[i].replace("\"", "");
            }
            df.getData().push_back(metadataRow);

            // Read the data part
            token = strtok(NULL, &DELIMITER);

            while (token != NULL)
            {
                dataRow[colIndex++] = atof(token);

                token = strtok(NULL, &DELIMITER);
            }
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
    ReadBody(fileName, df, matrix, numMetaCols);
}

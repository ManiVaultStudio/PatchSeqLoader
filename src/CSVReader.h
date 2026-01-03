#pragma once

#include <QString>

#include <sstream>

class CSVReader
{
public:
    void LoadCSV(QString filePath, std::vector<QString>& headers, std::vector<std::vector<QString>>& data);
    void LoadCSV(std::stringstream& sstream, std::vector<QString>& headers, std::vector<std::vector<QString>>& data);
private:

};

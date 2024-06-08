#include "NWBLoader.h"

#include "Electrophysiology/Recording.h"
#include "Electrophysiology/Experiment.h"

#include "H5Cpp.h"

#include <QDebug>
#include <iostream>
#include <string>
#include <fstream>

using namespace H5;

namespace
{
    void attr_op(H5::H5Location& loc, const std::string attr_name,
        void* operator_data) {
        std::cout << attr_name << std::endl;
    }

    herr_t op_func(hid_t obj, const char* name, const H5O_info2_t* info, void* op_data)
    {
        if (info->type == H5O_TYPE_GROUP)
        {
            auto vec = static_cast<std::vector<std::string>*>(op_data);
            vec->push_back(std::string(name));
        }
        if (info->type == H5O_TYPE_DATASET)
        {

        }

        return 0;
    }

    void ParseComments(const std::string& comments, Recording& recording)
    {
        QString str = QString::fromStdString(comments);

        QStringList tokens = str.split('\n');

        for (int i = 0; i < tokens.size(); i++)
        {
            if (!tokens[i].trimmed().isEmpty())
                recording.comments.push_back(tokens[i].toStdString());
        }
    }

    void ReadComments(const H5File& file, std::string groupName, Recording& recording)
    {
        Group group = file.openGroup(groupName);
        Attribute attribute = group.openAttribute("comments");

        std::cout << "Attribute obj name: " << attribute.getName() << std::endl;

        StrType type = attribute.getStrType();

        // Read the data
        std::string comments;
        attribute.read(type, comments);

        std::cout << "Comments: " << comments << std::endl;
        ParseComments(comments, recording);

        attribute.close();
    }

    void ReadTimeseries(const H5File& file, std::string groupName, Recording& recording)
    {
        ReadComments(file, groupName, recording);

        DataSet dataset = file.openDataSet(groupName + "/data");

        std::cout << "Dataset obj name: " << dataset.getObjName() << std::endl;

        DataSpace dataSpace = dataset.getSpace();
        int ndims = dataSpace.getSimpleExtentNdims();
        std::vector<hsize_t> dims(ndims);
        dataSpace.getSimpleExtentDims(dims.data());

        int totalDataSize = 1;
        for (hsize_t dimSize : dims)
            totalDataSize *= dimSize;

        FloatType type = dataset.getFloatType();

        recording.data.resize(totalDataSize);
        // Read the data
        dataset.read(recording.data.data(), type);

        dataset.close();
    }
}

void NWBLoader::LoadNWB(QString fileName)
{
    Experiment experiment;

    H5File file(fileName.toStdString(), H5F_ACC_RDONLY);

    std::vector<std::string> groups;
    herr_t status = H5Ovisit(file.getId(), H5_INDEX_NAME, H5_ITER_NATIVE, op_func, static_cast<void*>(&groups), H5O_INFO_ALL);

    bool written = false;
    for (int i = 0; i < groups.size(); i++)
    {
        QString group = QString::fromStdString(groups[i]);

        std::cout << i << ": " << groups[i] << std::endl;
        if (group.startsWith("acquisition/"))
        {
            Recording recording;
            ReadTimeseries(file, group.toStdString(), recording);

            experiment.addAcquisition(recording);
        }
        if (group.startsWith("stimulus/presentation/"))
        {
            Recording recording;
            ReadTimeseries(file, group.toStdString(), recording);

            experiment.addStimulus(recording);
        }
    }

    file.close();

    for (int i = 0; i < experiment.getAcquisitions()[0].comments.size(); i++)
    {
        std::cout << "Comment: " << experiment.getAcquisitions()[0].comments[i] << std::endl;
    }

    //dataset.iterateAttrs((H5::attr_operator_t) attr_op);
}

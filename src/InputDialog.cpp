#include "InputDialog.h"

#include "PatchSeqDataLoader.h"

using namespace mv::gui;

InputDialog::InputDialog(QWidget* parent, PatchSeqDataLoader& plugin) :
    QDialog(parent),
    _gexprFilePicker(this, "Transcriptomics File"),
    _ephysFilePicker(this, "Electrophysiology File"),
    _morphoFilePicker(this, "Morphology File"),
    _metadataFilePicker(this, "Metadata File"),
    _morphologiesDirPicker(this, "Morphologies Directory"),
    _loadAction(this, "Load"),

    _groupAction(this, "Settings")
{
    setWindowTitle(tr("Load Patch-seq Data"));

    //_numberOfDimensionsAction.setDefaultWidgetFlags(IntegralAction::WidgetFlag::SpinBox);

    //QStringList pointDataTypes;
    //for (const char* const typeName : PointData::getElementTypeNames())
    //{
    //    pointDataTypes.append(QString::fromLatin1(typeName));
    //}
    //_storeAsAction.setOptions(pointDataTypes);

    // Load some settings
    //_dataTypeAction.setCurrentIndex(plugin.getSetting("DataType").toInt());
    //_numberOfDimensionsAction.setValue(plugin.getSetting("NumberOfDimensions").toInt());
    //_storeAsAction.setCurrentIndex(plugin.getSetting("StoreAs").toInt());

    _groupAction.addAction(&_gexprFilePicker);
    _groupAction.addAction(&_ephysFilePicker);
    _groupAction.addAction(&_morphoFilePicker);
    _groupAction.addAction(&_metadataFilePicker);
    _groupAction.addAction(&_morphologiesDirPicker);

    //_groupAction.addAction(&_datasetNameAction);
    //_groupAction.addAction(&_dataTypeAction);
    //_groupAction.addAction(&_numberOfDimensionsAction);
    //_groupAction.addAction(&_storeAsAction);
    //_groupAction.addAction(&_isDerivedAction);
    //_groupAction.addAction(&_datasetPickerAction);
    _groupAction.addAction(&_loadAction);

    auto layout = new QVBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(_groupAction.createWidget(this));

    setLayout(layout);

    //// Update the state of the dataset picker
    //const auto updateDatasetPicker = [this]() -> void {
    //    if (_isDerivedAction.isChecked()) {

    //        // Get unique identifier and gui names from all point data sets in the core
    //        auto dataSets = mv::data().getAllDatasets(std::vector<mv::DataType> {PointType});

    //        // Assign found dataset(s)
    //        _datasetPickerAction.setDatasets(dataSets);
    //    }
    //    else {

    //        // Assign found dataset(s)
    //        _datasetPickerAction.setDatasets(mv::Datasets());
    //    }

    //    // Disable dataset picker when not marked as derived
    //    _datasetPickerAction.setEnabled(_isDerivedAction.isChecked());
    //    };

    //// Populate source datasets once the dataset is marked as derived
    //connect(&_isDerivedAction, &ToggleAction::toggled, this, updateDatasetPicker);

    //// Update dataset picker at startup
    //updateDatasetPicker();

    // Accept when the load action is triggered
    connect(&_loadAction, &TriggerAction::triggered, this, [this, &plugin]() {

        //// Save some settings
        //plugin.setSetting("DataType", _dataTypeAction.getCurrentIndex());
        //plugin.setSetting("NumberOfDimensions", _numberOfDimensionsAction.getValue());
        //plugin.setSetting("StoreAs", _storeAsAction.getCurrentIndex());

        accept();
    });
}

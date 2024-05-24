#pragma once

#include <Dataset.h>
#include <actions/GroupAction.h>
#include <actions/IntegralAction.h>
#include <actions/OptionAction.h>
#include <actions/StringAction.h>
#include <actions/TriggerAction.h>
#include <actions/ToggleAction.h>
#include <actions/DatasetPickerAction.h>

#include <QDialog>

class PatchSeqDataLoader;

class InputDialog : public QDialog
{
    Q_OBJECT

public:
    InputDialog(QWidget* parent, PatchSeqDataLoader& plugin, QString fileName);

    /** Get preferred size */
    QSize sizeHint() const override {
        return QSize(400, 50);
    }

    /** Get minimum size hint*/
    QSize minimumSizeHint() const override {
        return sizeHint();
    }

    /** Get the GUI name of the loaded dataset */
    QString getDatasetName() const {
        return _datasetNameAction.getString();
    }

    /** Get the number of dimensions */
    std::int32_t getNumberOfDimensions() const {
        return _numberOfDimensionsAction.getValue();
    }

    /** Get the desired storage type */
    QString getStoreAs() const {
        return _storeAsAction.getCurrentText();
    }

    /** Get whether the dataset will be marked as derived */
    bool getIsDerived() const {
        return _isDerivedAction.isChecked();
    }

    /** Get smart pointer to dataset (if any) */
    mv::Dataset<mv::DatasetImpl> getSourceDataset() {
        return _datasetPickerAction.getCurrentDataset();
    }

protected:
    mv::gui::StringAction            _datasetNameAction;             /** Dataset name action */
    mv::gui::OptionAction            _dataTypeAction;                /** Data type action */
    mv::gui::IntegralAction          _numberOfDimensionsAction;      /** Number of dimensions action */
    mv::gui::OptionAction            _storeAsAction;                 /** Store as action */
    mv::gui::ToggleAction            _isDerivedAction;               /** Mark dataset as derived action */
    mv::gui::DatasetPickerAction     _datasetPickerAction;           /** Dataset picker action for picking source datasets */
    mv::gui::TriggerAction           _loadAction;                    /** Load action */
    mv::gui::GroupAction             _groupAction;                   /** Group action */
};

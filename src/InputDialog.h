#pragma once

#include <Dataset.h>
#include <actions/GroupAction.h>
#include <actions/TriggerAction.h>
#include <actions/FilePickerAction.h>
#include <actions/DirectoryPickerAction.h>

#include <QDialog>

class PatchSeqDataLoader;

class InputDialog : public QDialog
{
    Q_OBJECT

public:
    InputDialog(QWidget* parent, PatchSeqDataLoader& plugin);

    /** Get preferred size */
    QSize sizeHint() const override {
        return QSize(400, 50);
    }

    /** Get minimum size hint*/
    QSize minimumSizeHint() const override {
        return sizeHint();
    }

    QString getTranscriptomicsFilePath() const
    {
        return _gexprFilePicker.getFilePath();
    }

    QString getElectrophysiologyFilePath() const
    {
        return _ephysFilePicker.getFilePath();
    }

    QString getMorphologyFilePath() const
    {
        return _morphoFilePicker.getFilePath();
    }

    QString getMetadataFilePath() const
    {
        return _metadataFilePicker.getFilePath();
    }

    QDir getMorphologiesDir() const
    {
        return QDir(_morphologiesDirPicker.getDirectory());
    }

protected:
    mv::gui::FilePickerAction       _gexprFilePicker;               /** File picker action */
    mv::gui::FilePickerAction       _ephysFilePicker;               /** File picker action */
    mv::gui::FilePickerAction       _morphoFilePicker;              /** File picker action */
    mv::gui::FilePickerAction       _metadataFilePicker;            /** File picker action */
    mv::gui::DirectoryPickerAction  _morphologiesDirPicker;         /** File picker action */

    mv::gui::TriggerAction          _loadAction;                    /** Load action */
    mv::gui::GroupAction            _groupAction;                   /** Group action */
};

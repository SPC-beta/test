// Copyright (c) 2011-2015 The BZX Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef BZX_QT_OPENURIDIALOG_H
#define BZX_QT_OPENURIDIALOG_H

#include <QDialog>

namespace Ui {
    class OpenURIDialog;
}

class OpenURIDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OpenURIDialog(QWidget *parent);
    ~OpenURIDialog();

    QString getURI();

protected Q_SLOTS:
    void accept() override;

private:
    Ui::OpenURIDialog *ui;
};

#endif // BZX_QT_OPENURIDIALOG_H

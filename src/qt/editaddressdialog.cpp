// Copyright (c) 2011-2016 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "editaddressdialog.h"
#include "ui_editaddressdialog.h"

#include "addresstablemodel.h"
#include "guiutil.h"
#include "bip47/paymentcode.h"

#include <QDataWidgetMapper>
#include <QMessageBox>

EditAddressDialog::EditAddressDialog(Mode _mode, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EditAddressDialog),
    mapper(0),
    mode(_mode),
    model(0)
{
    ui->setupUi(this);

    GUIUtil::setupAddressWidget(ui->addressEdit, this);

    ui->addressEdit->setEnabled(true);
    switch(mode)
    {
    case NewReceivingAddress:
        setWindowTitle(tr("New transparent receiving address"));
        ui->addressEdit->setEnabled(false);
        break;
    case NewSendingAddress:
        setWindowTitle(tr("New transparent sending address"));
        break;
    case EditReceivingAddress:
        setWindowTitle(tr("Edit transparent receiving address"));
        ui->addressEdit->setEnabled(false);
        break;
    case EditSendingAddress:
        setWindowTitle(tr("Edit transparent sending address"));
        break;
    case NewPcode:
        setWindowTitle(tr("New RAP payment code"));
        ui->label_2->setText(tr("RAP address"));
        break;
    case EditPcode:
        setWindowTitle(tr("Edit RAP payment code"));
        ui->label_2->setText(tr("RAP address"));
        break;
    case NewSparkSendingAddress:
        setWindowTitle(tr("New spark sending address"));
        ui->label_2->setText(tr("Spark address"));
        break;
    case EditSparkSendingAddress:
        setWindowTitle(tr("Edit spark sending address"));
        ui->label_2->setText(tr("Spark address"));
        break;
    case NewSparkReceivingAddress:
        setWindowTitle(tr("New spark receiving address"));
        ui->addressEdit->setEnabled(false);
        break;
    case EditSparkReceivingAddress:
        setWindowTitle(tr("Edit spark receiving address"));
        ui->addressEdit->setEnabled(false);
        break;
    }

    mapper = new QDataWidgetMapper(this);
    mapper->setSubmitPolicy(QDataWidgetMapper::ManualSubmit);
}

EditAddressDialog::~EditAddressDialog()
{
    delete ui;
}

void EditAddressDialog::setModel(AddressTableModel *_model)
{
    this->model = _model;
    if(!_model)
        return;

    mapper->setModel(_model);
    mapper->addMapping(ui->labelEdit, AddressTableModel::Label);
    mapper->addMapping(ui->addressEdit, AddressTableModel::Address);
}

void EditAddressDialog::loadRow(int row)
{
    mapper->setCurrentIndex(row);
}

bool EditAddressDialog::saveCurrentRow()
{
    if(!model)
        return false;

    switch(mode)
    {
    case NewReceivingAddress:
    case NewSendingAddress:
        address = model->addRow(
                mode == NewSendingAddress ? AddressTableModel::Send : AddressTableModel::Receive,
                ui->labelEdit->text(),
                ui->addressEdit->text(),
                AddressTableModel::Transparent);
        break;
    case EditReceivingAddress:
    case EditSendingAddress:
        if(mapper->submit())
        {
            address = ui->addressEdit->text();
        }
        break;
    case NewPcode:
    case EditPcode:
        address = model->addRow(AddressTableModel::Send, ui->labelEdit->text(), ui->addressEdit->text(), AddressTableModel::RAP);
        break;
    case NewSparkReceivingAddress:
    case NewSparkSendingAddress:
        address = model->addRow(
                mode == NewSparkSendingAddress ? AddressTableModel::Send : AddressTableModel::Receive,
                ui->labelEdit->text(),
                ui->addressEdit->text(),
                AddressTableModel::Spark);
        break;
    case EditSparkReceivingAddress:
    case EditSparkSendingAddress:
        if(mapper->submit())
        {
            address = ui->addressEdit->text();
        }
        break;
    }
    return !address.isEmpty();
}

void EditAddressDialog::accept()
{
    if(!model)
        return;

    if(!saveCurrentRow())
    {
        switch(model->getEditStatus())
        {
        case AddressTableModel::OK:
            // Failed with unknown reason. Just reject.
            break;
        case AddressTableModel::NO_CHANGES:
            // No changes were made during edit operation. Just reject.
            break;
        case AddressTableModel::INVALID_ADDRESS:
            QMessageBox::warning(this, windowTitle(),
                tr("The entered address \"%1\" is not a valid address.").arg(ui->addressEdit->text()),
                QMessageBox::Ok, QMessageBox::Ok);
            break;
        case AddressTableModel::INVALID_SPARK_ADDRESS:
            QMessageBox::warning(this, windowTitle(),
                tr("The entered address \"%1\" is not a valid spark address.").arg(ui->addressEdit->text()),
                QMessageBox::Ok, QMessageBox::Ok);
            break;
        case AddressTableModel::DUPLICATE_ADDRESS:
            QMessageBox::warning(this, windowTitle(),
                tr("The entered address \"%1\" is already in the address book.").arg(ui->addressEdit->text()),
                QMessageBox::Ok, QMessageBox::Ok);
            break;
        case AddressTableModel::WALLET_UNLOCK_FAILURE:
            QMessageBox::critical(this, windowTitle(),
                tr("Could not unlock wallet."),
                QMessageBox::Ok, QMessageBox::Ok);
            break;
        case AddressTableModel::KEY_GENERATION_FAILURE:
            QMessageBox::critical(this, windowTitle(),
                tr("New key generation failed."),
                QMessageBox::Ok, QMessageBox::Ok);
            break;
        case AddressTableModel::PCODE_VALIDATION_FAILURE:
            QMessageBox::critical(this, windowTitle(),
                tr("New RAP address validation failed."),
                QMessageBox::Ok, QMessageBox::Ok);
            break;
        case AddressTableModel::PCODE_CANNOT_BE_LABELED:
            QMessageBox::critical(this, windowTitle(),
                tr("Receiving RAP addresses cannot be relabeled."),
                QMessageBox::Ok, QMessageBox::Ok);
            break;

        }
        return;
    }
    QDialog::accept();
}

QString EditAddressDialog::getAddress() const
{
    return address;
}

void EditAddressDialog::setAddress(const QString &_address)
{
    this->address = _address;
    ui->addressEdit->setText(_address);
}

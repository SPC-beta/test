// Copyright (c) 2011-2014 The BZX Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef BZX_QT_BZXADDRESSVALIDATOR_H
#define BZX_QT_BZXADDRESSVALIDATOR_H

#include <QValidator>

#ifdef ENABLE_WALLET
#include "../spark/sparkwallet.h"
#endif

#include "../spark/state.h"

/** Base58 entry widget validator, checks for valid characters and
 * removes some whitespace.
 */
class BZXAddressEntryValidator : public QValidator
{
    Q_OBJECT

public:
    explicit BZXAddressEntryValidator(QObject *parent);

    State validate(QString &input, int &pos) const override;
};

/** BZX address widget validator, checks for a valid BZX address.
 */
class BZXAddressCheckValidator : public QValidator
{
    Q_OBJECT

public:
    explicit BZXAddressCheckValidator(QObject *parent);

    State validate(QString &input, int &pos) const override;

    bool validateSparkAddress(const std::string& address) const;
};

#endif // BZX_QT_BZXADDRESSVALIDATOR_H

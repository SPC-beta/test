// Copyright (c) 2011-2014 The BZX Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef BZX_QT_WALLETMODELTRANSACTION_H
#define BZX_QT_WALLETMODELTRANSACTION_H

#include "../hdmint/hdmint.h"
#include "../primitives/mint_spend.h"
#include "spark/state.h"

#include "walletmodel.h"

#include <QObject>

class SendCoinsRecipient;

class CReserveKey;
class CWallet;
class CWalletTx;

/** Data model for a walletmodel transaction. */
class WalletModelTransaction
{
public:
    explicit WalletModelTransaction(const QList<SendCoinsRecipient> &recipients);
    ~WalletModelTransaction();

    QList<SendCoinsRecipient> getRecipients();

    CWalletTx *getTransaction();
    unsigned int getTransactionSize();

    void setTransactionFee(const CAmount& newFee);
    CAmount getTransactionFee();

    CAmount getTotalTransactionAmount();

    void newPossibleKeyChange(CWallet *wallet);
    CReserveKey *getPossibleKeyChange();

    void reassignAmounts(int nChangePosRet); // needed for the subtract-fee-from-amount feature

    std::vector<CLelantusEntry>& getSpendCoins();
    std::vector<CSigmaEntry>& getSigmaSpendCoins();
    std::vector<CHDMint>& getMintCoins();

private:
    QList<SendCoinsRecipient> recipients;
    CWalletTx *walletTransaction;
    CReserveKey *keyChange;
    CAmount fee;

    // lelantus transaction
    std::vector<CLelantusEntry> spendCoins;
    std::vector<CSigmaEntry> sigmaSpendCoins;
    std::vector<CHDMint> mintCoins;
};

#endif // BZX_QT_WALLETMODELTRANSACTION_H

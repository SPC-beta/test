// Copyright (c) 2016 The BZX Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef BZX_WALLET_RPCWALLET_H
#define BZX_WALLET_RPCWALLET_H

class CRPCTable;
class JSONRPCRequest;

void RegisterWalletRPCCommands(CRPCTable &t);

/**
 * Figures out what wallet, if any, to use for a JSONRPCRequest.
 *
 * @param[in] request JSONRPCRequest that wishes to access a wallet
 * @return NULL if no wallet should be used, or a pointer to the CWallet
 */
CWallet *GetWalletForJSONRPCRequest(const JSONRPCRequest& request);

std::string HelpRequiringPassphrase(CWallet *);
void EnsureWalletIsUnlocked(CWallet *);
bool EnsureWalletIsAvailable(CWallet *, bool avoidException);

#endif //BZX_WALLET_RPCWALLET_H

// Copyright (c) 2014-2019 The Dash Core developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "masternode-utils.h"

#include "init.h"
#include "masternode-sync.h"
#ifdef ENABLE_WALLET
//#include "privatesend-client.h"
#endif
#include "validation.h"

struct CompareScoreMN
{
    bool operator()(const std::pair<arith_uint256, const CDeterministicMNCPtr&>& t1,
                    const std::pair<arith_uint256, const CDeterministicMNCPtr&>& t2) const
    {
        return (t1.first != t2.first) ? (t1.first < t2.first) : (t1.second->collateralOutpoint < t2.second->collateralOutpoint);
    }
};

bool CMasternodeUtils::GetMasternodeScores(const uint256& nBlockHash, score_pair_vec_t& vecMasternodeScoresRet)
{
    vecMasternodeScoresRet.clear();

    auto mnList = deterministicMNManager->GetListAtChainTip();
    auto scores = mnList.CalculateScores(nBlockHash);
    for (const auto& p : scores) {
        vecMasternodeScoresRet.emplace_back(p.first, p.second);
    }

    std::sort(vecMasternodeScoresRet.rbegin(), vecMasternodeScoresRet.rend(), CompareScoreMN());
    return !vecMasternodeScoresRet.empty();
}

bool CMasternodeUtils::GetMasternodeRank(const COutPoint& outpoint, int& nRankRet, uint256& blockHashRet, int nBlockHeight)
{
    nRankRet = -1;

    if (!masternodeSync.IsBlockchainSynced())
        return false;

    // make sure we know about this block
    blockHashRet = uint256();
    if (!GetBlockHash(blockHashRet, nBlockHeight)) {
        LogPrintf("CMasternodeUtils::%s -- ERROR: GetBlockHash() failed at nBlockHeight %d\n", __func__, nBlockHeight);
        return false;
    }

    score_pair_vec_t vecMasternodeScores;
    if (!GetMasternodeScores(blockHashRet, vecMasternodeScores))
        return false;

    int nRank = 0;
    for (const auto& scorePair : vecMasternodeScores) {
        nRank++;
        if(scorePair.second->collateralOutpoint == outpoint) {
            nRankRet = nRank;
            return true;
        }
    }

    return false;
}

void CMasternodeUtils::DoMaintenance(CConnman& connman)
{
    if(fLiteMode) return; // disable all Dash specific functionality

    if(!masternodeSync.IsBlockchainSynced() || ShutdownRequested())
        return;
}


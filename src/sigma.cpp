#include "validation.h"
#include "sigma.h"
#include "timedata.h"
#include "chainparams.h"
#include "util.h"
#include "base58.h"
#include "definition.h"
#include "txmempool.h"
#include "wallet/wallet.h"
#include "wallet/walletdb.h"
#include "crypto/sha256.h"
#include "sigma/coinspend.h"
#include "sigma/coin.h"
#include "primitives/mint_spend.h"
#include "batchproof_container.h"

#include <atomic>
#include <sstream>
#include <chrono>

#include <boost/foreach.hpp>
#include <boost/scope_exit.hpp>

#include <ios>

int64_t nMinimumInputValue = DUST_HARD_LIMIT;

namespace sigma {

static CSigmaState sigmaState;

bool CheckSigmaSpendSerial(
        CValidationState &state,
        CSigmaTxInfo *sigmaTxInfo,
        const Scalar &serial,
        int nHeight,
        bool fConnectTip) {

    return false;
}

bool IsSigmaAllowed()
{
    LOCK(cs_main);
    return IsSigmaAllowed(chainActive.Height());
}

bool IsSigmaAllowed(int height)
{
    return false;
}

bool IsRemintWindow(int height) {

    return false;
}

secp_primitives::GroupElement ParseSigmaMintScript(const CScript& script)
{
    if (script.size() < 1) {
        throw std::invalid_argument("Script is not a valid Sigma mint");
    }

    std::vector<unsigned char> serialized(script.begin() + 1, script.end());

    secp_primitives::GroupElement pub;
    if (serialized.size() < pub.memoryRequired()) {
        throw std::invalid_argument("Script is not a valid Sigma mint");
    }
    pub.deserialize(serialized.data());

    return pub;
}

std::pair<std::unique_ptr<sigma::CoinSpend>, uint32_t> ParseSigmaSpend(const CTxIn& in)
{
    uint32_t groupId = in.prevout.n;

    if (groupId < 1 || groupId >= INT_MAX || in.scriptSig.size() < 1) {
        throw CBadTxIn();
    }

    CDataStream serialized(
        std::vector<unsigned char>(in.scriptSig.begin() + 1, in.scriptSig.end()),
        SER_NETWORK,
        PROTOCOL_VERSION
    );

    std::unique_ptr<sigma::CoinSpend> spend(new sigma::CoinSpend(sigma::Params::get_default(), serialized));

    return std::make_pair(std::move(spend), groupId);
}

// This function will not report an error only if the transaction is sigma spend.
CAmount GetSpendAmount(const CTxIn& in) {
    if (in.IsSigmaSpend()) {
        std::unique_ptr<sigma::CoinSpend> spend;

        try {
            std::tie(spend, std::ignore) = ParseSigmaSpend(in);
        } catch (const std::ios_base::failure& e) {
            LogPrintf("GetSpendAmount: io error %s\n", e.what());
            return 0;
        } catch (const CBadTxIn& e) {
            LogPrintf("GetSpendAmount: %s\n", e.what());
            return 0;
        }

        return spend->getIntDenomination();
    }
    return 0;
}

CAmount GetSpendAmount(const CTransaction& tx) {
    CAmount sum(0);
    for (const auto& vin : tx.vin) {
        sum += GetSpendAmount(vin);
    }
    return sum;
}

bool CheckSigmaBlock(CValidationState &state, const CBlock& block) {

    return false;
}

// Will return false for V1, V1.5 and V2 spends.
// Mixing V2 and sigma spends into the same transaction will fail.
bool CheckSigmaSpendTransaction(
        const CTransaction &tx,
        const std::vector<sigma::CoinDenomination>& targetDenominations,
        CValidationState &state,
        uint256 hashTx,
        bool isVerifyDB,
        int nHeight,
        int nRealHeight,
        bool isCheckWallet,
        bool fStatefulSigmaCheck,
        CSigmaTxInfo *sigmaTxInfo) {

    return false;
}

bool CheckSigmaMintTransaction(
        const CTxOut &txout,
        CValidationState &state,
        uint256 hashTx,
        bool fStatefulSigmaCheck,
        CSigmaTxInfo *sigmaTxInfo) {

    return false;
}

bool CheckSigmaTransaction(
        const CTransaction &tx,
        CValidationState &state,
        uint256 hashTx,
        bool isVerifyDB,
        int nHeight,
        bool isCheckWallet,
        bool fStatefulSigmaCheck,
        CSigmaTxInfo *sigmaTxInfo)
{

    return false;
}

void RemoveSigmaSpendsReferencingBlock(CTxMemPool& pool, CBlockIndex* blockIndex) {
    LOCK2(cs_main, pool.cs);
    std::vector<CTransaction> txn_to_remove;
    for (CTxMemPool::txiter mi = pool.mapTx.begin(); mi != pool.mapTx.end(); ++mi) {
        const CTransaction& tx = mi->GetTx();
        if (tx.IsSigmaSpend()) {
            // Run over all the inputs, check if their Accumulator block hash is equal to
            // block removed. If any one is equal, remove txn from mempool.
            for (const CTxIn& txin : tx.vin) {
                if (txin.IsSigmaSpend()) {
                    std::unique_ptr<sigma::CoinSpend> spend;
                    uint32_t pubcoinId;
                    std::tie(spend, pubcoinId) = ParseSigmaSpend(txin);
                    uint256 accumulatorBlockHash = spend->getAccumulatorBlockHash();
                    if (accumulatorBlockHash == blockIndex->GetBlockHash()) {
                        // Do not remove transaction immediately, that will invalidate iterator mi.
                        txn_to_remove.push_back(tx);
                        break;
                    }
                }
            }
        }
    }
    for (const CTransaction& tx: txn_to_remove) {
        // Remove txn from mempool.
        pool.removeRecursive(tx);
        LogPrintf("DisconnectTipSigma: removed sigma spend which referenced a removed blockchain tip.");
    }
}

void DisconnectTipSigma(CBlock& block, CBlockIndex *pindexDelete) {
    sigmaState.RemoveBlock(pindexDelete);

    // Also remove from mempool sigma spends that reference given block hash.
    RemoveSigmaSpendsReferencingBlock(mempool, pindexDelete);
    RemoveSigmaSpendsReferencingBlock(txpools.getStemTxPool(), pindexDelete);
}

Scalar GetSigmaSpendSerialNumber(const CTransaction &tx, const CTxIn &txin) {
    if (!tx.IsSigmaSpend())
        return Scalar(uint64_t(0));

    try {
        // NOTE(martun): +1 on the next line stands for 1 byte in which the opcode of
        // OP_SIGMASPEND is written. In privcoin you will see +4 instead,
        // because the size of serialized spend is also written, probably in 3 bytes.
        CDataStream serializedCoinSpend(
                (const char *)&*(txin.scriptSig.begin() + 1),
                (const char *)&*txin.scriptSig.end(),
                SER_NETWORK, PROTOCOL_VERSION);
        sigma::CoinSpend spend(sigma::Params::get_default(), serializedCoinSpend);
        return spend.getCoinSerialNumber();
    }
    catch (const std::ios_base::failure &) {
        return Scalar(uint64_t(0));
    }
}

CAmount GetSigmaSpendInput(const CTransaction &tx) {
    if (!tx.IsSigmaSpend())
        return CAmount(0);

    try {
        CAmount sum(0);
        BOOST_FOREACH(const CTxIn& txin, tx.vin){
            // NOTE(martun): +1 on the next line stands for 1 byte in which the opcode of
            // OP_PRIVCOINSPENDV3 is written. In privcoin you will see +4 instead,
            // because the size of serialized spend is also written, probably in 3 bytes.
            CDataStream serializedCoinSpend(
                    (const char *)&*(txin.scriptSig.begin() + 1),
                    (const char *)&*txin.scriptSig.end(),
                    SER_NETWORK, PROTOCOL_VERSION);
            sigma::CoinSpend spend(sigma::Params::get_default(), serializedCoinSpend);
            sum += spend.getIntDenomination();
        }
        return sum;
    }
    catch (const std::runtime_error &) {
        return CAmount(0);
    }
}


/**
 * Connect a new ZCblock to chainActive. pblock is either NULL or a pointer to a CBlock
 * corresponding to pindexNew, to bypass loading it again from disk.
 */
bool ConnectBlockSigma(
        CValidationState &state,
        const CChainParams &chainparams,
        CBlockIndex *pindexNew,
        const CBlock *pblock,
        bool fJustCheck) {

    return false;
}

bool GetOutPointFromBlock(COutPoint& outPoint, const GroupElement &pubCoinValue, const CBlock &block){
    secp_primitives::GroupElement txPubCoinValue;
    // cycle transaction hashes, looking for this pubcoin.
    BOOST_FOREACH(CTransactionRef tx, block.vtx){
        uint32_t nIndex = 0;
        for (const CTxOut &txout: tx->vout) {
            if (txout.scriptPubKey.IsSigmaMint()){

                // If you wonder why +1, go to file wallet.cpp and read the comments in function
                // CWallet::CreatePrivcoinMintModelV3 around "scriptSerializedCoin << OP_PRIVCOINMINTV3";
                std::vector<unsigned char> coin_serialised(txout.scriptPubKey.begin() + 1,
                                                      txout.scriptPubKey.end());
                try {
                    txPubCoinValue.deserialize(&coin_serialised[0]);
                } catch (const std::exception &) {
                    return false;
                }
                if(pubCoinValue==txPubCoinValue){
                    outPoint = COutPoint(tx->GetHash(), nIndex);
                    return true;
                }
            }
            nIndex++;
        }
    }

    return false;
}

bool GetOutPoint(COutPoint& outPoint, const sigma::PublicCoin &pubCoin) {

    sigma::CSigmaState *sigmaState = sigma::CSigmaState::GetState();
    auto mintedCoinHeightAndId = sigmaState->GetMintedCoinHeightAndId(pubCoin);
    int mintHeight = mintedCoinHeightAndId.first;
    int coinId = mintedCoinHeightAndId.second;

    if(mintHeight==-1 && coinId==-1)
        return false;

    // get block containing mint
    CBlockIndex *mintBlock = chainActive[mintHeight];
    CBlock block;
    if(!ReadBlockFromDisk(block, mintBlock, ::Params().GetConsensus()))
        LogPrintf("can't read block from disk.\n");

    return GetOutPointFromBlock(outPoint, pubCoin.getValue(), block);
}

bool GetOutPoint(COutPoint& outPoint, const GroupElement &pubCoinValue) {
    int mintHeight = 0;
    int coinId = 0;

    sigma::CSigmaState *sigmaState = sigma::CSigmaState::GetState();
    std::vector<sigma::CoinDenomination> denominations;
    GetAllDenoms(denominations);
    BOOST_FOREACH(sigma::CoinDenomination denomination, denominations){
        sigma::PublicCoin pubCoin(pubCoinValue, denomination);
        auto mintedCoinHeightAndId = sigmaState->GetMintedCoinHeightAndId(pubCoin);
        mintHeight = mintedCoinHeightAndId.first;
        coinId = mintedCoinHeightAndId.second;
        if(mintHeight!=-1 && coinId!=-1)
            break;
    }

    if(mintHeight==-1 && coinId==-1)
        return false;

    // get block containing mint
    CBlockIndex *mintBlock = chainActive[mintHeight];
    CBlock block;
    if(!ReadBlockFromDisk(block, mintBlock, ::Params().GetConsensus()))
        LogPrintf("can't read block from disk.\n");

    return GetOutPointFromBlock(outPoint, pubCoinValue, block);
}

bool GetOutPoint(COutPoint& outPoint, const uint256 &pubCoinValueHash) {
    GroupElement pubCoinValue;
    sigma::CSigmaState *sigmaState = sigma::CSigmaState::GetState();
    if(!sigmaState->HasCoinHash(pubCoinValue, pubCoinValueHash)){
        return false;
    }

    return GetOutPoint(outPoint, pubCoinValue);
}

bool BuildSigmaStateFromIndex(CChain *chain) {
    for (CBlockIndex *blockIndex = chain->Genesis(); blockIndex; blockIndex=chain->Next(blockIndex))
    {
        sigmaState.AddBlock(blockIndex);
    }
    // DEBUG
    LogPrintf(
        "Latest IDs for sigma coin groups are %d, %d, %d, %d, %d\n",
        sigmaState.GetLatestCoinID(CoinDenomination::SIGMA_DENOM_0_1),
        sigmaState.GetLatestCoinID(CoinDenomination::SIGMA_DENOM_0_5),
        sigmaState.GetLatestCoinID(CoinDenomination::SIGMA_DENOM_1),
        sigmaState.GetLatestCoinID(CoinDenomination::SIGMA_DENOM_10),
        sigmaState.GetLatestCoinID(CoinDenomination::SIGMA_DENOM_100));
    return true;
}

// CSigmaTxInfo

void CSigmaTxInfo::Complete() {
    // We need to sort mints lexicographically by serialized value of pubCoin. That's the way old code
    // works, we need to stick to it. Denomination doesn't matter but we will sort by it as well
    sort(mints.begin(), mints.end(),
            [](decltype(mints)::const_reference m1, decltype(mints)::const_reference m2)->bool {
            CDataStream ds1(SER_DISK, CLIENT_VERSION), ds2(SER_DISK, CLIENT_VERSION);
            ds1 << m1;
            ds2 << m2;
            return ds1.str() < ds2.str();
            });

    // Mark this info as complete
    fInfoIsComplete = true;
}

/******************************************************************************/
// CSigmaState::Containers
/******************************************************************************/

CSigmaState::Containers::Containers(std::atomic<bool> & surgeCondition)
: surgeCondition(surgeCondition)
{}

void CSigmaState::Containers::AddMint(sigma::PublicCoin const & pubCoin, CMintedCoinInfo const & coinInfo) {
    mintedPubCoins.insert(std::make_pair(pubCoin, coinInfo));
    mintMetaInfo[coinInfo.coinGroupId][coinInfo.denomination] += 1;
    CheckSurgeCondition(coinInfo.coinGroupId, coinInfo.denomination);
}

void CSigmaState::Containers::RemoveMint(sigma::PublicCoin const & pubCoin) {
    mint_info_container::const_iterator iter = mintedPubCoins.find(pubCoin);
    if (iter != mintedPubCoins.end()) {
        mintMetaInfo[iter->second.coinGroupId][iter->second.denomination] -= 1;
        CMintedCoinInfo tmpMintInfo(iter->second);
        mintedPubCoins.erase(iter);
        CheckSurgeCondition(tmpMintInfo.coinGroupId, tmpMintInfo.denomination);
    }
}

void CSigmaState::Containers::AddSpend(Scalar const & serial, CSpendCoinInfo const & coinInfo) {
    usedCoinSerials[serial] = coinInfo;
    spendMetaInfo[coinInfo.coinGroupId][coinInfo.denomination] += 1;
    CheckSurgeCondition(coinInfo.coinGroupId, coinInfo.denomination);
}

void CSigmaState::Containers::RemoveSpend(Scalar const & serial) {
    spend_info_container::const_iterator iter = usedCoinSerials.find(serial);
    if (iter != usedCoinSerials.end()) {
        spendMetaInfo[iter->second.coinGroupId][iter->second.denomination] -= 1;
        CSpendCoinInfo tmpSpendInfo(iter->second);
        usedCoinSerials.erase(iter);
        CheckSurgeCondition(tmpSpendInfo.coinGroupId, tmpSpendInfo.denomination);
    }
}

mint_info_container const & CSigmaState::Containers::GetMints() const {
    return mintedPubCoins;
}

spend_info_container const & CSigmaState::Containers::GetSpends() const {
    return usedCoinSerials;
}

bool CSigmaState::Containers::IsSurgeCondition() const {
    return surgeCondition;
}

void CSigmaState::Containers::Reset() {
    mintedPubCoins.clear();
    usedCoinSerials.clear();
    mintMetaInfo.clear();
    spendMetaInfo.clear();
    surgeCondition = false;
}

void CSigmaState::Containers::CheckSurgeCondition(int groupId, CoinDenomination denom) {
    bool result = spendMetaInfo[groupId][denom] > mintMetaInfo[groupId][denom];
    if( result ) {
        std::ostringstream ostr;
        ostr << "Turning sigma surge protection ON: groupId: " << groupId << ", denomination: " << denom << '\n';
        error(ostr.str().c_str());
    }

    for(metainfo_container_t::const_iterator smi = spendMetaInfo.begin(); smi != spendMetaInfo.end() && !result; ++smi) {
        for(std::map<CoinDenomination, size_t>::const_iterator di = smi->second.begin(); di != smi->second.end() && !result; ++di) {
            if(di->second > mintMetaInfo[smi->first][di->first]) {
                result = true;
            }
        }
    }
    surgeCondition = result;
}

/******************************************************************************/
// CSigmaState
/******************************************************************************/

CSigmaState::CSigmaState()
:containers(surgeCondition)
{}

void CSigmaState::AddMintsToStateAndBlockIndex(
        CBlockIndex *index,
        const CBlock* pblock) {

    std::unordered_map<sigma::CoinDenomination, std::vector<sigma::PublicCoin>> blockDenomMints;
    for (const auto& mint : pblock->sigmaTxInfo->mints) {
        blockDenomMints[mint.getDenomination()].push_back(mint);
    }

    for (const auto& it : blockDenomMints) {
        const sigma::CoinDenomination denomination = it.first;
        const std::vector<sigma::PublicCoin>& mintsWithThisDenom = it.second;

        if (mintsWithThisDenom.empty())
            continue;

        if (latestCoinIds[denomination] < 1)
            latestCoinIds[denomination] = 1;
        auto mintCoinGroupId = latestCoinIds[denomination];


        SigmaCoinGroupInfo &coinGroup = coinGroups[std::make_pair(denomination, mintCoinGroupId)];

        if (coinGroup.nCoins + mintsWithThisDenom.size() <= 0) {
            if (coinGroup.nCoins == 0) {
                // first group of coins for given denomination
                assert(coinGroup.firstBlock == nullptr);
                assert(coinGroup.lastBlock == nullptr);

                coinGroup.firstBlock = coinGroup.lastBlock = index;
            } else {
                assert(coinGroup.firstBlock != nullptr);
                assert(coinGroup.lastBlock != nullptr);
                assert(coinGroup.lastBlock->nHeight <= index->nHeight);

                coinGroup.lastBlock = index;
            }
            coinGroup.nCoins += mintsWithThisDenom.size();
        }
        else {
            latestCoinIds[denomination] = ++mintCoinGroupId;
            SigmaCoinGroupInfo& newCoinGroup = coinGroups[std::make_pair(denomination, mintCoinGroupId)];
            newCoinGroup.firstBlock = newCoinGroup.lastBlock = index;
            newCoinGroup.nCoins = mintsWithThisDenom.size();
        }

        for (const auto& mint : mintsWithThisDenom) {
            containers.AddMint(mint, CMintedCoinInfo::make(denomination, mintCoinGroupId, index->nHeight));

            LogPrintf("AddMintsToStateAndBlockIndex: mint added denomination=%d, id=%d\n", denomination, mintCoinGroupId);
            index->sigmaMintedPubCoins[{denomination, mintCoinGroupId}].push_back(mint);
        }
    }
}

void CSigmaState::AddSpend(const Scalar &serial, CoinDenomination denom, int coinGroupId) {
    containers.AddSpend(serial, CSpendCoinInfo::make(denom, coinGroupId));
}

void CSigmaState::AddBlock(CBlockIndex *index) {
    BOOST_FOREACH(
        const PAIRTYPE(PAIRTYPE(sigma::CoinDenomination, int), std::vector<sigma::PublicCoin>) &pubCoins,
            index->sigmaMintedPubCoins) {

        if (pubCoins.second.empty())
            continue;

        SigmaCoinGroupInfo& coinGroup = coinGroups[pubCoins.first];

        if (coinGroup.firstBlock == NULL)
            coinGroup.firstBlock = index;
        coinGroup.lastBlock = index;
        coinGroup.nCoins += pubCoins.second.size();

        latestCoinIds[pubCoins.first.first] = pubCoins.first.second;
        BOOST_FOREACH(const sigma::PublicCoin &coin, pubCoins.second) {
            containers.AddMint(coin, CMintedCoinInfo::make(pubCoins.first.first, pubCoins.first.second, index->nHeight));
        }
    }

    BOOST_FOREACH(const spend_info_container::value_type &serial, index->sigmaSpentSerials) {
        AddSpend(serial.first, serial.second.denomination, serial.second.coinGroupId);
    }
}

void CSigmaState::RemoveBlock(CBlockIndex *index) {
    // roll back accumulator updates
    BOOST_FOREACH(
        const PAIRTYPE(PAIRTYPE(sigma::CoinDenomination, int),std::vector<sigma::PublicCoin>) &coin,
        index->sigmaMintedPubCoins)
    {
        SigmaCoinGroupInfo   &coinGroup = coinGroups[coin.first];
        int  nMintsToForget = coin.second.size();

        if (nMintsToForget == 0)
            continue;

        assert(coinGroup.nCoins >= nMintsToForget);

        if ((coinGroup.nCoins -= nMintsToForget) == 0) {
            // all the coins of this group have been erased, remove the group altogether
            coinGroups.erase(coin.first);
            // decrease pubcoin id for this denomination
            latestCoinIds[coin.first.first]--;
            if (0 == latestCoinIds[coin.first.first]) {
                latestCoinIds.erase(coin.first.first);
            }
        }
        else {
            // roll back lastBlock to previous position
            assert(coinGroup.lastBlock == index);

            do {
                assert(coinGroup.lastBlock != coinGroup.firstBlock);
                coinGroup.lastBlock = coinGroup.lastBlock->pprev;
            } while (coinGroup.lastBlock->sigmaMintedPubCoins.count(coin.first) == 0 ||
                        coinGroup.lastBlock->sigmaMintedPubCoins[coin.first].size() == 0);
        }
    }

    // roll back mints
    BOOST_FOREACH(const PAIRTYPE(PAIRTYPE(sigma::CoinDenomination, int),std::vector<sigma::PublicCoin>) &pubCoins,
                  index->sigmaMintedPubCoins) {
        BOOST_FOREACH(const sigma::PublicCoin &coin, pubCoins.second) {
            auto coins = containers.GetMints().equal_range(coin);
            auto coinIt = find_if(
                coins.first, coins.second,
                [&pubCoins](const mint_info_container::value_type &v) {
                    return v.second.denomination == pubCoins.first.first &&
                        v.second.coinGroupId == pubCoins.first.second;
                });
            assert(coinIt != coins.second);
            containers.RemoveMint(coinIt->first);
        }
    }

    // roll back spends
    BOOST_FOREACH(const spend_info_container::value_type &serial, index->sigmaSpentSerials) {
        containers.RemoveSpend(serial.first);
    }
}

bool CSigmaState::GetCoinGroupInfo(
        sigma::CoinDenomination denomination,
        int group_id,
        SigmaCoinGroupInfo& result) {
    std::pair<sigma::CoinDenomination, int> key =
        std::make_pair(denomination, group_id);
    if (coinGroups.count(key) == 0)
        return false;

    result = coinGroups[key];
    return true;
}

bool CSigmaState::IsUsedCoinSerial(const Scalar &coinSerial) {
    return containers.GetSpends().count(coinSerial) != 0;
}

bool CSigmaState::IsUsedCoinSerialHash(Scalar &coinSerial, const uint256 &coinSerialHash) {
    for ( auto it = GetSpends().begin(); it != GetSpends().end(); ++it ){
        if(primitives::GetSerialHash(it->first)==coinSerialHash){
            coinSerial = it->first;
            return true;
        }
    }
    return false;
}

bool CSigmaState::HasCoin(const sigma::PublicCoin& pubCoin) {
    return containers.GetMints().find(pubCoin) != containers.GetMints().end();
}

bool CSigmaState::HasCoinHash(GroupElement &pubCoinValue, const uint256 &pubCoinValueHash) {
    for ( auto it = GetMints().begin(); it != GetMints().end(); ++it ){
        const sigma::PublicCoin & pubCoin = (*it).first;
        if(pubCoin.getValueHash()==pubCoinValueHash){
            pubCoinValue = pubCoin.getValue();
            return true;
        }
    }
    return false;
}

int CSigmaState::GetCoinSetForSpend(
        CChain *chain,
        int maxHeight,
        sigma::CoinDenomination denomination,
        int coinGroupID,
        uint256& blockHash_out,
        std::vector<sigma::PublicCoin>& coins_out) {

    coins_out.clear();

    std::pair<sigma::CoinDenomination, int> denomAndId = std::make_pair(denomination, coinGroupID);

    if (coinGroups.count(denomAndId) == 0)
        return 0;

    SigmaCoinGroupInfo coinGroup = coinGroups[denomAndId];

    int numberOfCoins = 0;
    for (CBlockIndex *block = coinGroup.lastBlock;
            ;
            block = block->pprev) {
        if (block->sigmaMintedPubCoins.count(denomAndId) > 0 &&
                block->sigmaMintedPubCoins[denomAndId].size() > 0) {
            if (block->nHeight <= maxHeight) {
                if (numberOfCoins == 0) {
                    // latest block satisfying given conditions
                    // remember block hash
                    blockHash_out = block->GetBlockHash();
                }
                BOOST_FOREACH(const sigma::PublicCoin& pubCoinValue,
                        block->sigmaMintedPubCoins[denomAndId]) {
                    if (true) {
                        if (::Params().GetConsensus().sigmaBlacklist.count(pubCoinValue.getValue()) > 0) {
                            continue;
                        }
                    }
                    coins_out.push_back(pubCoinValue);
                    numberOfCoins++;
                }
            }
        }
        if (block == coinGroup.firstBlock) {
            break ;
        }
    }
    return numberOfCoins;
}

void CSigmaState::GetAnonymitySet(
        sigma::CoinDenomination denomination,
        int coinGroupID,
        bool fStartSigmaBlacklist,
        std::vector<GroupElement>& coins_out) {

    coins_out.clear();

    std::pair<sigma::CoinDenomination, int> denomAndId = std::make_pair(denomination, coinGroupID);

    if (coinGroups.count(denomAndId) == 0)
        return;

    SigmaCoinGroupInfo coinGroup = coinGroups[denomAndId];
    int maxHeight = (chainActive.Height() - (ZC_MINT_CONFIRMATIONS - 1));

    for (CBlockIndex *block = coinGroup.lastBlock;
            ;
            block = block->pprev) {
        if (block->sigmaMintedPubCoins.count(denomAndId) > 0 &&
                block->sigmaMintedPubCoins[denomAndId].size() > 0) {
            if (block->nHeight <= maxHeight) {
                BOOST_FOREACH(const sigma::PublicCoin& pubCoinValue,
                        block->sigmaMintedPubCoins[denomAndId]) {
                        if (::Params().GetConsensus().sigmaBlacklist.count(pubCoinValue.getValue()) > 0) {
                            continue;
                        }
                    coins_out.push_back(pubCoinValue.getValue());
                }
            }
        }
        if (block == coinGroup.firstBlock) {
            break ;
        }
    }
}

std::pair<int, int> CSigmaState::GetMintedCoinHeightAndId(
        const sigma::PublicCoin& pubCoin) {
    auto coinIt = containers.GetMints().find(pubCoin);

    if (coinIt != containers.GetMints().end()) {
        return std::make_pair(coinIt->second.nHeight, coinIt->second.coinGroupId);
    }
    return std::make_pair(-1, -1);
}

bool CSigmaState::AddSpendToMempool(const std::vector<Scalar> &coinSerials, uint256 txHash) {
    BOOST_FOREACH(Scalar coinSerial, coinSerials){
        if (IsUsedCoinSerial(coinSerial) || mempoolCoinSerials.count(coinSerial))
            return false;

        mempoolCoinSerials[coinSerial] = txHash;
    }

    return true;
}

bool CSigmaState::AddSpendToMempool(const Scalar &coinSerial, uint256 txHash) {
    if (IsUsedCoinSerial(coinSerial) || mempoolCoinSerials.count(coinSerial))
        return false;

    mempoolCoinSerials[coinSerial] = txHash;
    return true;
}

void CSigmaState::RemoveSpendFromMempool(const Scalar& coinSerial) {
    mempoolCoinSerials.erase(coinSerial);
}

void CSigmaState::AddMintsToMempool(const std::vector<GroupElement>& pubCoins){
    BOOST_FOREACH(const GroupElement& pubCoin, pubCoins){
        mempoolMints.insert(pubCoin);
    }
}

void CSigmaState::RemoveMintFromMempool(const GroupElement& pubCoin){
    mempoolMints.erase(pubCoin);
}

uint256 CSigmaState::GetMempoolConflictingTxHash(const Scalar& coinSerial) {
    if (mempoolCoinSerials.count(coinSerial) == 0)
        return uint256();

    return mempoolCoinSerials[coinSerial];
}

bool CSigmaState::CanAddSpendToMempool(const Scalar& coinSerial) {
    return !IsUsedCoinSerial(coinSerial) && mempoolCoinSerials.count(coinSerial) == 0;
}

bool CSigmaState::CanAddMintToMempool(const GroupElement& pubCoin){
    return mempoolMints.count(pubCoin) == 0;
}

void CSigmaState::Reset() {
    coinGroups.clear();
    latestCoinIds.clear();
    mempoolCoinSerials.clear();
    mempoolMints.clear();
    containers.Reset();
}

CSigmaState* CSigmaState::GetState() {
    return &sigmaState;
}

int CSigmaState::GetLatestCoinID(sigma::CoinDenomination denomination) const {
    auto iter = latestCoinIds.find(denomination);
    if (iter == latestCoinIds.end()) {
        // Do not throw here, if there was no sigma mint, that's fine.
        return 0;
    }
    return iter->second;
}

bool CSigmaState::IsSurgeConditionDetected() const {
    return surgeCondition;
}

mint_info_container const & CSigmaState::GetMints() const {
    return containers.GetMints();
}

spend_info_container const & CSigmaState::GetSpends() const {
    return containers.GetSpends();
}

std::unordered_map<std::pair<CoinDenomination, int>, CSigmaState::SigmaCoinGroupInfo, CSigmaState::pairhash> const & CSigmaState::GetCoinGroups() const {
    return coinGroups;
}

std::unordered_map<CoinDenomination, int> const & CSigmaState::GetLatestCoinIds() const {
    return latestCoinIds;
}

std::unordered_map<Scalar, uint256, CScalarHash> const & CSigmaState::GetMempoolCoinSerials() const {
    return mempoolCoinSerials;
}

} // end of namespace sigma.

// Copyright (c) 2011-2016 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "clientmodel.h"

#include "bantablemodel.h"
#include "guiconstants.h"
#include "guiutil.h"
#include "peertablemodel.h"

#include "chainparams.h"
#include "checkpoints.h"
#include "clientversion.h"
#include "validation.h"
#include "net.h"
#include "txmempool.h"
#include "ui_interface.h"
#include "util.h"

#include <stdint.h>

#include <QDebug>
#include <QTimer>

class CBlockIndex;

static const int64_t nClientStartupTime = GetTime();
static int64_t nLastHeaderTipUpdateNotification = 0;
static int64_t nLastBlockTipUpdateNotification = 0;

ClientModel::ClientModel(OptionsModel *_optionsModel, QObject *parent) :
    QObject(parent),
    optionsModel(_optionsModel),
    peerTableModel(0),
    banTableModel(0),
    pollTimer(0),
    lockedElysiumStateChanged(false),
    lockedElysiumBalanceChanged(false)
{
    cachedBestHeaderHeight = -1;
    cachedBestHeaderTime = -1;
    peerTableModel = new PeerTableModel(this);
    banTableModel = new BanTableModel(this);
    pollTimer = new QTimer(this);
    connect(pollTimer, SIGNAL(timeout()), this, SLOT(updateTimer()));
    pollTimer->start(MODEL_UPDATE_DELAY);

    subscribeToCoreSignals();
}

ClientModel::~ClientModel()
{
    unsubscribeFromCoreSignals();
}

int ClientModel::getNumConnections(unsigned int flags) const
{
    CConnman::NumConnections connections = CConnman::CONNECTIONS_NONE;

    if(flags == CONNECTIONS_IN)
        connections = CConnman::CONNECTIONS_IN;
    else if (flags == CONNECTIONS_OUT)
        connections = CConnman::CONNECTIONS_OUT;
    else if (flags == CONNECTIONS_ALL)
        connections = CConnman::CONNECTIONS_ALL;

    if(g_connman)
         return g_connman->GetNodeCount(connections);
    return 0;
}

void ClientModel::setMasternodeList(const CDeterministicMNList& mnList)
{
    LOCK(cs_mnlinst);
    if (mnListCached.GetBlockHash() == mnList.GetBlockHash()) {
        return;
    }
    mnListCached = mnList;
    Q_EMIT masternodeListChanged();
}

CDeterministicMNList ClientModel::getMasternodeList() const
{
    LOCK(cs_mnlinst);
    return mnListCached;
}

void ClientModel::refreshMasternodeList()
{
    LOCK(cs_mnlinst);
    setMasternodeList(deterministicMNManager->GetListAtChainTip());
}

int ClientModel::getNumBlocks() const
{
    LOCK(cs_main);
    return chainActive.Height();
}

int ClientModel::getHeaderTipHeight() const
{
    if (cachedBestHeaderHeight == -1) {
        // make sure we initially populate the cache via a cs_main lock
        // otherwise we need to wait for a tip update
        LOCK(cs_main);
        if (pindexBestHeader) {
            cachedBestHeaderHeight = pindexBestHeader->nHeight;
            cachedBestHeaderTime = pindexBestHeader->GetBlockTime();
        }
    }
    return cachedBestHeaderHeight;
}

int64_t ClientModel::getHeaderTipTime() const
{
    if (cachedBestHeaderTime == -1) {
        LOCK(cs_main);
        if (pindexBestHeader) {
            cachedBestHeaderHeight = pindexBestHeader->nHeight;
            cachedBestHeaderTime = pindexBestHeader->GetBlockTime();
        }
    }
    return cachedBestHeaderTime;
}

quint64 ClientModel::getTotalBytesRecv() const
{
    if(!g_connman)
        return 0;
    return g_connman->GetTotalBytesRecv();
}

quint64 ClientModel::getTotalBytesSent() const
{
    if(!g_connman)
        return 0;
    return g_connman->GetTotalBytesSent();
}

QDateTime ClientModel::getLastBlockDate() const
{
    LOCK(cs_main);

    if (chainActive.Tip())
        return QDateTime::fromTime_t(chainActive.Tip()->GetBlockTime());

    return QDateTime::fromTime_t(Params().GenesisBlock().GetBlockTime()); // Genesis block's time of current network
}

long ClientModel::getMempoolSize() const
{
    return mempool.size();
}

size_t ClientModel::getMempoolDynamicUsage() const
{
    return (mempool.DynamicMemoryUsage() / 2);
}

double ClientModel::getVerificationProgress(const CBlockIndex *tipIn) const
{
    CBlockIndex *tip = const_cast<CBlockIndex *>(tipIn);
    if (!tip)
    {
        LOCK(cs_main);
        tip = chainActive.Tip();
    }
    return Checkpoints::GuessVerificationProgress(tip);
}

void ClientModel::updateTimer()
{
    // no locking required at this point
    // the following calls will acquire the required lock
    Q_EMIT mempoolSizeChanged(getMempoolSize(), getMempoolDynamicUsage());
    Q_EMIT bytesChanged(getTotalBytesRecv(), getTotalBytesSent());
}

void ClientModel::updateNumConnections(int numConnections)
{
    Q_EMIT numConnectionsChanged(numConnections);
}


void ClientModel::invalidateElysiumState()
{
    Q_EMIT reinitElysiumState();
}

void ClientModel::updateElysiumState()
{
    lockedElysiumStateChanged = false;
    Q_EMIT refreshElysiumState();
}

bool ClientModel::tryLockElysiumStateChanged()
{
    // Try to avoid Elysium queuing too many messages for the UI
    if (lockedElysiumStateChanged) {
        return false;
    }

    lockedElysiumStateChanged = true;
    return true;
}

void ClientModel::updateElysiumBalance()
{
    lockedElysiumBalanceChanged = false;
    Q_EMIT refreshElysiumBalance();
}

bool ClientModel::tryLockElysiumBalanceChanged()
{
    // Try to avoid Elysium queuing too many messages for the UI
    if (lockedElysiumBalanceChanged) {
        return false;
    }

    lockedElysiumBalanceChanged = true;
    return true;
}

void ClientModel::updateElysiumPending(bool pending)
{
    Q_EMIT refreshElysiumPending(pending);
}

void ClientModel::updateNetworkActive(bool networkActive)
{
    Q_EMIT networkActiveChanged(networkActive);
}

void ClientModel::updateAlert()
{
    Q_EMIT alertsChanged(getStatusBarWarnings());
}

bool ClientModel::inInitialBlockDownload() const
{
    return IsInitialBlockDownload();
}

enum BlockSource ClientModel::getBlockSource() const
{
    if (fReindex)
        return BLOCK_SOURCE_REINDEX;
    else if (fImporting)
        return BLOCK_SOURCE_DISK;
    else if (getNumConnections() > 0)
        return BLOCK_SOURCE_NETWORK;

    return BLOCK_SOURCE_NONE;
}

void ClientModel::setNetworkActive(bool active)
{
    if (g_connman) {
         g_connman->SetNetworkActive(active);
    }
}

bool ClientModel::getNetworkActive() const
{
    if (g_connman) {
        return g_connman->GetNetworkActive();
    }
    return false;
}

QString ClientModel::getStatusBarWarnings() const
{
    return QString::fromStdString(GetWarnings("gui"));
}

OptionsModel *ClientModel::getOptionsModel()
{
    return optionsModel;
}

PeerTableModel *ClientModel::getPeerTableModel()
{
    return peerTableModel;
}

BanTableModel *ClientModel::getBanTableModel()
{
    return banTableModel;
}

QString ClientModel::formatFullVersion() const
{
    return QString::fromStdString(FormatFullVersion());
}

QString ClientModel::formatSubVersion() const
{
    return QString::fromStdString(strSubVersion);
}

bool ClientModel::isReleaseVersion() const
{
    return CLIENT_VERSION_IS_RELEASE;
}

QString ClientModel::formatClientStartupTime() const
{
    return QDateTime::fromTime_t(nClientStartupTime).toString();
}

QString ClientModel::dataDir() const
{
    return GUIUtil::boostPathToQString(GetDataDir());
}

void ClientModel::updateBanlist()
{
    banTableModel->refresh();
}

// Handlers for core signals
static void ElysiumStateInvalidated(ClientModel *clientmodel)
{
    // This will be triggered if a reorg invalidates the state
    QMetaObject::invokeMethod(clientmodel, "invalidateElysiumState", Qt::QueuedConnection);
}

static void ElysiumStateChanged(ClientModel *clientmodel)
{
    // This will be triggered for each block that contains Elysium layer transactions
    if (clientmodel->tryLockElysiumStateChanged()) {
        QMetaObject::invokeMethod(clientmodel, "updateElysiumState", Qt::QueuedConnection);
    }
}

static void ElysiumBalanceChanged(ClientModel *clientmodel)
{
    // Triggered when a balance for a wallet address changes
    if (clientmodel->tryLockElysiumBalanceChanged()) {
        QMetaObject::invokeMethod(clientmodel, "updateElysiumBalance", Qt::QueuedConnection);
    }
}

static void ElysiumPendingChanged(ClientModel *clientmodel, bool pending)
{
    // Triggered when Elysium pending map adds/removes transactions
    QMetaObject::invokeMethod(clientmodel, "updateElysiumPending", Qt::QueuedConnection, Q_ARG(bool, pending));
}

static void ShowProgress(ClientModel *clientmodel, const std::string &title, int nProgress)
{
    // emits signal "showProgress"
    QMetaObject::invokeMethod(clientmodel, "showProgress", Qt::QueuedConnection,
                              Q_ARG(QString, QString::fromStdString(title)),
                              Q_ARG(int, nProgress));
}

static void NotifyNumConnectionsChanged(ClientModel *clientmodel, int newNumConnections)
{
    // Too noisy: qDebug() << "NotifyNumConnectionsChanged: " + QString::number(newNumConnections);
    QMetaObject::invokeMethod(clientmodel, "updateNumConnections", Qt::QueuedConnection,
                              Q_ARG(int, newNumConnections));
}

static void NotifyMasternodeListChanged(ClientModel *clientmodel, const CDeterministicMNList& newList)
{
    clientmodel->setMasternodeList(newList);
}

static void NotifyAdditionalDataSyncProgressChanged(ClientModel *clientmodel, double nSyncProgress)
{
    QMetaObject::invokeMethod(clientmodel, "additionalDataSyncProgressChanged", Qt::QueuedConnection,
                              Q_ARG(double, nSyncProgress));
}

static void NotifyNetworkActiveChanged(ClientModel *clientmodel, bool networkActive)
{
    QMetaObject::invokeMethod(clientmodel, "updateNetworkActive", Qt::QueuedConnection,
                              Q_ARG(bool, networkActive));
}

static void NotifyAlertChanged(ClientModel *clientmodel)
{
    qDebug() << "NotifyAlertChanged";
    QMetaObject::invokeMethod(clientmodel, "updateAlert", Qt::QueuedConnection);
}

static void BannedListChanged(ClientModel *clientmodel)
{
    qDebug() << QString("%1: Requesting update for peer banlist").arg(__func__);
    QMetaObject::invokeMethod(clientmodel, "updateBanlist", Qt::QueuedConnection);
}

static void BlockTipChanged(ClientModel *clientmodel, bool initialSync, const CBlockIndex *pIndex, bool fHeader)
{
    // lock free async UI updates in case we have a new block tip
    // during initial sync, only update the UI if the last update
    // was > 250ms (MODEL_UPDATE_DELAY) ago
    int64_t now = 0;
    if (initialSync)
        now = GetTimeMillis();

    int64_t& nLastUpdateNotification = fHeader ? nLastHeaderTipUpdateNotification : nLastBlockTipUpdateNotification;

    if (fHeader) {
        // cache best headers time and height to reduce future cs_main locks
        clientmodel->cachedBestHeaderHeight = pIndex->nHeight;
        clientmodel->cachedBestHeaderTime = pIndex->GetBlockTime();
    }
    // if we are in-sync, update the UI regardless of last update time
    if (!initialSync || now - nLastUpdateNotification > MODEL_UPDATE_DELAY) {
        //pass a async signal to the UI thread
        QMetaObject::invokeMethod(clientmodel, "numBlocksChanged", Qt::QueuedConnection,
                                  Q_ARG(int, pIndex->nHeight),
                                  Q_ARG(QDateTime, QDateTime::fromTime_t(pIndex->GetBlockTime())),
                                  Q_ARG(double, clientmodel->getVerificationProgress(pIndex)),
                                  Q_ARG(bool, fHeader));
        nLastUpdateNotification = now;
    }
}

void ClientModel::subscribeToCoreSignals()
{
    // Connect signals to client
    uiInterface.ShowProgress.connect(boost::bind(ShowProgress, this, _1, _2));
    uiInterface.NotifyNumConnectionsChanged.connect(boost::bind(NotifyNumConnectionsChanged, this, _1));
    uiInterface.NotifyNetworkActiveChanged.connect(boost::bind(NotifyNetworkActiveChanged, this, _1));
    uiInterface.NotifyAlertChanged.connect(boost::bind(NotifyAlertChanged, this));
    uiInterface.BannedListChanged.connect(boost::bind(BannedListChanged, this));
    uiInterface.NotifyBlockTip.connect(boost::bind(BlockTipChanged, this, _1, _2, false));
    uiInterface.NotifyHeaderTip.connect(boost::bind(BlockTipChanged, this, _1, _2, true));
    uiInterface.NotifyAdditionalDataSyncProgressChanged.connect(boost::bind(NotifyAdditionalDataSyncProgressChanged, this, _1));
    uiInterface.NotifyMasternodeListChanged.connect(boost::bind(NotifyMasternodeListChanged, this, _1));

    // Connect Elysium signals
    uiInterface.ElysiumStateChanged.connect(boost::bind(ElysiumStateChanged, this));
    uiInterface.ElysiumPendingChanged.connect(boost::bind(ElysiumPendingChanged, this, _1));
    uiInterface.ElysiumBalanceChanged.connect(boost::bind(ElysiumBalanceChanged, this));
    uiInterface.ElysiumStateInvalidated.connect(boost::bind(ElysiumStateInvalidated, this));
}

void ClientModel::unsubscribeFromCoreSignals()
{
    // Disconnect signals from client
    uiInterface.ShowProgress.disconnect(boost::bind(ShowProgress, this, _1, _2));
    uiInterface.NotifyNumConnectionsChanged.disconnect(boost::bind(NotifyNumConnectionsChanged, this, _1));
    uiInterface.NotifyNetworkActiveChanged.disconnect(boost::bind(NotifyNetworkActiveChanged, this, _1));
    uiInterface.NotifyAlertChanged.disconnect(boost::bind(NotifyAlertChanged, this));
    uiInterface.BannedListChanged.disconnect(boost::bind(BannedListChanged, this));
    uiInterface.NotifyBlockTip.disconnect(boost::bind(BlockTipChanged, this, _1, _2, false));
    uiInterface.NotifyHeaderTip.disconnect(boost::bind(BlockTipChanged, this, _1, _2, true));
    uiInterface.NotifyAdditionalDataSyncProgressChanged.disconnect(boost::bind(NotifyAdditionalDataSyncProgressChanged, this, _1));
    uiInterface.NotifyMasternodeListChanged.disconnect(boost::bind(NotifyMasternodeListChanged, this, _1));

    // Disconnect Elysium signals
    uiInterface.ElysiumStateChanged.disconnect(boost::bind(ElysiumStateChanged, this));
    uiInterface.ElysiumPendingChanged.disconnect(boost::bind(ElysiumPendingChanged, this, _1));
    uiInterface.ElysiumBalanceChanged.disconnect(boost::bind(ElysiumBalanceChanged, this));
	uiInterface.ElysiumStateInvalidated.disconnect(boost::bind(ElysiumStateInvalidated, this));
}

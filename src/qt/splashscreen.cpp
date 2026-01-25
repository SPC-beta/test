// Copyright (c) 2011-2016 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#if defined(HAVE_CONFIG_H)
#include "config/bitcoin-config.h"
#endif

#include "splashscreen.h"

#include "networkstyle.h"

#include "clientversion.h"
#include "init.h"
#include "util.h"
#include "version.h"

#ifdef ENABLE_WALLET
#include "wallet/wallet.h"
#endif

#include <QApplication>
#include <QCloseEvent>
#include <QScreen>
#include <QPainter>
#include <QRadialGradient>
#include <QTimer>
#include <QtMath>

SplashScreen::SplashScreen(const QPixmap &pixmap, Qt::WindowFlags f) : QSplashScreen(pixmap, f)
{
    float fontFactor            = 1.0;

    // define text to place
    QString titleText       = QString(QApplication::applicationName()).replace(QString("-testnet"), QString(""), Qt::CaseSensitive); // cut of testnet, place it as single object further down
    QString font            = "Arial";

    // load the bitmap for writing some text over it
    QPixmap newPixmap;
    newPixmap     = QPixmap(":/images/splash");

    QPainter pixPaint(&newPixmap);
    pixPaint.setPen(QColor(70,70,70));

    pixPaint.setFont(QFont(font, 9*fontFactor));

    // draw copyright stuff
    pixPaint.setFont(QFont(font, 9*fontFactor));
    pixPaint.end();

    this->setPixmap(newPixmap);
    subscribeToCoreSignals();
}

void SplashScreen::slotFinish(QWidget *mainWin)
{
    Q_UNUSED(mainWin);

    /* If the window is minimized, hide() will be ignored. */
    /* Make sure we de-minimize the splashscreen window before hiding */
    if (isMinimized())
        showNormal();
    hide();
    unsubscribeFromCoreSignals();
    deleteLater(); // No more need for this
}

static void InitMessage(SplashScreen *splash, const std::string &message)
{
    QMetaObject::invokeMethod(splash, "showMessage",
        Qt::QueuedConnection,
        Q_ARG(QString, QString::fromStdString(message)),
        Q_ARG(int, Qt::AlignBottom|Qt::AlignHCenter),
        Q_ARG(QColor, QColor(255,255,255)));
}

static void ShowProgress(SplashScreen *splash, const std::string &title, int nProgress)
{
    InitMessage(splash, title + strprintf("%d", nProgress) + "%");
}

#ifdef ENABLE_WALLET
void SplashScreen::ConnectWallet(CWallet* wallet)
{
    wallet->ShowProgress.connect(boost::bind(ShowProgress, this, _1, _2));
    connectedWallets.push_back(wallet);
}
#endif

void SplashScreen::subscribeToCoreSignals()
{
    // Connect signals to client
    uiInterface.InitMessage.connect(boost::bind(InitMessage, this, _1));
    uiInterface.ShowProgress.connect(boost::bind(ShowProgress, this, _1, _2));
#ifdef ENABLE_WALLET
    uiInterface.LoadWallet.connect(boost::bind(&SplashScreen::ConnectWallet, this, _1));
#endif
}

void SplashScreen::unsubscribeFromCoreSignals()
{
    // Disconnect signals from client
    uiInterface.InitMessage.disconnect(boost::bind(InitMessage, this, _1));
    uiInterface.ShowProgress.disconnect(boost::bind(ShowProgress, this, _1, _2));
#ifdef ENABLE_WALLET
    for (CWallet* const & pwallet : connectedWallets) {
        pwallet->ShowProgress.disconnect(boost::bind(ShowProgress, this, _1, _2));
    }
#endif
}

void SplashScreen::showMessage(const QString &message, int alignment, const QColor &color)
{
    curMessage = message;
    curAlignment = alignment;
    curColor = color;
    update();
}

void SplashScreen::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    painter.drawPixmap(0, 0, this->pixmap);
    painter.setPen(curColor);

    QRect textRect = rect().adjusted(5, 5, -5, -60);
    painter.drawText(textRect, Qt::AlignHCenter | Qt::AlignBottom, curMessage);

    static int rotation = 0;
    rotation = (rotation + 10) % 360;

    painter.setRenderHint(QPainter::Antialiasing);
    int x = width() / 2;
    int y = height() - 120;
    int radius = 5;
    int circleRadius = 18;

    for (int i = 0; i < 12; ++i) {
        int alpha = 255 - i * 20;
        painter.setBrush(QColor(255, 255, 255, alpha));
        painter.setPen(Qt::NoPen);

        qreal angleRad = qDegreesToRadians(qreal(rotation + i * 30));
        int dx = int(x + qCos(angleRad) * circleRadius);
        int dy = int(y + qSin(angleRad) * circleRadius);

        painter.drawEllipse(QPoint(dx, dy), radius, radius);
    }

    QTimer::singleShot(50, this, SLOT(update()));
}

void SplashScreen::closeEvent(QCloseEvent *event)
{
    StartShutdown(); // allows an "emergency" shutdown during startup
    event->ignore();
}

void SplashScreen::mousePressEvent(QMouseEvent* event)
{
    event->ignore();
}

# Bitcoinzero [BZX] (Lelantus) Core 2025

## Bitcoinzero

- Coin Suffix: BZX
- Algorithm: Lyra2Z
- Target Spacing: 150 Seconds
- Retarget: every block
- Confirmation: 6 Blocks
- Maturity: 120 Blocks
- Max Coins: N/A
- Min TX Fee: 0.001 BZX
- Block Size: 4MB

## Net Parameters

- P2P Port=29301
- RPC Port=29201
- Client core= based on Firo 14.14.X
- Client name=bitcoinzero.qt
- Conf file=bitcoinzero.conf

## Installation folder

- Windows: C:\Users\Username\AppData\Roaming\bitcoinzero
- Mac: /Library/Application Support/bitcoinzero
- Unix: /.bitcoinzero

# 1. Debian/Ubuntu Linux Daemon Build Instructions

    install dependencies:
    Build a node or qt:

    if you need a swap memory:
    free
    dd if=/dev/zero of=/var/swap.img bs=2048 count=1048576
    mkswap /var/swap.img
    swapon /var/swap.img
    free


    sudo apt-get update
    sudo apt-get upgrade

    sudo apt-get install ???

    git clone https://github.com/BitcoinZeroOfficial/bitcoinzero

    cd bitcoinzero
    gmake -C depends HOST=x86_64-w64-mingw32  # Append "-j N" for N parallel jobs.
    cmake -B build --toolchain depends/x86_64-w64-mingw32/toolchain.cmake
    cmake --build build     # Append "-j N" for N parallel jobs.
	
	Hosts:
	x86_64-pc-linux-gnu for Linux x86 64 bit
	x86_64-w64-mingw32 for Win64
	x86_64-apple-darwin for macOS

    cd src
    strip bitcoinzerod
    strip bitcoinzero-cli
    or:
    cd src
    cd qt
    strip bitcoinzero-qt

    files are:
    bitcoinzerod
    bitcoinzero-cli
    bitcoinzero-qt
    bitcoinzero.conf

    data folder:
    bitcoinzero

    port 29301
    rpc port 29201



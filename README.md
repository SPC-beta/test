# BZX [BZX] (Lelantus) Core 2024

[![Build Status](https://travis-ci.org/BZXOfficial/BZX.svg?branch=master)](https://travis-ci.org/BZXOfficial/BZX)

## BZX

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
- Client core= based on Firo 14.12.X
- Client name=BZX.qt
- Conf file=BZX.conf

## Installation folder

- Windows: C:\Users\Username\AppData\Roaming\BZX
- Mac: /Library/Application Support/BZX
- Unix: /.BZX

# Debian/Ubuntu Linux Daemon Build Instructions

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

    sudo apt-get install make automake cmake curl g++-multilib libtool binutils-gold bsdmainutils pkg-config python3 patch bison

    git clone https://github.com/BZXOfficial/BZX

    cd BZX
    for vps:
    cd depends
    make -j4   (-j is optional, number of your cores, -j4)
    cd ..
    ./autogen.sh
    CONFIG_SITE=$PWD/depends/x86_64-pc-linux-gnu/share/config.site ./configure --disable-option-checking --prefix=$PWD/depends/x86_64-pc-linux-gnu --disable-dependency-tracking --enable-zmq --with-gui=no --enable-glibc-back-compat --enable-reduce-exports --disable-shared --with-pic --enable-module-recovery --disable-jni
    cd src/bls-signatures/build && make chiabls_la && cd -  (optional, only if bls libs fail to build)
    cd src/bls-signatures/build && make && cd - (optional, only if bls libs fail to build)
    make -j4   (-j is optional, number of your cores, -j4)

    for qt:
    cd depends
    make
    cd ..
    ./autogen.sh
    CONFIG_SITE=$PWD/depends/x86_64-pc-linux-gnu/share/config.site ./configure --disable-option-checking --prefix=$PWD/depends/x86_64-pc-linux-gnu --disable-dependency-tracking --enable-zmq --with-gui --enable-glibc-back-compat --enable-reduce-exports --disable-shared --with-pic --enable-module-recovery --disable-jni
    cd src/bls-signatures/build && make chiabls_la && cd -   (optional, only if bls libs fail to build)
    cd src/bls-signatures/build && make && cd - (optional, only if bls libs fail to build)
    make -j4   (-j is optional, number of your cores, -j4)

    cd src
    strip BZXd
    strip BZX-cli
    or:
    cd src
    cd qt
    strip BZX-qt

    files are:
    BZXd
    BZX-cli
    BZX-qt
    BZX.conf

    data folder:
    BZX

    port 29301
    rpc port 29201


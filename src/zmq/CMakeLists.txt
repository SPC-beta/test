# Copyright (c) 2023-present The BZX Core developers
# Distributed under the MIT software license, see the accompanying
# file COPYING or https://opensource.org/license/mit/.

add_library(BZX_zmq STATIC EXCLUDE_FROM_ALL
  zmqabstractnotifier.cpp
  zmqnotificationinterface.cpp
  zmqpublishnotifier.cpp
)
target_compile_definitions(BZX_zmq
  INTERFACE
    ENABLE_ZMQ=1
  PRIVATE
    $<$<BOOL:${WIN32}>:ZMQ_STATIC>
)
target_link_libraries(BZX_zmq
  PRIVATE
    core_interface
    univalue
    zeromq
    secp256k1
    leveldb
)

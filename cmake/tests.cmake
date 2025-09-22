# Copyright (c) 2023-present The BZX Core developers
# Distributed under the MIT software license, see the accompanying
# file COPYING or https://opensource.org/license/mit/.

if(TARGET BZX-util AND TARGET BZX-tx AND PYTHON_COMMAND)
  add_test(NAME util_test_runner
    COMMAND ${CMAKE_COMMAND} -E env BZXUTIL=$<TARGET_FILE:BZX-util> BZXTX=$<TARGET_FILE:BZX-tx> ${PYTHON_COMMAND} ${PROJECT_BINARY_DIR}/test/util/test_runner.py
  )
endif()
# Copyright (c) 2023-present The BZX Core developers
# Distributed under the MIT software license, see the accompanying
# file COPYING or https://opensource.org/license/mit/.

function(generate_setup_nsi)
  set(abs_top_srcdir ${PROJECT_SOURCE_DIR})
  set(abs_top_builddir ${PROJECT_BINARY_DIR})
  set(CLIENT_URL ${PROJECT_HOMEPAGE_URL})
  set(PACKAGE_URL ${PROJECT_HOMEPAGE_URL})
  set(CLIENT_TARNAME "BZX")
  set(BZX_GUI_NAME "BZX-qt")
  set(BZX_DAEMON_NAME "BZXd")
  set(BZX_CLI_NAME "BZX-cli")
  set(BZX_TX_NAME "BZX-tx")
  set(BZX_WALLET_TOOL_NAME "BZX-wallet")
  set(BZX_TEST_NAME "test_BZX")
  set(EXEEXT ${CMAKE_EXECUTABLE_SUFFIX})
  # Set executable directory based on build system
  # CMake/Guix builds put executables in bin/, autotools in release/
  set(EXECUTABLE_DIR "bin")
  configure_file(${PROJECT_SOURCE_DIR}/share/setup.nsi.in ${PROJECT_BINARY_DIR}/BZX-win64-setup.nsi USE_SOURCE_PERMISSIONS @ONLY)
endfunction()

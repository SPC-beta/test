# Copyright (c) 2023-present The Bitcoin Core developers
# Distributed under the MIT software license, see the accompanying
# file COPYING or https://opensource.org/license/mit/.

include_guard(GLOBAL)

function(setup_split_debug_script)
  if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
    set(OBJCOPY ${CMAKE_OBJCOPY})
    set(STRIP ${CMAKE_STRIP})
    configure_file(
      contrib/devtools/split-debug.sh.in split-debug.sh
      FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE
                       GROUP_READ GROUP_EXECUTE
                       WORLD_READ
      @ONLY
    )
  endif()
endfunction()

function(add_maintenance_targets)
  if(NOT PYTHON_COMMAND)
    return()
  endif()

  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(exe_format MACHO)
  elseif(WIN32)
    set(exe_format PE)
  else()
    set(exe_format ELF)
  endif()

  # In CMake, the components of the compiler invocation are separated into distinct variables:
  #  - CMAKE_CXX_COMPILER: the full path to the compiler binary itself (e.g., /usr/bin/clang++).
  #  - CMAKE_CXX_COMPILER_ARG1: a string containing initial compiler options (e.g., --target=x86_64-apple-darwin -nostdlibinc).
  # By concatenating these variables, we form the complete command line to be passed to a Python script via the CXX environment variable.
  string(STRIP "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}" cxx_compiler_command)
  add_custom_target(test-security-check
    COMMAND ${CMAKE_COMMAND} -E env CXX=${cxx_compiler_command} CXXFLAGS=${CMAKE_CXX_FLAGS} LDFLAGS=${CMAKE_EXE_LINKER_FLAGS} ${PYTHON_COMMAND} ${PROJECT_SOURCE_DIR}/contrib/devtools/test-security-check.py TestSecurityChecks.test_${exe_format}
    COMMAND ${CMAKE_COMMAND} -E env CXX=${cxx_compiler_command} CXXFLAGS=${CMAKE_CXX_FLAGS} LDFLAGS=${CMAKE_EXE_LINKER_FLAGS} ${PYTHON_COMMAND} ${PROJECT_SOURCE_DIR}/contrib/devtools/test-symbol-check.py TestSymbolChecks.test_${exe_format}
    VERBATIM
  )

  foreach(target IN ITEMS BZXd BZX-qt BZX-cli BZX-tx BZX-util BZX-wallet test_BZX bench_BZX)
    if(TARGET ${target})
      list(APPEND executables $<TARGET_FILE:${target}>)
    endif()
  endforeach()

  add_custom_target(check-symbols
    COMMAND ${CMAKE_COMMAND} -E echo "Running symbol and dynamic library checks..."
    COMMAND ${PYTHON_COMMAND} ${PROJECT_SOURCE_DIR}/contrib/devtools/symbol-check.py ${executables}
    VERBATIM
  )

  if(ENABLE_HARDENING)
    add_custom_target(check-security
      COMMAND ${CMAKE_COMMAND} -E echo "Checking binary security..."
      COMMAND ${PYTHON_COMMAND} ${PROJECT_SOURCE_DIR}/contrib/devtools/security-check.py ${executables}
      VERBATIM
    )
  else()
    add_custom_target(check-security)
  endif()
endfunction()

function(add_windows_deploy_target)
  if(MINGW AND TARGET BZX-qt AND TARGET BZXd AND TARGET BZX-cli AND TARGET BZX-tx AND TARGET BZX-wallet AND TARGET BZX-util AND TARGET test_BZX)
    # TODO: Consider replacing this code with the CPack NSIS Generator.
    #       See https://cmake.org/cmake/help/latest/cpack_gen/nsis.html
    include(GenerateSetupNsi)
    generate_setup_nsi()
    add_custom_command(
      OUTPUT ${PROJECT_BINARY_DIR}/BZX-win64-setup.exe
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/release
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:BZX-qt> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:BZX-qt>
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:BZXd> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:BZXd>
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:BZX-cli> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:BZX-cli>
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:BZX-tx> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:BZX-tx>
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:BZX-wallet> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:BZX-wallet>
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:BZX-util> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:BZX-util>
      COMMAND ${CMAKE_STRIP} $<TARGET_FILE:test_BZX> -o ${PROJECT_BINARY_DIR}/release/$<TARGET_FILE_NAME:test_BZX>
      COMMAND makensis -V2 ${PROJECT_BINARY_DIR}/BZX-win64-setup.nsi
      VERBATIM
    )
    add_custom_target(deploy DEPENDS ${PROJECT_BINARY_DIR}/BZX-win64-setup.exe)
  endif()
endfunction()

function(add_macos_deploy_target)
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND TARGET BZX-qt)
    set(macos_app "BZX-Qt.app")
    # Populate Contents subdirectory.
    configure_file(${PROJECT_SOURCE_DIR}/share/qt/Info.plist.in ${macos_app}/Contents/Info.plist NO_SOURCE_PERMISSIONS)
    file(CONFIGURE OUTPUT ${macos_app}/Contents/PkgInfo CONTENT "APPL????")
    # Populate Contents/Resources subdirectory.
    file(CONFIGURE OUTPUT ${macos_app}/Contents/Resources/empty.lproj CONTENT "")
    configure_file(${PROJECT_SOURCE_DIR}/src/qt/res/icons/BZX.icns ${macos_app}/Contents/Resources/BZX.icns NO_SOURCE_PERMISSIONS COPYONLY)
    file(CONFIGURE OUTPUT ${macos_app}/Contents/Resources/Base.lproj/InfoPlist.strings
      CONTENT "{ CFBundleDisplayName = \"@CLIENT_NAME@\"; CFBundleName = \"@CLIENT_NAME@\"; }"
    )

    add_custom_command(
      OUTPUT ${PROJECT_BINARY_DIR}/${macos_app}/Contents/MacOS/BZX-Qt
      COMMAND ${CMAKE_COMMAND} --install ${PROJECT_BINARY_DIR} --config $<CONFIG> --component GUI --prefix ${macos_app}/Contents/MacOS --strip
      COMMAND ${CMAKE_COMMAND} -E rename ${macos_app}/Contents/MacOS/bin/$<TARGET_FILE_NAME:BZX-qt> ${macos_app}/Contents/MacOS/BZX-Qt
      COMMAND ${CMAKE_COMMAND} -E rm -rf ${macos_app}/Contents/MacOS/bin
      VERBATIM
    )

    string(REPLACE " " "-" osx_volname ${CLIENT_NAME})
    if(CMAKE_HOST_APPLE)
      add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/${osx_volname}.zip
        COMMAND ${PYTHON_COMMAND} ${PROJECT_SOURCE_DIR}/contrib/macdeploy/macdeployqtplus ${macos_app} ${osx_volname} -translations-dir=${QT_TRANSLATIONS_DIR} -zip
        DEPENDS ${PROJECT_BINARY_DIR}/${macos_app}/Contents/MacOS/BZX-Qt
        VERBATIM
      )
      add_custom_target(deploydir
        DEPENDS ${PROJECT_BINARY_DIR}/${osx_volname}.zip
      )
      add_custom_target(deploy
        DEPENDS ${PROJECT_BINARY_DIR}/${osx_volname}.zip
      )
    else()
      add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/dist/${macos_app}/Contents/MacOS/BZX-Qt
        COMMAND OBJDUMP=${CMAKE_OBJDUMP} ${PYTHON_COMMAND} ${PROJECT_SOURCE_DIR}/contrib/macdeploy/macdeployqtplus ${macos_app} ${osx_volname} -translations-dir=${QT_TRANSLATIONS_DIR}
        DEPENDS ${PROJECT_BINARY_DIR}/${macos_app}/Contents/MacOS/BZX-Qt
        VERBATIM
      )
      add_custom_target(deploydir
        DEPENDS ${PROJECT_BINARY_DIR}/dist/${macos_app}/Contents/MacOS/BZX-Qt
      )

      find_program(ZIP_COMMAND zip REQUIRED)
      add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/dist/${osx_volname}.zip
        WORKING_DIRECTORY dist
        COMMAND ${PROJECT_SOURCE_DIR}/cmake/script/macos_zip.sh ${ZIP_COMMAND} ${osx_volname}.zip
        VERBATIM
      )
      add_custom_target(deploy
        DEPENDS ${PROJECT_BINARY_DIR}/dist/${osx_volname}.zip
      )
    endif()
    add_dependencies(deploydir BZX-qt)
    add_dependencies(deploy deploydir)
  endif()
endfunction()

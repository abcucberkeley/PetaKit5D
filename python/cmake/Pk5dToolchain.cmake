# The toolchain of a source build on Windows (docs/platforms.md, "Windows: the toolchain of a source
# build"; docs/decisions/0011-windows-toolchain.md): MinGW-builds' GCC 12.2.0 with the UCRT, the
# toolchain of the release wheels and of the cpp-zarr and cpp-tiff wheels, where it is found, else
# MSVC. CMakeLists.txt calls pk5d_windows_toolchain() before project(), on Windows. The exact
# add-on, whose source tree stands alone, has a copy (cpp/fft/exact/cmake/Pk5dToolchain.cmake),
# which tests/release/test_windows_toolchain.py keeps the same.
#
# pyproject.toml's overrides set PK5D_WINDOWS_TOOLCHAIN for pip's builds on Windows when the
# environment names no toolchain (CC, CXX, CMAKE_GENERATOR, CMAKE_ARGS, CMAKE_TOOLCHAIN_FILE):
#   default   the first attempt, with Ninja: the toolchain if the first gcc on PATH is it, else if
#             tools/get_mingw.sh's folder has it; an error without it, after which scikit-build-core
#             builds again with its default there, Visual Studio's generator and MSVC (if.failed)
#   fallback  that second attempt
#   mingw     PETAKIT5D_WINDOWS_TOOLCHAIN=mingw: as default, but no second attempt
# Unset (a CMake build tree, a toolchain the environment names, PETAKIT5D_WINDOWS_TOOLCHAIN=msvc):
# the compilers the generator, the environment or the command line give. A compiler, toolchain
# file or generator on CMake's command line (-C cmake.define.CMAKE_CXX_COMPILER=...) wins over
# default and mingw too. PK5D_TOOLCHAIN_WHY (cache) holds the reason, which Pk5dExtension.cmake
# reports with the compiler.

# The folder tools/get_mingw.sh unpacks the toolchain into without a folder argument
# (%LOCALAPPDATA%\petakit5d\mingw-12.2.0, its compilers in mingw64/bin), in `out`; empty without
# LOCALAPPDATA
function(pk5d_get_mingw_dir out)
  set(_d "")
  if(NOT "$ENV{LOCALAPPDATA}" STREQUAL "")
    file(TO_CMAKE_PATH "$ENV{LOCALAPPDATA}/petakit5d/mingw-12.2.0" _d)
  endif()
  set(${out} "${_d}" PARENT_SCOPE)
endfunction()

# Whether the compiler <gcc> is MinGW-builds' GCC 12.2.0 with the UCRT, in <ok_var> (TRUE or FALSE),
# and what it is otherwise, in <why_var>:
# - the version: its banner, that of MinGW-builds' x86_64 build of 12.2.0 with posix threads and SEH
#   (rev2, as tools/get_mingw.sh fetches it): "(x86_64-posix-seh-rev2, Built by MinGW-W64 project)
#   12.2.0";
# - the runtime: the C library its headers build for, the UCRT (__MSVCRT_VERSION__ 0xE00), not
#   msvcrt.dll (0x700), which MinGW-builds' msvcrt builds of 12.2.0 take under the same banner.
function(pk5d_mingw_toolchain_check gcc ok_var why_var)
  set(${ok_var} FALSE PARENT_SCOPE)
  execute_process(COMMAND "${gcc}" --version RESULT_VARIABLE _rc OUTPUT_VARIABLE _out
                  ERROR_QUIET)
  string(REGEX MATCH "^[^\r\n]+" _banner "${_out}")
  if(NOT _rc EQUAL 0 OR _banner STREQUAL "")
    set(${why_var} "${gcc} does not run" PARENT_SCOPE)
    return()
  endif()
  if(NOT _banner MATCHES "\\(x86_64-posix-seh-rev2, Built by MinGW-W64 project\\) 12\\.2\\.0$")
    set(${why_var} "${gcc} is ${_banner}" PARENT_SCOPE)
    return()
  endif()
  set(_probe "${CMAKE_BINARY_DIR}/CMakeFiles/pk5d_crt_probe.c")
  file(WRITE "${_probe}" "#include <stdio.h>\n")
  execute_process(COMMAND "${gcc}" -x c -dM -E "${_probe}" RESULT_VARIABLE _rc
                  OUTPUT_VARIABLE _macros ERROR_VARIABLE _err)
  if(NOT _rc EQUAL 0)
    # (as when its programs can't load their DLLs)
    string(REGEX MATCH "^[^\r\n]*" _err "${_err}")
    set(${why_var} "${gcc} (12.2.0) fails to preprocess a C file (${_rc}: ${_err})" PARENT_SCOPE)
    return()
  endif()
  set(_crt "")
  if(_macros MATCHES "#define __MSVCRT_VERSION__ (0[xX][0-9A-Fa-f]+|[0-9]+)")
    math(EXPR _crt "${CMAKE_MATCH_1}")
  endif()
  if(_crt STREQUAL "")
    set(${why_var} "${gcc} (12.2.0) names no C library (no __MSVCRT_VERSION__)" PARENT_SCOPE)
  elseif(_crt LESS 3584)
    set(${why_var} "${gcc} is MinGW-builds' 12.2.0 for msvcrt.dll, not the UCRT" PARENT_SCOPE)
  else()
    set(${ok_var} TRUE PARENT_SCOPE)
    set(${why_var} "MinGW-builds' GCC 12.2.0 with the UCRT" PARENT_SCOPE)
  endif()
endfunction()

# Before project(): the compilers of a source build on Windows, as above
function(pk5d_windows_toolchain)
  set(_mode "${PK5D_WINDOWS_TOOLCHAIN}")
  set(_switch "$ENV{PETAKIT5D_WINDOWS_TOOLCHAIN}")
  if(NOT _switch MATCHES "^(mingw|msvc|)$")
    message(FATAL_ERROR "PETAKIT5D_WINDOWS_TOOLCHAIN is '${_switch}': mingw, msvc or unset "
                        "(docs/platforms.md)")
  endif()
  if(NOT _mode MATCHES "^(default|fallback|mingw|)$")
    message(FATAL_ERROR "PK5D_WINDOWS_TOOLCHAIN is '${_mode}': default, fallback, mingw or unset")
  endif()
  # a configure again of a build tree: the first one chose
  if(DEFINED CACHE{PK5D_TOOLCHAIN_WHY} AND DEFINED CMAKE_CXX_COMPILER)
    return()
  endif()
  set(_why "")
  if(_mode STREQUAL "fallback")
    set(_why "the fallback: the attempt with MinGW-w64 failed (its log is above)")
  elseif(_mode STREQUAL "" AND _switch STREQUAL "msvc")
    set(_why "PETAKIT5D_WINDOWS_TOOLCHAIN=msvc")
  elseif(_mode STREQUAL "")
    set(_why "as CC, CXX, CMAKE_GENERATOR, CMAKE_ARGS, CMAKE_TOOLCHAIN_FILE or CMake's command "
             "line name it")
  elseif(DEFINED CMAKE_C_COMPILER OR DEFINED CMAKE_CXX_COMPILER OR CMAKE_TOOLCHAIN_FILE
         OR NOT "$ENV{CMAKE_TOOLCHAIN_FILE}" STREQUAL "" OR NOT CMAKE_GENERATOR STREQUAL "Ninja")
    set(_why "as CMake's command line names it")
  else()
    # default or mingw: the first gcc on PATH, then the one in tools/get_mingw.sh's folder
    set(_gccs "")
    set(_wheres "")
    set(_tried "")
    unset(_pk5d_tc_path_gcc)
    find_program(_pk5d_tc_path_gcc NAMES gcc NO_CACHE NO_PACKAGE_ROOT_PATH NO_CMAKE_PATH
                 NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_INSTALL_PREFIX)
    if(_pk5d_tc_path_gcc)
      list(APPEND _gccs "${_pk5d_tc_path_gcc}")
      list(APPEND _wheres "on PATH")
    else()
      list(APPEND _tried "no gcc on PATH")
    endif()
    pk5d_get_mingw_dir(_dir)
    if(_dir)
      unset(_pk5d_tc_dir_gcc)
      find_program(_pk5d_tc_dir_gcc NAMES gcc PATHS "${_dir}/mingw64/bin" NO_DEFAULT_PATH NO_CACHE)
      if(_pk5d_tc_dir_gcc AND NOT _pk5d_tc_dir_gcc STREQUAL _pk5d_tc_path_gcc)
        list(APPEND _gccs "${_pk5d_tc_dir_gcc}")
        list(APPEND _wheres "in tools/get_mingw.sh's folder")
      elseif(NOT _pk5d_tc_dir_gcc)
        list(APPEND _tried "none in tools/get_mingw.sh's folder ${_dir}")
      endif()
    endif()
    foreach(_gcc _where IN ZIP_LISTS _gccs _wheres)
      pk5d_mingw_toolchain_check("${_gcc}" _ok _what)
      if(_ok)
        get_filename_component(_bin "${_gcc}" DIRECTORY)
        unset(_pk5d_tc_gxx)
        find_program(_pk5d_tc_gxx NAMES g++ PATHS "${_bin}" NO_DEFAULT_PATH NO_CACHE)
        if(NOT _pk5d_tc_gxx)
          set(_ok FALSE)
          set(_what "${_gcc} has no g++ next to it")
        endif()
      endif()
      if(_ok)
        set(CMAKE_C_COMPILER "${_gcc}" CACHE FILEPATH "C compiler")
        set(CMAKE_CXX_COMPILER "${_pk5d_tc_gxx}" CACHE FILEPATH "C++ compiler")
        if(_mode STREQUAL "mingw")
          set(_why "${_what} ${_where}: PETAKIT5D_WINDOWS_TOOLCHAIN=mingw")
        else()
          set(_why "${_what} ${_where}: the default of source builds")
        endif()
        break()
      endif()
      list(APPEND _tried "gcc ${_where}: ${_what}")
    endforeach()
    if(_why STREQUAL "")
      list(JOIN _tried "; " _tried)
      if(_mode STREQUAL "mingw")
        message(FATAL_ERROR "petakit5d: PETAKIT5D_WINDOWS_TOOLCHAIN=mingw, but no MinGW-builds GCC "
                            "12.2.0 with the UCRT: ${_tried} (bash tools/get_mingw.sh fetches it)")
      endif()
      message(STATUS "petakit5d: no MinGW-builds GCC 12.2.0 with the UCRT, the default toolchain "
                     "of source builds on Windows: ${_tried} (bash tools/get_mingw.sh fetches it)")
      message(FATAL_ERROR "petakit5d: MinGW-w64 not found (above): scikit-build-core builds again "
                          "with MSVC")
    endif()
  endif()
  set(PK5D_TOOLCHAIN_WHY "${_why}" CACHE INTERNAL "Windows: why the build takes its compilers")
endfunction()

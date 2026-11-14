# OpenMP for the kernels: the interface target pk5d_openmp, and
#   PK5D_OPENMP_FOUND    TRUE when the kernels are built with OpenMP
#   PK5D_OPENMP_DESC     what was used (reported by petakit5d._ext.core.build_info())
#   PK5D_OPENMP_RUNTIME  a runtime DLL to put next to the modules (Windows), or empty
#   PK5D_OPENMP_RPATH    the runtime's folder when it is not a system folder (Linux, macOS), or empty
#   PK5D_OPENMP_SHARED   the SONAME of the I/O packages' libgomp when the modules share it (Linux),
#                        "off" when PK5D_OPENMP_SHARE is OFF, else empty
#
# Linux and GCC anywhere: -fopenmp, libgomp (FindOpenMP); on Linux the cpp-zarr/cpp-tiff wheels'
#   libgomp when they are installed (PK5D_OPENMP_SHARE, below).
# macOS with Apple Clang: -Xpreprocessor -fopenmp and libomp (Homebrew's, or OpenMP_ROOT).
# Windows with MSVC: /openmp:llvm (OpenMP 3.1 and more, as the kernels use: collapse,
#   min/max reductions, unsigned loop indices), runtime libomp140.<arch>.dll.
#
# One OpenMP runtime per process on Linux (docs/decisions/0008-openmp-runtime.md). The cpp-zarr
# and cpp-tiff wheels bundle a libgomp that auditwheel renamed (libgomp-<hash>.so.1.0.0, also its
# SONAME). A process that loads it and modules linked against the toolchain's libgomp.so.1 has two
# runtimes, two pools of threads that spin on the cores after their regions while the other's
# run. The loader maps a library once per SONAME, so modules linked against that same file, with
# a copy of it next to them (found through $ORIGIN), use the wheels' runtime whichever is loaded
# first.
#   AUTO (default)  the libgomp that the build's Python loads for cppzarr (else cpptiff), when it
#                   has a SONAME of its own and links the kernels' constructs
#                   (pk5d_openmp_probe.cpp); else the toolchain's runtime
#   OFF             the toolchain's runtime
#   <file>          that libgomp (an error if it can't be used)
set(PK5D_OPENMP_SHARE "AUTO" CACHE STRING
    "Linux: the cpp-zarr/cpp-tiff wheels' libgomp (AUTO), the toolchain's (OFF), or a libgomp file")

add_library(pk5d_openmp INTERFACE)
set(PK5D_OPENMP_FOUND FALSE)
set(PK5D_OPENMP_DESC "none")
set(PK5D_OPENMP_RUNTIME "")
set(PK5D_OPENMP_RPATH "")
set(PK5D_OPENMP_SHARED "")
set(PK5D_READELF "")
set(_pk5d_openmp_dir "${CMAKE_CURRENT_LIST_DIR}")

# The libgomp to share (PK5D_OPENMP_SHARE), copied to <PK5D_EXT_OUTPUT_DIR>/<its SONAME>: sets
# PK5D_OPENMP_SHARED to the SONAME ("off", or empty for the toolchain's runtime), PK5D_READELF
# and _pk5d_shared_from (where the file came from) in the caller.
function(pk5d_find_shared_openmp toolchain_lib other_libs)
  string(TOUPPER "${PK5D_OPENMP_SHARE}" _mode)
  if(_mode MATCHES "^(OFF|NO|FALSE|N|0|)$")
    set(PK5D_OPENMP_SHARED "off" PARENT_SCOPE)
    return()
  endif()
  set(_auto FALSE)
  if(_mode STREQUAL "AUTO")
    set(_auto TRUE)
  endif()
  macro(_pk5d_not_shared why)
    if(NOT _auto)
      message(FATAL_ERROR "PK5D_OPENMP_SHARE=${PK5D_OPENMP_SHARE}: ${why}")
    endif()
    message(STATUS "petakit5d OpenMP: the toolchain's runtime (${why})")
    return()
  endmacro()

  if(_auto)
    execute_process(COMMAND "${Python_EXECUTABLE}" "${_pk5d_openmp_dir}/pk5d_io_openmp.py"
                    OUTPUT_VARIABLE _out RESULT_VARIABLE _rc ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _rc EQUAL 0 OR NOT _out)
      _pk5d_not_shared("no libgomp of cpp-zarr or cpp-tiff found for ${Python_EXECUTABLE}")
    endif()
    string(REPLACE "\n" ";" _lines "${_out}")
    list(GET _lines 0 _first)
    string(REPLACE "\t" ";" _first "${_first}")
    list(GET _first 0 _from)
    list(GET _first 1 _file)
    list(LENGTH _lines _n)
    if(_n GREATER 1)
      message(WARNING "petakit5d OpenMP: cpp-zarr and cpp-tiff load different libgomp files "
                      "(${_lines}); the modules share ${_from}'s")
    endif()
  else()
    set(_from "PK5D_OPENMP_SHARE")
    set(_file "${PK5D_OPENMP_SHARE}")
    if(NOT EXISTS "${_file}")
      _pk5d_not_shared("no file ${_file}")
    endif()
  endif()
  get_filename_component(_real "${_file}" REALPATH)
  get_filename_component(_own "${toolchain_lib}" REALPATH)
  if(_real STREQUAL _own)
    _pk5d_not_shared("${_from} loads it too")
  endif()

  set(_readelf "${CMAKE_READELF}")
  if(NOT _readelf)
    find_program(PK5D_READELF_EXE NAMES readelf llvm-readelf eu-readelf)
    set(_readelf "${PK5D_READELF_EXE}")
  endif()
  if(NOT _readelf)
    _pk5d_not_shared("no readelf to read ${_file}")
  endif()
  execute_process(COMMAND "${_readelf}" -d -W "${_file}" OUTPUT_VARIABLE _dyn ERROR_QUIET)
  if(NOT _dyn MATCHES "\\(SONAME\\)[^\n]*\\[([^]\n]+)\\]")
    _pk5d_not_shared("${_file} has no SONAME")
  endif()
  set(_soname "${CMAKE_MATCH_1}")
  if(_soname STREQUAL "libgomp.so.1")
    # the toolchain runtime's own name: the loader maps one of the two for both already
    _pk5d_not_shared("${_from}'s libgomp is named libgomp.so.1")
  endif()

  # the runtime must define every entry point the compiler emits for the kernels' constructs
  separate_arguments(_cxxflags NATIVE_COMMAND "${OpenMP_CXX_FLAGS}")
  try_compile(_ok SOURCES "${_pk5d_openmp_dir}/pk5d_openmp_probe.cpp"
              COMPILE_DEFINITIONS ${_cxxflags}
              LINK_LIBRARIES "${_file}" ${other_libs}
              CXX_STANDARD 17
              OUTPUT_VARIABLE _log NO_CACHE)
  if(NOT _ok)
    string(REGEX MATCHALL "undefined (reference to|symbol:) [^\n]*" _undef "${_log}")
    _pk5d_not_shared("${_file} can't link the kernels' OpenMP constructs: ${_undef}")
  endif()

  file(MAKE_DIRECTORY "${PK5D_EXT_OUTPUT_DIR}")
  file(COPY_FILE "${_file}" "${PK5D_EXT_OUTPUT_DIR}/${_soname}" ONLY_IF_DIFFERENT)
  message(STATUS "petakit5d OpenMP: the modules share ${_from}'s runtime ${_soname} (${_file})")
  set(PK5D_OPENMP_SHARED "${_soname}" PARENT_SCOPE)
  set(PK5D_READELF "${_readelf}" PARENT_SCOPE)
  set(_pk5d_shared_from "${_from}" PARENT_SCOPE)
endfunction()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.28)
    message(WARNING "MSVC ${CMAKE_CXX_COMPILER_VERSION} has no /openmp:llvm (needs VS 2019 16.9 or newer)")
  else()
    target_compile_options(pk5d_openmp INTERFACE /openmp:llvm)
    set(PK5D_OPENMP_FOUND TRUE)
    set(PK5D_OPENMP_DESC "MSVC /openmp:llvm")
    # the runtime DLL ships with the compiler, not with Windows: find it in the toolset's
    # redistributable folders and put it next to the modules
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(ARM64|arm64|aarch64)$")
      set(_pk5d_arch arm64)
    else()
      set(_pk5d_arch x64)
    endif()
    get_filename_component(_pk5d_cl_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_pk5d_vc_dir "${_pk5d_cl_dir}/../../../../../.." ABSOLUTE)
    set(_pk5d_omp_globs
        "$ENV{VCToolsRedistDir}/${_pk5d_arch}/Microsoft.VC*.OpenMP.LLVM/libomp140.*.dll"
        "${_pk5d_vc_dir}/Redist/MSVC/*/${_pk5d_arch}/Microsoft.VC*.OpenMP.LLVM/libomp140.*.dll"
        "${_pk5d_cl_dir}/libomp140.*.dll")
    set(_pk5d_omp_dlls "")
    foreach(_g IN LISTS _pk5d_omp_globs)
      file(GLOB _found "${_g}")
      list(APPEND _pk5d_omp_dlls ${_found})
    endforeach()
    list(FILTER _pk5d_omp_dlls EXCLUDE REGEX "[Dd]ebug")
    if(NOT _pk5d_omp_dlls)
      # other layouts (Visual Studio 2026): anywhere under the redistributable and tool folders
      foreach(_d "${_pk5d_vc_dir}/Redist/MSVC" "${_pk5d_vc_dir}/Tools/MSVC")
        file(GLOB_RECURSE _found "${_d}/libomp140.*.dll")
        list(APPEND _pk5d_omp_dlls ${_found})
      endforeach()
      list(FILTER _pk5d_omp_dlls EXCLUDE REGEX "[Dd]ebug|libomp140d")
      if(_pk5d_arch STREQUAL "x64")
        list(FILTER _pk5d_omp_dlls INCLUDE REGEX "x86_64\\.dll$")
      else()
        list(FILTER _pk5d_omp_dlls INCLUDE REGEX "aarch64\\.dll$|arm64\\.dll$")
      endif()
    endif()
    if(_pk5d_omp_dlls)
      list(SORT _pk5d_omp_dlls COMPARE NATURAL ORDER DESCENDING)
      list(GET _pk5d_omp_dlls 0 PK5D_OPENMP_RUNTIME)
    else()
      message(WARNING "libomp140 DLL not found under ${_pk5d_vc_dir}/Redist; "
                      "the modules need it on the DLL search path")
    endif()
  endif()
elseif(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(_pk5d_omp_hints "${OpenMP_ROOT}" "$ENV{OpenMP_ROOT}")
  find_program(PK5D_BREW brew)
  if(PK5D_BREW)
    execute_process(COMMAND "${PK5D_BREW}" --prefix libomp OUTPUT_VARIABLE _pk5d_brew_omp
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    list(APPEND _pk5d_omp_hints "${_pk5d_brew_omp}")
  endif()
  list(APPEND _pk5d_omp_hints /opt/homebrew/opt/libomp /usr/local/opt/libomp)
  find_path(PK5D_LIBOMP_INCLUDE_DIR omp.h HINTS ${_pk5d_omp_hints} PATH_SUFFIXES include)
  find_library(PK5D_LIBOMP_LIBRARY omp HINTS ${_pk5d_omp_hints} PATH_SUFFIXES lib)
  if(PK5D_LIBOMP_INCLUDE_DIR AND PK5D_LIBOMP_LIBRARY)
    target_compile_options(pk5d_openmp INTERFACE "SHELL:-Xpreprocessor -fopenmp")
    target_include_directories(pk5d_openmp INTERFACE "${PK5D_LIBOMP_INCLUDE_DIR}")
    target_link_libraries(pk5d_openmp INTERFACE "${PK5D_LIBOMP_LIBRARY}")
    set(PK5D_OPENMP_FOUND TRUE)
    set(PK5D_OPENMP_DESC "AppleClang -Xpreprocessor -fopenmp, ${PK5D_LIBOMP_LIBRARY}")
    # a libomp named @rpath/libomp.dylib (conda's, one built from LLVM's sources) is found
    # through its folder on the modules' RPATH; Homebrew's has an absolute name
    get_filename_component(PK5D_OPENMP_RPATH "${PK5D_LIBOMP_LIBRARY}" DIRECTORY)
  endif()
else()
  find_package(OpenMP COMPONENTS C CXX)
  if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
    set(PK5D_OPENMP_FOUND TRUE)
    set(PK5D_OPENMP_DESC "${CMAKE_CXX_COMPILER_ID} ${OpenMP_CXX_FLAGS} (OpenMP ${OpenMP_CXX_VERSION})")
    set(_pk5d_omp_lib "")
    set(_pk5d_omp_other_libs "")
    foreach(_lib IN LISTS OpenMP_CXX_LIBRARIES)
      if(NOT _pk5d_omp_lib AND _lib MATCHES "(gomp|omp|iomp5)[^/]*$" AND EXISTS "${_lib}")
        set(_pk5d_omp_lib "${_lib}")
      else()
        list(APPEND _pk5d_omp_other_libs "${_lib}")
      endif()
    endforeach()
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND _pk5d_omp_lib MATCHES "gomp[^/]*$")
      pk5d_find_shared_openmp("${_pk5d_omp_lib}" "${_pk5d_omp_other_libs}")
    endif()
    if(PK5D_OPENMP_SHARED AND NOT PK5D_OPENMP_SHARED STREQUAL "off")
      # the compile flags of FindOpenMP's targets, the shared runtime's copy instead of its library
      foreach(_lang C CXX)
        separate_arguments(_flags NATIVE_COMMAND "${OpenMP_${_lang}_FLAGS}")
        foreach(_f IN LISTS _flags)
          target_compile_options(pk5d_openmp INTERFACE "$<$<COMPILE_LANGUAGE:${_lang}>:${_f}>")
        endforeach()
        target_include_directories(pk5d_openmp INTERFACE ${OpenMP_${_lang}_INCLUDE_DIRS})
      endforeach()
      add_library(pk5d_shared_libgomp SHARED IMPORTED)
      set_target_properties(pk5d_shared_libgomp PROPERTIES
        IMPORTED_LOCATION "${PK5D_EXT_OUTPUT_DIR}/${PK5D_OPENMP_SHARED}"
        IMPORTED_SONAME "${PK5D_OPENMP_SHARED}")
      target_link_libraries(pk5d_openmp INTERFACE pk5d_shared_libgomp ${_pk5d_omp_other_libs})
      install(FILES "${PK5D_EXT_OUTPUT_DIR}/${PK5D_OPENMP_SHARED}" DESTINATION petakit5d/_ext)
      set(PK5D_OPENMP_RPATH "$ORIGIN")
      string(APPEND PK5D_OPENMP_DESC ", runtime ${PK5D_OPENMP_SHARED} (${_pk5d_shared_from}'s)")
    else()
      target_link_libraries(pk5d_openmp INTERFACE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
      # a compiler newer than the system's needs its own runtime (newer GOMP symbol versions):
      # its folder goes on the modules' RPATH, also when installed
      if(_pk5d_omp_lib)
        get_filename_component(_dir "${_pk5d_omp_lib}" REALPATH)
        get_filename_component(_dir "${_dir}" DIRECTORY)
        if(NOT _dir MATCHES "^/(usr/)?lib(32|64)?(/[a-z0-9_]+-linux-gnu)?$")
          set(PK5D_OPENMP_RPATH "${_dir}")
        endif()
      endif()
    endif()
  endif()
endif()

if(NOT PK5D_OPENMP_FOUND)
  set(_pk5d_omp_hint "")
  if(APPLE)
    set(_pk5d_omp_hint " On macOS: brew install libomp (or set OpenMP_ROOT).")
  endif()
  if(PK5D_REQUIRE_OPENMP)
    message(FATAL_ERROR "OpenMP not found for ${CMAKE_CXX_COMPILER_ID}.${_pk5d_omp_hint} "
                        "-DPK5D_REQUIRE_OPENMP=OFF builds serial kernels.")
  endif()
  message(WARNING "OpenMP not found: the kernels are built serial.${_pk5d_omp_hint}")
endif()
message(STATUS "petakit5d OpenMP: ${PK5D_OPENMP_DESC}")

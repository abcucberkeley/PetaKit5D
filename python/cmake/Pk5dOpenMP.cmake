# OpenMP for the kernels: the interface target pk5d_openmp, and
#   PK5D_OPENMP_FOUND    TRUE when the kernels are built with OpenMP
#   PK5D_OPENMP_DESC     what was used (reported by petakit5d._ext.core.build_info())
#   PK5D_OPENMP_RUNTIME  a runtime DLL to put next to the modules (Windows), or empty
#   PK5D_OPENMP_RPATH    the runtime's folder when it is not a system folder (Linux), or empty
#
# Linux and GCC anywhere: -fopenmp, libgomp (FindOpenMP).
# macOS with Apple Clang: -Xpreprocessor -fopenmp and libomp (Homebrew's, or OpenMP_ROOT).
# Windows with MSVC: /openmp:llvm (OpenMP 3.1 and more, as the kernels use: collapse,
#   min/max reductions, unsigned loop indices), runtime libomp140.<arch>.dll.

add_library(pk5d_openmp INTERFACE)
set(PK5D_OPENMP_FOUND FALSE)
set(PK5D_OPENMP_DESC "none")
set(PK5D_OPENMP_RUNTIME "")
set(PK5D_OPENMP_RPATH "")

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
  endif()
else()
  find_package(OpenMP COMPONENTS C CXX)
  if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
    target_link_libraries(pk5d_openmp INTERFACE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
    set(PK5D_OPENMP_FOUND TRUE)
    set(PK5D_OPENMP_DESC "${CMAKE_CXX_COMPILER_ID} ${OpenMP_CXX_FLAGS} (OpenMP ${OpenMP_CXX_VERSION})")
    # a compiler newer than the system's needs its own runtime (newer GOMP symbol versions):
    # its folder goes on the modules' RPATH, also when installed
    foreach(_lib IN LISTS OpenMP_CXX_LIBRARIES)
      if(_lib MATCHES "(gomp|omp|iomp5)[^/]*$" AND EXISTS "${_lib}")
        get_filename_component(_dir "${_lib}" REALPATH)
        get_filename_component(_dir "${_dir}" DIRECTORY)
        if(NOT _dir MATCHES "^/(usr/)?lib(32|64)?(/[a-z0-9_]+-linux-gnu)?$")
          set(PK5D_OPENMP_RPATH "${_dir}")
        endif()
        break()
      endif()
    endforeach()
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

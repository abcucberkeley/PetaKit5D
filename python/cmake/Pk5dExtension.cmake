# pk5d_add_extension(<group> <sources>...
#                    [AVX2] [AVX2_SOURCES <sources>...] [NO_OPENMP]
#                    [LIBRARIES <targets or libs>...] [INCLUDE_DIRS <dirs>...]
#                    [DEFINITIONS <defs>...])
#
# Builds the nanobind module petakit5d._ext.<group> (target pk5d_ext_<group>) from the
# sources of python/cpp/<group>/, with the flags every kernel gets:
#   -O3, OpenMP, -ffp-contract=off, no -march=native, no fast-math
#   (MSVC: /O2, /fp:precise, /openmp:llvm; VS 2022 does not contract under /fp:precise;
#   MinGW-w64 GCC on Windows: GCC's, docs/decisions/0011-windows-toolchain.md).
# AVX2 adds -mavx2 (/arch:AVX2) on x86-64 to every source of the group, AVX2_SOURCES only to
# those; neither does anything on other processors. The sources write NB_MODULE(<group>, m)
# (or NB_MODULE(PK5D_MODULE_NAME, m)); python/cpp/common is on the include path.

set(PK5D_X86_64 FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|x64|X64)$")
  set(PK5D_X86_64 TRUE)
endif()
# MinGW-w64's GCC on Windows (CC=gcc CXX=g++, the Ninja generator): the compiler of the lab's
# cpp-zarr and cpp-tiff, an alternative to MSVC
set(PK5D_MINGW FALSE)
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(PK5D_MINGW TRUE)
endif()

# flags of every kernel
add_library(pk5d_flags INTERFACE)
if(MSVC)
  set(PK5D_POLICY_FLAGS "/O2;/fp:precise")
  set(PK5D_AVX2_FLAG "/arch:AVX2")
  target_compile_options(pk5d_flags INTERFACE $<$<NOT:$<CONFIG:Debug>>:/O2> /fp:precise /utf-8)
else()
  set(PK5D_POLICY_FLAGS "-O3;-ffp-contract=off")
  set(PK5D_AVX2_FLAG "-mavx2")
  target_compile_options(pk5d_flags INTERFACE $<$<NOT:$<CONFIG:Debug>>:-O3> -ffp-contract=off)
endif()
if(PK5D_MINGW AND PK5D_X86_64)
  # - the instruction set and tuning of Linux's GCC (x86-64, generic), unless CFLAGS/CXXFLAGS
  #   name their own: MinGW-w64's builds default to nocona and core2
  # - unaligned vector moves: GCC can't align the stack beyond 16 bytes with Windows' SEH
  #   (MAX_STACK_ALIGNMENT), but spills 32-byte AVX values with aligned moves (GCC bug 54412),
  #   which fault on a stack slot that isn't 32-byte aligned; binutils 2.38 and newer assemble
  #   every aligned vector move as an unaligned one, as fast on aligned data
  foreach(_f "-march=x86-64" "-mtune=generic")
    string(REGEX REPLACE "=.*" "=" _key "${_f}")
    if(NOT "${CMAKE_C_FLAGS} ${CMAKE_CXX_FLAGS}" MATCHES "${_key}")
      list(APPEND PK5D_POLICY_FLAGS "${_f}")
      target_compile_options(pk5d_flags INTERFACE "${_f}")
    endif()
  endforeach()
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_FLAGS "-mavx2 -Wa,-muse-unaligned-vector-move")
  check_cxx_source_compiles("int main() { return 0; }" PK5D_HAVE_UNALIGNED_VECTOR_MOVE)
  unset(CMAKE_REQUIRED_FLAGS)
  if(NOT PK5D_HAVE_UNALIGNED_VECTOR_MOVE)
    message(FATAL_ERROR "MinGW-w64: the assembler doesn't take -muse-unaligned-vector-move "
                        "(binutils 2.38 or newer), without which GCC's AVX code faults on Windows' "
                        "stack (GCC bug 54412)")
  endif()
  list(APPEND PK5D_POLICY_FLAGS "-Wa,-muse-unaligned-vector-move")
  target_compile_options(pk5d_flags INTERFACE "-Wa,-muse-unaligned-vector-move")
endif()

# Linux: the C++ runtime linked statically and hidden, so that a module built with a newer
# GCC than the system's loads also after NumPy or the I/O packages loaded the system's
# libstdc++ (the loader reuses a loaded libstdc++.so.6 whatever a module's RPATH says).
# MinGW-w64: libstdc++ and libgcc linked statically, as in the cpp-zarr and cpp-tiff
# libraries; winpthreads stays the toolchain's DLL, which libgomp-1.dll loads too
set(_pk5d_static_default OFF)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(_pk5d_static_default ON)
endif()
if(PK5D_MINGW)
  set(_pk5d_static_default ON)
endif()
option(PK5D_STATIC_CXX_RUNTIME "Link libstdc++ and libgcc statically (Linux, MinGW-w64)"
       ${_pk5d_static_default})
if(PK5D_STATIC_CXX_RUNTIME AND PK5D_MINGW)
  target_link_options(pk5d_flags INTERFACE -static-libstdc++ -static-libgcc)
elseif(PK5D_STATIC_CXX_RUNTIME)
  target_link_options(pk5d_flags INTERFACE -static-libstdc++ -static-libgcc "LINKER:--exclude-libs,ALL")
endif()

# MinGW-w64: the C math functions the modules call, made independent of the calling thread's x87
# precision and right (cpp/common/mingw_libm.c): a wrapper per function in pk5d_mingw_libm (a module
# carries the ones it calls), the modules linked with --wrap=<function>. The functions left to
# MinGW-w64 are exact at any precision (floor, ceil, trunc, nearbyint, rint, ldexp, frexp, modf, fmod,
# fabs, copysign, sqrt: GCC computes it with the instruction).
if(PK5D_MINGW)
  set(PK5D_MINGW_LIBM
      D_D:exp D_D:exp2 D_D:expm1 D_D:log D_D:log2 D_D:log10 D_D:log1p
      D_D:sin D_D:cos D_D:tan D_D:asin D_D:acos D_D:atan D_D:sinh D_D:cosh D_D:tanh
      D_D:asinh D_D:acosh D_D:atanh D_D:erf D_D:erfc D_D:cbrt D_D:tgamma D_D:lgamma
      D_DD:pow D_DD:atan2 D_DD:hypot
      F_F:expf F_F:exp2f F_F:expm1f F_F:logf F_F:log2f F_F:log10f F_F:log1pf
      F_F:sinf F_F:cosf F_F:tanf F_F:asinf F_F:acosf F_F:atanf F_F:sinhf F_F:coshf F_F:tanhf
      F_F:asinhf F_F:acoshf F_F:atanhf F_F:erff F_F:erfcf F_F:cbrtf F_F:tgammaf F_F:lgammaf
      F_FF:powf F_FF:atan2f F_FF:hypotf
      SINCOS:sincos SINCOSF:sincosf UCRT_FMA:fma UCRT_FMAF:fmaf
      ROUND:round LROUND:lround LLROUND:llround ROUNDF:roundf LROUNDF:lroundf LLROUNDF:llroundf)
  set(_pk5d_libm_src "${CMAKE_CURRENT_LIST_DIR}/../cpp/common/mingw_libm.c")
  get_filename_component(_pk5d_libm_src "${_pk5d_libm_src}" ABSOLUTE)
  set(_pk5d_libm_sources "")
  foreach(_e IN LISTS PK5D_MINGW_LIBM)
    string(REPLACE ":" ";" _kv "${_e}")
    list(GET _kv 0 _kind)
    list(GET _kv 1 _name)
    set(_f "${CMAKE_BINARY_DIR}/pk5d_mingw_libm/${_name}.c")
    file(CONFIGURE OUTPUT "${_f}" CONTENT
         "#define PK5D_LIBM_${_kind} 1\n#define PK5D_LIBM_NAME ${_name}\n#include \"${_pk5d_libm_src}\"\n"
         @ONLY)
    list(APPEND _pk5d_libm_sources "${_f}")
    target_link_options(pk5d_flags INTERFACE "LINKER:--wrap=${_name}")
  endforeach()
  add_library(pk5d_mingw_libm STATIC ${_pk5d_libm_sources})
  target_compile_options(pk5d_mingw_libm PRIVATE $<$<NOT:$<CONFIG:Debug>>:-O2> -ffp-contract=off)
  set_target_properties(pk5d_mingw_libm PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# MinGW-w64: the toolchain's notices of what the modules link in and the wheels carry (GCC's
# runtime: libgomp-1.dll, libgcc_s_seh-1.dll, libstdc++ and libgcc linked in; MinGW-w64's runtime,
# linked in; winpthreads: libwinpthread-1.dll), copied as the toolchain ships them, from its
# licenses/ folder (the MinGW-builds layout) next to the modules (petakit5d/_ext/licenses/), as its
# DLLs are. Built with MinGW-builds' GCC 12.2.0, the toolchain of tools/get_mingw.sh, their SHA-256
# must be that toolchain's; another toolchain's are copied as they are, their SHA-256 shown.
#   PK5D_MINGW_NOTICES  "<file>|<path under licenses/>" of each, for pk5d_finalize
set(PK5D_MINGW_NOTICES "")
if(PK5D_MINGW)
  get_filename_component(_pk5d_tc "${CMAKE_CXX_COMPILER}" DIRECTORY)
  get_filename_component(_pk5d_tc "${_pk5d_tc}/../licenses" ABSOLUTE)
  set(PK5D_MINGW_NOTICES_DIR "${_pk5d_tc}" CACHE PATH
      "MinGW-w64: the toolchain's licenses folder (gcc/, mingw-w64/, winpthreads/)")
  set(_pk5d_notices
      "gcc/COPYING3=8ceb4b9ee5adedde47b31e975c1d90c73ad27b6b165a1dcd80c7c545eb65b903"
      "gcc/COPYING.RUNTIME=9d6b43ce4d8de0c878bf16b54d8e7a10d9bd42b75178153e3af6a815bdc90f74"
      "mingw-w64/COPYING.MinGW-w64-runtime.txt=36a4857b6557a97baa82762444ffd6a74c7575015595dc26ff02215c299f81f6"
      "winpthreads/COPYING=63263614cdd29f2f93cba85e992f041b31f9fc7b4033692f31269489a8a1b177")
  execute_process(COMMAND "${CMAKE_CXX_COMPILER}" --version OUTPUT_VARIABLE _pk5d_banner ERROR_QUIET)
  set(_pk5d_pinned FALSE)
  if(_pk5d_banner MATCHES "x86_64-posix-seh-rev2, Built by MinGW-W64 project\\) 12\\.2\\.0")
    set(_pk5d_pinned TRUE)
  endif()
  foreach(_e IN LISTS _pk5d_notices)
    string(REPLACE "=" ";" _kv "${_e}")
    list(GET _kv 0 _rel)
    list(GET _kv 1 _sha)
    set(_f "${PK5D_MINGW_NOTICES_DIR}/${_rel}")
    if(NOT EXISTS "${_f}")
      message(WARNING "MinGW-w64: no ${_f}: the wheel would lack that notice of the toolchain "
                      "(PK5D_MINGW_NOTICES_DIR; tools/check_wheel.py fails such a wheel)")
      continue()
    endif()
    file(SHA256 "${_f}" _got)
    if(_pk5d_pinned AND NOT _got STREQUAL _sha)
      message(FATAL_ERROR "MinGW-w64: ${_f} has SHA-256 ${_got}, not MinGW-builds GCC 12.2.0's ${_sha}")
    elseif(NOT _pk5d_pinned)
      message(STATUS "petakit5d MinGW-w64 notice ${_rel}: SHA-256 ${_got} (not MinGW-builds 12.2.0's toolchain)")
    endif()
    list(APPEND PK5D_MINGW_NOTICES "${_f}|${_rel}")
  endforeach()
  if(_pk5d_pinned)
    list(LENGTH PK5D_MINGW_NOTICES _n)
    message(STATUS "petakit5d MinGW-w64: ${_n} notices of the toolchain, SHA-256 those of MinGW-builds GCC 12.2.0")
  endif()
endif()

# flags from the environment (CFLAGS, CXXFLAGS, toolchain files) that change results
set(_pk5d_bad_flags
    "-march=native" "-mtune=native" "-mcpu=native" "-ffast-math" "-Ofast" "-mfma"
    "-funsafe-math-optimizations" "-ffp-contract=fast" "-ffp-contract=on" "-ffp-model=fast"
    "/fp:fast" "/fp:contract" "-fp-model=fast")
foreach(_lang C CXX)
  foreach(_var CMAKE_${_lang}_FLAGS CMAKE_${_lang}_FLAGS_RELEASE CMAKE_${_lang}_FLAGS_RELWITHDEBINFO
               CMAKE_${_lang}_FLAGS_MINSIZEREL CMAKE_${_lang}_FLAGS_DEBUG)
    separate_arguments(_flags NATIVE_COMMAND "${${_var}}")
    foreach(_f IN LISTS _flags)
      string(TOLOWER "${_f}" _fl)
      foreach(_bad IN LISTS _pk5d_bad_flags)
        string(TOLOWER "${_bad}" _badl)
        if(_fl STREQUAL _badl OR (_badl STREQUAL "-march=native" AND _fl MATCHES "^-march=native"))
          message(FATAL_ERROR "${_var} has ${_f}, which changes floating-point results; "
                              "remove it (docs/decisions/0001-build.md)")
        endif()
      endforeach()
    endforeach()
  endforeach()
endforeach()

# the stable ABI (abi3) when scikit-build-core builds a wheel for it (wheel.py-api = "cp312" and
# Python 3.12 or newer: SKBUILD_SABI_COMPONENT); build trees and older Pythons build for their
# Python only
set(PK5D_STABLE_ABI "")
if(SKBUILD_SABI_COMPONENT)
  set(PK5D_STABLE_ABI STABLE_ABI)
endif()

define_property(GLOBAL PROPERTY PK5D_GROUPS
  BRIEF_DOCS "extension module groups" FULL_DOCS "groups added by pk5d_add_extension")

function(pk5d_add_extension group)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "AVX2;NO_OPENMP" ""
                        "AVX2_SOURCES;LIBRARIES;INCLUDE_DIRS;DEFINITIONS")
  if(NOT group MATCHES "^[a-z][a-z0-9_]*$")
    message(FATAL_ERROR "pk5d_add_extension: group '${group}' must be a lower-case identifier")
  endif()
  get_property(_groups GLOBAL PROPERTY PK5D_GROUPS)
  if(group IN_LIST _groups)
    message(FATAL_ERROR "pk5d_add_extension: group '${group}' is defined twice")
  endif()
  set(_sources ${ARG_UNPARSED_ARGUMENTS} ${ARG_AVX2_SOURCES})
  if(NOT _sources)
    message(FATAL_ERROR "pk5d_add_extension(${group}): no sources")
  endif()

  set(_target pk5d_ext_${group})
  # NOMINSIZE: nanobind would compile the module with -Os otherwise
  nanobind_add_module(${_target} NOMINSIZE ${PK5D_STABLE_ABI} ${_sources})
  set_target_properties(${_target} PROPERTIES
    OUTPUT_NAME ${group}
    LIBRARY_OUTPUT_DIRECTORY "${PK5D_EXT_OUTPUT_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${PK5D_EXT_OUTPUT_DIR}")
  foreach(_cfg IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER "${_cfg}" _cfg)
    set_target_properties(${_target} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY_${_cfg} "${PK5D_EXT_OUTPUT_DIR}"
      RUNTIME_OUTPUT_DIRECTORY_${_cfg} "${PK5D_EXT_OUTPUT_DIR}")
  endforeach()
  target_compile_definitions(${_target} PRIVATE PK5D_MODULE_NAME=${group} ${ARG_DEFINITIONS})
  target_include_directories(${_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${ARG_INCLUDE_DIRS})
  target_link_libraries(${_target} PRIVATE pk5d_common pk5d_flags ${ARG_LIBRARIES})
  if(NOT ARG_NO_OPENMP)
    target_link_libraries(${_target} PRIVATE pk5d_openmp)
    if(PK5D_OPENMP_RPATH)
      set_property(TARGET ${_target} APPEND PROPERTY BUILD_RPATH "${PK5D_OPENMP_RPATH}")
      set_property(TARGET ${_target} APPEND PROPERTY INSTALL_RPATH "${PK5D_OPENMP_RPATH}")
    endif()
    if(PK5D_OPENMP_SHARED AND NOT PK5D_OPENMP_SHARED STREQUAL "off")
      # the module needs the shared runtime only, and nothing it lacks (Pk5dCheckOpenMP.cmake)
      add_custom_command(TARGET ${_target} POST_BUILD
        COMMAND "${CMAKE_COMMAND}" "-DPK5D_MODULE=$<TARGET_FILE:${_target}>"
                "-DPK5D_SONAME=${PK5D_OPENMP_SHARED}" "-DPK5D_READELF=${PK5D_READELF}"
                -P "${PROJECT_SOURCE_DIR}/cmake/Pk5dCheckOpenMP.cmake"
        VERBATIM)
    endif()
  endif()
  if(PK5D_MINGW)
    # (last, so that the wrappers resolve the calls of every library before them; in one group with
    # MinGW-w64's library, whose functions call wrapped ones too: erf calls exp)
    target_link_libraries(${_target} PRIVATE "$<LINK_GROUP:RESCAN,pk5d_mingw_libm,mingwex>")
  endif()
  if(PK5D_X86_64)
    if(ARG_AVX2)
      target_compile_options(${_target} PRIVATE ${PK5D_AVX2_FLAG})
    elseif(ARG_AVX2_SOURCES)
      set_source_files_properties(${ARG_AVX2_SOURCES} TARGET_DIRECTORY ${_target}
                                  PROPERTIES COMPILE_OPTIONS "${PK5D_AVX2_FLAG}")
    endif()
  endif()
  install(TARGETS ${_target} LIBRARY DESTINATION petakit5d/_ext RUNTIME DESTINATION petakit5d/_ext)
  set_property(GLOBAL APPEND PROPERTY PK5D_GROUPS ${group})
endfunction()

# after all groups: OpenMP runtime and CUDA sources next to the modules, summary
function(pk5d_finalize)
  foreach(_dll IN LISTS PK5D_OPENMP_RUNTIME)
    file(COPY "${_dll}" DESTINATION "${PK5D_EXT_OUTPUT_DIR}")
    install(FILES "${_dll}" DESTINATION petakit5d/_ext)
  endforeach()
  # MinGW-w64: the toolchain's notices, as it ships them
  foreach(_e IN LISTS PK5D_MINGW_NOTICES)
    string(REPLACE "|" ";" _kv "${_e}")
    list(GET _kv 0 _file)
    list(GET _kv 1 _rel)
    get_filename_component(_sub "${_rel}" DIRECTORY)
    file(COPY "${_file}" DESTINATION "${PK5D_EXT_OUTPUT_DIR}/licenses/${_sub}")
    install(FILES "${_file}" DESTINATION "petakit5d/_ext/licenses/${_sub}")
  endforeach()
  # GPU kernels are compiled at run time from their sources (docs/decisions/0003-gpu.md)
  set(_cuda_dir "${PROJECT_SOURCE_DIR}/cpp/cuda")
  if(EXISTS "${_cuda_dir}")
    file(GLOB_RECURSE _cuda_files RELATIVE "${_cuda_dir}" CONFIGURE_DEPENDS
         "${_cuda_dir}/*.cu" "${_cuda_dir}/*.cuh" "${_cuda_dir}/*.h" "${_cuda_dir}/*.hpp")
    foreach(_f IN LISTS _cuda_files)
      configure_file("${_cuda_dir}/${_f}" "${PK5D_EXT_OUTPUT_DIR}/cuda/${_f}" COPYONLY)
      get_filename_component(_sub "${_f}" DIRECTORY)
      install(FILES "${_cuda_dir}/${_f}" DESTINATION "petakit5d/_ext/cuda/${_sub}")
    endforeach()
  endif()
  get_property(_groups GLOBAL PROPERTY PK5D_GROUPS)
  message(STATUS "petakit5d extension groups: ${_groups}")
  message(STATUS "petakit5d kernel flags: ${PK5D_POLICY_FLAGS}; AVX2 groups add ${PK5D_AVX2_FLAG} on x86-64")
endfunction()

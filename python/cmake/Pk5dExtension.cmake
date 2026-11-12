# pk5d_add_extension(<group> <sources>...
#                    [AVX2] [AVX2_SOURCES <sources>...] [NO_OPENMP]
#                    [LIBRARIES <targets or libs>...] [INCLUDE_DIRS <dirs>...]
#                    [DEFINITIONS <defs>...])
#
# Builds the nanobind module petakit5d._ext.<group> (target pk5d_ext_<group>) from the
# sources of python/cpp/<group>/, with the flags every kernel gets:
#   -O3, OpenMP, -ffp-contract=off, no -march=native, no fast-math
#   (MSVC: /O2, /fp:precise, /openmp:llvm; VS 2022 does not contract under /fp:precise).
# AVX2 adds -mavx2 (/arch:AVX2) on x86-64 to every source of the group, AVX2_SOURCES only to
# those; neither does anything on other processors. The sources write NB_MODULE(<group>, m)
# (or NB_MODULE(PK5D_MODULE_NAME, m)); python/cpp/common is on the include path.

set(PK5D_X86_64 FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|x64|X64)$")
  set(PK5D_X86_64 TRUE)
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

# Linux: the C++ runtime linked statically and hidden, so that a module built with a newer
# GCC than the system's loads also after NumPy or the I/O packages loaded the system's
# libstdc++ (the loader reuses a loaded libstdc++.so.6 whatever a module's RPATH says)
set(_pk5d_static_default OFF)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(_pk5d_static_default ON)
endif()
option(PK5D_STATIC_CXX_RUNTIME "Link libstdc++ and libgcc statically (Linux)" ${_pk5d_static_default})
if(PK5D_STATIC_CXX_RUNTIME)
  target_link_options(pk5d_flags INTERFACE -static-libstdc++ -static-libgcc "LINKER:--exclude-libs,ALL")
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
  nanobind_add_module(${_target} NOMINSIZE ${_sources})
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
  if(PK5D_OPENMP_RUNTIME)
    file(COPY "${PK5D_OPENMP_RUNTIME}" DESTINATION "${PK5D_EXT_OUTPUT_DIR}")
    install(FILES "${PK5D_OPENMP_RUNTIME}" DESTINATION petakit5d/_ext)
  endif()
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

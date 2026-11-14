# Run after the link of each extension module that shares an OpenMP runtime (Pk5dOpenMP.cmake):
#   cmake -DPK5D_MODULE=<module> -DPK5D_SONAME=<the runtime's SONAME> -DPK5D_READELF=<readelf>
#         -P Pk5dCheckOpenMP.cmake
# Fails when the module needs another OpenMP runtime as well, or imports an OpenMP symbol that the
# shared runtime doesn't define (the linker leaves such a symbol without a version, and the module
# wouldn't load).
execute_process(COMMAND "${PK5D_READELF}" -d -W "${PK5D_MODULE}"
                OUTPUT_VARIABLE _dyn RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "${PK5D_READELF} -d ${PK5D_MODULE} failed")
endif()
string(REGEX MATCHALL "\\(NEEDED\\)[^\n]*" _needed "${_dyn}")
set(_runtimes "")
foreach(_l IN LISTS _needed)
  if(_l MATCHES "\\[(lib(gomp|omp|iomp5)[^]]*)\\]")
    list(APPEND _runtimes "${CMAKE_MATCH_1}")
  endif()
endforeach()
if(NOT _runtimes STREQUAL PK5D_SONAME)
  message(FATAL_ERROR "${PK5D_MODULE} needs the OpenMP runtimes '${_runtimes}' instead of the "
                      "shared ${PK5D_SONAME} alone (OpenMP linked twice?); configure with "
                      "-DPK5D_OPENMP_SHARE=OFF")
endif()
execute_process(COMMAND "${PK5D_READELF}" --dyn-syms -W "${PK5D_MODULE}"
                OUTPUT_VARIABLE _syms RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "${PK5D_READELF} --dyn-syms ${PK5D_MODULE} failed")
endif()
string(REGEX MATCHALL " UND (GOMP_|GOACC_|omp_|acc_)[^\n]*" _und "${_syms}")
set(_missing "")
foreach(_l IN LISTS _und)
  if(NOT _l MATCHES "@")
    string(REGEX REPLACE "^ UND ([A-Za-z0-9_]+).*$" "\\1" _s "${_l}")
    list(APPEND _missing "${_s}")
  endif()
endforeach()
if(_missing)
  message(FATAL_ERROR "${PK5D_MODULE} calls OpenMP entry points that the shared runtime "
                      "${PK5D_SONAME} lacks: ${_missing}; configure with -DPK5D_OPENMP_SHARE=OFF "
                      "(the toolchain's runtime)")
endif()

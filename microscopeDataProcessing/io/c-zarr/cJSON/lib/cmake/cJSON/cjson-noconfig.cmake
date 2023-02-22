#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cjson" for configuration ""
set_property(TARGET cjson APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(cjson PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "m"
  IMPORTED_LOCATION_NOCONFIG "/home/xruan/Projects/XR_Repository/LLSM5DTools/microscopeDataProcessing/io/c-zarr/cJSON/lib/libcjson.so.1.7.15"
  IMPORTED_SONAME_NOCONFIG "libcjson.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS cjson )
list(APPEND _IMPORT_CHECK_FILES_FOR_cjson "/home/xruan/Projects/XR_Repository/LLSM5DTools/microscopeDataProcessing/io/c-zarr/cJSON/lib/libcjson.so.1.7.15" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

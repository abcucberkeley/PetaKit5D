# Whether the utils lib was build.
set(CJSON_UTILS_FOUND OFF)

# The include directories used by cJSON
set(CJSON_INCLUDE_DIRS "/home/xruan/Projects/XR_Repository/LLSM3DTools/microscopeDataProcessing/io/Parallel_Zarr_Reader/parallelReadZarr/cJSON/include")
set(CJSON_INCLUDE_DIR "/home/xruan/Projects/XR_Repository/LLSM3DTools/microscopeDataProcessing/io/Parallel_Zarr_Reader/parallelReadZarr/cJSON/include")

get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)

# The cJSON library
set(CJSON_LIBRARY "cjson")
if(ON)
  # Include the target
  include("${_dir}/cjson.cmake")
endif()

if(CJSON_UTILS_FOUND)
  # The cJSON utils library
  set(CJSON_UTILS_LIBRARY )
  # All cJSON libraries
  set(CJSON_LIBRARIES "" "cjson")
  if(ON)
    # Include the target
    include("${_dir}/cjson_utils.cmake")
  endif()
else()
  # All cJSON libraries
  set(CJSON_LIBRARIES "cjson")
endif()

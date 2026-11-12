# Header-only helpers shared by every group (python/cpp/common/pk5d/*.hpp)
add_library(pk5d_common INTERFACE)
target_include_directories(pk5d_common INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/common")

# Every folder with a CMakeLists.txt is a kernel group (pk5d_add_extension); a new folder is
# picked up without editing this file. common/ holds headers, cuda/ the GPU kernel sources
# that are compiled at run time.
file(GLOB _pk5d_group_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt")
list(SORT _pk5d_group_files)
foreach(_f IN LISTS _pk5d_group_files)
  get_filename_component(_dir "${_f}" DIRECTORY)
  if(NOT _dir STREQUAL "common" AND NOT _dir STREQUAL "cuda")
    add_subdirectory("${_dir}")
  endif()
endforeach()

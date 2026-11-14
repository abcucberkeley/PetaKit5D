# petakit5d._ext.core: build information and run-time checks (module.cpp), and the threaded
# element-wise primitives of ops/ (K08), which register with PK5D_BINDINGS
file(GLOB _pk5d_core_ops CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/ops/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/ops/*.c")
pk5d_add_extension(core module.cpp ${_pk5d_core_ops})

string(REPLACE ";" " " _pk5d_policy "${PK5D_POLICY_FLAGS}")
string(STRIP "${CMAKE_CXX_FLAGS} ${_pk5d_policy}" _pk5d_flags)
target_compile_definitions(pk5d_ext_core PRIVATE
  "PK5D_BUILD_FLAGS=\"${_pk5d_flags}\""
  "PK5D_BUILD_CONFIG=\"$<CONFIG>\""
  "PK5D_OPENMP_DESC=\"${PK5D_OPENMP_DESC}\""
  "PK5D_OPENMP_SHARED=\"${PK5D_OPENMP_SHARED}\"")

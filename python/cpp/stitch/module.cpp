// petakit5d._ext.stitch: the stitching kernels. Each source registers its bindings with
// PK5D_BINDINGS; the Python side is petakit5d.kernels.stitch.
#include <pk5d/common.hpp>

NB_MODULE(stitch, m) {
    m.doc() = "Stitching kernels (feather and tile blending, distance map resize, integral image)";
    pk5d::run_bindings(m);
}

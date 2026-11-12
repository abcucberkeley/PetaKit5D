// petakit5d._ext.morph: MATLAB's binary morphology, labeling, distance transform, convex hull
// images and median filters, bit for bit (petakit5d.matlabcompat.morphology, label, distance).
#include <pk5d/common.hpp>

NB_MODULE(morph, m) {
    m.doc() = "Binary morphology, connected components, distance transform, convex hulls and median "
              "filters as MATLAB computes them";
    pk5d::run_bindings(m);
}

// NaN replaced by a value: replace_nan_with_value.c (pointDetection/mex_functions), used by the
// stitching blending and the point detection. The Python side passes the value already in the
// input's class (the mex's mxGetScalar, then (float) for single input).
#include "stitch.hpp"

#include <cmath>

namespace {

using namespace pk5d_stitch;

template <typename T>
void bind_replace(nb::module_ &m) {
    m.def(
        "replace_nan_with_value",
        [](in1d<T> a, out1d<T> out, T val, int num_threads) {
            pk5d::check_threads(num_threads);
            const int64_t n = (int64_t)a.shape(0);
            need_size(out, (uint64_t)n, "output");
            const T *in = a.data();
            T *o = out.data();
            nb::gil_scoped_release release;
#pragma omp parallel for num_threads(num_threads) if (num_threads > 1)
            for (int64_t i = 0; i < n; i++) {
                if (std::isnan(in[i]))
                    o[i] = val;
                else
                    o[i] = in[i];
            }
        },
        nb::arg("a").noconvert(), nb::arg("out").noconvert(), nb::arg("value"),
        nb::arg("num_threads"), "out = a with NaN replaced by value (out may be a).");
}

}  // namespace

PK5D_BINDINGS(bind_replace_nan) {
    bind_replace<float>(m);
    bind_replace<double>(m);
}

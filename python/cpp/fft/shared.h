// For an extension module other than petakit5d._ext.fft that links pk5d_fft_driver (its kernels use the
// backends, the planner lock or the cached plans of backend.h): the module uses the process's one FFT state,
// petakit5d._ext.fft's, instead of a copy of its own. Call attach_to_fft_module() first thing in NB_MODULE (the
// interpreter lock is held there), then def_shared_state_hooks(m) for the tests (tests/fft/test_shared_state.py).
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstdint>
#include <exception>
#include <string>
#include <vector>

#include "backend.h"

namespace pk5d::fft {

inline void attach_to_fft_module() {
    namespace nb = nanobind;
    nb::object cap = nb::module_::import_("petakit5d._ext.fft").attr("_shared_state");
    void* p = PyCapsule_GetPointer(cap.ptr(), kSharedCapsule);
    if (p == nullptr) throw nb::python_error();
    attach_shared(static_cast<const SharedApi*>(p));
}

// _fft_state_id(): the identity of the FFT state the module uses; _fft_cached_plan(be, prec, kind, n,
// howmany): the cached plan (as an integer); _fft_plan_lengths(be, lengths): MATLAB's four plans of each length,
// single and double precision, 4 lines per execute, made without the interpreter lock (the errors, if any)
inline void def_shared_state_hooks(nanobind::module_& m) {
    namespace nb = nanobind;
    m.def("_fft_state_id", [] { return reinterpret_cast<std::uintptr_t>(&state()); });
    m.def("_fft_cached_plan", [](int be, int prec, int kind, std::int64_t n, int howmany) {
        nb::gil_scoped_release g;
        return reinterpret_cast<std::uintptr_t>(cached_plan(be, prec, kind, n, howmany));
    });
    m.def("_fft_plan_lengths", [](int be, const std::vector<std::int64_t>& lengths) {
        std::vector<std::string> errors;
        nb::gil_scoped_release g;
        for (std::int64_t n : lengths) {
            try {
                for (int prec = 0; prec < 2; prec++) {
                    for (int kind : {kR2C, kC2CF, kC2CB, kC2R}) {
                        if (cached_plan(be, prec, kind, n, 4) == nullptr) {
                            errors.push_back("no plan of length " + std::to_string(n));
                        }
                    }
                }
            } catch (const std::exception& e) {
                errors.push_back(e.what());
            }
        }
        return errors;
    });
}

}  // namespace pk5d::fft

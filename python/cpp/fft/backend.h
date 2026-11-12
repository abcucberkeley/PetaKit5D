// FFTW entry points used by the MATLAB fft driver, for one precision: our statically linked FFTW
// (built with the codelet sets of MATLAB R2024b's libmwfftw3f / libmwfftw3) or, for tests, a
// library loaded at run time (MATLAB's own, through PETAKIT5D_TEST_FFTW).
#pragma once

#include <cstddef>
#include <mutex>
#include <string>

namespace pk5d::fft {

// FFTW's fftw_iodim64
struct Iodim {
    std::ptrdiff_t n, is, os;
};

struct FftwApi {
    void* (*plan_dft)(int, const Iodim*, int, const Iodim*, void*, void*, int, unsigned) = nullptr;
    void* (*plan_r2c)(int, const Iodim*, int, const Iodim*, void*, void*, unsigned) = nullptr;
    void* (*plan_c2r)(int, const Iodim*, int, const Iodim*, void*, void*, unsigned) = nullptr;
    void (*exec_dft)(void*, void*, void*) = nullptr;
    void (*exec_r2c)(void*, void*, void*) = nullptr;
    void (*exec_c2r)(void*, void*, void*) = nullptr;
    void (*destroy)(void*) = nullptr;
    char* (*sprint)(void*) = nullptr;
    void (*free_string)(void*) = nullptr;
    const char* version = "";
};

// the two precisions of one library (index 0: single, 1: double)
struct Backend {
    FftwApi api[2];
    std::string name;  // "builtin" or the loaded files
    void* handles[2] = {nullptr, nullptr};
};

// the backend in use (the built-in one unless load_backend was called)
const Backend& backend();
// loads the single and double precision libraries at libf and libd (throws std::runtime_error);
// an empty libf returns to the built-in FFTW. Cached plans of the old backend are dropped by the caller.
void set_backend(const std::string& libf, const std::string& libd);
// FFTW's planner is not thread-safe: every plan is made and destroyed under this lock
std::mutex& planner_lock();

}  // namespace pk5d::fft

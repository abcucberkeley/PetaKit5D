// The CPU FFT backends behind MATLAB's fft driver, each a table of FFTW's guru64 entry points (rank-1 plans
// of several lines, new-array execute) for single and double precision:
//   permissive  pocketfft (BSD-3), built in: cpp/fft/pocketfft_api.cpp
//   exact       FFTW 3.3.8 with MATLAB R2024b's codelet sets, an optional GPL add-on of two shared libraries
//               (cpp/fft/exact) loaded at run time; with it the driver gives MATLAB's bits
// No file of the default build includes fftw3.h: the exact backend is reached only through this table.
#pragma once

#include <cstddef>
#include <mutex>
#include <string>

namespace pk5d::fft {

enum BackendId { kPermissive = 0, kExact = 1 };

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

// the two precisions of one backend (index 0: single, 1: double)
struct Backend {
    FftwApi api[2];
    std::string name;   // "permissive" or "exact"
    std::string where;  // "built-in" or the loaded files
    void* handles[2] = {nullptr, nullptr};
    bool ready = false;
};

// the backend id (kPermissive, kExact); throws std::runtime_error for the exact one when it is not loaded
const Backend& backend(int id);
bool exact_loaded();
// loads the exact backend's single and double precision libraries (throws std::runtime_error); the
// caller drops the exact backend's cached plans first
void load_exact(const std::string& libf, const std::string& libd);
void unload_exact();
// the backend the calls use when they name none (kPermissive unless set)
int default_backend();
void set_default_backend(int id);
// FFTW's planner is not thread-safe: every plan is made and destroyed under this lock
std::mutex& planner_lock();

// pocketfft's table (pocketfft_api.cpp)
FftwApi pocketfft_api(int prec);

}  // namespace pk5d::fft

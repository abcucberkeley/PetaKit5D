// The CPU FFT backends behind MATLAB's fft driver, each a table of FFTW's guru64 entry points (rank-1 plans
// of several lines, new-array execute) for single and double precision:
//   permissive  pocketfft (BSD-3), built in: cpp/fft/pocketfft_api.cpp
//   exact       FFTW 3.3.8 with MATLAB R2024b's codelet sets, an optional GPL add-on of two shared libraries
//               (cpp/fft/exact) loaded at run time; with it the driver gives MATLAB's bits
// No file of the default build includes fftw3.h: the exact backend is reached only through this table.
#pragma once

#include <cstddef>
#include <cstdint>
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
// loads the exact backend's single and double precision libraries (throws std::runtime_error), after
// destroying the exact backend's cached plans
void load_exact(const std::string& libf, const std::string& libd);
void unload_exact();
// the backend the calls use when they name none (kPermissive unless set)
int default_backend();
void set_default_backend(int id);
// FFTW's planner is not thread-safe: every plan is made and destroyed under this lock
std::mutex& planner_lock();

// MATLAB's plans: kind kR2C, kC2CF (forward), kC2CB (backward) or kC2R
enum PlanKind { kR2C = 0, kC2CF = 1, kC2CB = 2, kC2R = 3 };
// a rank-1 plan of `howmany` contiguous lines of length n, as MATLAB makes them (FFTW_ESTIMATE, out of place,
// line distances n, or n / 2 + 1 on the complex side of r2c and c2r, scratch aligned as MATLAB's), on backend
// be in precision prec (0 single, 1 double): made once under the planner lock and kept for the process until
// the backend is unloaded
void* cached_plan(int be, int prec, int kind, std::int64_t n, int howmany);
// destroys backend be's cached plans
void clear_plans(int be);

// --- one state per process -------------------------------------------------------------------------------
// The backends, the default backend, the planner lock and the plan cache are one per process. They live in
// petakit5d._ext.fft, which exports this table as the capsule petakit5d._ext.fft._shared_state; another
// extension module that links this library attaches to it at import (shared.h), so that FFTW is loaded once,
// every plan is made under one lock and plans are cached and destroyed in one place.
inline constexpr std::uint32_t kSharedAbi = 1;
inline constexpr const char* kSharedCapsule = "petakit5d._ext.fft._shared_state";

struct SharedApi {
    std::uint32_t abi;           // kSharedAbi
    std::uint32_t size;          // sizeof(SharedApi)
    std::uint32_t backend_size;  // sizeof(Backend), whose fields both sides read
    const Backend& (*backend)(int id);
    bool (*exact_loaded)();
    void (*load_exact)(const std::string& libf, const std::string& libd);
    void (*unload_exact)();
    int (*default_backend)();
    void (*set_default_backend)(int id);
    std::mutex& (*planner_lock)();
    void* (*cached_plan)(int be, int prec, int kind, std::int64_t n, int howmany);
    void (*clear_plans)(int be);
};

// this module's own state (petakit5d._ext.fft exports it)
const SharedApi& own_state();
// the state the functions above use: the attached one, else this module's own
const SharedApi& state();
// makes this module use another module's state; throws std::runtime_error when the two modules were built from
// other sources (the ABI, the table's size or the Backend layout differ)
void attach_shared(const SharedApi* other);

// pocketfft's table (pocketfft_api.cpp)
FftwApi pocketfft_api(int prec);

}  // namespace pk5d::fft

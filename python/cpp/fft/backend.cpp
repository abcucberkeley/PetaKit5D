// The two backends: pocketfft built in, FFTW (the exact add-on) loaded at run time by path.
#include "backend.h"

#include <atomic>
#include <cstdlib>
#include <stdexcept>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace pk5d::fft {
namespace {

Backend make_permissive() {
    Backend b;
    b.api[0] = pocketfft_api(0);
    b.api[1] = pocketfft_api(1);
    b.name = "permissive";
    b.where = "built-in";
    b.ready = true;
    return b;
}

Backend& slot(int id) {
    static Backend b[2] = {make_permissive(), Backend()};
    if (id != kPermissive && id != kExact) throw std::invalid_argument("unknown FFT backend");
    return b[id];
}

std::atomic<int>& default_id() {
    static std::atomic<int> d(kPermissive);
    return d;
}

#ifdef _WIN32
void* open_library(const std::string& path) {
    const int n = MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
    std::wstring w(static_cast<size_t>(n > 0 ? n : 1), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, &w[0], n);
    return reinterpret_cast<void*>(LoadLibraryW(w.c_str()));
}
void* symbol(void* h, const char* name) {
    return reinterpret_cast<void*>(GetProcAddress(reinterpret_cast<HMODULE>(h), name));
}
void close_library(void* h) { FreeLibrary(reinterpret_cast<HMODULE>(h)); }
#else
void* open_library(const std::string& path) {
    int fl = RTLD_NOW | RTLD_LOCAL;
#ifdef RTLD_DEEPBIND
    // the loaded library resolves its own symbols first
    fl |= RTLD_DEEPBIND;
#endif
    return dlopen(path.c_str(), fl);
}
void* symbol(void* h, const char* name) { return dlsym(h, name); }
void close_library(void* h) { dlclose(h); }
#endif

// sprint_plan's string comes from the C library's malloc (FFTW 3.3.8: "free it with free()")
void free_string(void* s) { std::free(s); }

template <typename F>
void bind(F& f, void* h, const std::string& prefix, const char* name, const std::string& path) {
    void* s = symbol(h, (prefix + name).c_str());
    if (s == nullptr) {
        throw std::runtime_error("FFTW library " + path + " has no " + prefix + name);
    }
    f = reinterpret_cast<F>(s);
}

void load_api(FftwApi& a, void*& handle, const std::string& path, const std::string& prefix) {
    handle = open_library(path);
    if (handle == nullptr) {
        throw std::runtime_error("cannot load the FFTW library " + path);
    }
    bind(a.plan_dft, handle, prefix, "plan_guru64_dft", path);
    bind(a.plan_r2c, handle, prefix, "plan_guru64_dft_r2c", path);
    bind(a.plan_c2r, handle, prefix, "plan_guru64_dft_c2r", path);
    bind(a.exec_dft, handle, prefix, "execute_dft", path);
    bind(a.exec_r2c, handle, prefix, "execute_dft_r2c", path);
    bind(a.exec_c2r, handle, prefix, "execute_dft_c2r", path);
    bind(a.destroy, handle, prefix, "destroy_plan", path);
    bind(a.sprint, handle, prefix, "sprint_plan", path);
#ifdef _WIN32
    // the loaded DLL may use another C runtime: its plan strings (test-only) are not freed
    a.free_string = [](void*) {};
#else
    a.free_string = free_string;
#endif
    const void* v = symbol(handle, (prefix + "version").c_str());
    a.version = v != nullptr ? static_cast<const char*>(v) : "";
}

}  // namespace

const Backend& backend(int id) {
    const Backend& b = slot(id);
    if (!b.ready) {
        throw std::runtime_error(
            "the exact FFT backend is not loaded (FFTW 3.3.8 with MATLAB's codelets, the optional add-on)");
    }
    return b;
}

bool exact_loaded() { return slot(kExact).ready; }

int default_backend() { return default_id().load(); }

void set_default_backend(int id) {
    if (id == kExact && !exact_loaded()) backend(kExact);  // throws
    slot(id);
    default_id() = id;
}

std::mutex& planner_lock() {
    static std::mutex m;
    return m;
}

void unload_exact() {
    Backend& b = slot(kExact);
    for (void*& h : b.handles) {
        if (h != nullptr) close_library(h);
        h = nullptr;
    }
    b = Backend();
    if (default_id() == kExact) default_id() = kPermissive;
}

void load_exact(const std::string& libf, const std::string& libd) {
    unload_exact();
    Backend nb;
    try {
        load_api(nb.api[0], nb.handles[0], libf, "fftwf_");
        load_api(nb.api[1], nb.handles[1], libd, "fftw_");
    } catch (...) {
        for (void* h : nb.handles) {
            if (h != nullptr) close_library(h);
        }
        throw;
    }
    nb.name = "exact";
    nb.where = libf + ";" + libd;
    nb.ready = true;
    slot(kExact) = nb;
}

}  // namespace pk5d::fft

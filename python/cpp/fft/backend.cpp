// The built-in FFTW (static, symbols hidden in this module) and libraries loaded at run time.
#include "backend.h"

#include <cstdlib>
#include <stdexcept>

#include "fftw3.h"

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

static_assert(sizeof(Iodim) == sizeof(fftwf_iodim64), "iodim layout");

// typed wrappers of the built-in API, single precision
void* f_plan_dft(int r, const Iodim* d, int hr, const Iodim* h, void* in, void* out, int sign, unsigned fl) {
    return fftwf_plan_guru64_dft(r, reinterpret_cast<const fftwf_iodim64*>(d), hr, reinterpret_cast<const fftwf_iodim64*>(h),
        static_cast<fftwf_complex*>(in), static_cast<fftwf_complex*>(out), sign, fl);
}
void* f_plan_r2c(int r, const Iodim* d, int hr, const Iodim* h, void* in, void* out, unsigned fl) {
    return fftwf_plan_guru64_dft_r2c(r, reinterpret_cast<const fftwf_iodim64*>(d), hr,
        reinterpret_cast<const fftwf_iodim64*>(h), static_cast<float*>(in), static_cast<fftwf_complex*>(out), fl);
}
void* f_plan_c2r(int r, const Iodim* d, int hr, const Iodim* h, void* in, void* out, unsigned fl) {
    return fftwf_plan_guru64_dft_c2r(r, reinterpret_cast<const fftwf_iodim64*>(d), hr,
        reinterpret_cast<const fftwf_iodim64*>(h), static_cast<fftwf_complex*>(in), static_cast<float*>(out), fl);
}
void f_exec_dft(void* p, void* in, void* out) {
    fftwf_execute_dft(static_cast<fftwf_plan>(p), static_cast<fftwf_complex*>(in), static_cast<fftwf_complex*>(out));
}
void f_exec_r2c(void* p, void* in, void* out) {
    fftwf_execute_dft_r2c(static_cast<fftwf_plan>(p), static_cast<float*>(in), static_cast<fftwf_complex*>(out));
}
void f_exec_c2r(void* p, void* in, void* out) {
    fftwf_execute_dft_c2r(static_cast<fftwf_plan>(p), static_cast<fftwf_complex*>(in), static_cast<float*>(out));
}
void f_destroy(void* p) { fftwf_destroy_plan(static_cast<fftwf_plan>(p)); }
char* f_sprint(void* p) { return fftwf_sprint_plan(static_cast<fftwf_plan>(p)); }

// double precision
void* d_plan_dft(int r, const Iodim* d, int hr, const Iodim* h, void* in, void* out, int sign, unsigned fl) {
    return fftw_plan_guru64_dft(r, reinterpret_cast<const fftw_iodim64*>(d), hr, reinterpret_cast<const fftw_iodim64*>(h),
        static_cast<fftw_complex*>(in), static_cast<fftw_complex*>(out), sign, fl);
}
void* d_plan_r2c(int r, const Iodim* d, int hr, const Iodim* h, void* in, void* out, unsigned fl) {
    return fftw_plan_guru64_dft_r2c(r, reinterpret_cast<const fftw_iodim64*>(d), hr,
        reinterpret_cast<const fftw_iodim64*>(h), static_cast<double*>(in), static_cast<fftw_complex*>(out), fl);
}
void* d_plan_c2r(int r, const Iodim* d, int hr, const Iodim* h, void* in, void* out, unsigned fl) {
    return fftw_plan_guru64_dft_c2r(r, reinterpret_cast<const fftw_iodim64*>(d), hr,
        reinterpret_cast<const fftw_iodim64*>(h), static_cast<fftw_complex*>(in), static_cast<double*>(out), fl);
}
void d_exec_dft(void* p, void* in, void* out) {
    fftw_execute_dft(static_cast<fftw_plan>(p), static_cast<fftw_complex*>(in), static_cast<fftw_complex*>(out));
}
void d_exec_r2c(void* p, void* in, void* out) {
    fftw_execute_dft_r2c(static_cast<fftw_plan>(p), static_cast<double*>(in), static_cast<fftw_complex*>(out));
}
void d_exec_c2r(void* p, void* in, void* out) {
    fftw_execute_dft_c2r(static_cast<fftw_plan>(p), static_cast<fftw_complex*>(in), static_cast<double*>(out));
}
void d_destroy(void* p) { fftw_destroy_plan(static_cast<fftw_plan>(p)); }
char* d_sprint(void* p) { return fftw_sprint_plan(static_cast<fftw_plan>(p)); }

// sprint_plan's string comes from the C library's malloc (FFTW 3.3.8: "free it with free()")
void free_string(void* s) { std::free(s); }

Backend make_builtin() {
    Backend b;
    FftwApi& f = b.api[0];
    f.plan_dft = f_plan_dft;
    f.plan_r2c = f_plan_r2c;
    f.plan_c2r = f_plan_c2r;
    f.exec_dft = f_exec_dft;
    f.exec_r2c = f_exec_r2c;
    f.exec_c2r = f_exec_c2r;
    f.destroy = f_destroy;
    f.sprint = f_sprint;
    f.free_string = free_string;
    f.version = fftwf_version;
    FftwApi& d = b.api[1];
    d.plan_dft = d_plan_dft;
    d.plan_r2c = d_plan_r2c;
    d.plan_c2r = d_plan_c2r;
    d.exec_dft = d_exec_dft;
    d.exec_r2c = d_exec_r2c;
    d.exec_c2r = d_exec_c2r;
    d.destroy = d_destroy;
    d.sprint = d_sprint;
    d.free_string = free_string;
    d.version = fftw_version;
    b.name = "builtin";
    return b;
}

Backend& current() {
    static Backend b = make_builtin();
    return b;
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

const Backend& backend() { return current(); }

std::mutex& planner_lock() {
    static std::mutex m;
    return m;
}

void set_backend(const std::string& libf, const std::string& libd) {
    Backend& b = current();
    for (void*& h : b.handles) {
        if (h != nullptr) close_library(h);
        h = nullptr;
    }
    if (libf.empty()) {
        b = make_builtin();
        return;
    }
    Backend nb;
    try {
        load_api(nb.api[0], nb.handles[0], libf, "fftwf_");
        load_api(nb.api[1], nb.handles[1], libd, "fftw_");
    } catch (...) {
        for (void* h : nb.handles) {
            if (h != nullptr) close_library(h);
        }
        b = make_builtin();
        throw;
    }
    nb.name = libf + ";" + libd;
    b = nb;
}

}  // namespace pk5d::fft

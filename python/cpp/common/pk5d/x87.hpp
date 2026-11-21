// The x87 unit's precision on Windows builds with MinGW-w64's GCC (docs/decisions/0011-windows-toolchain.md).
//
// GCC computes long double (80 bits there, as on Linux) on the x87 unit, whose results depend on the precision
// in the thread's x87 control word: 64 bits on Linux and in MinGW-w64's own threads (winpthreads sets it for
// every thread it starts, OpenMP's among them), 53 bits in the threads of a Python built with MSVC (Windows'
// default). pk5d::x87_extended sets 64 bits for its scope on the thread that runs it and restores the thread's
// own afterwards, so that such code gives the same bits in every thread, and Linux's. MinGW-w64's own C math
// functions are made independent of the thread the same way (cpp/common/mingw_libm.c). Elsewhere (Linux,
// macOS, MSVC: no x87 code, or the setting is already 64 bits) the guard is empty and compiles to nothing.
#pragma once

namespace pk5d {

#if defined(__MINGW32__) && (defined(__x86_64__) || defined(__i386__))
class x87_extended {
public:
    // ("memory": the code in the scope stays between the two)
    x87_extended() {
        __asm__ __volatile__("fnstcw %0" : "=m"(saved_) : : "memory");
        unsigned short cw = static_cast<unsigned short>(saved_ | 0x0300);  // precision control: 64 bits
        __asm__ __volatile__("fldcw %0" : : "m"(cw) : "memory");
    }
    ~x87_extended() { __asm__ __volatile__("fldcw %0" : : "m"(saved_) : "memory"); }
    x87_extended(const x87_extended&) = delete;
    x87_extended& operator=(const x87_extended&) = delete;

private:
    unsigned short saved_;
};
#else
class x87_extended {
public:
    x87_extended() = default;
    x87_extended(const x87_extended&) = delete;
    x87_extended& operator=(const x87_extended&) = delete;
};
#endif

}  // namespace pk5d

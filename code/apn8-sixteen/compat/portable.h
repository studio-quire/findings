// portable.h: lets the sources (written for MSVC) build with gcc or clang. Force-included by the gcc/clang build
// line (-include compat/portable.h); MSVC builds don't use it. It changes no logic: each shim reproduces the MSVC
// intrinsic's result exactly, including at 0.
//   __popcnt(x), __lzcnt(x), _tzcnt_u32(x)   MSVC bit intrinsics on unsigned int (lzcnt/tzcnt of 0 = 32)
//   __popcnt64(x), _tzcnt_u64(x)             the same on 64 bits (tzcnt of 0 = 64)
//   _BitScanReverse(&i, m)                   index of the highest set bit; returns 0 (i untouched) if m = 0
//   _BitScanReverse64(&i, m)                 the same on 64 bits (six_check.cpp)
//   _lzcnt_u32(x)                            as __lzcnt (cover.cpp)
//   _aligned_malloc(n, a), _aligned_free(p)  MSVC's aligned allocation, by posix_memalign (six_check.cpp)
//   _pdep_u32(src, mask)                     BMI2 parallel deposit, in software (no -mbmi2 needed)
//   _pext_u64(src, mask)                     BMI2 parallel extract on 64 bits, in software (seven_check.cpp)
//   __argc, __argv                          MSVC CRT globals, read by review_engine.cpp and review_inv2.cpp
//                                            from a static initialiser (so before main)
#pragma once
#if !defined(_MSC_VER)

// On x86, pull in the compiler's own intrinsic headers first, so their declarations (some compilers declare
// __lzcnt or _tzcnt_u32 there) come before the macros below and aren't rewritten by them.
#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#endif

#include <stdlib.h>
static inline void* quire_aligned_malloc(size_t n, size_t a) { void* p = 0; return posix_memalign(&p, a, n) ? 0 : p; }
#define _aligned_malloc(n, a) quire_aligned_malloc((n), (a))
#define _aligned_free(p) free(p)

static inline unsigned quire_popcnt(unsigned x) { return (unsigned)__builtin_popcount(x); }
static inline unsigned quire_lzcnt(unsigned x) { return x ? (unsigned)__builtin_clz(x) : 32u; }
static inline unsigned quire_tzcnt(unsigned x) { return x ? (unsigned)__builtin_ctz(x) : 32u; }
#undef __popcnt
#undef __lzcnt
#undef _tzcnt_u32
#define __popcnt(x) quire_popcnt((unsigned)(x))
#define __lzcnt(x) quire_lzcnt((unsigned)(x))
#define _tzcnt_u32(x) quire_tzcnt((unsigned)(x))

static inline unsigned long long quire_popcnt64(unsigned long long x) { return (unsigned long long)__builtin_popcountll(x); }
static inline unsigned long long quire_tzcnt64(unsigned long long x) { return x ? (unsigned long long)__builtin_ctzll(x) : 64ull; }
static inline unsigned char quire_bsr(unsigned long* i, unsigned long m) {
    if (!m) return 0;
    *i = (unsigned long)(8 * sizeof(unsigned long) - 1 - __builtin_clzl(m));
    return 1;
}
static inline unsigned char quire_bsr64(unsigned long* i, unsigned long long m) {
    if (!m) return 0;
    *i = (unsigned long)(63 - __builtin_clzll(m));
    return 1;
}
static inline unsigned quire_pdep32(unsigned src, unsigned mask) {
    unsigned r = 0;
    for (unsigned bit = 1; mask; bit <<= 1, mask &= mask - 1)
        if (src & bit) r |= mask & (0u - mask);
    return r;
}
static inline unsigned long long quire_pext64(unsigned long long src, unsigned long long mask) {
    unsigned long long r = 0;
    for (unsigned long long bit = 1; mask; bit <<= 1, mask &= mask - 1)
        if (src & mask & (0ull - mask)) r |= bit;
    return r;
}
#undef __popcnt64
#undef _tzcnt_u64
#undef _BitScanReverse
#undef _BitScanReverse64
#undef _lzcnt_u32
#undef _pdep_u32
#undef _pext_u64
#define __popcnt64(x) quire_popcnt64((unsigned long long)(x))
#define _tzcnt_u64(x) quire_tzcnt64((unsigned long long)(x))
#define _BitScanReverse(i, m) quire_bsr((i), (unsigned long)(m))
#define _BitScanReverse64(i, m) quire_bsr64((i), (unsigned long long)(m))
#define _lzcnt_u32(x) quire_lzcnt((unsigned)(x))
#define _pdep_u32(s, m) quire_pdep32((unsigned)(s), (unsigned)(m))
#define _pext_u64(s, m) quire_pext64((unsigned long long)(s), (unsigned long long)(m))

#if defined(__APPLE__)
#include <crt_externs.h>
#define __argc (*_NSGetArgc())
#define __argv (*_NSGetArgv())
#elif defined(__linux__)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
// The command line from /proc/self/cmdline, read on first use (works during static initialisation).
struct QuireArgs {
    int argc = 0; char** argv = nullptr;
    QuireArgs() {
        static std::vector<char> buf; static std::vector<char*> ptrs;
        FILE* f = fopen("/proc/self/cmdline", "rb"); if (!f) return;
        char tmp[4096]; size_t n;
        while ((n = fread(tmp, 1, sizeof tmp, f)) > 0) buf.insert(buf.end(), tmp, tmp + n);
        fclose(f);
        for (size_t i = 0; i < buf.size(); i += strlen(&buf[i]) + 1) ptrs.push_back(&buf[i]);
        ptrs.push_back(nullptr);
        argc = (int)ptrs.size() - 1; argv = ptrs.data();
    }
};
static inline QuireArgs& quire_args() { static QuireArgs a; return a; }
#define __argc (quire_args().argc)
#define __argv (quire_args().argv)
#endif

#endif // !_MSC_VER

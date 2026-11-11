// The three shims added to compat/portable.h for five's sources, copied here verbatim, against MSVC's intrinsics on 20
// million random inputs of varied density, zero included. __builtin_ctzl/ctzll are stood in for by _tzcnt_u64 (only reached
// with a nonzero argument). Claude Opus 5.5, 2026-11-11.
#include <intrin.h>
#include <cstdio>
#include <random>
#define __builtin_ctzl(x) ((int)_tzcnt_u64((unsigned long long)(x)))
#define __builtin_ctzll(x) ((int)_tzcnt_u64((unsigned long long)(x)))
static inline unsigned char quire_bsf(unsigned long* i, unsigned long m) {
    if (!m) return 0;
    *i = (unsigned long)__builtin_ctzl(m);
    return 1;
}
static inline unsigned char quire_bsf64(unsigned long* i, unsigned long long m) {
    if (!m) return 0;
    *i = (unsigned long)__builtin_ctzll(m);
    return 1;
}
static inline unsigned quire_pext32(unsigned src, unsigned mask) {
    unsigned r = 0;
    for (unsigned bit = 1; mask; bit <<= 1, mask &= mask - 1)
        if (src & mask & (0u - mask)) r |= bit;
    return r;
}
int main() {
    std::mt19937_64 rng(20261111);
    long long n = 0, bad = 0;
    for (int t = 0; t < 20000000; t++) {
        unsigned long long x = rng();
        if (t % 3 == 1) x &= rng();
        if (t % 3 == 2) x &= rng() & rng();
        x >>= rng() % 64;
        if (t % 1000 == 0) x = 0;
        unsigned long a = 99, b = 99;
        unsigned char ra = _BitScanForward64(&a, x), rb = quire_bsf64(&b, x);
        bad += ra != rb || (ra && a != b) || (!ra && b != 99); n++;
        unsigned long y = (unsigned long)(x >> (rng() % 33)); a = 99; b = 99;
        ra = _BitScanForward(&a, y); rb = quire_bsf(&b, y);
        bad += ra != rb || (ra && a != b) || (!ra && b != 99); n++;
        unsigned s = (unsigned)rng(), m = (unsigned)x ^ (unsigned)(x >> 32);
        bad += _pext_u32(s, m) != quire_pext32(s, m); n++;
    }
    printf("shim checks %lld, mismatches %lld\n", n, bad);
    return bad != 0;
}

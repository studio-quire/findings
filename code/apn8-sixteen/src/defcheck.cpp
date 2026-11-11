// defcheck: the tree below a span of alternating forms on F_2^8, straight from the definition.
// Claude Fable 5.2, 2026-11-11, for the review of five's page. Shares no code with cover.cpp.
//
// A node is a span S in which no nonzero point lies in the radicals of two nonzero forms. Its children: take the
// lowest point p not in any radical, run over ALL 2^21 - 1 nonzero forms B with B(p,.) = 0, and keep B when, for every
// s in S, the radical of B + s (found by comparing B(x,.) with s(x,.) at all 255 points) is trivial or a line, and
// those radicals are disjoint from each other and from the radicals of S. No filter, no reduction mod F_y, no rank
// tables: equality of two functionals at a point is the only test. Rank-4 or lower forms in a new coset are rejected,
// as search A does at five (a sixth rank-4 form is six's case).
//
// input lines: tag c0 c1 ... (generators as 28-bit coefficient words, pair (i<j) number k <-> B(e_i, e_j))
// output: tag, nodes by span dimension from the start dimension to 8, solutions, full-below-8, dim-8 non-APN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <immintrin.h>
#include <intrin.h>
#include <omp.h>

typedef uint64_t u64;
typedef unsigned char u8;

struct Tab { alignas(32) u8 t[256]; };   // t[x] = B(x, .) as a byte
struct Mask { u64 m[4]; };

static inline int par(unsigned x) { return __popcnt(x) & 1; }

static void tab_from_coef(unsigned c, Tab& T) {
    u8 rows[8] = {0}; int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) if (c >> k & 1) { rows[i] |= 1 << j; rows[j] |= 1 << i; }
    for (int x = 0; x < 256; x++) { u8 v = 0; for (int a = 0; a < 8; a++) if (x >> a & 1) v ^= rows[a]; T.t[x] = v; }
}
static unsigned coef_from_tab(const Tab& T) {
    unsigned c = 0; int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) if (T.t[1 << i] >> j & 1) c |= 1u << k;
    return c;
}
// points x != 0 where the two functionals agree
static inline Mask eq(const Tab& A, const Tab& B) {
    Mask r;
    for (int w = 0; w < 4; w++) {
        __m256i a0 = _mm256_load_si256((const __m256i*)(A.t + 64 * w)), b0 = _mm256_load_si256((const __m256i*)(B.t + 64 * w));
        __m256i a1 = _mm256_load_si256((const __m256i*)(A.t + 64 * w + 32)), b1 = _mm256_load_si256((const __m256i*)(B.t + 64 * w + 32));
        u64 lo = (unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a0, b0)), hi = (unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a1, b1));
        r.m[w] = lo | hi << 32;
    }
    r.m[0] &= ~(u64)1;
    return r;
}
static inline int cnt(const Mask& a) { return (int)(__popcnt64(a.m[0]) + __popcnt64(a.m[1]) + __popcnt64(a.m[2]) + __popcnt64(a.m[3])); }
static inline bool meets(const Mask& a, const Mask& b) { return (a.m[0] & b.m[0]) | (a.m[1] & b.m[1]) | (a.m[2] & b.m[2]) | (a.m[3] & b.m[3]); }
static inline void join(Mask& a, const Mask& b) { for (int w = 0; w < 4; w++) a.m[w] |= b.m[w]; }

struct Counts { long long bydim[9] = {0}, sols = 0, fulls = 0, nonapn = 0; std::vector<std::string> found; };

static void explore(const std::vector<Tab>& el, const Mask& cov, int dim, Counts& C) {
    C.bydim[dim]++;
    if (cnt(cov) == 255) {
        if (dim == 8) {
            C.sols++;
            std::string s; char hx[16];
            for (int k = 0; k < 8; k++) { snprintf(hx, sizeof hx, " %07x", coef_from_tab(el[(size_t)1 << k])); s += hx; }
            C.found.push_back(s);
        } else C.fulls++;
        return;
    }
    if (dim == 8) { C.nonapn++; return; }
    int p = -1;
    for (int x = 1; x < 256; x++) if (!(cov.m[x >> 6] >> (x & 63) & 1)) { p = x; break; }
    // seven independent functionals vanishing at p
    int u[7], n = 0; bool in[256] = {false}; in[0] = true;
    for (int c = 1; c < 256 && n < 7; c++) {
        if (par(c & p) || in[c]) continue;
        u[n++] = c;
        for (int x = 255; x >= 0; x--) if (in[x]) in[x ^ c] = true;
    }
    if (n != 7) { fprintf(stderr, "functionals\n"); exit(1); }
    std::vector<Tab> Eloc(21);
    unsigned piv[28] = {0}; int rk = 0, k = 0;
    for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++, k++) {
        for (int x = 0; x < 256; x++) Eloc[k].t[x] = (u8)((par(u[i] & x) ? u[j] : 0) ^ (par(u[j] & x) ? u[i] : 0));
        if (Eloc[k].t[p]) { fprintf(stderr, "p not in radical\n"); exit(1); }
        unsigned v = coef_from_tab(Eloc[k]);
        Tab back; tab_from_coef(v, back);
        if (memcmp(back.t, Eloc[k].t, 256)) { fprintf(stderr, "not an alternating form\n"); exit(1); }
        for (int b = 27; b >= 0 && v; b--) if (v >> b & 1) { if (piv[b]) v ^= piv[b]; else { piv[b] = v; rk++; break; } }
    }
    if (rk != 21) { fprintf(stderr, "basis rank %d\n", rk); exit(1); }
    Tab R; memset(R.t, 0, 256);
    const size_t nel = el.size();
    for (unsigned m = 1; m < (1u << 21); m++) {
        const Tab& e = Eloc[_tzcnt_u32(m)];
        for (int w = 0; w < 8; w++) {
            __m256i a = _mm256_load_si256((const __m256i*)(R.t + 32 * w)), b = _mm256_load_si256((const __m256i*)(e.t + 32 * w));
            _mm256_store_si256((__m256i*)(R.t + 32 * w), _mm256_xor_si256(a, b));
        }
        Mask acc = cov; bool ok = true;
        for (size_t i = 0; i < nel; i++) {
            Mask r = eq(R, el[i]);            // radical of B + el[i], minus 0
            int c = cnt(r);
            if (c == 0) continue;
            if (c != 3 || meets(r, acc)) { ok = false; break; }
            join(acc, r);
        }
        if (!ok) continue;
        std::vector<Tab> ch(2 * nel);
        for (size_t i = 0; i < nel; i++) { ch[i] = el[i]; for (int x = 0; x < 256; x++) ch[nel + i].t[x] = el[i].t[x] ^ R.t[x]; }
        explore(ch, acc, dim + 1, C);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: defcheck <input> [maxdepth-ignored]\n"); return 2; }
    FILE* f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    struct Item { std::string tag; std::vector<unsigned> g; };
    std::vector<Item> items; char line[4096];
    while (fgets(line, sizeof line, f)) {
        Item it; char tag[256]; int off;
        if (sscanf(line, "%255s%n", tag, &off) != 1) continue;
        it.tag = tag; const char* q = line + off;
        while (true) { char* e; unsigned long v = strtoul(q, &e, 16); if (e == q) break; it.g.push_back((unsigned)v); q = e; }
        items.push_back(it);
    }
    fclose(f);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int t = 0; t < (int)items.size(); t++) {
        const Item& it = items[t];
        std::vector<Tab> el(1); memset(el[0].t, 0, 256);
        for (unsigned c : it.g) { Tab T; tab_from_coef(c, T); size_t n = el.size(); for (size_t i = 0; i < n; i++) { Tab X; for (int x = 0; x < 256; x++) X.t[x] = el[i].t[x] ^ T.t[x]; el.push_back(X); } }
        // the start span, from the definition: radicals of all nonzero elements, pairwise disjoint
        Tab Z; memset(Z.t, 0, 256);
        Mask cov{{0, 0, 0, 0}}; int r4 = 0, r6 = 0, bad = 0;
        for (size_t i = 1; i < el.size(); i++) {
            Mask r = eq(el[i], Z); int c = cnt(r);
            if (c == 255) bad = 1;                    // dependent generators
            if (c == 15) r4++; else if (c == 3) r6++; else if (c) bad = 1;
            if (meets(r, cov)) bad = 1;
            join(cov, r);
        }
        Counts C; int dim = (int)it.g.size();
        if (!bad) explore(el, cov, dim, C);
        std::string out; char buf[256];
        snprintf(buf, sizeof buf, "ITEM %s gens", it.tag.c_str()); out = buf;
        for (unsigned c : it.g) { snprintf(buf, sizeof buf, " %07x", c); out += buf; }
        snprintf(buf, sizeof buf, " start r4 %d r6 %d%s nodes", r4, r6, bad ? " ILLEGAL-START" : ""); out += buf;
        for (int k = 1; k < 9; k++) { snprintf(buf, sizeof buf, " %lld", C.bydim[k]); out += buf; }
        snprintf(buf, sizeof buf, " sols %lld full %lld nonapn %lld", C.sols, C.fulls, C.nonapn); out += buf;
        #pragma omp critical
        {
            printf("%s\n", out.c_str());
            for (auto& s : C.found) printf("SOLUTION %s basis%s\n", it.tag.c_str(), s.c_str());
            fflush(stdout);
        }
    }
    printf("done: items %zu\n", items.size());
    return 0;
}

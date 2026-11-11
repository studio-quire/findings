// Reviewer's independent reference for search B (five-check). Shares no code with fivecheck_*.cpp.
// The only thing taken from search B is forms.txt, to learn the *index* search B gave each form, so that slices can be compared;
// the forms themselves are re-derived here by a nullspace and the two sets asserted equal.
//
// A form is a 28-bit vector over pairs i<j. The sieve here is not B's: for a point a, Z_a = {B : a in rad B} is written down
// explicitly (21 generators), pushed to the quotient Alt(V)/S, and every element of its image is marked dead in a bitset.
// No K_a, no packed words, no P/Q split, no Gray-coded coset walk. Assembly is plain membership over a sorted list.
//
//   count                       all 22 * 28^5 tuples, admissibility by plain radicals, no prefix test; per (class, c1) counts
//   slice <cls> <q0> <q1>       every tuple with first index in [q0,q1): f0, f1, planes, T's; compare with runs/check5/r5
//   list <file> <nou?>          lines "cls d0 d1 d2 d3 d4": same, per tuple
//   known <file>                known S8s in raw coordinates: S = rank-4 forms + one coset; must return the S8
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <omp.h>
#include <intrin.h>
#include <immintrin.h>
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8;

static int PI[8][8]; // pair index
static int PA[28], PB[28];
static void initPairs() { int n = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { PI[i][j] = PI[j][i] = n; PA[n] = i; PB[n] = j; n++; } }

struct Rows { u8 r[8]; };
static inline Rows rowsOf(u32 f) { Rows R; memset(&R, 0, 8); while (f) { unsigned long b; _BitScanForward(&b, f); f &= f - 1; R.r[PA[b]] |= (u8)(1 << PB[b]); R.r[PB[b]] |= (u8)(1 << PA[b]); } return R; }
static u32 fromU64(u64 x) { u32 f = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) if (x >> (8 * i + j) & 1) { f |= 1u << PI[i][j]; if (!(x >> (8 * j + i) & 1)) { puts("not symmetric"); exit(1); } } return f; }

struct P256 {
    u64 w[4] = {0, 0, 0, 0};
    inline void set(int v) { w[v >> 6] |= 1ull << (v & 63); }
    inline bool get(int v) const { return w[v >> 6] >> (v & 63) & 1; }
    inline bool meets(const P256& o) const { return (w[0] & o.w[0]) | (w[1] & o.w[1]) | (w[2] & o.w[2]) | (w[3] & o.w[3]); }
    inline void add(const P256& o) { for (int i = 0; i < 4; i++) w[i] |= o.w[i]; }
    inline int count() const { return (int)(__popcnt64(w[0]) + __popcnt64(w[1]) + __popcnt64(w[2]) + __popcnt64(w[3])); }
    inline bool operator==(const P256& o) const { return !memcmp(w, o.w, 32); }
};
// radical as a point set, by combining rows over all 255 vectors
static inline P256 radical(u32 f) {
    Rows R = rowsOf(f); u8 comb[256]; comb[0] = 0; P256 m;
    for (int v = 1; v < 256; v++) {
        unsigned long b; _BitScanForward(&b, v);
        comb[v] = comb[v & (v - 1)] ^ R.r[b];
        if (!comb[v]) m.set(v);
    }
    return m;
}
static inline bool nondegenerate(u32 f) { // quick path only
    Rows R = rowsOf(f); u8 piv[8]; int n = 0;
    for (int i = 0; i < 8; i++) { u8 x = R.r[i]; for (int t = 0; t < n; t++) { u8 hi = piv[t]; if ((x ^ hi) < x) x ^= hi; } if (!x) return false; piv[n++] = x; }
    return true;
}

// ---- solids, forms
struct Cls { int T4, T5; P256 sol[5], all; vector<u32> forms[5]; };
static int matv(int T, int x) { int y = 0; for (int i = 0; i < 4; i++) if (__popcnt((T >> (4 * i)) & 15 & x) & 1) y |= 1 << i; return y; }
static vector<u32> formsWithRadical(const P256& W) {
    // constraints B(w, e_j) = 0 for w in W: functionals on the 28 coordinates
    vector<u32> rows;
    for (int w = 1; w < 256; w++) if (W.get(w)) for (int j = 0; j < 8; j++) {
        u32 fn = 0;
        for (int i = 0; i < 8; i++) if ((w >> i & 1) && i != j) fn ^= 1u << PI[i][j];
        rows.push_back(fn);
    }
    // RREF
    vector<u32> rr; vector<int> pv;
    for (u32 x : rows) {
        for (size_t t = 0; t < rr.size(); t++) if (x >> pv[t] & 1) x ^= rr[t];
        if (!x) continue;
        unsigned long b; _BitScanForward(&b, x);
        for (size_t t = 0; t < rr.size(); t++) if (rr[t] >> b & 1) rr[t] ^= x;
        rr.push_back(x); pv.push_back((int)b);
    }
    u32 pm = 0; for (int p : pv) pm |= 1u << p;
    vector<u32> basis;
    for (int q = 0; q < 28; q++) if (!(pm >> q & 1)) {
        u32 v = 1u << q;
        for (size_t t = 0; t < rr.size(); t++) if (rr[t] >> q & 1) v |= 1u << pv[t];
        basis.push_back(v);
    }
    if (basis.size() != 6) { puts("nullspace dim"); exit(1); }
    vector<u32> out;
    for (int m = 1; m < 64; m++) { u32 f = 0; for (int i = 0; i < 6; i++) if (m >> i & 1) f ^= basis[i]; if (radical(f) == W) out.push_back(f); }
    return out;
}
static vector<Cls> C;
static void load() {
    FILE* f = fopen("runs/check5/configs.txt", "r"); if (!f) { puts("configs"); exit(1); }
    int r; unsigned t, u; long long a, g;
    while (fscanf(f, "%d %x %x %lld %lld", &r, &t, &u, &a, &g) == 5) if (r == 5) {
        Cls c; c.T4 = t; c.T5 = u;
        int Ts[3] = {0x8421, (int)t, (int)u};
        for (int x = 1; x < 16; x++) { c.sol[0].set(x); c.sol[1].set(x << 4); for (int i = 0; i < 3; i++) c.sol[2 + i].set(x | (matv(Ts[i], x) << 4)); }
        for (int i = 0; i < 5; i++) { if (c.all.meets(c.sol[i])) { puts("solids meet"); exit(1); } c.all.add(c.sol[i]); }
        C.push_back(c);
    }
    fclose(f);
    if (C.size() != 22) { puts("not 22"); exit(1); }
    for (auto& c : C) for (int i = 0; i < 5; i++) c.forms[i].assign(28, 0);
    FILE* g2 = fopen("runs/check5/review/forms.txt", "r"); if (!g2) { puts("forms.txt"); exit(1); }
    int ci, i, d; unsigned long long x; int n = 0;
    while (fscanf(g2, "%d %d %d %llx", &ci, &i, &d, &x) == 4) { C[ci].forms[i][d] = fromU64(x); n++; }
    fclose(g2);
    if (n != 22 * 5 * 28) { puts("forms count"); exit(1); }
    for (auto& c : C) for (int i2 = 0; i2 < 5; i2++) {
        vector<u32> mine = formsWithRadical(c.sol[i2]), theirs = c.forms[i2];
        sort(mine.begin(), mine.end()); sort(theirs.begin(), theirs.end());
        if (mine.size() != 28 || mine != theirs) { puts("forms differ from search B's"); exit(1); }
    }
}

// admissible? returns dim of span, or 0 if not. Plain definition: every point in at most one radical of a nonzero element of
// the span; a point of a solid only in the radical of that solid's own form.
static int admissible(const Cls& c, const u32* B) {
    u32 seen[32]; int ns = 0; P256 cover = c.all;
    for (int pass = 0; pass < 2; pass++) for (int m = 1; m < 32; m++) {
        int pc = (int)__popcnt(m);
        if (pc == 1) continue;
        if ((pass == 0) != (pc == 2)) continue; // pairs first: most tuples die there
        u32 x = 0; for (int i = 0; i < 5; i++) if (m >> i & 1) x ^= B[i];
        if (!x) continue;
        bool isB = false; for (int i = 0; i < 5; i++) if (x == B[i]) isB = true;
        if (isB) continue;
        bool dup = false; for (int i = 0; i < ns; i++) if (seen[i] == x) dup = true;
        if (dup) continue;
        seen[ns++] = x;
        if (nondegenerate(x)) continue;
        P256 r = radical(x);
        if (r.meets(cover)) return 0;
        cover.add(r);
    }
    // rank
    u32 piv[5]; int n = 0;
    for (int i = 0; i < 5; i++) { u32 x = B[i]; for (int t = 0; t < n; t++) if ((x ^ piv[t]) < x) x ^= piv[t]; if (x) piv[n++] = x; }
    return n;
}

// ---- the sieve by marking
struct Quot {
    int k, nq; u32 rr[8]; int pv[8]; u32 pm, npm;
    void build(const u32* b, int kk) {
        k = 0;
        for (int i = 0; i < kk; i++) {
            u32 x = b[i];
            for (int t = 0; t < k; t++) if (x >> pv[t] & 1) x ^= rr[t];
            if (!x) continue;
            unsigned long q; _BitScanForward(&q, x);
            for (int t = 0; t < k; t++) if (rr[t] >> q & 1) rr[t] ^= x;
            rr[k] = x; pv[k] = (int)q; k++;
        }
        pm = 0; for (int t = 0; t < k; t++) pm |= 1u << pv[t];
        npm = ((1u << 28) - 1) & ~pm; nq = 28 - k;
    }
    inline u32 phi(u32 f) const { for (int t = 0; t < k; t++) if (f >> pv[t] & 1) f ^= rr[t]; return _pext_u32(f, npm); }
    inline u32 lift(u32 x) const { return _pdep_u32(x, npm); }
};
// 21 generators of Z_a: forms pulled back along pi, the projection of V onto <e_i : i != p> killing a
static void zGens(int a, u32* out) {
    unsigned long p; _BitScanForward(&p, (unsigned)a);
    u8 pr[8]; // pi(e_s)
    for (int s = 0; s < 8; s++) pr[s] = (s == (int)p) ? (u8)((1 << s) ^ a) : (u8)(1 << s);
    int n = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) if (i != (int)p && j != (int)p) {
        u32 f = 0;
        for (int s = 0; s < 8; s++) for (int t = s + 1; t < 8; t++) {
            int v = ((pr[s] >> i & 1) & (pr[t] >> j & 1)) ^ ((pr[s] >> j & 1) & (pr[t] >> i & 1));
            if (v) f |= 1u << PI[s][t];
        }
        out[n++] = f;
    }
}
static void markPoint(const Quot& Q, int a, vector<u64>& dead) {
    u32 z[21]; zGens(a, z);
    u32 g[21]; int n = 0; int pv[21];
    for (int i = 0; i < 21; i++) {
        u32 x = Q.phi(z[i]);
        for (int t = 0; t < n; t++) if (x >> pv[t] & 1) x ^= g[t];
        if (!x) continue;
        unsigned long b; _BitScanForward(&b, x);
        g[n] = x; pv[n] = (int)b; n++;
    }
    if (n != 20) { printf("image dim %d at %d\n", n, a); exit(1); }
    u32 cur = 0; dead[0] |= 1;
    for (u32 it = 1; it < (1u << n); it++) {
        unsigned long b; _BitScanForward(&b, it);
        cur ^= g[b];
        dead[cur >> 6] |= 1ull << (cur & 63);
    }
}
struct Res { long long f0 = 0, f1 = 0, planes = 0, spaces = 0, full = 0, third = 0; int us = 0; };

static inline bool inList(const vector<u64>& bits, u32 x) { return bits[x >> 6] >> (x & 63) & 1; }

// subspaces of dimension e in L, all elements fitting pairwise; returns how many, with planes and 3-spaces counted on the way
static void assemble(const vector<u32>& L, const vector<P256>& U, const vector<u64>& bits, const map<u32, int>* unused, int e, bool useU,
                     const vector<int>& idx, Res& R, vector<vector<u32>>* outT) {
    auto fits = [&](u32 x, u32 y) { return !useU || !U[idx[x]].meets(U[idx[y]]); };
    vector<array<u32, 3>> planes;
    int n = (int)L.size();
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        u32 x = L[i], y = L[j], z = x ^ y;
        if (z <= y || !inList(bits, z)) continue;
        if (!fits(x, y) || !fits(x, z) || !fits(y, z)) continue;
        planes.push_back({x, y, z});
    }
    R.planes += (long long)planes.size();
    set<vector<u32>> cur;
    for (auto& p : planes) cur.insert(vector<u32>{p[0], p[1], p[2]});
    for (int d = 2; d < e; d++) {
        set<vector<u32>> nxt;
        for (auto& sp : cur) for (u32 t : L) {
            if (binary_search(sp.begin(), sp.end(), t)) continue;
            bool ok = true;
            for (u32 s : sp) if (!inList(bits, s ^ t)) { ok = false; break; }
            if (!ok) continue;
            vector<u32> el = sp; el.push_back(t);
            for (u32 s : sp) el.push_back(s ^ t);
            for (size_t i = 0; i < el.size() && ok; i++) for (size_t j = i + 1; j < el.size() && ok; j++) if (!fits(el[i], el[j])) ok = false;
            if (ok) { sort(el.begin(), el.end()); nxt.insert(el); }
        }
        if (d == 2) R.spaces += (long long)nxt.size();
        cur.swap(nxt);
    }
    R.full += (long long)cur.size();
    if (outT) for (auto& s : cur) outT->push_back(s);
}

struct Wk { vector<u64> dead, bits; vector<int> idx; };
// mode: useU = false stops at the solid points
static void analyse(Wk& w, const u32* b, int kk, const P256& solids, bool useU, int maxN, Res& R, vector<vector<u32>>* outT, Quot* qout = nullptr) {
    Quot Q; Q.build(b, kk);
    int k = Q.k, e = 8 - k; size_t nw = (size_t)1 << (Q.nq - 6);
    w.dead.assign(nw, 0);
    vector<u32> el{0};
    for (int i = 0; i < k; i++) { size_t m = el.size(); for (size_t j = 0; j < m; j++) el.push_back(el[j] ^ Q.rr[i]); }
    P256 US; for (size_t i = 1; i < el.size(); i++) US.add(radical(el[i]));
    R.us = US.count();
    for (int a = 1; a < 256; a++) if (solids.get(a)) markPoint(Q, a, w.dead);
    long long f0 = 0; for (size_t i = 0; i < nw; i++) f0 += __popcnt64(~w.dead[i]);
    R.f0 += f0;
    if (useU) for (int a = 1; a < 256; a++) if (US.get(a) && !solids.get(a)) markPoint(Q, a, w.dead);
    vector<u32> L;
    for (size_t i = 0; i < nw; i++) { u64 x = ~w.dead[i]; while (x) { unsigned long q; _BitScanForward64(&q, x); x &= x - 1; L.push_back((u32)(i * 64 + q)); } }
    R.f1 += (long long)L.size();
    (void)maxN;
    vector<P256> U(L.size());
    if (w.idx.size() != ((size_t)1 << Q.nq)) w.idx.assign((size_t)1 << Q.nq, -1);
    for (size_t i = 0; i < L.size(); i++) {
        u32 f = Q.lift(L[i]);
        P256 u;
        for (u32 s : el) {
            P256 r = radical(f ^ s);
            if (useU && (r.meets(u) || r.meets(US))) { puts("ref: survivor with clashing radicals"); exit(1); }
            u.add(r);
        }
        U[i] = u; w.idx[L[i]] = (int)i;
    }
    w.bits.assign(nw, 0);
    for (u32 x : L) w.bits[x >> 6] |= 1ull << (x & 63);
    assemble(L, U, w.bits, nullptr, e, useU, w.idx, R, outT);
    if (qout) *qout = Q;
}

int main(int argc, char** argv) {
    initPairs();
    string mode = argv[1];
    if (mode == "count") {
        load();
        FILE* o = fopen("runs/check5/review/ref_count.txt", "w");
        FILE* s = fopen("runs/check5/review/ref_sample.txt", "w");
        FILE* k4 = fopen("runs/check5/review/ref_k4.txt", "w");
        for (int ci = 0; ci < 22; ci++) {
            const Cls& c = C[ci]; long long tot = 0, tk4 = 0;
            long long sets[28] = {0}, k4s[28] = {0};
            #pragma omp parallel for schedule(dynamic, 1)
            for (int t = 0; t < 28 * 28; t++) {
                int d0 = t / 28, d1 = t % 28; u32 B[5]; B[0] = c.forms[0][d0]; B[1] = c.forms[1][d1];
                long long ls = 0, lk = 0;
                for (int d2 = 0; d2 < 28; d2++) for (int d3 = 0; d3 < 28; d3++) for (int d4 = 0; d4 < 28; d4++) {
                    B[2] = c.forms[2][d2]; B[3] = c.forms[3][d3]; B[4] = c.forms[4][d4];
                    int k = admissible(c, B);
                    if (!k) continue;
                    if (k != 5 && k != 4) { printf("dim %d\n", k); exit(1); }
                    ls++;
                    if (k == 4) lk++;
                    #pragma omp critical
                    {
                        if (k == 4) fprintf(k4, "%d %d %d %d %d %d\n", ci, d0, d1, d2, d3, d4);
                        else { u32 h = (u32)(((((ci * 28 + d0) * 28 + d1) * 28 + d2) * 28 + d3) * 28 + d4) * 2654435761u; if ((h >> 8) % 4999 == 0) fprintf(s, "%d %d %d %d %d %d\n", ci, d0, d1, d2, d3, d4); }
                    }
                }
                #pragma omp critical
                { sets[d0] += ls; k4s[d0] += lk; }
            }
            for (int q = 0; q < 28; q++) { fprintf(o, "%d %d %lld %lld\n", ci, q, sets[q], k4s[q]); tot += sets[q]; tk4 += k4s[q]; }
            fflush(o); fflush(s); fflush(k4);
            printf("class %d sets %lld k4 %lld\n", ci, tot, tk4); fflush(stdout);
        }
        fclose(o); fclose(s); fclose(k4);
        puts("REFCOUNTDONE");
        return 0;
    }
    if (mode == "slice") {
        load();
        int ci = atoi(argv[2]), q0 = atoi(argv[3]), q1 = atoi(argv[4]); bool nou = argc > 5 && !strcmp(argv[5], "nou");
        const Cls& c = C[ci];
        for (int q = q0; q < q1; q++) {
            Res T; long long nsets = 0, nk4 = 0;
            #pragma omp parallel
            {
                Wk w; Res loc; long long ls = 0, lk = 0;
                #pragma omp for schedule(dynamic, 1)
                for (int t = 0; t < 28 * 28; t++) {
                    u32 B[5]; B[0] = c.forms[0][q]; B[1] = c.forms[1][t / 28]; B[2] = c.forms[2][t % 28];
                    for (int d3 = 0; d3 < 28; d3++) for (int d4 = 0; d4 < 28; d4++) {
                        B[3] = c.forms[3][d3]; B[4] = c.forms[4][d4];
                        int k = admissible(c, B);
                        if (!k) continue;
                        ls++;
                        if (k == 4) { lk++; continue; }
                        analyse(w, B, 5, c.all, !nou, 0, loc, nullptr);
                    }
                }
                #pragma omp critical
                { T.f0 += loc.f0; T.f1 += loc.f1; T.planes += loc.planes; T.spaces += loc.spaces; T.full += loc.full; nsets += ls; nk4 += lk; }
            }
            printf("class %d q %d sets %lld k4 %lld f0 %lld f1 %lld planes %lld T %lld%s\n", ci, q, nsets, nk4, T.f0, T.f1, T.planes, T.full, nou ? " nou" : ""); fflush(stdout);
        }
        puts("REFSLICEDONE");
        return 0;
    }
    if (mode == "list") {
        load();
        bool nou = argc > 3 && !strcmp(argv[3], "nou");
        vector<array<int, 6>> tl;
        { FILE* f = fopen(argv[2], "r"); array<int, 6> x; while (fscanf(f, "%d %d %d %d %d %d", &x[0], &x[1], &x[2], &x[3], &x[4], &x[5]) == 6) tl.push_back(x); fclose(f); }
        vector<string> out(tl.size());
        #pragma omp parallel
        {
            Wk w;
            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < (int)tl.size(); i++) {
                auto& x = tl[i]; const Cls& c = C[x[0]];
                u32 B[5]; for (int j = 0; j < 5; j++) B[j] = c.forms[j][x[1 + j]];
                int k = admissible(c, B); Res R; char buf[256];
                if (!k) { snprintf(buf, sizeof buf, "%d %d %d %d %d %d inadmissible", x[0], x[1], x[2], x[3], x[4], x[5]); out[i] = buf; continue; }
                if (k == 4) { // survivors only; the assembly is hours here
                    Quot Q; Q.build(B, 5);
                    size_t nw = (size_t)1 << (Q.nq - 6); w.dead.assign(nw, 0);
                    vector<u32> el{0};
                    for (int t = 0; t < Q.k; t++) { size_t m = el.size(); for (size_t j = 0; j < m; j++) el.push_back(el[j] ^ Q.rr[t]); }
                    P256 US; for (size_t t = 1; t < el.size(); t++) US.add(radical(el[t]));
                    for (int a = 1; a < 256; a++) if (c.all.get(a)) markPoint(Q, a, w.dead);
                    long long f0 = 0; for (size_t t = 0; t < nw; t++) f0 += __popcnt64(~w.dead[t]);
                    for (int a = 1; a < 256; a++) if (US.get(a) && !c.all.get(a)) markPoint(Q, a, w.dead);
                    long long f1 = 0; for (size_t t = 0; t < nw; t++) f1 += __popcnt64(~w.dead[t]);
                    snprintf(buf, sizeof buf, "%d %d %d %d %d %d k4 f0 %lld f1 %lld", x[0], x[1], x[2], x[3], x[4], x[5], f0, f1);
                } else {
                    analyse(w, B, 5, c.all, !nou, 0, R, nullptr);
                    snprintf(buf, sizeof buf, "%d %d %d %d %d %d k5 f0 %lld f1 %lld planes %lld T %lld us %d", x[0], x[1], x[2], x[3], x[4], x[5], R.f0, R.f1, R.planes, R.full, R.us);
                }
                out[i] = buf;
            }
        }
        for (auto& s : out) puts(s.c_str());
        puts("REFLISTDONE");
        return 0;
    }
    if (mode == "known") {
        // raw coordinates, no normal form. Lines of known4.txt: n4 n2 n0 then eight hex forms.
        FILE* f = fopen(argv[2], "r"); char line[1024]; int nf = 0, ok = 0, tried = 0;
        Wk w;
        while (fgets(line, sizeof line, f)) {
            int a1, a2, a3; unsigned long long h[8]; int off = 0;
            if (sscanf(line, "%d %d %d %n", &a1, &a2, &a3, &off) != 3) continue;
            char* p = line + off; bool good = true;
            for (int i = 0; i < 8; i++) { char* e; h[i] = strtoull(p, &e, 16); if (e == p) good = false; p = e; }
            if (!good) continue;
            u32 b[8]; for (int i = 0; i < 8; i++) b[i] = fromU64(h[i]);
            vector<u32> all; for (int m = 1; m < 256; m++) { u32 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= b[i]; all.push_back(x); }
            vector<u32> S; P256 sol;
            for (u32 x : all) { P256 r = radical(x); if (r.count() == 15) { S.push_back(x); sol.add(r); } }
            nf++;
            Quot Q0; Q0.build(S.data(), (int)S.size());
            set<u32> cos; for (u32 x : all) { u32 y = Q0.phi(x); if (y) cos.insert(y); }
            int done = 0;
            for (u32 y : cos) {
                if (done == 3) break;  // three cosets per function: this is a check of the reference
                vector<u32> B = S; B.push_back(Q0.lift(y));
                Res R; vector<vector<u32>> Ts; Quot Q;
                analyse(w, B.data(), (int)B.size(), sol, true, 0, R, &Ts, &Q);
                tried++; done++;
                // S8 mod S must be among the T's
                set<u32> want; for (u32 x : all) { u32 z = Q.phi(x); if (z) want.insert(z); }
                bool hit = false;
                for (auto& t : Ts) { set<u32> g(t.begin(), t.end()); if (g == want) hit = true; }
                if (hit) ok++; else printf("function %d coset %x missed: f1 %lld T %zu\n", nf, y, R.f1, Ts.size());
            }
        }
        printf("known: %d functions, %d cosets tried, S8 returned %d\n", nf, tried, ok);
        return 0;
    }
    return 1;
}

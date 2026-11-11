// Five-check (search B): is there a quadratic APN function on F_2^8 with >= r rank-4 components? r = 5 is the question,
// r = 4 the control. See five-check.md.
//
// Decomposition. For one representative of each class of r solids and every tuple of forms (28^r, nothing normalised
// away), S = span(B_1..B_r) is known. S8 = S + T. No growth by points: every coset t + S of Alt(V)/S is sieved at once,
// block by block. V = R1 + R2 and t = (P, M, Q):
//   a in R1 sees only (P, M), a in R2 only (Q, M), so for a fixed M the P- and Q-lifts are sieved separately,
//   then the surviving (P, Q) pairs meet the other solids.
// A coset survives iff no form of it has a radical point on a solid (necessary: each point lies in one radical only).
// T is then assembled as a subspace of the survivors. Nothing is assumed about ranks outside the r forms, so a function
// with more than r rank-4 components is reached as well.
//
// Forms are 8 row bytes packed in a u64. Vector byte: low nibble x1, high nibble x2.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <chrono>
#include <omp.h>
#include <intrin.h>
#include <immintrin.h>
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16; typedef uint8_t u8;

static inline int par(u32 x) { return (int)(__popcnt(x) & 1); }
static inline int ctz64(u64 x) { unsigned long i; _BitScanForward64(&i, x); return (int)i; }

struct Mask {
    u64 w[4];
    Mask() { w[0] = w[1] = w[2] = w[3] = 0; }
    inline void set(int i) { w[i >> 6] |= 1ull << (i & 63); }
    inline bool get(int i) const { return w[i >> 6] >> (i & 63) & 1; }
    inline bool meets(const Mask& o) const { return (w[0] & o.w[0]) | (w[1] & o.w[1]) | (w[2] & o.w[2]) | (w[3] & o.w[3]); }
    inline void add(const Mask& o) { for (int i = 0; i < 4; i++) w[i] |= o.w[i]; }
    inline int count() const { return (int)(__popcnt64(w[0]) + __popcnt64(w[1]) + __popcnt64(w[2]) + __popcnt64(w[3])); }
    inline bool operator==(const Mask& o) const { return !memcmp(w, o.w, 32); }
};

static inline u8 rowOf(u64 f, int i) { return (u8)(f >> (8 * i)); }
static inline u8 apply(u64 f, int a) {  // the byte of B(a, .)
    u8 r = 0; int x = a;
    while (x) { int i = ctz64((u64)x); x &= x - 1; r ^= rowOf(f, i); }
    return r;
}
static inline int evalForm(u64 f, int a, int b) { return par(apply(f, a) & b); }

// radical of f: basis in out, returns its dimension. Row combinations that vanish are kernel vectors (f is symmetric).
static inline int radical(u64 f, u8* out) {
    u8 pr[8], pt[8]; u8 have = 0; int n = 0;
    for (int i = 0; i < 8; i++) {
        u8 x = rowOf(f, i), tg = (u8)(1 << i);
        while (x) {
            int b = ctz64(x);
            if (have >> b & 1) { x ^= pr[b]; tg ^= pt[b]; }
            else { have |= (u8)(1 << b); pr[b] = x; pt[b] = tg; break; }
        }
        if (!x) out[n++] = tg;
    }
    return n;
}
static inline Mask spanMask(const u8* b, int n) {
    Mask m; u8 el[256]; int c = 1; el[0] = 0;
    for (int i = 0; i < n; i++) { for (int j = 0; j < c; j++) el[c + j] = el[j] ^ b[i]; c *= 2; }
    for (int j = 1; j < c; j++) m.set(el[j]);
    return m;
}

// ---- geometry
struct Config {
    int r; int nT; u16 T[3]; // T[0] = I, then T4, T5 for graphs after R1, R2
    u8 pts[5][15];
    Mask solid[5], solidsAll;
    u64 forms[5][28];
};
static inline int matVec(u16 T, int v) {
    int r = 0;
    for (int i = 0; i < 4; i++) r |= par((T >> (4 * i)) & 15 & v) << i;
    return r;
}
static vector<array<u8, 4>>* alt28p;
static void buildAlt() {
    alt28p = new vector<array<u8, 4>>();
    for (int m = 0; m < 64; m++) {
        u8 a[4] = {0, 0, 0, 0}; int bit = 0;
        for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) { if (m >> bit & 1) { a[i] |= 1 << j; a[j] |= 1 << i; } bit++; }
        // Pfaffian
        int pf = (a[0] >> 1 & 1) * (a[2] >> 3 & 1) ^ (a[0] >> 2 & 1) * (a[1] >> 3 & 1) ^ (a[0] >> 3 & 1) * (a[1] >> 2 & 1);
        if (pf) alt28p->push_back({a[0], a[1], a[2], a[3]});
    }
    if (alt28p->size() != 28) { puts("alt28 != 28"); exit(1); }
}
// L: V -> F_2^4 with kernel solid i. i=0: x2; i=1: x1; else T x1 + x2.
static inline int Lmap(const Config& c, int i, int v) {
    int x1 = v & 15, x2 = v >> 4;
    if (i == 0) return x2;
    if (i == 1) return x1;
    return matVec(c.T[i - 2], x1) ^ x2;
}
static void buildConfig(Config& c, int r, u16 T4, u16 T5) {
    c.r = r; c.T[0] = 0x8421; c.T[1] = T4; c.T[2] = T5;
    memset(c.pts, 0, sizeof c.pts);
    c.solidsAll = Mask();
    for (int i = 0; i < r; i++) {
        int n = 0; c.solid[i] = Mask();
        for (int v = 1; v < 256; v++) if (Lmap(c, i, v) == 0) { c.pts[i][n++] = (u8)v; c.solid[i].set(v); }
        if (n != 15) { puts("solid size"); exit(1); }
        c.solidsAll.add(c.solid[i]);
        for (int d = 0; d < 28; d++) {
            const auto& a = (*alt28p)[d];
            u64 f = 0;
            for (int p = 0; p < 8; p++) {
                u8 row = 0;
                for (int q = 0; q < 8; q++) {
                    int lp = Lmap(c, i, 1 << p), lq = Lmap(c, i, 1 << q), val = 0;
                    for (int k = 0; k < 4; k++) if (lp >> k & 1) val ^= par(a[k] & lq);
                    row |= (u8)(val << q);
                }
                f |= (u64)row << (8 * p);
            }
            c.forms[i][d] = f;
            u8 b[8]; int n2 = radical(f, b);
            if (n2 != 4 || !(spanMask(b, 4) == c.solid[i])) { puts("form radical wrong"); exit(1); }
        }
    }
    if (c.solidsAll.count() != 15 * r) { puts("solids meet"); exit(1); }
}

// ---- echelon on forms
struct Ech {
    u64 row[28]; int piv[28]; int n = 0;
    inline u64 reduce(u64 f) const { for (int i = 0; i < n; i++) if (f >> piv[i] & 1) f ^= row[i]; return f; }
    inline bool add(u64 f) { f = reduce(f); if (!f) return false; piv[n] = ctz64(f); row[n] = f; n++; return true; }
};
static u64 gensM[16], gensP[6], gensQ[6];
static void buildGens() {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        int a = i, b = 4 + j;
        gensM[4 * i + j] = (1ull << (8 * a + b)) | (1ull << (8 * b + a));
    }
    int k = 0;
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) {
        gensP[k] = (1ull << (8 * i + j)) | (1ull << (8 * j + i));
        gensQ[k] = (1ull << (8 * (4 + i) + 4 + j)) | (1ull << (8 * (4 + j) + 4 + i));
        k++;
    }
}

struct Stats {
    long long tuples = 0, sets = 0, setsK[9] = {0}, skippedK3 = 0;
    long long f0 = 0, f1 = 0, f0max = 0, f1max = 0;
    long long nodes[5] = {0}, found = 0;
    long long hist[32] = {0}; // log2 bucket of f1
    long long ranks6 = 0;
    void add(const Stats& o) {
        tuples += o.tuples; sets += o.sets; skippedK3 += o.skippedK3; for (int i = 0; i < 9; i++) setsK[i] += o.setsK[i];
        f0 += o.f0; f1 += o.f1; f0max = max(f0max, o.f0max); f1max = max(f1max, o.f1max);
        for (int i = 0; i < 5; i++) nodes[i] += o.nodes[i];
        found += o.found; for (int i = 0; i < 32; i++) hist[i] += o.hist[i];
        ranks6 += o.ranks6;
    }
};

struct Item { u64 form; u32 key; int w; Mask U; };
static inline bool allNib(u64 x) {  // every 4-bit slot nonzero; unused slots are set by the caller
    u64 y = x | (x >> 1) | (x >> 2) | (x >> 3);
    return (y & 0x1111111111111111ull) == 0x1111111111111111ull;
}

struct Flags { bool useU = true, pivot = false, keepF0 = false; };

struct Worker {
    const Config* c; Flags fl; Stats st;
    vector<vector<u64>> results; // bases of S8 found
    vector<Item> items; vector<int> hashT; int hmask = 0;
    vector<u32> hashK;
    vector<u64> F0forms; // for flat comparison
    // current
    int k, e; u64 Sb[8]; Mask US; u64 Sel[256]; int nSel;
    set<vector<u32>> seen;
    bool skipAssembly = false;
    u64 prof[4] = {0, 0, 0, 0}, prof2[5] = {0, 0, 0, 0, 0};

    void hashBuild() {
        size_t sz = 16; while (sz < items.size() * 3) sz *= 2;
        hashT.assign(sz, -1); hashK.resize(sz); hmask = (int)sz - 1;
        for (int i = 0; i < (int)items.size(); i++) {
            u32 h = (items[i].key * 2654435761u) & hmask;
            while (hashT[h] >= 0) h = (h + 1) & hmask;
            hashT[h] = i; hashK[h] = items[i].key;
        }
    }
    inline int lookup(u32 key) const {
        u32 h = (key * 2654435761u) & hmask;
        while (hashT[h] >= 0) { if (hashK[h] == key) return hashT[h]; h = (h + 1) & hmask; }
        return -1;
    }
    void record(const vector<u32>& span) {
        // canonical id: sorted keys
        vector<u32> s(span); sort(s.begin(), s.end());
        if (!seen.insert(s).second) return;
        st.found++;
        vector<u64> basis(Sb, Sb + k);
        Ech ec; for (int i = 0; i < k; i++) ec.add(Sb[i]);
        for (u32 key : span) if (key) { int i = lookup(key); if (ec.add(items[i].form)) basis.push_back(items[i].form); }
        results.push_back(basis);
    }
    // Second level on: cand holds items that fit the span so far (one of each pair u, u + c is enough).
    void rec(int d, vector<u32>& span, Mask Uspan, const vector<int>& cand) {
        if (d == e) { record(span); return; }
        for (size_t ci = 0; ci < cand.size(); ci++) {
            const Item& cI = items[cand[ci]];
            Mask Un = Uspan; bool ok = true;
            for (u32 s : span) {
                int id = lookup(s ^ cI.key);
                if (id < 0) { ok = false; break; }
                if (fl.useU) { if (items[id].U.meets(Un)) { ok = false; break; } Un.add(items[id].U); }
            }
            if (!ok) continue;
            vector<u32> ns(span); for (u32 s : span) ns.push_back(s ^ cI.key);
            st.nodes[d]++;
            vector<int> nc;
            if (d + 1 < e) {
                for (size_t cj = ci + 1; cj < cand.size(); cj++) {
                    u32 uk = items[cand[cj]].key; bool good = true;
                    if (fl.useU && items[cand[cj]].U.meets(Un)) continue;
                    for (u32 s : ns) {
                        if (s == 0) continue;
                        if (s == uk) { good = false; break; }
                        int id = lookup(s ^ uk);
                        if (id < 0) { good = false; break; }
                        if (fl.useU && items[id].U.meets(Un)) { good = false; break; }
                    }
                    if (good) nc.push_back(cand[cj]);
                }
            }
            rec(d + 1, ns, Un, nc);
        }
    }

    // First level. Items in order (by weight if pivoting); T is reached from its first element c, so everything else
    // in T comes after c. With the pivot, c has weight <= rem/(2^e - 1) because T's radicals partition what S leaves.
    vector<u64> pointBits, present;
    void assemble() {
        int n = (int)items.size();
        size_t words = (n + 63) / 64;
        int bound = fl.pivot ? (255 - US.count()) / ((1 << e) - 1) : 1 << 30;
        if (present.empty()) present.assign(1 << 19, 0); // 2^25 bits
        for (int i = 0; i < n; i++) present[items[i].key >> 6] |= 1ull << (items[i].key & 63);
        if (fl.useU) {
            pointBits.assign(256 * words, 0);
            for (int i = 0; i < n; i++) for (int q = 0; q < 4; q++) {
                u64 x = items[i].U.w[q];
                while (x) { int b = ctz64(x); x &= x - 1; pointBits[(size_t)(64 * q + b) * words + (i >> 6)] |= 1ull << (i & 63); }
            }
        }
        vector<u64> dis(words);
        vector<int> nc;
        for (int ci = 0; ci < n; ci++) {
            const Item& cI = items[ci];
            if (cI.w > bound) break;
            st.nodes[0]++;
            fill(dis.begin(), dis.end(), 0);
            if (fl.useU) for (int q = 0; q < 4; q++) {
                u64 x = cI.U.w[q];
                while (x) {
                    int b = ctz64(x); x &= x - 1;
                    const u64* pb = &pointBits[(size_t)(64 * q + b) * words];
                    for (size_t t = 0; t < words; t++) dis[t] |= pb[t];
                }
            }
            nc.clear();
            for (size_t t = (ci + 1) >> 6; t < words; t++) {
                u64 x = ~dis[t];
                if (t == (size_t)((ci + 1) >> 6)) x &= ~0ull << ((ci + 1) & 63);
                while (x) {
                    int b = ctz64(x); x &= x - 1;
                    int u = (int)(t * 64 + b);
                    if (u >= n) break;
                    u32 sum = items[u].key ^ cI.key;
                    if (!(present[sum >> 6] >> (sum & 63) & 1)) continue;
                    int id = lookup(sum);
                    if (id <= u) continue; // the pair is kept once, from its earlier member; id < ci can't be in T with c first
                    if (fl.useU && (items[id].U.meets(cI.U) || items[id].U.meets(items[u].U))) continue;
                    nc.push_back(u);
                }
            }
            vector<u32> span{0, cI.key};
            rec(1, span, cI.U, nc);
        }
        for (int i = 0; i < n; i++) present[items[i].key >> 6] = 0;
    }

    // process the set S (basis b, dimension kk), which already passed the checks on S itself
    void process(const u64* b, int kk) {
        k = kk; e = 8 - k;
        for (int i = 0; i < k; i++) Sb[i] = b[i];
        // elements of S and the radical union
        nSel = 1 << k; Sel[0] = 0;
        for (int i = 0; i < k; i++) for (int j = 0; j < (1 << i); j++) Sel[(1 << i) + j] = Sel[j] ^ b[i];
        US = Mask();
        for (int i = 1; i < nSel; i++) { u8 bb[8]; int n = radical(Sel[i], bb); if (n) US.add(spanMask(bb, n)); }
        u64 tA = __rdtsc();
        // complement generators: M, then P, then Q
        Ech ec; for (int i = 0; i < k; i++) ec.add(b[i]);
        // (P and Q claimed first, so that M runs over M_4 / M(S): the outer loop is the small one)
        u64 gP[6], gQ[6], gM[16]; int md = 0, pd = 0, qd = 0;
        for (int i = 0; i < 6; i++) if (ec.add(gensP[i])) gP[pd++] = gensP[i];
        for (int i = 0; i < 6; i++) if (ec.add(gensQ[i])) gQ[qd++] = gensQ[i];
        for (int i = 0; i < 16; i++) if (ec.add(gensM[i])) gM[md++] = gensM[i];
        u64 gen[28]; int ng = 0;
        for (int i = 0; i < md; i++) gen[ng++] = gM[i];
        for (int i = 0; i < pd; i++) gen[ng++] = gP[i];
        for (int i = 0; i < qd; i++) gen[ng++] = gQ[i];
        if (ng != 28 - k) { puts("complement size"); exit(1); }
        // K_a per solid point, and the words
        // Points of S's own rank-6 radicals act like solid points: each is in one radical only, and two forms of a
        // coset sharing a point differ by an element of S with that point in its radical, so this one test also
        // gives disjointness inside a coset.
        int r = c->r;
        u8 pts[12][15]; int npts[12]; memset(npts, 0, sizeof npts);
        for (int i = 0; i < r; i++) { memcpy(pts[i], c->pts[i], 15); npts[i] = 15; }
        int nw = r;
        if (fl.useU) for (int v = 1; v < 256; v++) if (US.get(v) && !c->solidsAll.get(v)) {
            if (npts[nw] == 15) nw++;
            pts[nw][npts[nw]++] = (u8)v;
        }
        if (fl.useU && npts[nw] > 0) nw++;
        if (nw > 12) { puts("too many words"); exit(1); }
        u64 um[12];
        for (int i = 0; i < nw; i++) { um[i] = 0; for (int p = npts[i]; p < 16; p++) um[i] |= 1ull << (4 * p); }
        u64 W[28][12]; memset(W, 0, sizeof W);
        for (int i = 0; i < nw; i++) for (int p = 0; p < npts[i]; p++) {
            int a = pts[i][p];
            u8 w[8]; for (int j = 0; j < k; j++) w[j] = apply(b[j], a);
            // nullspace of the bytes w
            u8 piv[8], prow[8]; int np = 0; u8 used = 0;
            for (int j = 0; j < k; j++) {
                u8 x = w[j];
                for (int t = 0; t < np; t++) if (x >> piv[t] & 1) x ^= prow[t];
                if (x) {
                    int q = ctz64(x);
                    for (int t = 0; t < np; t++) if (prow[t] >> q & 1) prow[t] ^= x;
                    piv[np] = (u8)q; prow[np] = x; np++; used |= (u8)(1 << q);
                }
            }
            if (np != k - 1) { puts("rank at solid point"); exit(1); }
            u8 kb[8]; int nk = 0;
            for (int q = 0; q < 8; q++) if (!(used >> q & 1)) {
                u8 v = (u8)(1 << q);
                for (int t = 0; t < np; t++) if (prow[t] >> q & 1) v |= (u8)(1 << piv[t]);
                kb[nk++] = v;
            }
            // pick e vectors completing a in K_a
            u8 red[8]; int nr = 0, ne = 0; u8 vs[8];
            auto addv = [&](u8 v) -> bool {
                for (int t = 0; t < nr; t++) { int hb = ctz64(red[t]); if (v >> hb & 1) v ^= red[t]; }
                if (!v) return false;
                red[nr++] = v; return true;
            };
            addv((u8)a);
            for (int t = 0; t < nk; t++) if (addv(kb[t])) vs[ne++] = kb[t];
            if (ne != e) { puts("K_a size"); exit(1); }
            for (int g = 0; g < ng; g++) {
                u8 ga = apply(gen[g], a);
                for (int t = 0; t < e; t++) if (par(ga & vs[t])) W[g][i] |= 1ull << (4 * p + t);
            }
        }
        // tables for P and Q
        static thread_local vector<u64> XP, XQ;
        int np2 = 1 << pd, nq2 = 1 << qd;
        XP.assign((size_t)np2 * nw, 0); XQ.assign((size_t)nq2 * nw, 0);
        for (int x = 1; x < np2; x++) { int g = ctz64(x), y = x & (x - 1); for (int i = 0; i < nw; i++) XP[(size_t)x * nw + i] = XP[(size_t)y * nw + i] ^ W[md + g][i]; }
        for (int x = 1; x < nq2; x++) { int g = ctz64(x), y = x & (x - 1); for (int i = 0; i < nw; i++) XQ[(size_t)x * nw + i] = XQ[(size_t)y * nw + i] ^ W[md + pd + g][i]; }
        // the structural fact the split rests on
        for (int g = 0; g < pd; g++) if (W[md + g][1]) { puts("P generator seen on R2"); exit(1); }
        for (int g = 0; g < qd; g++) if (W[md + pd + g][0]) { puts("Q generator seen on R1"); exit(1); }
        u64 tB = __rdtsc(); u64 tRad = 0;
        items.clear();
        const u64* const XPp = XP.data(); const u64* const XQp = XQ.data();
        u64 XM[12]; memset(XM, 0, sizeof XM); u64 fm = 0; int mc = 0;
        int pl[64], ql[64]; alignas(32) u64 q2[68];
        long long nM = 1ll << md;
        F0forms.clear();
        long long f0 = 0;
        for (long long it = 0; it < nM; it++) {
            if (it) {
                int g = ctz64((u64)it);
                mc ^= 1 << g; fm ^= gen[g];
                for (int i = 0; i < nw; i++) XM[i] ^= W[g][i];
            }
            int npl = 0, nql = 0;
            for (int x = 0; x < np2; x++) if (allNib((XM[0] ^ XPp[(size_t)x * nw]) | um[0])) pl[npl++] = x;
            if (!npl) continue;
            for (int x = 0; x < nq2; x++) if (allNib((XM[1] ^ XQp[(size_t)x * nw + 1]) | um[1])) ql[nql++] = x;
            // word 2 for all Q at once, four lanes at a time; the rest per surviving pair
            prof2[0]++; prof2[1] += npl; prof2[2] += nql; prof2[3] += (u64)npl * nql;
            for (int bq = 0; bq < nql; bq++) q2[bq] = XM[2] ^ XQp[(size_t)ql[bq] * nw + 2];
            for (int a = 0; a < npl; a++) for (int b4 = 0; b4 < nql; b4 += 4) {
                const u64* xp = &XPp[(size_t)pl[a] * nw];
                __m256i x = _mm256_or_si256(_mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(q2 + b4)), _mm256_set1_epi64x((long long)xp[2])), _mm256_set1_epi64x((long long)um[2]));
                __m256i y = _mm256_or_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 1)), _mm256_or_si256(_mm256_srli_epi64(x, 2), _mm256_srli_epi64(x, 3)));
                const __m256i ones = _mm256_set1_epi64x(0x1111111111111111ll);
                int lanes = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(_mm256_and_si256(y, ones), ones)));
                if (nql - b4 < 4) lanes &= (1 << (nql - b4)) - 1;
              while (lanes) {
                int bq = b4 + ctz64((u64)lanes); lanes &= lanes - 1;
                prof2[4]++;
                const u64* xq = &XQp[(size_t)ql[bq] * nw];
                bool ok = true;
                for (int i = 3; i < r; i++) if (!allNib((XM[i] ^ xp[i] ^ xq[i]) | um[i])) { ok = false; break; }
                if (!ok) continue;
                f0++;
                for (int i = r; i < nw; i++) if (!allNib((XM[i] ^ xp[i] ^ xq[i]) | um[i])) { ok = false; break; }
                if (!ok) continue;
                u64 f = fm;
                for (int g = 0; g < pd; g++) if (pl[a] >> g & 1) f ^= gen[md + g];
                for (int g = 0; g < qd; g++) if (ql[bq] >> g & 1) f ^= gen[md + pd + g];
                u32 key = (u32)mc | ((u32)pl[a] << md) | ((u32)ql[bq] << (md + pd));
                if (fl.keepF0) F0forms.push_back(f);
                // The radicals of the coset, by direct computation; the word test above must have made them clean.
                Mask U;
                u64 tr0 = __rdtsc();
                for (int s = 0; s < nSel; s++) {
                    u8 bb[8]; int n = radical(f ^ Sel[s], bb);
                    if (n) {
                        Mask m = spanMask(bb, n);
                        if (m.meets(c->solidsAll) || (fl.useU && (m.meets(U) || m.meets(US)))) { puts("word test and radicals disagree"); exit(1); }
                        U.add(m);
                    }
                }
                Item I; I.form = f; I.key = key; I.U = U; I.w = U.count();
                items.push_back(I);
                tRad += __rdtsc() - tr0;
              }
            }
        }
        u64 tC = __rdtsc();
        prof[0] += tB - tA; prof[1] += tC - tB - tRad; prof[2] += tRad;
        st.f0 += f0; st.f1 += items.size();
        st.f0max = max(st.f0max, f0); st.f1max = max(st.f1max, (long long)items.size());
        { int h = 0; size_t n = items.size(); while (n) { h++; n >>= 1; } st.hist[h]++; }
        if (k == 4 && getenv("FC_DEBUG")) {
            map<int, int> wh; for (auto& it : items) wh[it.w]++;
            printf("k4 set: md %d pd %d qd %d f0 %lld f1 %zu US %d w:", md, pd, qd, f0, items.size(), US.count());
            for (auto& kv : wh) printf(" %d:%d", kv.first, kv.second);
            puts(""); fflush(stdout);
        }
        // assembly
        if (!skipAssembly && (int)items.size() >= (1 << e) - 1) {
            if (fl.pivot) sort(items.begin(), items.end(), [](const Item& x, const Item& y) { return x.w < y.w; });
            hashBuild(); seen.clear();
            assemble();
        }
        prof[3] += __rdtsc() - tC;
    }
};

// S-level search over all tuples with first index c1 and second c2, third c3; on each complete tuple check S and call proc.
struct Level {
    const Config* c;
    template <class F> void run(int c1, int c2, int c3, Stats& st, F proc) {
        u64 B[5]; int idx[5] = {c1, c2, c3, 0, 0};
        go(0, idx, B, st, proc);
    }
    // Span bookkeeping is recomputed per level: small.
    template <class F> void go(int lvl, int* idx, u64* B, Stats& st, F& proc) {
        int r = c->r;
        if (lvl == r) { leaf(B, st, proc); return; }
        int lo = lvl < 3 ? idx[lvl] : 0, hi = lvl < 3 ? idx[lvl] + 1 : 28;
        for (int x = lo; x < hi; x++) {
            B[lvl] = c->forms[lvl][x];
            if (!partial(B, lvl + 1)) { if (lvl == r - 1) st.tuples++; else { long long m = 1; for (int q = lvl + 1; q < r; q++) m *= 28; st.tuples += m; } continue; }
            go(lvl + 1, idx, B, st, proc);
        }
    }
    // Cheap necessary test on a prefix: every sum's radical is off the solids unless it is exactly a solid not yet placed,
    // and radicals are pairwise disjoint.
    bool partial(const u64* B, int n) {
        Mask placed; for (int i = 0; i < n; i++) placed.add(c->solid[i]);
        Mask seen = placed;
        u64 el[32]; int ne = 1; el[0] = 0; Ech ec;
        for (int i = 0; i < n; i++) {
            if (!ec.add(B[i])) continue;
            for (int j = 0; j < ne; j++) el[ne + j] = el[j] ^ B[i];
            ne *= 2;
        }
        for (int j = 1; j < ne; j++) {
            bool isB = false; for (int i = 0; i < n; i++) if (el[j] == B[i]) isB = true;
            if (isB) continue;
            u8 bb[8]; int d = radical(el[j], bb);
            if (!d) continue;
            Mask m = spanMask(bb, d);
            if (m.meets(seen)) return false;
            if (m.meets(c->solidsAll)) {
                bool whole = false;
                for (int i = n; i < c->r; i++) if (m == c->solid[i]) whole = true;
                if (!whole) return false;
                continue; // a form with the radical of an unplaced solid: it has to turn out to be that B
            }
            seen.add(m);
        }
        return true;
    }
    template <class F> void leaf(u64* B, Stats& st, F& proc) {
        st.tuples++;
        int r = c->r;
        // full check on S
        u64 basis[5]; int k = 0; Ech ec;
        for (int i = 0; i < r; i++) if (ec.add(B[i])) basis[k++] = B[i];
        u64 el[32]; int ne = 1; el[0] = 0;
        for (int i = 0; i < k; i++) { for (int j = 0; j < ne; j++) el[ne + j] = el[j] ^ basis[i]; ne *= 2; }
        Mask seen = c->solidsAll;
        for (int j = 1; j < ne; j++) {
            bool isB = false; for (int i = 0; i < r; i++) if (el[j] == B[i]) isB = true;
            if (isB) continue;
            u8 bb[8]; int d = radical(el[j], bb);
            if (!d) continue;
            Mask m = spanMask(bb, d);
            if (m.meets(seen)) return;
            seen.add(m);
        }
        st.sets++; st.setsK[k]++;
        if (8 - k > 4) { st.skippedK3++; return; }
        proc(basis, k);
    }
};

// from-scratch check of a candidate S8: returns true if the radicals partition V \ 0.
static bool verifyS8(const vector<u64>& b, int* cnt) {
    if (b.size() != 8) return false;
    int cover[256] = {0}; memset(cnt, 0, 9 * sizeof(int));
    for (int m = 1; m < 256; m++) {
        u64 f = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) f ^= b[i];
        if (!f) return false;
        int d = 0;
        for (int v = 1; v < 256; v++) { // plain evaluation, not radical()
            bool in = true;
            for (int i = 0; i < 8 && in; i++) if (evalForm(f, v, 1 << i)) in = false;
            if (in) { cover[v]++; d++; }
        }
        int dim = 0; while ((1 << dim) - 1 < d) dim++;
        cnt[dim]++;
    }
    for (int v = 1; v < 256; v++) if (cover[v] != 1) return false;
    return true;
}

static vector<Config> loadConfigs(int r) {
    vector<Config> out;
    FILE* f = fopen("runs/check5/configs.txt", "r");
    if (!f) { puts("no configs.txt"); exit(1); }
    int rr; unsigned t, u; long long a, g;
    while (fscanf(f, "%d %x %x %lld %lld", &rr, &t, &u, &a, &g) == 5) if (rr == r) { out.emplace_back(); buildConfig(out.back(), r, (u16)t, (u16)u); }
    fclose(f);
    return out;
}

static void printStats(FILE* o, const Stats& s) {
    fprintf(o, "tuples %lld sets %lld k4 %lld k5 %lld skippedK3 %lld f0 %lld f1 %lld f0max %lld f1max %lld nodes %lld %lld %lld %lld found %lld\n",
        s.tuples, s.sets, s.setsK[4], s.setsK[5], s.skippedK3, s.f0, s.f1, s.f0max, s.f1max, s.nodes[0], s.nodes[1], s.nodes[2], s.nodes[3], s.found);
    fprintf(o, "hist");
    for (int i = 0; i < 32; i++) if (s.hist[i]) fprintf(o, " %d:%lld", i, s.hist[i]);
    fprintf(o, "\n");
}

// flat <r> <cfg> <nsamples> <seed>: the block sieve against the definition. All cosets of Alt(V)/S, taken from a
// complement built in a different order, Gray-coded; a coset passes iff for every solid point a, no form of it kills a.
static int flatCheck(int r, int cfg, int ns, unsigned seed) {
    vector<Config> cs = loadConfigs(r);
    const Config& c = cs[cfg];
    srand(seed);
    Level L; L.c = &c;
    int done = 0, bad = 0;
    while (done < ns) {
        u64 B[5]; for (int i = 0; i < r; i++) B[i] = c.forms[i][(rand() * 32768 + rand()) % 28];
        // must pass the same check as the search
        Stats ls; bool got = false; u64 basis[5]; int k = 0;
        auto proc = [&](const u64* b, int kk) { got = true; k = kk; memcpy(basis, b, kk * 8); };
        L.leaf(B, ls, proc);
        if (!got) continue;
        Worker w; w.c = &c; w.fl.useU = false; w.fl.keepF0 = true;
        w.fl.pivot = true; // assembly is not what's measured here
        w.skipAssembly = true;
        w.process(basis, k);
        Ech es; for (int i = 0; i < k; i++) es.add(basis[i]);
        vector<u64> a; for (u64 f : w.F0forms) a.push_back(es.reduce(f));
        sort(a.begin(), a.end());
        // flat
        Ech ec; for (int i = 0; i < k; i++) ec.add(basis[i]);
        u64 gen[28]; int ng = 0;
        for (int i = 5; i >= 0; i--) if (ec.add(gensQ[i])) gen[ng++] = gensQ[i];
        for (int i = 15; i >= 0; i--) if (ec.add(gensM[i])) gen[ng++] = gensM[i];
        for (int i = 5; i >= 0; i--) if (ec.add(gensP[i])) gen[ng++] = gensP[i];
        vector<array<u64, 4>> vals; vector<int> ptsv;
        for (int i = 0; i < r; i++) for (int p = 0; p < 15; p++) {
            int pt = c.pts[i][p]; array<u64, 4> m = {0, 0, 0, 0};
            for (int sI = 0; sI < (1 << k); sI++) {
                u64 sf = 0; for (int j = 0; j < k; j++) if (sI >> j & 1) sf ^= basis[j];
                int v = apply(sf, pt); m[v >> 6] |= 1ull << (v & 63);
            }
            vals.push_back(m); ptsv.push_back(pt);
        }
        vector<u64> b; u64 f = 0;
        for (long long it = 0; it < (1ll << ng); it++) {
            if (it) f ^= gen[ctz64((u64)it)];
            bool ok = true;
            for (size_t q = 0; q < ptsv.size(); q++) { int v = apply(f, ptsv[q]); if (vals[q][v >> 6] >> (v & 63) & 1) { ok = false; break; } }
            if (ok) b.push_back(es.reduce(f));
        }
        sort(b.begin(), b.end());
        bool same = a == b;
        printf("cfg %d k %d: sieve %zu flat %zu %s%c", cfg, k, a.size(), b.size(), same ? "same" : "DIFFERENT", 10); fflush(stdout);
        if (!same) bad++;
        done++;
    }
    return bad;
}

// run <r> <cfg> <c1lo> <c1hi> <outdir> [nou] [pivot] [every=N] [onlyk=K]
int main(int argc, char** argv) {
    buildAlt(); buildGens();
    if (argc >= 6 && !strcmp(argv[1], "flat")) return flatCheck(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    if (argc < 6) { puts("usage: run r cfg c1lo c1hi outdir [nou] [pivot] [every=N]"); return 1; }
    string mode = argv[1];
    int r = atoi(argv[2]), cfg = atoi(argv[3]), c1lo = atoi(argv[4]), c1hi = atoi(argv[5]);
    string outdir = argc > 6 ? argv[6] : "runs/check5/t";
    Flags fl; int every = 1, onlyk = 0;
    for (int i = 7; i < argc; i++) {
        if (!strcmp(argv[i], "nou")) fl.useU = false;
        if (!strcmp(argv[i], "pivot")) fl.pivot = true;
        if (!strncmp(argv[i], "every=", 6)) every = atoi(argv[i] + 6);
        if (!strncmp(argv[i], "onlyk=", 6)) onlyk = atoi(argv[i] + 6);
    }
    vector<Config> cs = loadConfigs(r);
    if (cfg >= (int)cs.size()) { puts("cfg out of range"); return 1; }
    const Config& c = cs[cfg];
    for (int c1 = c1lo; c1 < c1hi; c1++) {
        char path[256]; snprintf(path, sizeof path, "%s/r%d_c%02d_q%02d.txt", outdir.c_str(), r, cfg, c1);
        { // finished slices end in DONE
            FILE* f = fopen(path, "r");
            if (f) {
                char line[4096], last[4096] = ""; while (fgets(line, sizeof line, f)) if (line[0] != '\n') strcpy(last, line);
                fclose(f);
                if (!strncmp(last, "DONE", 4)) continue;
            }
        }
        auto t0 = chrono::steady_clock::now();
        Stats tot; vector<vector<u64>> res;
        #pragma omp parallel
        {
            Worker w; w.c = &c; w.fl = fl;
            Level L; L.c = &c;
            long long ctr = 0;
            #pragma omp for schedule(dynamic, 1)
            for (int c23 = 0; c23 < 784; c23++) {
                int c2 = c23 / 28, c3 = c23 % 28;
                Stats ls;
                auto proc = [&](const u64* b, int k) {
                    if (onlyk && k != onlyk) return;
                    if (every == 0 || (every > 1 && (ctr++ % every))) return;
                    w.process(b, k);
                };
                L.run(c1, c2, c3, ls, proc);
                #pragma omp critical
                { tot.add(ls); }
            }
            #pragma omp critical
            { tot.add(w.st); for (auto& v : w.results) res.push_back(v);
              if (getenv("FC_PROF")) printf("prof Mcycles: setup %.1f sieve %.1f radicals %.1f assembly %.1f%c", w.prof[0] / 1e6, w.prof[1] / 1e6, w.prof[2] / 1e6, w.prof[3] / 1e6, 10), printf("Ms %llu npl %llu nql %llu pairs %llu word2pass %llu%c", w.prof2[0], w.prof2[1], w.prof2[2], w.prof2[3], w.prof2[4], 10); }
        }
        double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        string tmp = path; tmp += ".tmp";
        FILE* o = fopen(tmp.c_str(), "w");
        fprintf(o, "r %d cfg %d c1 %d every %d useU %d pivot %d\n", r, cfg, c1, every, (int)fl.useU, (int)fl.pivot);
        printStats(o, tot);
        int good = 0;
        for (auto& b : res) {
            int cnt[9]; bool ok = verifyS8(b, cnt);
            good += ok;
            fprintf(o, "S8 %s rad0 %d rad2 %d rad4 %d rad6 %d :", ok ? "APN" : "notAPN", cnt[0], cnt[2], cnt[4], cnt[6]);
            for (u64 x : b) fprintf(o, " %016llx", (unsigned long long)x);
            fprintf(o, "\n");
        }
        fprintf(o, "APN %d of %zu\nseconds %.1f\nDONE\n", good, res.size(), sec);
        fclose(o);
        remove(path); rename(tmp.c_str(), path);
        printf("c%d q%d: ", cfg, c1); printStats(stdout, tot); printf("  APN %d of %zu, %.1f s\n", good, res.size(), sec); fflush(stdout);
    }
    return 0;
}

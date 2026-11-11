// cover: every quadratic APN space S8 on F_2^8 whose rank-4 radicals are exactly a given set of solids (pairwise
// disjoint 4-spaces), with all rank-6 radicals inside the holes (the nonzero points outside the solids).
// Rowen (Claude Opus 5.5), 2026-10-24. For i = 4 (13 rank-4 components): README "Thirteen"; inputs from
// cover_prep.py (one cover per line: tag, k, 4k hex basis bytes, optional "seed" forms as 28-bit coefficients).
//
// Point branching, complete by the APN condition (every nonzero point lies in exactly one radical of S8): take the
// lowest-index solid not yet a radical, and branch over its 28 forms (d(pi x, pi y), pi: V -> V/R, d nondegenerate);
// once all solids are radicals, take the lowest hole point not yet covered and branch over every rank-6 form whose
// radical is a line inside the holes through it (13,888 forms per line). Each candidate is checked against the whole
// span: every element must be bent, rank 6 with its radical a line of holes, or rank 4 with radical one of the
// solids (each at most once), and all radicals pairwise disjoint. dim 8 with all 255 points covered = SOLUTION
// (printed with its rank counts). All points covered below dim 8 = FULL (the rest would be bent; printed, extend
// separately). No symmetry is used by default: all 28 forms on the first solid are tried. --reps (five only) keeps one core per
// Aut(P)-orbit, so with it the search finds at least one image of every S8, not every S8.
//
// usage: cover <input> [maxlines|-1] [splitdim]
// splitdim D (2026-10-31, for seven's p7-2222, whose dim-4 core held one thread for most of a day): phase 1 runs the
// tasks but stops at every node with all solids radicals, holes left and span dimension D, keeping its span as a
// frontier item; phase 2 runs the frontier items in parallel. Same tree, same counts (a frontier node is counted
// once, in phase 2).
//
// Hole stage by a linear filter (2026-11-04, for six, whose hole-stage nodes cost ~0.7 s each, ~0.7 us per candidate
// in push): for a covered point y (in the radical of some s0 in S), a legal new B must have B(y,.) outside
// F_y = {s(y,.) : s in S}, else y lies in the radicals of both s0 and some B + s. So instead of the lines through p
// times their 13,888 forms, Gray-code over the 21-dim space {B : B(p,.) = 0}, keep B(y,.) for the covered y up to
// date, and send only what passes F_y for every covered y to push, unchanged. Every form push would accept is in that
// space (radical a line through p) and passes the filter (a necessary condition), so the children are the same set
// and the counts identical; only their order differs. In six only ~1 of 2M candidates per node passes, so the
// default (holes_tree) solves the filter as a system of "psi_y(m) != 0" constraints over affine subspaces instead
// of enumerating it. --gray runs the enumeration, --old the original lines x forms.
//
// --reps FILE (five_orbits' output, with --only): a core not listed is dropped with everything below it; --p1all: --frange and
// --fsample keep the hole stage in phase 1, without which cores below the split dimension never reach the frontier.
// Options after the positional ones: --old | --gray; --only TAG (just that cover); --frange A B (profiling: phase 1
// explores nothing below the solid stage except to collect the frontier, phase 2 runs frontier items [A, B) only);
// --fsample N SEED (the same, with N items per cover at random, and an ITEM line per item to compare the modes on).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <atomic>
#include <mutex>
#include <algorithm>
#include <random>
#include <array>
#include <intrin.h>
#include <omp.h>
static std::vector<std::array<unsigned, 5>> g_reps;   // --reps: searched cores, sorted; empty means every core

typedef uint64_t u64;
typedef uint32_t u32;

static inline int par(u32 x) { return __popcnt(x) & 1; }

struct PS { u64 m[4]; };
static inline bool has(const PS& s, int x) { return s.m[x >> 6] >> (x & 63) & 1; }
static inline void add(PS& s, int x) { s.m[x >> 6] |= (u64)1 << (x & 63); }

struct Rad { u64 m[4]; int size; };
static inline void radical(u64 B, Rad& r) {
    r.m[0] = 1; r.m[1] = r.m[2] = r.m[3] = 0; r.size = 1;
    int v = 0, x = 0;
    for (unsigned g = 1; g < 256; g++) {
        int a = _tzcnt_u32(g);
        x ^= 1 << a;
        v ^= (B >> (8 * a)) & 255;
        if (!v) { r.m[x >> 6] |= (u64)1 << (x & 63); r.size++; }
    }
}
static u64 to_coef(u64 B) {
    u64 c = 0; int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) if (B >> (8 * i + j) & 1) c |= (u64)1 << k;
    return c;
}
static u64 from_coef(u64 c) {
    u64 B = 0; int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) if (c >> k & 1) { B |= (u64)1 << (8 * i + j); B |= (u64)1 << (8 * j + i); }
    return B;
}

// nondegenerate alternating forms on F_2^r as r rows of r bits
static std::vector<std::vector<int>> alt_nondeg(int r) {
    std::vector<std::pair<int, int>> pr;
    for (int i = 0; i < r; i++) for (int j = i + 1; j < r; j++) pr.push_back({i, j});
    std::vector<std::vector<int>> out;
    for (u32 bits = 0; bits < (1u << pr.size()); bits++) {
        std::vector<int> rows(r, 0);
        for (size_t k = 0; k < pr.size(); k++) if (bits >> k & 1) { rows[pr[k].first] |= 1 << pr[k].second; rows[pr[k].second] |= 1 << pr[k].first; }
        // rank by elimination
        std::vector<int> t = rows; int rk = 0;
        for (int c = 0; c < r; c++) {
            int p = -1; for (int i = rk; i < r; i++) if (t[i] >> c & 1) { p = i; break; }
            if (p < 0) continue;
            std::swap(t[p], t[rk]);
            for (int i = 0; i < r; i++) if (i != rk && (t[i] >> c & 1)) t[i] ^= t[rk];
            rk++;
        }
        if (rk == r) out.push_back(rows);
    }
    return out;
}
static std::vector<std::vector<int>> D4, D6;

// functionals u[0..8-k) spanning the annihilator of span(basis)
static void annihilator(const int* basis, int k, int* u) {
    bool inspan[256] = {false}; inspan[0] = true;
    int n = 0;
    for (int c = 1; c < 256 && n < 8 - k; c++) {
        bool ok = true;
        for (int i = 0; i < k; i++) if (par(c & basis[i])) { ok = false; break; }
        if (!ok || inspan[c]) continue;
        u[n++] = c;
        bool nx[256]; memcpy(nx, inspan, sizeof nx);
        for (int x = 0; x < 256; x++) if (inspan[x]) nx[x ^ c] = true;
        memcpy(inspan, nx, sizeof nx);
    }
    if (n != 8 - k) { fprintf(stderr, "annihilator: dependent basis\n"); exit(1); }
}
// B(x, y) = d(pi x, pi y), pi_i(x) = u_i . x
static u64 form_on(const int* u, int r, const std::vector<int>& d) {
    int p[8];
    for (int a = 0; a < 8; a++) { p[a] = 0; for (int i = 0; i < r; i++) p[a] |= ((u[i] >> a) & 1) << i; }
    u64 B = 0;
    for (int a = 0; a < 8; a++) {
        int v = 0;
        for (int i = 0; i < r; i++) if (p[a] >> i & 1) v ^= d[i];
        for (int b = 0; b < 8; b++) if (a != b && par(v & p[b])) B |= (u64)1 << (8 * a + b);
    }
    return B;
}

struct Cover {
    std::string tag;
    int k = 0;
    std::vector<std::vector<int>> basis;   // per solid, 4 points
    std::vector<PS> solmask;               // per solid, its 16 points (0 included)
    signed char pt2sol[256];
    PS holes;
    std::vector<std::vector<u64>> sforms;  // per solid, its 28 forms
    std::vector<u64> seeds;
};

static bool parse(const char* line, Cover& c) {
    char tag[256]; int k, off;
    if (sscanf(line, "%255s %d%n", tag, &k, &off) != 2) return false;
    c.tag = tag; c.k = k;
    const char* p = line + off;
    memset(c.pt2sol, -1, sizeof c.pt2sol);
    c.holes = PS{{0, 0, 0, 0}};
    for (int s = 0; s < k; s++) {
        std::vector<int> b(4);
        for (int i = 0; i < 4; i++) { char* e; b[i] = (int)strtoul(p, &e, 16); if (e == p) return false; p = e; }
        PS m{{0, 0, 0, 0}};
        std::vector<int> pts = {0};
        for (int i = 0; i < 4; i++) { size_t n = pts.size(); for (size_t j = 0; j < n; j++) pts.push_back(pts[j] ^ b[i]); }
        for (int x : pts) {
            if (has(m, x)) { fprintf(stderr, "%s: solid %d not 4-dim\n", tag, s); exit(1); }
            add(m, x);
            if (x) { if (c.pt2sol[x] >= 0) { fprintf(stderr, "%s: solids meet\n", tag); exit(1); } c.pt2sol[x] = (signed char)s; }
        }
        c.basis.push_back(b); c.solmask.push_back(m);
        int u[4]; annihilator(b.data(), 4, u);
        std::vector<u64> f;
        for (auto& d : D4) f.push_back(form_on(u, 4, d));
        c.sforms.push_back(f);
    }
    for (int x = 1; x < 256; x++) if (c.pt2sol[x] < 0) add(c.holes, x);
    const char* sd = strstr(p, "seed");
    if (sd) {
        p = sd + 4;
        while (true) { char* e; u64 v = strtoull(p, &e, 16); if (e == p) break; c.seeds.push_back(from_coef(v)); p = e; }
    }
    return true;
}

static std::mutex g_mu;

struct Search {
    const Cover* C;
    u64 elems[256]; int nel = 1, dim = 0;
    PS cov{{1, 0, 0, 0}};
    int ncov = 1;          // points in cov, 0 included
    int covered = 0;       // solids that are radicals
    long long nodes = 0, sols = 0, fulls = 0, nonapn = 0, nears = 0, bydim[9] = {0};
    int split = -1;                       // phase 1: stop at hole-level nodes of this dim
    std::vector<std::vector<u64>>* frontier = nullptr;
    int hmode = 0;                        // hole stage: 0 affine tree, 1 Gray code + filter, 2 old (lines x forms)
    bool p1only = false;                  // p1only: phase 1 skips hole-stage nodes off the split dimension
    long long cand = 0, passed = 0;       // hole-stage candidates enumerated, and those reaching push

    // -2 reject, -1 ok (bent or rank 6 in holes), s >= 0 rank 4 with radical solid s
    int classify(u64 B, Rad& r) const {
        radical(B, r);
        if (r.size == 1) return -1;
        if (r.size == 4) {
            if (((r.m[0] & ~(u64)1) & ~C->holes.m[0]) | (r.m[1] & ~C->holes.m[1]) | (r.m[2] & ~C->holes.m[2]) | (r.m[3] & ~C->holes.m[3])) return -2;
            return -1;
        }
        if (r.size == 16) {
            int x = r.m[0] & ~(u64)1 ? (int)_tzcnt_u64(r.m[0] & ~(u64)1) : -1;
            for (int w = 1; x < 0 && w < 4; w++) if (r.m[w]) x = 64 * w + (int)_tzcnt_u64(r.m[w]);
            int s = C->pt2sol[x];
            if (s < 0) return -2;
            const PS& m = C->solmask[s];
            if (m.m[0] != r.m[0] || m.m[1] != r.m[1] || m.m[2] != r.m[2] || m.m[3] != r.m[3]) return -2;
            return s;
        }
        return -2;
    }

    // add B to the span if legal; returns false (state unchanged) if not
    bool push(u64 B, PS& svcov, int& svncov, int& svcovered) {
        PS ncv = cov; int nc = ncov, ncd = covered;
        Rad r;
        for (int i = 0; i < nel; i++) {
            int s = classify(B ^ elems[i], r);
            if (s == -2) return false;
            if (s >= 0) { if (ncd >> s & 1) return false; ncd |= 1 << s; }
            u64 a0 = r.m[0] & ~(u64)1;
            if ((a0 & ncv.m[0]) | (r.m[1] & ncv.m[1]) | (r.m[2] & ncv.m[2]) | (r.m[3] & ncv.m[3])) return false;
            ncv.m[0] |= a0; ncv.m[1] |= r.m[1]; ncv.m[2] |= r.m[2]; ncv.m[3] |= r.m[3];
            nc += r.size - 1;
        }
        svcov = cov; svncov = ncov; svcovered = covered;
        for (int i = 0; i < nel; i++) elems[nel + i] = B ^ elems[i];
        nel *= 2; dim++; cov = ncv; ncov = nc; covered = ncd;
        return true;
    }
    void pop(const PS& svcov, int svncov, int svcovered) { nel /= 2; dim--; cov = svcov; ncov = svncov; covered = svcovered; }

    // --reps: this core's span in five_orbits' key (reduced echelon, descending), looked up in the list of orbit representatives
    bool is_rep() {
        unsigned b[5]; int r = 0;
        for (int i = 1; i < nel && r < dim; i++) {
            unsigned x = (unsigned)to_coef(elems[i]);
            for (int k = 0; k < r; k++) { unsigned long p; _BitScanReverse(&p, b[k]); if (x >> p & 1) x ^= b[k]; }
            if (!x) continue;
            unsigned long h; _BitScanReverse(&h, x);
            for (int k = 0; k < r; k++) if (b[k] >> h & 1) b[k] ^= x;
            b[r++] = x;
        }
        std::sort(b, b + r, [](unsigned a, unsigned c) { return a > c; });
        std::array<unsigned, 5> key{};
        for (int i = 0; i < r; i++) key[i] = b[i];
        return std::binary_search(g_reps.begin(), g_reps.end(), key);
    }

    void report(const char* what) {
        u64 piv[28] = {0}; std::vector<u64> basis;
        for (int i = 1; i < nel && (int)basis.size() < dim; i++) {
            u64 v = to_coef(elems[i]);
            for (int b = 27; b >= 0 && v; b--) if (v >> b & 1) {
                if (piv[b]) v ^= piv[b];
                else { piv[b] = v; basis.push_back(elems[i]); break; }
            }
        }
        int cnt[257] = {0}; Rad r;
        for (int i = 1; i < nel; i++) { radical(elems[i], r); cnt[r.size]++; }
        char buf[512]; int o = snprintf(buf, sizeof buf, "%s %s dim %d ranks 8:%d 6:%d 4:%d basis", what, C->tag.c_str(), dim, cnt[1], cnt[4], cnt[16]);
        for (u64 b : basis) o += snprintf(buf + o, sizeof buf - o, " %07llx", (unsigned long long)to_coef(b));
        std::lock_guard<std::mutex> g(g_mu);
        printf("%s\n", buf); fflush(stdout);
    }

    void dfs() {
        if (split >= 0 && dim == split && covered == (1 << C->k) - 1 && ncov < 256) {
            std::vector<u64> g;
            for (int k = 0; k < dim; k++) g.push_back(elems[1 << k]);
            frontier->push_back(std::move(g));
            return;
        }
        nodes++; bydim[dim]++;
        if (ncov == 256) {
            if (dim == 8) { sols++; report("SOLUTION"); } else { fulls++; report("FULL"); }
            return;
        }
        if (dim == 8) { nonapn++; return; }
        PS sc; int sn, sd;
        const int all = (1 << C->k) - 1;
        if (covered != all) {
            int s = _tzcnt_u32(~covered & all);
            for (u64 B : C->sforms[s]) if (push(B, sc, sn, sd)) {
                bool go = true;
                if (covered == all && ncov < 256) {   // every solid a radical, holes left
                    if (g_reps.empty() || is_rep()) { nears++; report("NEAR"); } else go = false;
                }
                if (go) dfs();
                pop(sc, sn, sd);
            }
            return;
        }
        if (p1only) return;
        int p = -1;
        for (int x = 1; x < 256; x++) if (!has(cov, x)) { p = x; break; }
        if (hmode == 0) { holes_tree(p); return; }
        if (hmode == 1) { holes_lin(p); return; }
        // every uncovered point is a hole now (all solids are radicals)
        for (int q = 1; q < 256; q++) {
            int t = p ^ q;
            if (q == p || q > t || !has(C->holes, q) || !has(C->holes, t) || has(cov, q) || has(cov, t)) continue;
            int lb[2] = {p, q}, u[6];
            annihilator(lb, 2, u);
            for (auto& d : D6) {
                u64 B = form_on(u, 6, d);
                cand++; passed++;
                if (push(B, sc, sn, sd)) { dfs(); pop(sc, sn, sd); }
            }
        }
    }

    // rows of B at every point: t[y] = B(y, .) as a byte
    static void rowtab(u64 B, unsigned char* t) {
        t[0] = 0; int x = 0, v = 0;
        for (unsigned g = 1; g < 256; g++) { int a = _tzcnt_u32(g); x ^= 1 << a; v ^= (B >> (8 * a)) & 255; t[x] = (unsigned char)v; }
    }

    void holes_lin(int p) {
        alignas(32) unsigned char ys[256];
        int L = 0;
        for (int y = 1; y < 256; y++) if (has(cov, y)) ys[L++] = (unsigned char)y;
        int Lp = (L + 31) & ~31;
        std::vector<PS> F(L, PS{{0, 0, 0, 0}});
        unsigned char t[256];
        for (int i = 0; i < nel; i++) {
            rowtab(elems[i], t);
            for (int j = 0; j < L; j++) add(F[j], t[ys[j]]);
        }
        int u[7]; int pb[1] = {p};
        annihilator(pb, 1, u);
        u64 E[21];
        std::vector<unsigned char> RK(21 * Lp, 0);
        int k = 0;
        for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++, k++) {
            u64 B = 0;
            for (int a = 0; a < 8; a++) {
                int row = (par(u[i] & (1 << a)) ? u[j] : 0) ^ (par(u[j] & (1 << a)) ? u[i] : 0);
                B |= (u64)(row & 255) << (8 * a);
            }
            E[k] = B;
            unsigned char* rk = &RK[k * Lp];
            for (int jj = 0; jj < L; jj++) { int y = ys[jj]; rk[jj] = (unsigned char)((par(u[i] & y) ? u[j] : 0) ^ (par(u[j] & y) ? u[i] : 0)); }
        }
        alignas(32) unsigned char R[256] = {0};
        u64 B = 0;
        PS sc; int sn, sd;
        for (u32 m = 1; m < (1u << 21); m++) {
            int b = _tzcnt_u32(m);
            B ^= E[b];
            const unsigned char* rk = &RK[b * Lp];
            for (int jj = 0; jj < Lp; jj++) R[jj] ^= rk[jj];
            cand++;
            int jj = 0;
            for (; jj < L; jj++) if (has(F[jj], R[jj])) break;
            if (jj < L) continue;
            passed++;
            if (push(B, sc, sn, sd)) { dfs(); pop(sc, sn, sd); }
        }
    }

    // The same filter solved instead of enumerated. Coordinates m in F_2^21 over the basis E_k = u_i ^ u_j of
    // {B : B(p,.) = 0}. Covered y gives psi_y(m) = B(y,.) reduced mod F_y (a linear map; bad iff 0). A tree node is
    // an affine space a + <K> and the constraints still live on it. Per constraint: v = psi(a), I = psi(<K>);
    // v not in I: satisfied on the whole space, dropped; I = 0 (so v = 0): nothing here passes; else rank rho of I.
    // Branch on a least-rho constraint: <K> = <P> + ker, with psi injective on <P>; the cosets a + c + ker with
    // psi(a + c) != 0 (2^rho - 1 of them) partition the passing part of a + <K>, and on each that constraint is
    // satisfied. Leaves (dim 0, every constraint checked) are exactly the m the Gray loop would pass, in another order.
    struct Tree {
        const unsigned char* T; int L;   // T[j * 384 + 128 c + x]: psi_j on bits 7c..7c+6 of m
        u64 E[21];
        long long tnodes = 0;
        unsigned psi(int j, u32 m) const { const unsigned char* t = T + j * 384; return t[m & 127] ^ t[128 + (m >> 7 & 127)] ^ t[256 + (m >> 14)]; }
    };

    void tree_leaf(const Tree& tr, u32 m) {
        u64 B = 0;
        for (u32 x = m; x; x &= x - 1) B ^= tr.E[_tzcnt_u32(x)];
        passed++;
        PS sc; int sn, sd;
        if (push(B, sc, sn, sd)) { dfs(); pop(sc, sn, sd); }
    }

    void tree_solve(Tree& tr, u32 a, const u32* K, int k, const unsigned char* live, int n) {
        tr.tnodes++;
        if (k == 0) {
            for (int i = 0; i < n; i++) if (!tr.psi(live[i], a)) return;
            tree_leaf(tr, a);
            return;
        }
        unsigned char nl[256]; int nn = 0;
        int best = -1, brho = 99;
        unsigned w[21];
        for (int i = 0; i < n; i++) {
            int j = live[i];
            if (brho == 1) { nl[nn++] = (unsigned char)j; continue; }   // a forced move is found: keep the rest unchecked
            unsigned v = tr.psi(j, a);
            unsigned piv[8] = {0}; int rho = 0;
            for (int t = 0; t < k; t++) {
                unsigned x = tr.psi(j, K[t]);
                for (int b = 7; b >= 0 && x; b--) if (x >> b & 1) { if (piv[b]) x ^= piv[b]; else { piv[b] = x; rho++; break; } }
            }
            for (int b = 7; b >= 0 && v; b--) if ((v >> b & 1) && piv[b]) v ^= piv[b];
            if (v) continue;                 // psi != 0 on all of a + <K>
            if (rho == 0) return;            // psi = 0 on all of it
            nl[nn++] = (unsigned char)j;
            if (rho < brho) { brho = rho; best = nn - 1; }
        }
        if (nn == 0) {                       // nothing left to fail: every point passes
            for (u32 g = 0; g < (1u << k); g++) {
                u32 m = a;
                for (int t = 0; t < k; t++) if (g >> t & 1) m ^= K[t];
                tree_leaf(tr, m);
            }
            return;
        }
        int j = nl[best];
        nl[best] = nl[--nn];
        // split <K> into P (independent images under psi_j) and ker
        unsigned pimg[8]; u32 pvec[8]; int np = 0;
        u32 ker[21]; int nk = 0;
        for (int t = 0; t < k; t++) {
            unsigned x = tr.psi(j, K[t]); u32 vec = K[t];
            for (int q = 0; q < np && x; q++) {
                int hb = 31 - (int)_lzcnt_u32(pimg[q]);
                if (x >> hb & 1) { x ^= pimg[q]; vec ^= pvec[q]; }
            }
            if (x) { pimg[np] = x; pvec[np] = vec; np++;
                     // keep pivots with distinct top bits, reduced against the new one
                     int hb = 31 - (int)_lzcnt_u32(x);
                     for (int q = 0; q < np - 1; q++) if (pimg[q] >> hb & 1) { pimg[q] ^= x; pvec[q] ^= vec; }
            } else ker[nk++] = vec;
        }
        for (u32 g = 0; g < (1u << np); g++) {
            u32 m = a;
            for (int q = 0; q < np; q++) if (g >> q & 1) m ^= pvec[q];
            if (!tr.psi(j, m)) continue;
            tree_solve(tr, m, ker, nk, nl, nn);
        }
    }

    void holes_tree(int p) {
        unsigned char ys[256];
        int L = 0;
        for (int y = 1; y < 256; y++) if (has(cov, y)) ys[L++] = (unsigned char)y;
        int u[7]; int pb[1] = {p};
        annihilator(pb, 1, u);
        Tree tr; tr.L = L;
        int k = 0;
        for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++, k++) {
            u64 B = 0;
            for (int a = 0; a < 8; a++) {
                int row = (par(u[i] & (1 << a)) ? u[j] : 0) ^ (par(u[j] & (1 << a)) ? u[i] : 0);
                B |= (u64)(row & 255) << (8 * a);
            }
            tr.E[k] = B;
        }
        std::vector<unsigned char> T((size_t)L * 384);
        for (int jj = 0; jj < L; jj++) {
            int y = ys[jj];
            unsigned piv[8] = {0};             // F_y = span of the generators' rows at y, echelon by top bit
            for (int g = 0; g < dim; g++) {
                u64 G = elems[1 << g]; unsigned x = 0;
                for (int a = 0; a < 8; a++) if (y >> a & 1) x ^= (G >> (8 * a)) & 255;
                for (int b = 7; b >= 0 && x; b--) if (x >> b & 1) { if (piv[b]) x ^= piv[b]; else { piv[b] = x; break; } }
            }
            unsigned char M[21];
            for (int kk = 0; kk < 21; kk++) {
                unsigned x = 0;
                for (int a = 0; a < 8; a++) if (y >> a & 1) x ^= (tr.E[kk] >> (8 * a)) & 255;
                for (int b = 7; b >= 0; b--) if ((x >> b & 1) && piv[b]) x ^= piv[b];
                M[kk] = (unsigned char)x;
            }
            unsigned char* t = &T[(size_t)jj * 384];
            for (int c = 0; c < 3; c++) {
                t[128 * c] = 0;
                for (int x = 1; x < 128; x++) t[128 * c + x] = t[128 * c + (x & (x - 1))] ^ M[7 * c + _tzcnt_u32(x)];
            }
        }
        tr.T = T.data();
        u32 K[21]; for (int t = 0; t < 21; t++) K[t] = 1u << t;
        unsigned char live[256]; for (int jj = 0; jj < L; jj++) live[jj] = (unsigned char)jj;
        tree_solve(tr, 0, K, 21, live, L);
        cand += tr.tnodes;
    }
};

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: cover <input> [maxlines|-1] [splitdim] [--old|--gray] [--p1all] [--only TAG] [--frange A B] [--fsample N SEED]\n"); return 2; }
    std::vector<char*> pos; int hm = 0; bool p1all = false; const char* only = nullptr; long fa = -1, fb = -1, fsample = 0; unsigned long long fseed = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--old")) hm = 2;
        else if (!strcmp(argv[i], "--reps") && i + 1 < argc) {   // five_orbits' output for the one configuration given by --only
            FILE* rf = fopen(argv[++i], "r"); if (!rf) { perror(argv[i]); return 1; }
            int d; std::array<unsigned, 5> k;
            while (fscanf(rf, "%d %x %x %x %x %x", &d, &k[0], &k[1], &k[2], &k[3], &k[4]) == 6) g_reps.push_back(k);
            fclose(rf);
            std::sort(g_reps.begin(), g_reps.end());
            if (g_reps.empty()) { fprintf(stderr, "empty reps file\n"); return 1; }
        }
        else if (!strcmp(argv[i], "--p1all")) p1all = true;   // --frange/--fsample keep the hole stage in phase 1: without it, cores below the split dimension never reach the frontier
        else if (!strcmp(argv[i], "--gray")) hm = 1;
        else if (!strcmp(argv[i], "--only") && i + 1 < argc) only = argv[++i];
        else if (!strcmp(argv[i], "--fsample") && i + 2 < argc) { fsample = atol(argv[i + 1]); fseed = strtoull(argv[i + 2], nullptr, 10); fa = 0; fb = 2000000000L; i += 2; }
        else if (!strcmp(argv[i], "--frange") && i + 2 < argc) { fa = atol(argv[i + 1]); fb = atol(argv[i + 2]); i += 2; }
        else pos.push_back(argv[i]);
    }
    argc = (int)pos.size() + 1;
    for (int i = 0; i < (int)pos.size(); i++) argv[i + 1] = pos[i];
    long maxl = argc > 2 ? atol(argv[2]) : -1;
    D4 = alt_nondeg(4); D6 = alt_nondeg(6);
    printf("nondegenerate alternating forms: dim 4 %zu (expect 28), dim 6 %zu (expect 13888)\n", D4.size(), D6.size());
    FILE* f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    std::vector<Cover> cv;
    static char line[1 << 16];
    while (fgets(line, sizeof line, f)) {
        if (maxl >= 0 && (long)cv.size() >= maxl) break;
        Cover c;
        if (parse(line, c) && (!only || c.tag == only)) cv.push_back(std::move(c));
    }
    fclose(f);
    // tasks: (cover, first form on solid 0) without seeds; (cover, -1) with seeds
    std::vector<std::pair<int, int>> tasks;
    for (int i = 0; i < (int)cv.size(); i++) {
        if (cv[i].seeds.empty()) for (int j = 0; j < 28; j++) tasks.push_back({i, j});
        else tasks.push_back({i, -1});
    }
    int splitdim = argc > 3 ? atoi(argv[3]) : -1;
    bool profile = fa >= 0;
    if (profile && splitdim < 0) { fprintf(stderr, "--frange needs a split dimension\n"); return 2; }
    printf("covers %zu, tasks %zu, threads %d, split %d, hole stage %s%s\n", cv.size(), tasks.size(), omp_get_max_threads(), splitdim,
           hm == 2 ? "old (lines x forms)" : hm == 1 ? "Gray code + filter" : "affine tree", profile ? ", profiling a frontier range (counts partial)" : ""); fflush(stdout);
    std::atomic<long long> nodes(0), sols(0), fulls(0), nonapn(0), badseed(0), done(0), nears(0), bydim[9], cand(0), passed(0);
    for (auto& b : bydim) b = 0;
    double t0 = omp_get_wtime();
    std::vector<std::pair<int, std::vector<u64>>> front;   // (cover, span generators)
    auto add_counts = [&](const Search& se) {
        nodes += se.nodes; sols += se.sols; fulls += se.fulls; nonapn += se.nonapn; nears += se.nears; cand += se.cand; passed += se.passed;
        for (int k = 0; k < 9; k++) bydim[k] += se.bydim[k];
    };
    #pragma omp parallel for schedule(dynamic, 1)
    for (long t = 0; t < (long)tasks.size(); t++) {
        Search se; se.C = &cv[tasks[t].first]; se.elems[0] = 0;
        std::vector<std::vector<u64>> fr;
        se.split = splitdim; se.frontier = &fr; se.hmode = hm; se.p1only = profile && !p1all;
        bool ok = true;
        PS sc; int sn, sd;
        for (u64 B : se.C->seeds) if (!se.push(B, sc, sn, sd)) { ok = false; break; }
        if (!ok) {
            badseed++;
            std::lock_guard<std::mutex> g(g_mu);
            printf("BADSEED %s\n", se.C->tag.c_str()); fflush(stdout);
        } else if (tasks[t].second >= 0) {
            if (se.push(se.C->sforms[0][tasks[t].second], sc, sn, sd)) se.dfs();
        } else se.dfs();
        add_counts(se);
        long long d = ++done;
        std::lock_guard<std::mutex> g(g_mu);
        for (auto& v : fr) front.push_back({tasks[t].first, std::move(v)});
        if (d % 2000 == 0) {
            printf("tasks %lld/%zu  %.0f s  nodes %lld  solutions %lld  full %lld\n", d, tasks.size(), omp_get_wtime() - t0, nodes.load(), sols.load(), fulls.load());
            fflush(stdout);
        }
    }
    if (splitdim >= 0) {
        std::sort(front.begin(), front.end());   // deterministic order
        std::vector<int> per(cv.size(), 0);
        for (auto& f : front) per[f.first]++;
        printf("phase 1 done %.0f s: frontier %zu items at dim %d\n", omp_get_wtime() - t0, front.size(), splitdim);
        if (profile) {
            fb = std::min<long>(fb, (long)front.size()); fa = std::min(fa, fb);
            printf("phase 1 counts: nodes %lld, by dimension", nodes.load());   // the reset below drops them from the totals
            for (int k = 1; k < 9; k++) printf(" %lld", bydim[k].load());
            printf(", near misses %lld, candidates %lld, passed %lld\n", nears.load(), cand.load(), passed.load());
            front = std::vector<std::pair<int, std::vector<u64>>>(front.begin() + fa, front.begin() + fb);
            nodes = 0; for (auto& b : bydim) b = 0; cand = 0; passed = 0;
            printf("profiling frontier items [%ld, %ld)\n", fa, fb);
        }
        if (fsample > 0) {   // fsample items per cover, chosen from the sorted frontier by a fixed seed, so every mode runs the same ones
            std::vector<std::pair<int, std::vector<u64>>> pick;
            std::mt19937_64 rng(fseed);
            for (size_t i = 0, j; i < front.size(); i = j) {
                for (j = i; j < front.size() && front[j].first == front[i].first; j++) {}
                std::vector<size_t> idx;
                for (size_t q = i; q < j; q++) idx.push_back(q);
                for (size_t q = idx.size(); q > 1; q--) std::swap(idx[q - 1], idx[rng() % q]);
                idx.resize(std::min<size_t>(idx.size(), fsample));
                std::sort(idx.begin(), idx.end());
                for (size_t q : idx) pick.push_back(front[q]);
            }
            front = std::move(pick);
            nodes = 0; for (auto& b : bydim) b = 0; cand = 0; passed = 0;
            printf("sampling %ld frontier items per cover (seed %llu): %zu items\n", fsample, (unsigned long long)fseed, front.size());
        }
        double t1 = omp_get_wtime();
        for (size_t i = 0; i < cv.size(); i++) if (per[i] >= 100) printf("  frontier %s %d\n", cv[i].tag.c_str(), per[i]);
        fflush(stdout);
        std::atomic<long long> fdone(0);
        long long step = std::max<long long>(1, (long long)front.size() / 200);
        #pragma omp parallel for schedule(dynamic, 1)
        for (long t = 0; t < (long)front.size(); t++) {
            Search se; se.C = &cv[front[t].first]; se.elems[0] = 0; se.hmode = hm;
            PS sc; int sn, sd;
            for (u64 B : front[t].second) if (!se.push(B, sc, sn, sd)) { fprintf(stderr, "frontier push failed\n"); exit(1); }
            se.dfs();
            add_counts(se);
            if (fsample > 0) {   // per item, for comparing modes
                std::string g; char hx[32];
                for (u64 b : front[t].second) { snprintf(hx, sizeof hx, " %07llx", (unsigned long long)to_coef(b)); g += hx; }
                std::string by; for (int k = 1; k < 9; k++) { snprintf(hx, sizeof hx, " %lld", se.bydim[k]); by += hx; }
                std::lock_guard<std::mutex> lk(g_mu);
                printf("ITEM %s gens%s nodes%s sols %lld full %lld nonapn %lld\n", se.C->tag.c_str(), g.c_str(), by.c_str(), se.sols, se.fulls, se.nonapn);
            }
            long long d = ++fdone;
            if (d % step == 0) {
                std::lock_guard<std::mutex> g(g_mu);
                printf("frontier %lld/%zu  %.0f s  nodes %lld  solutions %lld  full %lld\n", d, front.size(), omp_get_wtime() - t0, nodes.load(), sols.load(), fulls.load());
                fflush(stdout);
            }
        }
        if (profile) printf("phase 2: %zu items, %.2f s, %.3f s per item per thread\n", front.size(), omp_get_wtime() - t1,
                            front.empty() ? 0.0 : (omp_get_wtime() - t1) * omp_get_max_threads() / front.size());
    }
    printf("done: covers %zu, tasks %zu, nodes %lld, solutions %lld, full (dim < 8) %lld, dim 8 non-APN %lld, bad seeds %lld, %.0f s\n",
           cv.size(), tasks.size(), nodes.load(), sols.load(), fulls.load(), nonapn.load(), badseed.load(), omp_get_wtime() - t0);
    printf("nodes by span dimension 1..8:");
    for (int k = 1; k < 9; k++) printf(" %lld", bydim[k].load());
    printf("\nnear misses (every solid a radical, holes uncovered, below dim 8 or not): %lld\n", nears.load());
    printf("hole-stage candidates %lld, passing the filter (sent to push) %lld\n", cand.load(), passed.load());
    return 0;
}

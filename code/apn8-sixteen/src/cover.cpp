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
// separately). No symmetry is used: all 28 forms on the first solid are tried.
//
// usage: cover <input> [maxlines|-1] [splitdim]
// splitdim D (2026-10-31, for seven's p7-2222, whose dim-4 core held one thread for most of a day): phase 1 runs the
// tasks but stops at every node with all solids radicals, holes left and span dimension D, keeping its span as a
// frontier item; phase 2 runs the frontier items in parallel. Same tree, same counts (a frontier node is counted
// once, in phase 2).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <atomic>
#include <mutex>
#include <algorithm>
#include <intrin.h>
#include <omp.h>

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
                if (covered == all && ncov < 256) { nears++; report("NEAR"); }  // every solid a radical, holes left
                dfs(); pop(sc, sn, sd);
            }
            return;
        }
        int p = -1;
        for (int x = 1; x < 256; x++) if (!has(cov, x)) { p = x; break; }
        // every uncovered point is a hole now (all solids are radicals)
        for (int q = 1; q < 256; q++) {
            int t = p ^ q;
            if (q == p || q > t || !has(C->holes, q) || !has(C->holes, t) || has(cov, q) || has(cov, t)) continue;
            int lb[2] = {p, q}, u[6];
            annihilator(lb, 2, u);
            for (auto& d : D6) {
                u64 B = form_on(u, 6, d);
                if (push(B, sc, sn, sd)) { dfs(); pop(sc, sn, sd); }
            }
        }
    }
};

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: cover <input> [maxlines|-1] [splitdim]\n"); return 2; }
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
        if (parse(line, c)) cv.push_back(std::move(c));
    }
    fclose(f);
    // tasks: (cover, first form on solid 0) without seeds; (cover, -1) with seeds
    std::vector<std::pair<int, int>> tasks;
    for (int i = 0; i < (int)cv.size(); i++) {
        if (cv[i].seeds.empty()) for (int j = 0; j < 28; j++) tasks.push_back({i, j});
        else tasks.push_back({i, -1});
    }
    int splitdim = argc > 3 ? atoi(argv[3]) : -1;
    printf("covers %zu, tasks %zu, threads %d, split %d\n", cv.size(), tasks.size(), omp_get_max_threads(), splitdim); fflush(stdout);
    std::atomic<long long> nodes(0), sols(0), fulls(0), nonapn(0), badseed(0), done(0), nears(0), bydim[9];
    for (auto& b : bydim) b = 0;
    double t0 = omp_get_wtime();
    std::vector<std::pair<int, std::vector<u64>>> front;   // (cover, span generators)
    auto add_counts = [&](const Search& se) {
        nodes += se.nodes; sols += se.sols; fulls += se.fulls; nonapn += se.nonapn; nears += se.nears;
        for (int k = 0; k < 9; k++) bydim[k] += se.bydim[k];
    };
    #pragma omp parallel for schedule(dynamic, 1)
    for (long t = 0; t < (long)tasks.size(); t++) {
        Search se; se.C = &cv[tasks[t].first]; se.elems[0] = 0;
        std::vector<std::vector<u64>> fr;
        se.split = splitdim; se.frontier = &fr;
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
        for (size_t i = 0; i < cv.size(); i++) if (per[i] >= 100) printf("  frontier %s %d\n", cv[i].tag.c_str(), per[i]);
        fflush(stdout);
        std::atomic<long long> fdone(0);
        long long step = std::max<long long>(1, (long long)front.size() / 200);
        #pragma omp parallel for schedule(dynamic, 1)
        for (long t = 0; t < (long)front.size(); t++) {
            Search se; se.C = &cv[front[t].first]; se.elems[0] = 0;
            PS sc; int sn, sd;
            for (u64 B : front[t].second) if (!se.push(B, sc, sn, sd)) { fprintf(stderr, "frontier push failed\n"); exit(1); }
            se.dfs();
            add_counts(se);
            long long d = ++fdone;
            if (d % step == 0) {
                std::lock_guard<std::mutex> g(g_mu);
                printf("frontier %lld/%zu  %.0f s  nodes %lld  solutions %lld  full %lld\n", d, front.size(), omp_get_wtime() - t0, nodes.load(), sols.load(), fulls.load());
                fflush(stdout);
            }
        }
    }
    printf("done: covers %zu, tasks %zu, nodes %lld, solutions %lld, full (dim < 8) %lld, dim 8 non-APN %lld, bad seeds %lld, %.0f s\n",
           cv.size(), tasks.size(), nodes.load(), sols.load(), fulls.load(), nonapn.load(), badseed.load(), omp_get_wtime() - t0);
    printf("nodes by span dimension 1..8:");
    for (int k = 1; k < 9; k++) printf(" %lld", bydim[k].load());
    printf("\nnear misses (every solid a radical, holes uncovered, below dim 8 or not): %lld\n", nears.load());
    return 0;
}

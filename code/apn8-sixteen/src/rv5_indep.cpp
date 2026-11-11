// Reviewer's independent recomputation for apn/review-five-orbits.md (fresh Claude, 2026-11-06). Shares no code with cover.cpp,
// five_aut.py or five_orbits.cpp, and uses different methods on purpose:
//  - forms with radical R: nullspace of the linear system B(r, e_j) = 0 over the 28 coefficients, not annihilators and alt_nondeg;
//  - radicals: Gaussian elimination with combination tracking, not Gray-code scans;
//  - cores: every tuple (f_0..f_{k-1}), f_s one of the forms on solid s, whole span rechecked from scratch at each depth. A core is
//    exactly such a tuple, so there's no "lowest uncovered solid" order and no incremental state;
//  - key: order-free 64-bit hash of the span's nonzero elements, not reduced echelon form;
//  - group: all elements by BFS from aut.txt's generators, each checked to permute the solids;
//  - action: B'(i, j) = B(g e_i, g e_j) evaluated pointwise on rows, not a table on unit forms.
// Per configuration: enumerated cores == near file as sets; union of the images of the reps under the whole group == the cores.
// usage: rv5_indep <p5.txt> <aut.txt> <near dir> <reps dir> [only tag]
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <intrin.h>
#include <omp.h>
typedef unsigned long long u64;
typedef unsigned char u8;

static inline u64 mix(u64 x) { x += 0x9E3779B97F4A7C15ull; x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull; x = (x ^ (x >> 27)) * 0x94D049BB133111EBull; return x ^ (x >> 31); }
static inline int parity(unsigned x) { return __popcnt(x) & 1; }

static u64 form_from_coef(unsigned c) {   // 28-bit coefficient word, pairs i < j in lexicographic order -> row bytes
    u64 B = 0; int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) if (c >> k & 1) { B |= 1ull << (8 * i + j); B |= 1ull << (8 * j + i); }
    return B;
}

// nonzero radical points into pts (size 256), returns dim of the radical
static int rad_pts(u64 B, u8* pts, int* n) {
    u8 row[8], cmb[8];
    for (int i = 0; i < 8; i++) { row[i] = (u8)(B >> (8 * i)); cmb[i] = (u8)(1 << i); }
    int r = 0;
    for (int c = 0; c < 8; c++) {
        int p = -1;
        for (int i = r; i < 8; i++) if (row[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(row[p], row[r]); std::swap(cmb[p], cmb[r]);
        for (int i = r + 1; i < 8; i++) if (row[i] >> c & 1) { row[i] ^= row[r]; cmb[i] ^= cmb[r]; }
        r++;
    }
    int d = 8 - r, m = 1; u8 all[256]; all[0] = 0;
    for (int i = r; i < 8; i++) { for (int j = 0; j < m; j++) all[m + j] = all[j] ^ cmb[i]; m *= 2; }
    *n = m - 1;
    for (int i = 1; i < m; i++) pts[i - 1] = all[i];
    return d;
}

struct Cfg {
    std::string tag; int k = 0; long long aut_tag = 0;
    std::vector<std::vector<u8>> sol;       // each solid's 15 nonzero points
    signed char where[256];
    std::vector<std::vector<u64>> forms;    // per solid
    std::vector<u64> cores;                 // sorted span hashes
    long long bydim[9] = {0};
};

static std::vector<u64> forms_on(const std::vector<u8>& basis, const std::vector<u8>& pts) {
    // equations over 28 unknowns: for each basis r and each j, sum_i r_i c_{ij} = 0
    int idx[8][8], k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) idx[i][j] = idx[j][i] = k;
    std::vector<unsigned> eq;
    for (u8 r : basis) for (int j = 0; j < 8; j++) {
        unsigned e = 0;
        for (int i = 0; i < 8; i++) if (i != j && (r >> i & 1)) e ^= 1u << idx[i][j];
        eq.push_back(e);
    }
    std::vector<unsigned> piv; std::vector<int> pc;   // reduced rows, pivot columns
    for (unsigned e : eq) {
        for (size_t t = 0; t < piv.size(); t++) if (e >> pc[t] & 1) e ^= piv[t];
        if (!e) continue;
        int c = (int)_tzcnt_u32(e);
        for (size_t t = 0; t < piv.size(); t++) if (piv[t] >> c & 1) piv[t] ^= e;
        piv.push_back(e); pc.push_back(c);
    }
    std::vector<unsigned> ns;
    for (int f = 0; f < 28; f++) {
        if (std::find(pc.begin(), pc.end(), f) != pc.end()) continue;
        unsigned v = 1u << f;
        for (size_t t = 0; t < piv.size(); t++) if (piv[t] >> f & 1) v |= 1u << pc[t];
        ns.push_back(v);
    }
    if (ns.size() != 6) { fprintf(stderr, "nullspace dim %zu\n", ns.size()); exit(1); }
    std::vector<u64> out;
    for (unsigned m = 1; m < 64; m++) {
        unsigned c = 0;
        for (int i = 0; i < 6; i++) if (m >> i & 1) c ^= ns[i];
        u64 B = form_from_coef(c);
        u8 p[256]; int n; int d = rad_pts(B, p, &n);
        if (d != 4) continue;
        std::vector<u8> a(p, p + n), b = pts; std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
        if (a != b) { fprintf(stderr, "rank-4 form on a solid with another radical\n"); exit(1); }
        out.push_back(B);
    }
    if (out.size() != 28) { fprintf(stderr, "%zu forms\n", out.size()); exit(1); }
    return out;
}

// the whole span of chosen, from scratch: legal? fills span (0 first)
static bool legal(const Cfg& C, const u64* ch, int n, u64* span, int* sn) {
    int m = 1; span[0] = 0;
    for (int i = 0; i < n; i++) {
        bool in = false;
        for (int j = 0; j < m; j++) if (span[j] == ch[i]) in = true;
        if (in) continue;
        for (int j = 0; j < m; j++) span[m + j] = span[j] ^ ch[i];
        m *= 2;
    }
    *sn = m;
    bool used[256] = {false}; unsigned seen = 0;
    for (int e = 1; e < m; e++) {
        u8 p[256]; int c; int d = rad_pts(span[e], p, &c);
        if (d == 0) continue;
        if (d == 2) { for (int i = 0; i < c; i++) if (C.where[p[i]] >= 0) return false; }
        else if (d == 4) {
            int s = C.where[p[0]];
            if (s < 0 || (seen >> s & 1)) return false;
            for (int i = 0; i < c; i++) if (C.where[p[i]] != s) return false;
            seen |= 1u << s;
        } else return false;
        for (int i = 0; i < c; i++) { if (used[p[i]]) return false; used[p[i]] = true; }
    }
    return true;
}

static u64 span_hash(const u64* span, int m) { u64 h = 0; for (int i = 1; i < m; i++) h += mix(span[i]); return h; }

static void rec(Cfg& C, u64* ch, int depth) {
    u64 span[64]; int m;
    for (u64 f : C.forms[depth]) {
        ch[depth] = f;
        if (!legal(C, ch, depth + 1, span, &m)) continue;
        if (depth + 1 == C.k) { C.cores.push_back(span_hash(span, m)); int d = 0; while ((1 << d) < m) d++; C.bydim[d]++; }
        else rec(C, ch, depth + 1);
    }
}

static inline u64 act(u64 B, const std::array<u8, 8>& g) {
    u8 w[8];
    for (int i = 0; i < 8; i++) { u8 v = 0; for (int a = 0; a < 8; a++) if (g[i] >> a & 1) v ^= (u8)(B >> (8 * a)); w[i] = v; }
    u64 o = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) if (parity(w[i] & g[j])) { o |= 1ull << (8 * i + j); o |= 1ull << (8 * j + i); }
    return o;
}
static inline u8 app(const std::array<u8, 8>& g, int x) { u8 y = 0; for (int i = 0; i < 8; i++) if (x >> i & 1) y ^= g[i]; return y; }

int main(int argc, char** argv) {
    if (argc < 5) return 2;
    const char* only = argc > 5 ? argv[5] : nullptr;
    std::vector<Cfg> cf;
    {
        FILE* f = fopen(argv[1], "r"); static char line[1 << 16];
        while (fgets(line, sizeof line, f)) {
            char* t = strtok(line, " \r\n"); if (!t) continue;
            Cfg C; C.tag = t; C.k = atoi(strtok(nullptr, " \r\n"));
            const char* a = strstr(t, "-a"); C.aut_tag = a ? atoll(a + 2) : 0;
            memset(C.where, -1, sizeof C.where);
            for (int s = 0; s < C.k; s++) {
                std::vector<u8> b, pts;
                for (int i = 0; i < 4; i++) b.push_back((u8)strtoul(strtok(nullptr, " \r\n"), nullptr, 16));
                for (int m = 1; m < 16; m++) { u8 v = 0; for (int i = 0; i < 4; i++) if (m >> i & 1) v ^= b[i]; if (!v) { fprintf(stderr, "dep\n"); return 1; } pts.push_back(v); }
                for (u8 v : pts) { if (C.where[v] >= 0) { fprintf(stderr, "meet\n"); return 1; } C.where[v] = (signed char)s; }
                { std::vector<u8> q = pts; std::sort(q.begin(), q.end()); if (std::unique(q.begin(), q.end()) != q.end()) { fprintf(stderr, "dep\n"); return 1; } }
                C.sol.push_back(pts); C.forms.push_back(forms_on(b, pts));
            }
            if (!only || C.tag == only) cf.push_back(C);
        }
        fclose(f);
    }
    long long tcores = 0, treps = 0; int fails = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int ci = 0; ci < (int)cf.size(); ci++) {
        Cfg& C = cf[ci];
        u64 ch[8]; rec(C, ch, 0);
        std::sort(C.cores.begin(), C.cores.end());
        bool dupes = std::adjacent_find(C.cores.begin(), C.cores.end()) != C.cores.end();
        // near file
        std::vector<u64> near; near.reserve(C.cores.size());
        {
            std::string p = std::string(argv[3]) + "/" + C.tag + ".txt";
            FILE* f = fopen(p.c_str(), "r"); if (!f) { fprintf(stderr, "no %s\n", p.c_str()); exit(1); }
            char line[1024];
            while (fgets(line, sizeof line, f)) {
                if (strncmp(line, "NEAR ", 5)) continue;
                char* b = strstr(line, " basis ") + 7; u64 sp[64]; int m = 1; sp[0] = 0; char* q;
                while (true) { unsigned long c = strtoul(b, &q, 16); if (q == b) break; b = q; u64 B = form_from_coef((unsigned)c); for (int j = 0; j < m; j++) sp[m + j] = sp[j] ^ B; m *= 2; }
                near.push_back(span_hash(sp, m));
            }
            fclose(f);
        }
        std::sort(near.begin(), near.end());
        bool near_eq = near == C.cores;
        // group
        std::vector<std::array<u8, 8>> G; long long autline = -1; bool stab = true;
        {
            FILE* f = fopen(argv[2], "r"); char* line = (char*)malloc(1 << 16); std::vector<std::array<u8, 8>> gens; char* sv;
            while (fgets(line, 1 << 16, f)) {
                char* t = strtok_s(line, " \r\n", &sv); if (!t || C.tag != t) continue;
                autline = atoll(strtok_s(nullptr, " \r\n", &sv));
                while ((t = strtok_s(nullptr, " \r\n", &sv))) { std::array<u8, 8> g; for (int i = 0; i < 8; i++) { char h[3] = {t[2 * i], t[2 * i + 1], 0}; g[i] = (u8)strtol(h, nullptr, 16); } gens.push_back(g); }
            }
            fclose(f); free(line);
            std::array<u8, 8> id; for (int i = 0; i < 8; i++) id[i] = (u8)(1 << i);
            std::set<std::array<u8, 8>> seen{id}; G.push_back(id);
            for (size_t i = 0; i < G.size(); i++) for (auto& h : gens) {
                std::array<u8, 8> y; for (int j = 0; j < 8; j++) y[j] = app(h, G[i][j]);
                if (seen.insert(y).second) G.push_back(y);
            }
            for (auto& g : G) {
                bool inj[256] = {false}; for (int x = 0; x < 256; x++) inj[app(g, x)] = true;
                for (int x = 0; x < 256; x++) if (!inj[x]) stab = false;
                std::set<int> img;
                for (auto& s : C.sol) { int t = C.where[app(g, s[0])]; if (t < 0) { stab = false; continue; } img.insert(t); for (u8 v : s) if (C.where[app(g, v)] != t) stab = false; }
                if ((int)img.size() != C.k) stab = false;
            }
        }
        // reps: images under the whole group
        std::vector<char> mark(C.cores.size(), 0); long long nreps = 0, missing = 0, repeat_orbit = 0, repnotcore = 0, sum = 0; std::map<int, long long> repdim;
        {
            std::string p = std::string(argv[4]) + "/" + C.tag + ".txt";
            FILE* f = fopen(p.c_str(), "r"); if (!f) { fprintf(stderr, "no %s\n", p.c_str()); exit(1); }
            int d; unsigned w[5];
            while (fscanf(f, "%d %x %x %x %x %x", &d, &w[0], &w[1], &w[2], &w[3], &w[4]) == 6) {
                nreps++; repdim[d]++;
                u64 B[5]; for (int i = 0; i < d; i++) B[i] = form_from_coef(w[i]);
                { u64 sp[64]; int m = 1; sp[0] = 0; for (int i = 0; i < d; i++) { for (int j = 0; j < m; j++) sp[m + j] = sp[j] ^ B[i]; m *= 2; } if (!std::binary_search(C.cores.begin(), C.cores.end(), span_hash(sp, m))) repnotcore++; }
                bool first = true;
                for (auto& g : G) {
                    u64 sp[64]; int m = 1; sp[0] = 0;
                    for (int i = 0; i < d; i++) { u64 I = act(B[i], g); for (int j = 0; j < m; j++) sp[m + j] = sp[j] ^ I; m *= 2; }
                    u64 h = span_hash(sp, m);
                    auto it = std::lower_bound(C.cores.begin(), C.cores.end(), h);
                    if (it == C.cores.end() || *it != h) { missing++; continue; }
                    size_t i = it - C.cores.begin();
                    if (first) { if (mark[i]) repeat_orbit++; first = false; }
                    if (!mark[i]) { mark[i] = 1; sum++; }
                }
            }
            fclose(f);
        }
        bool ok = !dupes && near_eq && stab && (long long)G.size() == autline && autline == C.aut_tag && !missing && !repeat_orbit && !repnotcore && sum == (long long)C.cores.size();
        #pragma omp critical
        {
            printf("%s %s: enumerated cores %zu (dim<=3 %lld, dim 4 %lld, dim 5 %lld), hash dupes %d; near set equal %d (%zu); |G| %zu, line %lld, tag %lld, all permute solids %d; reps %lld (dim 4 %lld, dim 5 %lld), reps not cores %lld, images off the list %lld, reps in one orbit %lld, covered %lld\n",
                   ok ? "OK" : "FAIL", C.tag.c_str(), C.cores.size(), C.bydim[0] + C.bydim[1] + C.bydim[2] + C.bydim[3], C.bydim[4], C.bydim[5], (int)dupes, (int)near_eq, near.size(), G.size(), autline, C.aut_tag, (int)stab, nreps, repdim[4], repdim[5], repnotcore, missing, repeat_orbit, sum);
            fflush(stdout);
            tcores += C.cores.size(); treps += nreps; if (!ok) fails++;
        }
    }
    printf("total cores %lld reps %lld configurations failing %d\n", tcores, treps, fails);
    return fails ? 1 : 0;
}

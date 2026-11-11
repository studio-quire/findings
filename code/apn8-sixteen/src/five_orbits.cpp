// five_orbits <near file of one configuration> <aut.txt> <tag> <reps out>
// Orbits of Aut(P) on a configuration's cores (the span of its rank-4 forms: every solid a radical). Each S8 over P has exactly one
// core, g in Aut(P) sends S8 to an S8 over P and its core to the image core, and so only one core per orbit needs searching.
//
// Checks, any failure aborting: the near file has no core twice; the image of every core under every generator is in the file (a
// wrong action or a wrong group would fail this, and so would a list missing part of an orbit, though not one missing a whole orbit:
// that the list is complete is shown separately, apn/review-five-orbits.md); every orbit's size divides |Aut|; orbit sizes sum to the
// number of cores. Out: one line per orbit, the least key: `<dim> <reduced basis, 5 hex words>`.
// B acts as (g B)(x, y) = B(g x, g y); over a group that gives the same orbits as the inverse action.
// Rowen (Claude Opus 5.5), 2026-11-06.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <intrin.h>
typedef uint32_t u32;
typedef std::array<u32, 5> Key;
static inline u32 high_bit(u32 x) { unsigned long p; _BitScanReverse(&p, x); return 1u << p; }
static inline int low_index(u32 x) { unsigned long p; _BitScanForward(&p, x); return (int)p; }

static int PI[28], PJ[28], IDX[8][8];

static Key reduce(const u32* v, int n) {   // reduced echelon form, descending, zero padded
    u32 b[5]; int r = 0;
    for (int i = 0; i < n; i++) {
        u32 x = v[i];
        for (int k = 0; k < r; k++) if (x & high_bit(b[k])) x ^= b[k];
        if (!x) continue;
        u32 hb = high_bit(x);
        for (int k = 0; k < r; k++) if (b[k] & hb) b[k] ^= x;
        b[r++] = x;
    }
    std::sort(b, b + r, [](u32 a, u32 c) { return a > c; });
    Key o{}; for (int i = 0; i < r; i++) o[i] = b[i];
    return o;
}

struct Gen { u32 t[28]; };   // image of each unit form

static Gen make_gen(const unsigned char* g) {
    Gen G;
    for (int k = 0; k < 28; k++) {   // unit form B(e_p, e_q) = 1
        int p = PI[k], q = PJ[k]; u32 c = 0;
        for (int m = 0; m < 28; m++) {   // (gB)(e_i, e_j) = B(g e_i, g e_j)
            int gi = g[PI[m]], gj = g[PJ[m]];
            int v = ((gi >> p & 1) & (gj >> q & 1)) ^ ((gi >> q & 1) & (gj >> p & 1));
            if (v) c |= 1u << m;
        }
        G.t[k] = c;
    }
    return G;
}

static inline u32 act(const Gen& G, u32 c) { u32 o = 0; while (c) { int k = low_index(c); c &= c - 1; o ^= G.t[k]; } return o; }

int main(int argc, char** argv) {
    if (argc != 5) { fprintf(stderr, "usage: five_orbits <near> <aut.txt> <tag> <reps out>\n"); return 2; }
    int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) { PI[k] = i; PJ[k] = j; IDX[i][j] = k; }
    const char* tag = argv[3];
    std::vector<Gen> gens; long long aut = 0;
    {
        FILE* f = fopen(argv[2], "r"); if (!f) { perror(argv[2]); return 1; }
        static char line[1 << 16];
        while (fgets(line, sizeof line, f)) {
            char* p = strtok(line, " \r\n");
            if (!p || strcmp(p, tag)) continue;
            aut = atoll(strtok(nullptr, " \r\n"));
            while ((p = strtok(nullptr, " \r\n"))) {
                if (strlen(p) != 16) { fprintf(stderr, "bad generator\n"); return 1; }
                unsigned char g[8];
                for (int i = 0; i < 8; i++) { char h[3] = {p[2 * i], p[2 * i + 1], 0}; g[i] = (unsigned char)strtol(h, nullptr, 16); }
                gens.push_back(make_gen(g));
            }
        }
        fclose(f);
    }
    if (!aut) { fprintf(stderr, "no group for %s\n", tag); return 1; }
    std::vector<Key> keys; long long bydim[6] = {0};
    {
        FILE* f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
        static char line[1024];
        while (fgets(line, sizeof line, f)) {
            if (strncmp(line, "NEAR ", 5)) continue;
            char* b = strstr(line, " basis "); if (!b) { fprintf(stderr, "bad line\n"); return 1; }
            int d = atoi(strstr(line, " dim ") + 5);
            if (d > 5) { fprintf(stderr, "core of dim %d\n", d); return 1; }
            u32 v[5]; int n = 0; char* e = b + 7;
            while (n < 5) { char* q; unsigned long x = strtoul(e, &q, 16); if (q == e) break; v[n++] = (u32)x; e = q; }
            if (n != d) { fprintf(stderr, "dim %d, %d words\n", d, n); return 1; }
            keys.push_back(reduce(v, n)); bydim[d]++;
        }
        fclose(f);
    }
    std::sort(keys.begin(), keys.end());
    size_t n0 = keys.size();
    if (std::adjacent_find(keys.begin(), keys.end()) != keys.end()) { fprintf(stderr, "%s: a core listed twice\n", tag); return 1; }
    std::vector<char> seen(n0, 0);
    std::vector<size_t> st;
    std::map<long long, long long> hist;
    long long orbits = 0, total = 0; std::vector<long long> repdim(6, 0);
    FILE* o = fopen(argv[4], "w"); if (!o) { perror(argv[4]); return 1; }
    for (size_t i = 0; i < n0; i++) {
        if (seen[i]) continue;
        long long size = 0; seen[i] = 1; st.push_back(i);
        while (!st.empty()) {
            size_t c = st.back(); st.pop_back(); size++;
            for (auto& G : gens) {
                u32 w[5]; int n = 0;
                for (u32 x : keys[c]) if (x) w[n++] = act(G, x);
                Key y = reduce(w, n);
                auto it = std::lower_bound(keys.begin(), keys.end(), y);
                if (it == keys.end() || *it != y) { fprintf(stderr, "%s: image of a core is not a core\n", tag); return 1; }
                size_t j = it - keys.begin();
                if (!seen[j]) { seen[j] = 1; st.push_back(j); }
            }
        }
        if (aut % size) { fprintf(stderr, "%s: orbit of %lld, |Aut| %lld\n", tag, size, aut); return 1; }
        hist[size]++; orbits++; total += size;
        int d = 0; for (u32 x : keys[i]) if (x) d++;
        repdim[d]++;
        fprintf(o, "%d %07x %07x %07x %07x %07x\n", d, keys[i][0], keys[i][1], keys[i][2], keys[i][3], keys[i][4]);
    }
    fclose(o);
    if (total != (long long)n0) { fprintf(stderr, "%s: orbit sum\n", tag); return 1; }
    printf("%s |Aut| %lld generators %zu: cores %zu (dim 4 %lld, dim 5 %lld), orbits %lld (dim 4 %lld, dim 5 %lld), sizes", tag, aut, gens.size(), n0,
           bydim[4], bydim[5], orbits, repdim[4], repdim[5]);
    for (auto& h : hist) printf(" %lld:%lld", h.first, h.second);
    printf("; closed, sum OK\n");
    return 0;
}

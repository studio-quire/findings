// Reviewer's group checks, sharing no code with fivecheck_classes.cpp.
//   gamma            for every class of 5 and of 4 solids: |Gamma| by my own enumeration, against configs.txt
//   orbits           the dim-4 tuples: orbits under my Gamma by union-find; compare with k4all/k4orbit/k4 lists and my own list
//   closure <r> <cfg> <outfile...>   S8s in search B's output files: the set must be closed under the stabiliser of the tuple
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <string>
#include <numeric>
#include <algorithm>
#include <intrin.h>
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8;

static int par(u32 x) { return __popcnt(x) & 1; }
static int matv(int T, int x) { int y = 0; for (int i = 0; i < 4; i++) if (par((T >> (4 * i)) & 15 & x)) y |= 1 << i; return y; }
typedef array<u8, 256> Set; // membership
struct Lin { u8 c[8]; int ap(int v) const { int r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= c[i]; return r; } };

static vector<int> GL4;
static void gl4() {
    for (int m = 0; m < 65536; m++) {
        // columns b_i = bits 4i..: invertible iff the 15 combinations are nonzero
        bool ok = true;
        for (int s = 1; s < 16 && ok; s++) { int x = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) x ^= (m >> (4 * i)) & 15; if (!x) ok = false; }
        if (ok) GL4.push_back(m);
    }
    if (GL4.size() != 20160) { puts("GL4"); exit(1); }
}
struct Solids { int r; vector<vector<int>> pts; vector<Set> mem; };
static Solids mk(int r, int T4, int T5) {
    Solids S; S.r = r; int Ts[3] = {0x8421, T4, T5};
    for (int i = 0; i < r; i++) {
        vector<int> p; Set m; m.fill(0);
        for (int x = 1; x < 16; x++) { int v = i == 0 ? x : i == 1 ? x << 4 : x | (matv(Ts[i - 2], x) << 4); p.push_back(v); m[v] = 1; }
        S.pts.push_back(p); S.mem.push_back(m);
    }
    return S;
}
// All invertible A permuting the solids. A sends R1 onto W_i through a basis choice, R2 into W_j so that graph(I) lands in W_k.
static vector<Lin> gammaOf(const Solids& S) {
    vector<Lin> out; int r = S.r;
    vector<int> basisVecs;
    for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) for (int k = 0; k < r; k++) {
        if (i == j || i == k || j == k) continue;
        // fixed basis of W_i
        int b[4], n = 0; { set<int> sp{0}; for (int v : S.pts[i]) if (!sp.count(v)) { b[n++] = v; set<int> t = sp; for (int x : sp) t.insert(x ^ v); sp = t; } }
        // for each u in W_i, the v in W_j with u + v in W_k
        int part[256]; memset(part, 0, sizeof part);
        for (int u : S.pts[i]) { int cnt = 0; for (int v : S.pts[j]) if (S.mem[k][u ^ v]) { part[u] = v; cnt++; } if (cnt != 1) { puts("not complementary"); exit(1); } }
        for (int g : GL4) {
            int im[4]; bool ok = true;
            Lin A;
            for (int q = 0; q < 4; q++) {
                int col = (g >> (4 * q)) & 15, u = 0;
                for (int t = 0; t < 4; t++) if (col >> t & 1) u ^= b[t];
                im[q] = u; A.c[q] = (u8)u; A.c[4 + q] = (u8)part[u];
            }
            (void)im;
            // the other solids must go to other solids
            int used = (1 << i) | (1 << j) | (1 << k);
            for (int t = 3; t < r && ok; t++) {
                int hit = -1;
                for (int s = 0; s < r; s++) if (!(used >> s & 1)) { bool all = true; for (int v : S.pts[t]) if (!S.mem[s][A.ap(v)]) { all = false; break; } if (all) hit = s; }
                if (hit < 0) ok = false; else used |= 1 << hit;
            }
            if (ok) out.push_back(A);
        }
    }
    return out;
}
// forms as u64 row bytes, as in the output files
static u64 pull(u64 f, const Lin& A) { // (x, y) -> f(Ax, Ay)
    u8 rowA[8];
    for (int i = 0; i < 8; i++) { int v = A.c[i], r = 0; for (int t = 0; t < 8; t++) if (v >> t & 1) r ^= (int)((f >> (8 * t)) & 255); rowA[i] = (u8)r; }
    u64 o = 0;
    for (int i = 0; i < 8; i++) { int row = 0; for (int j = 0; j < 8; j++) if (par(rowA[i] & A.c[j])) row |= 1 << j; o |= (u64)row << (8 * i); }
    return o;
}
static vector<u64> canon(vector<u64> b) {
    vector<u64> rr; vector<int> pv;
    for (u64 x : b) {
        for (size_t t = 0; t < rr.size(); t++) if (x >> pv[t] & 1) x ^= rr[t];
        if (!x) continue;
        unsigned long q; _BitScanForward64(&q, x);
        for (size_t t = 0; t < rr.size(); t++) if (rr[t] >> q & 1) rr[t] ^= x;
        rr.push_back(x); pv.push_back((int)q);
    }
    sort(rr.begin(), rr.end()); return rr;
}
struct Cfg { int r, T4, T5; long long om, gam; };
static vector<Cfg> readCfg() {
    vector<Cfg> v; FILE* f = fopen("runs/check5/configs.txt", "r"); Cfg c; unsigned a, b;
    while (fscanf(f, "%d %x %x %lld %lld", &c.r, &a, &b, &c.om, &c.gam) == 5) { c.T4 = a; c.T5 = b; v.push_back(c); }
    fclose(f); return v;
}
struct UF { vector<int> p; UF(int n) : p(n) { iota(p.begin(), p.end(), 0); } int f(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; } void u(int a, int b) { p[f(a)] = f(b); } };

int main(int argc, char** argv) {
    gl4();
    string mode = argv[1];
    auto cf = readCfg();
    if (mode == "gamma") {
        long long sum = 0; int bad = 0;
        for (auto& c : cf) {
            auto G = gammaOf(mk(c.r, c.T4, c.T5));
            printf("r %d T %04x %04x |Gamma| mine %zu configs %lld %s\n", c.r, c.T4, c.T5, G.size(), c.gam, (long long)G.size() == c.gam ? "same" : "DIFFERENT");
            if ((long long)G.size() != c.gam) bad++;
            if (c.r == 5) sum += c.om;
        }
        printf("sum of |Omega| over r=5 classes %lld, different %d\n", sum, bad);
        return 0;
    }
    if (mode == "orbits") {
        // forms table for indexing
        map<u64, array<int, 2>> id[22]; vector<array<u64, 28>> tbl[22];
        for (int i = 0; i < 22; i++) tbl[i].resize(5);
        { FILE* f = fopen("runs/check5/review/forms.txt", "r"); int ci, i, d; unsigned long long x; while (fscanf(f, "%d %d %d %llx", &ci, &i, &d, &x) == 4) { id[ci][x] = {i, d}; tbl[ci][i][d] = x; } fclose(f); }
        map<int, vector<array<int, 5>>> mine;
        { FILE* f = fopen("runs/check5/review/ref_k4.txt", "r"); int c; array<int, 5> t; while (fscanf(f, "%d %d %d %d %d %d", &c, &t[0], &t[1], &t[2], &t[3], &t[4]) == 6) mine[c].push_back(t); fclose(f); }
        vector<Cfg> c5; for (auto& c : cf) if (c.r == 5) c5.push_back(c);
        int totT = 0, totO = 0, disagree = 0;
        for (auto& kv : mine) {
            int ci = kv.first; auto& T = kv.second; sort(T.begin(), T.end());
            map<array<int, 5>, int> ix; for (size_t i = 0; i < T.size(); i++) ix[T[i]] = (int)i;
            auto G = gammaOf(mk(5, c5[ci].T4, c5[ci].T5));
            UF uf((int)T.size());
            for (auto& A : G) for (size_t i = 0; i < T.size(); i++) {
                array<int, 5> im{-1, -1, -1, -1, -1};
                for (int s = 0; s < 5; s++) { auto it = id[ci].find(pull(tbl[ci][s][T[i][s]], A)); if (it == id[ci].end()) { puts("image not a form of the rep"); return 1; } im[it->second[0]] = it->second[1]; }
                auto f = ix.find(im); if (f == ix.end()) { puts("image tuple not in the list: not closed"); return 1; }
                uf.u((int)i, f->second);
            }
            // B's lists
            char p[128]; snprintf(p, sizeof p, "runs/check5/k4all_c%02d.txt", ci);
            vector<array<int, 5>> theirs; { FILE* f = fopen(p, "r"); array<int, 5> t; while (f && fscanf(f, "%d %d %d %d %d", &t[0], &t[1], &t[2], &t[3], &t[4]) == 5) theirs.push_back(t); if (f) fclose(f); }
            snprintf(p, sizeof p, "runs/check5/k4orbit_c%02d.txt", ci);
            vector<int> orb; { FILE* f = fopen(p, "r"); int x; while (f && fscanf(f, "%d", &x) == 1) orb.push_back(x); if (f) fclose(f); }
            snprintf(p, sizeof p, "runs/check5/k4_c%02d.txt", ci);
            vector<array<int, 5>> reps; { FILE* f = fopen(p, "r"); array<int, 5> t; while (f && fscanf(f, "%d %d %d %d %d", &t[0], &t[1], &t[2], &t[3], &t[4]) == 5) reps.push_back(t); if (f) fclose(f); }
            bool sameList = theirs == T;
            // same partition?
            map<int, int> a2b; map<int, int> b2a; bool samePart = orb.size() == T.size();
            for (size_t i = 0; i < T.size() && samePart; i++) {
                int a = uf.f((int)i), b = orb[i];
                if (a2b.count(a) && a2b[a] != b) samePart = false;
                if (b2a.count(b) && b2a[b] != a) samePart = false;
                a2b[a] = b; b2a[a == a ? b : b] = a;
            }
            set<int> roots; for (size_t i = 0; i < T.size(); i++) roots.insert(uf.f((int)i));
            set<int> repRoots; bool repsIn = true; for (auto& t : reps) { auto f = ix.find(t); if (f == ix.end()) repsIn = false; else repRoots.insert(uf.f(f->second)); }
            bool oneEach = repsIn && repRoots == roots && reps.size() == roots.size();
            printf("class %d |Gamma| %zu tuples %zu orbits %zu; list same %d, partition same %d, one rep per orbit %d\n", ci, G.size(), T.size(), roots.size(), (int)sameList, (int)samePart, (int)oneEach);
            totT += (int)T.size(); totO += (int)roots.size(); if (!sameList || !samePart || !oneEach) disagree++;
        }
        printf("tuples %d orbits %d classes disagreeing %d\n", totT, totO, disagree);
        return 0;
    }
    if (mode == "closure") {
        int r = atoi(argv[2]), cfg = atoi(argv[3]);
        vector<Cfg> cr; for (auto& c : cf) if (c.r == r) cr.push_back(c);
        auto S = mk(r, cr[cfg].T4, cr[cfg].T5);
        auto G = gammaOf(S);
        printf("|Gamma| %zu\n", G.size());
        long long files = 0, withS8 = 0, total = 0, checked = 0, broken = 0;
        for (int a = 4; a < argc; a++) {
            FILE* f = fopen(argv[a], "r"); char line[2048];
            vector<vector<u64>> s8;
            while (fgets(line, sizeof line, f)) {
                if (strncmp(line, "S8 ", 3)) continue;
                char* p = strchr(line, ':') + 1; vector<u64> b; for (int i = 0; i < 8; i++) { char* e; b.push_back(strtoull(p, &e, 16)); p = e; }
                s8.push_back(canon(b));
            }
            fclose(f); files++;
            if (s8.empty()) continue;
            withS8++; total += (long long)s8.size();
            // the r rank-4 forms of the first S8 (the tuple is common to all S8 of the file)
            auto rank4 = [&](const vector<u64>& b) {
                vector<u64> out;
                for (int m = 1; m < 256; m++) {
                    u64 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= b[i];
                    u8 rows[8]; for (int i = 0; i < 8; i++) rows[i] = (u8)(x >> (8 * i));
                    int cnt = 0; for (int v = 1; v < 256; v++) { int c = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) c ^= rows[i]; if (!c) cnt++; }
                    if (cnt == 15) out.push_back(x);
                }
                sort(out.begin(), out.end()); return out;
            };
            vector<u64> tup = rank4(s8[0]);
            set<vector<u64>> all(s8.begin(), s8.end());
            for (auto& b : s8) if (rank4(b) != tup) printf("file %s: differing rank-4 forms\n", argv[a]);
            long long stab = 0, misses = 0;
            for (auto& A : G) {
                vector<u64> im; for (u64 x : tup) im.push_back(pull(x, A)); sort(im.begin(), im.end());
                if (im != tup) continue;
                stab++;
                for (auto& b : s8) {
                    vector<u64> ib; for (u64 x : b) ib.push_back(pull(x, A));
                    checked++;
                    if (!all.count(canon(ib))) { misses++; broken++; }
                }
            }
            printf("%s: %zu S8, tuple stabiliser %lld, images outside the found set %lld\n", argv[a], s8.size(), stab, misses);
        }
        printf("files %lld with S8 %lld total S8 %lld images checked %lld broken %lld\n", files, withS8, total, checked, broken);
        return 0;
    }
    return 1;
}

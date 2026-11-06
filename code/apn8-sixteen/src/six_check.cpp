// six_check.cpp -- independent check of Open Problem VI.1 (arXiv 2510.12008) for 6 rank-4 components:
// [0^194, 2^55, 4^6]. Subagent (Claude Opus 5.5) for Rowen, 2026-11-03 ("search B" for i = 11, six rank-4). Written
// without reading pspread*.cpp, cover*.cpp, cover_prep.py, cover_chunks*.cmd, cover_tally.py, holes_tally.py, README's
// "Seven"/"Six" parts, runs/i10/, runs/i11/, notes/ or the git log.
// Built on seven_check.cpp (search B for 7 rank-4), generalised to 6 members. New here: the linear search (`lin`, its own
// off-D phase and D phase, sharing no primitive with the point branching: forms of a solid from a complement basis and a
// 4x4 Pfaffian, legality as "g a not in S a" for every point a covered by S, K_a by Gray code with byte lanes), and
// controls planted so that the planted space is NOT reached in the off-D phase (see ctlplant).
//
// Coordinates: v = (x, y), x = low nibble, y = high nibble. X = {y = 0}, Y = {x = 0}, Z = graph(I) = {y = x}.
// A 4-space complementary to Y is graph(N) = {(x, Nx)}; complementary to X iff N invertible, to Z iff N + I is.
// V5824 = {N : N, N + I invertible}; two graphs are disjoint iff N + N' is invertible ("adjacent").
//
// parta   : every 3-clique of V5824 (= 6-partial spread with X, Y, Z among its members) with minimal class c, up to
//           conjugation by the minimal-class trick, each canonicalised.
// merge   : the per-class files merged; per-orbit identity; per-orbit extension invariants; the link to seven-check's
//           3731 orbits (adding a solid), orbit by orbit.
// configs : hole sets of the orbits up to GL(8,2) (Kurz Table 9.3, r = 6: 341 and 0), exact covers of each,
//           independent sources; writes the configuration list.
// search  : per configuration, point-branching over the forms (radicals fixed); D phase by point branching or sweep.
// lin     : the linear search (own off-D and D phases).
// ctlbatch / ctlrot / ctlplant : known functions through the same search paths.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <mutex>
#include <set>
#include <map>
#include <chrono>
#include <functional>
#include <array>
#include <unordered_map>
#include <climits>
#include <random>
#include <atomic>
#include <memory>
#include <intrin.h>
#include <immintrin.h>
#include <omp.h>
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16; typedef uint8_t u8;

static inline int par(u32 v) { return __popcnt(v) & 1; }
static inline u8 rowof(u64 f, int i) { return (u8)(f >> (8 * i)); }
static inline u8 applyf(u64 f, u8 v) { u8 r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= rowof(f, i); return r; }
static inline int evalf(u64 f, u8 u, u8 w) { return par(u & applyf(f, w)); }

static int rank8(u64 f) {
    u8 r[8]; for (int i = 0; i < 8; i++) r[i] = rowof(f, i);
    int rk = 0;
    for (int c = 0; c < 8; c++) {
        int p = -1; for (int i = rk; i < 8; i++) if (r[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        u8 t = r[p]; r[p] = r[rk]; r[rk] = t;
        for (int i = rk + 1; i < 8; i++) if (r[i] >> c & 1) r[i] ^= t;
        rk++;
    }
    return rk;
}

struct Mask { u64 w[4]; bool operator<(const Mask& o) const { return memcmp(w, o.w, sizeof w) < 0; } bool operator==(const Mask& o) const { return !memcmp(w, o.w, sizeof w); } };
static inline void mclear(Mask& m) { m.w[0] = m.w[1] = m.w[2] = m.w[3] = 0; }
static inline void mset(Mask& m, int v) { m.w[v >> 6] |= 1ULL << (v & 63); }
static inline bool mget(const Mask& m, int v) { return m.w[v >> 6] >> (v & 63) & 1; }
static inline bool mdisj(const Mask& a, const Mask& b) { return ((a.w[0] & b.w[0]) | (a.w[1] & b.w[1]) | (a.w[2] & b.w[2]) | (a.w[3] & b.w[3])) == 0; }
static inline bool msub(const Mask& a, const Mask& b) { return ((a.w[0] & ~b.w[0]) | (a.w[1] & ~b.w[1]) | (a.w[2] & ~b.w[2]) | (a.w[3] & ~b.w[3])) == 0; }
static inline void mor(Mask& a, const Mask& b) { for (int i = 0; i < 4; i++) a.w[i] |= b.w[i]; }
static inline int mcount(const Mask& m) { return (int)(__popcnt64(m.w[0]) + __popcnt64(m.w[1]) + __popcnt64(m.w[2]) + __popcnt64(m.w[3])); }
static inline bool meq(const Mask& a, const Mask& b) { return !memcmp(&a, &b, sizeof(Mask)); }
static Mask mcompl(const Mask& m) { Mask r; for (int i = 0; i < 4; i++) r.w[i] = ~m.w[i]; r.w[0] &= ~1ULL; return r; } // nonzero points not in m
static vector<u8> mpoints(const Mask& m) { vector<u8> p; for (int v = 1; v < 256; v++) if (mget(m, v)) p.push_back((u8)v); return p; }

// radical by the Gray walk over all 255 vectors (eight-check's; kept for the selftest's comparison)
static int radical_walk(u64 f, Mask& m) {
    mclear(m); int cnt = 0; u8 cur = 0, v = 0;
    for (u32 i = 1; i < 256; i++) { int b = _tzcnt_u32(i); v ^= (u8)(1u << b); cur ^= rowof(f, b); if (!cur) { mset(m, v); cnt++; } }
    return cnt;
}
// radical as the null space of the rows: reduce the rows, carrying for each the combination of original rows that
// gives it; the combinations of rows that reduce to zero span the radical (f symmetric). Seven-check: the Gray walk
// above dominated the search's time.
static int radical(u64 f, Mask& m) {
    u8 row[8], cmb[8]; for (int i = 0; i < 8; i++) { row[i] = rowof(f, i); cmb[i] = (u8)(1 << i); }
    int rk = 0;
    for (int c = 0; c < 8 && rk < 8; c++) {
        int p = -1; for (int i = rk; i < 8; i++) if (row[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        swap(row[p], row[rk]); swap(cmb[p], cmb[rk]);
        for (int i = rk + 1; i < 8; i++) if (row[i] >> c & 1) { row[i] ^= row[rk]; cmb[i] ^= cmb[rk]; }
        rk++;
    }
    mclear(m); int k = 8 - rk; u8 pts[256]; pts[0] = 0; int n = 1;
    for (int i = rk; i < 8; i++) { for (int j = 0; j < n; j++) pts[n + j] = pts[j] ^ cmb[i]; n *= 2; }
    for (int j = 1; j < n; j++) mset(m, pts[j]);
    return (1 << k) - 1;
}

// ---- forms with a prescribed radical: pullbacks of nondegenerate forms along functionals u_1..u_m
static u64 pairform(u8 u, u8 w) {
    u64 f = 0;
    for (int i = 0; i < 8; i++) { u8 r = 0; if (u >> i & 1) r ^= w; if (w >> i & 1) r ^= u; f |= (u64)r << (8 * i); }
    return f;
}
static vector<u32> nondeg_codes(int m) {
    int np = m * (m - 1) / 2; vector<u32> out;
    for (u32 c = 0; c < (1u << np); c++) {
        u8 rows[8] = { 0 }; int p = 0;
        for (int k = 0; k < m; k++) for (int l = k + 1; l < m; l++) { if (c >> p & 1) { rows[k] |= 1 << l; rows[l] |= 1 << k; } p++; }
        int rk = 0;
        for (int col = 0; col < m; col++) {
            int piv = -1; for (int i = rk; i < m; i++) if (rows[i] >> col & 1) { piv = i; break; }
            if (piv < 0) continue;
            swap(rows[piv], rows[rk]);
            for (int i = 0; i < m; i++) if (i != rk && (rows[i] >> col & 1)) rows[i] ^= rows[rk];
            rk++;
        }
        if (rk == m) out.push_back(c);
    }
    return out;
}
static vector<u32> ND4, ND6;
static void forms_from(const u8* u, int m, const vector<u32>& codes, vector<u64>& out) {
    u64 Q[15]; int p = 0;
    for (int k = 0; k < m; k++) for (int l = k + 1; l < m; l++) Q[p++] = pairform(u[k], u[l]);
    for (u32 c : codes) { u64 f = 0; u32 cc = c; while (cc) { int b = _tzcnt_u32(cc); cc &= cc - 1; f ^= Q[b]; } out.push_back(f); }
}

// ---- 4x4 matrices: u16, row r in bits 4r..4r+3
static inline u8 mv4(u16 N, u8 x) { u8 r = 0; for (int i = 0; i < 4; i++) r |= par(((N >> (4 * i)) & 15) & x) << i; return r; }
static inline u16 mm4(u16 A, u16 B) {
    u16 C = 0;
    for (int r = 0; r < 4; r++) { u8 row = 0; for (int k = 0; k < 4; k++) if (A >> (4 * r + k) & 1) row ^= (B >> (4 * k)) & 15; C |= (u16)row << (4 * r); }
    return C;
}
static bool inv4(u16 N) { for (u8 x = 1; x < 16; x++) if (!mv4(N, x)) return false; return true; }
static u16 fromcols4(const u8 col[4]) { u16 N = 0; for (int i = 0; i < 4; i++) for (int r = 0; r < 4; r++) if (col[i] >> r & 1) N |= 1 << (4 * r + i); return N; }
static const u16 ID4 = 0x8421;

static vector<u16> GL4; static vector<u16> GL4inv;
static int gl4idx[65536];
static u16 inverse4(u16 N) { return GL4inv[gl4idx[N]]; }
// ---- V5824 and its graph
static const int NW = 91; // 5824 / 64
struct BS { u64 w[NW]; };
static vector<u16> VN;              // vertex -> matrix
static vector<int> cls;             // vertex -> conjugacy class index
static vector<int> classSize, classRep;
static int vid[65536];
static vector<BS> adj;
static inline void bset(BS& b, int i) { b.w[i >> 6] |= 1ULL << (i & 63); }
static inline bool bget(const BS& b, int i) { return b.w[i >> 6] >> (i & 63) & 1; }
static inline int bcount(const BS& b) { int c = 0; for (int i = 0; i < NW; i++) c += (int)__popcnt64(b.w[i]); return c; }
static inline void band(BS& o, const BS& a, const BS& b) { for (int i = 0; i < NW; i++) o.w[i] = a.w[i] & b.w[i]; }
static inline bool bempty(const BS& b) { for (int i = 0; i < NW; i++) if (b.w[i]) return false; return true; }

static void init_V() {
    memset(gl4idx, -1, sizeof gl4idx);
    for (u32 N = 0; N < 65536; N++) if (inv4((u16)N)) GL4.push_back((u16)N);
    GL4inv.resize(GL4.size()); for (size_t i = 0; i < GL4.size(); i++) gl4idx[GL4[i]] = (int)i;
    for (size_t i = 0; i < GL4.size(); i++) { u16 g = GL4[i], p = g, prev = ID4; while (p != ID4) { prev = p; p = mm4(p, g); } GL4inv[i] = (g == ID4) ? ID4 : prev; if (mm4(GL4inv[i], g) != ID4) { fprintf(stderr, "inverse\n"); exit(1); } }
    vector<u16> raw; for (u16 N : GL4) if (inv4(N ^ ID4)) raw.push_back(N);
    unordered_map<u16, int> c0; int nc = 0;
    for (u16 N : raw) if (!c0.count(N)) { for (size_t i = 0; i < GL4.size(); i++) c0[mm4(mm4(GL4[i], N), GL4inv[i])] = nc; nc++; }
    vector<pair<int, u16>> srt; for (u16 N : raw) srt.push_back({ c0[N], N }); sort(srt.begin(), srt.end());
    memset(vid, -1, sizeof vid);
    classSize.assign(nc, 0); classRep.assign(nc, -1);
    for (auto& p : srt) { int v = (int)VN.size(); VN.push_back(p.second); cls.push_back(p.first); vid[p.second] = v; classSize[p.first]++; if (classRep[p.first] < 0) classRep[p.first] = v; }
    adj.assign(VN.size(), BS{});
    for (size_t i = 0; i < VN.size(); i++) for (size_t j = 0; j < VN.size(); j++) if (i != j && inv4(VN[i] ^ VN[j])) bset(adj[i], (int)j);
}

// ---- conjugation tables: conjTab[h][v] = vertex of h N_v h^-1; toRep[v] = some h with h N_v h^-1 = class rep;
// cent[c] = the centraliser of the class-c rep (as GL4 indices)
static vector<vector<u16>> conjTab; static vector<int> toRep; static vector<vector<int>> cent;
static void init_conj() {
    int n = (int)VN.size();
    conjTab.assign(GL4.size(), vector<u16>(n));
#pragma omp parallel for
    for (int h = 0; h < (int)GL4.size(); h++) for (int v = 0; v < n; v++) conjTab[h][v] = (u16)vid[mm4(mm4(GL4[h], VN[v]), GL4inv[h])];
    toRep.assign(n, -1); cent.assign(classSize.size(), {});
    for (int h = 0; h < (int)GL4.size(); h++) for (int v = 0; v < n; v++) if (toRep[v] < 0 && conjTab[h][v] == classRep[cls[v]]) toRep[v] = h;
    for (int c = 0; c < (int)classSize.size(); c++) for (int h = 0; h < (int)GL4.size(); h++) if (conjTab[h][classRep[c]] == classRep[c]) cent[c].push_back(h);
    for (int v = 0; v < n; v++) if (toRep[v] < 0) { fprintf(stderr, "toRep\n"); exit(1); }
    for (int c = 0; c < (int)classSize.size(); c++) if ((long long)cent[c].size() * classSize[c] != (long long)GL4.size()) { fprintf(stderr, "centraliser size\n"); exit(1); }
}

// ---- clique enumeration (each clique once: increasing vertex index)
struct CliqueStats { long long nodes[16] = { 0 }; };
template <class Leaf> static void cliques(vector<int>& cur, const BS& cand, int need, CliqueStats& st, Leaf& leaf) {
    st.nodes[cur.size()]++;
    if (need == 0) { leaf(cur); return; }
    if (bcount(cand) < need) return;
    for (int wi = 0; wi < NW; wi++) {
        u64 w = cand.w[wi];
        while (w) {
            int b = (int)_tzcnt_u64(w); w &= w - 1; int v = wi * 64 + b;
            BS nc; for (int i = 0; i < NW; i++) nc.w[i] = 0;
            nc.w[wi] = cand.w[wi] & adj[v].w[wi] & (b == 63 ? 0 : (~0ULL << (b + 1)));
            for (int i = wi + 1; i < NW; i++) nc.w[i] = cand.w[i] & adj[v].w[i];
            cur.push_back(v);
            cliques(cur, nc, need - 1, st, leaf);
            cur.pop_back();
        }
    }
}
static inline void above(BS& o, const BS& a, const BS& b, int v) { // o = a & b & {> v}
    int wi = v >> 6, bb = v & 63;
    for (int i = 0; i < wi; i++) o.w[i] = 0;
    o.w[wi] = a.w[wi] & b.w[wi] & (bb == 63 ? 0 : (~0ULL << (bb + 1)));
    for (int i = wi + 1; i < NW; i++) o.w[i] = a.w[i] & b.w[i];
}
static long long count4(const BS& Q) {
    long long n = 0;
    for (int wi = 0; wi < NW; wi++) { u64 w = Q.w[wi]; while (w) { int u = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
        BS a; above(a, Q, adj[u], u);
        for (int wj = 0; wj < NW; wj++) { u64 x = a.w[wj]; while (x) { int t = wj * 64 + (int)_tzcnt_u64(x); x &= x - 1;
            BS b; above(b, a, adj[t], t);
            for (int wk = 0; wk < NW; wk++) { u64 y = b.w[wk]; while (y) { int s = wk * 64 + (int)_tzcnt_u64(y); y &= y - 1;
                BS c; above(c, b, adj[s], s); n += bcount(c); } }
        } }
    } }
    return n;
}
static Mask graph_mask(u16 N) { Mask m; mclear(m); for (u8 x = 1; x < 16; x++) mset(m, x | (mv4(N, x) << 4)); return m; }
static Mask XMASK() { Mask m; mclear(m); for (int v = 1; v < 16; v++) mset(m, v); return m; }
static Mask YMASK() { Mask m; mclear(m); for (int v = 1; v < 16; v++) mset(m, v << 4); return m; }
static vector<Mask> members_of(const vector<u16>& Ns) { vector<Mask> M = { XMASK(), YMASK(), graph_mask(ID4) }; for (u16 N : Ns) M.push_back(graph_mask(N)); return M; }

// ---- hole-set invariants
static void pair_degrees(const Mask& H, int deg[256]) { // deg[v] = #{p in H : p ^ v in H}, for every v
    vector<u8> p = mpoints(H); for (int v = 0; v < 256; v++) deg[v] = 0;
    for (size_t i = 0; i < p.size(); i++) for (size_t j = i + 1; j < p.size(); j++) deg[p[i] ^ p[j]] += 2;
}
// pc[v] = number of planes inside H through v (v in H), or number of planes inside H + v through v (v not in H)
static int plane_counts(const Mask& H, int pc[256]) {
    vector<u8> p = mpoints(H); set<Mask> planes; vector<array<u8, 2>> lines;
    for (int v = 0; v < 256; v++) pc[v] = 0;
    for (size_t i = 0; i < p.size(); i++) for (size_t j = i + 1; j < p.size(); j++) {
        u8 a = p[i], b = p[j]; if (!mget(H, a ^ b)) continue;
        if (a < (a ^ b) && b < (a ^ b)) lines.push_back({ a, b });
        for (u8 c : p) { if (c == a || c == b || c == (a ^ b)) continue;
            if (mget(H, c ^ a) && mget(H, c ^ b) && mget(H, c ^ a ^ b)) {
                Mask m; mclear(m); u8 q[7] = { a, b, (u8)(a ^ b), c, (u8)(c ^ a), (u8)(c ^ b), (u8)(c ^ a ^ b) }; for (u8 x : q) mset(m, x); planes.insert(m); } }
    }
    for (auto& m : planes) for (int v = 1; v < 256; v++) if (mget(m, v)) pc[v]++;
    // v not in H: lines L of H with the other 4 points of <L, v> in H (planes meeting H in all but v)
    for (int v = 1; v < 256; v++) if (!mget(H, v)) {
        int c = 0;
        for (auto& l : lines) { u8 a = l[0], b = l[1]; if (mget(H, a ^ v) && mget(H, b ^ v) && mget(H, a ^ b ^ v)) c++; }
        pc[v] = c;
    }
    return (int)planes.size();
}
static string hole_key(const Mask& H) {
    vector<u8> p = mpoints(H); int deg[256]; pair_degrees(H, deg);
    int pc[256]; int np = plane_counts(H, pc);
    int lines = 0; map<int, int> dh, ho, hy, ph, po;
    for (u8 v : p) { lines += deg[v]; dh[deg[v] / 2]++; ph[pc[v]]++; }
    for (int v = 1; v < 256; v++) if (!mget(H, v)) po[pc[v]]++;
    for (int v = 1; v < 256; v++) if (!mget(H, v)) ho[deg[v] / 2]++;
    for (int u = 1; u < 256; u++) { int c = 0; for (u8 v : p) c += !par(u & v); hy[c]++; }
    char buf[64]; string s; sprintf(buf, "n%zu L%d deg", p.size(), lines / 6); s += buf;
    for (auto& kv : dh) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    s += " out"; for (auto& kv : ho) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    s += " hyp"; for (auto& kv : hy) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    sprintf(buf, " planes %d pdeg", np); s += buf; for (auto& kv : ph) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    s += " pout"; for (auto& kv : po) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    return s;
}

// ---- linear maps of V by column images
typedef array<u8, 8> Lin;
static inline u8 lin_apply(const Lin& A, u8 v) { u8 r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= A[i]; return r; }
static Lin lin_comp(const Lin& A, const Lin& B) { Lin C; for (int i = 0; i < 8; i++) C[i] = lin_apply(A, B[i]); return C; } // A after B
static bool lin_inv_ok(const Lin& A, Lin& out) {
    u8 vec[8], comb[8]; for (int i = 0; i < 8; i++) { vec[i] = A[i]; comb[i] = (u8)(1 << i); }
    for (int bit = 0; bit < 8; bit++) {
        int p = -1; for (int i = bit; i < 8; i++) if (vec[i] >> bit & 1) { p = i; break; }
        if (p < 0) return false;
        swap(vec[p], vec[bit]); swap(comb[p], comb[bit]);
        for (int i = 0; i < 8; i++) if (i != bit && (vec[i] >> bit & 1)) { vec[i] ^= vec[bit]; comb[i] ^= comb[bit]; }
    }
    for (int k = 0; k < 8; k++) out[k] = comb[k];
    return true;
}
static Lin lin_inv(const Lin& A) { Lin C; if (!lin_inv_ok(A, C)) { fprintf(stderr, "singular\n"); exit(1); } return C; }
static Lin lin_id() { Lin I; for (int i = 0; i < 8; i++) I[i] = (u8)(1 << i); return I; }
static Mask mask_image(const Lin& A, const Mask& m) { Mask r; mclear(r); for (int v = 1; v < 256; v++) if (mget(m, v)) mset(r, lin_apply(A, (u8)v)); return r; }
static Lin diag_hh(u16 h) { Lin D; for (int i = 0; i < 4; i++) { u8 he = mv4(h, (u8)(1 << i)); D[i] = he; D[4 + i] = (u8)(he << 4); } return D; }
static Lin diag_ab(u16 a, u16 b) { Lin D; for (int i = 0; i < 4; i++) { D[i] = mv4(a, (u8)(1 << i)); D[4 + i] = (u8)(mv4(b, (u8)(1 << i)) << 4); } return D; }
static Lin random_gl8(mt19937& g) { while (true) { Lin A; for (int i = 0; i < 8; i++) A[i] = (u8)(g() & 255); Lin B; if (lin_inv_ok(A, B)) return A; } }

static void pick_independent(const vector<u8>& cand, int need, u8* out) {
    bool inspan[256] = { 0 }; inspan[0] = true; int k = 0;
    for (u8 u : cand) {
        if (k == need) break;
        if (inspan[u]) continue;
        out[k++] = u;
        bool nw[256]; memcpy(nw, inspan, 256);
        for (int v = 0; v < 256; v++) if (inspan[v]) nw[v ^ u] = true;
        memcpy(inspan, nw, 256);
    }
    if (k != need) { fprintf(stderr, "pick_independent failed\n"); exit(1); }
}
static void basis_of_mask(const Mask& m, int dim, u8* out) { pick_independent(mpoints(m), dim, out); }
static void annihilator(const Mask& m, int need, u8* out) {
    vector<u8> ann; for (int u = 1; u < 256; u++) { bool ok = true; for (int v = 1; v < 256 && ok; v++) if (mget(m, v) && par(u & v)) ok = false; if (ok) ann.push_back((u8)u); }
    pick_independent(ann, need, out);
}
// N with graph(N) = m (m complementary to Y), or -1
static int matrix_of_graph(const Mask& m) {
    u8 col[4]; for (int i = 0; i < 4; i++) { col[i] = 0xff; for (int v = 1; v < 256; v++) if (mget(m, v) && (v & 15) == (1 << i)) col[i] = (u8)(v >> 4); if (col[i] == 0xff) return -1; }
    u16 N = fromcols4(col); return meq(graph_mask(N), m) ? N : -1;
}

// ================================================================================ isomorphism of point sets
// Linear maps g with g(H1) = H2, by backtracking over images of a basis inside H1. Pruning: on every element of
// the span built so far, v in H1 <=> g(v) in H2 and deg1(v) = deg2(g v) (deg = number of H-pairs summing to v).
struct IsoCtx { Mask H1, H2; int d1[256], d2[256]; u8 b[8]; u8 vs[256], im[256]; long long count; bool first; Lin found; vector<u8> cand2; };
static void iso_rec(IsoCtx& X, int k) {
    if (k == 8) { X.count++; if (X.count == 1) for (int i = 0; i < 8; i++) X.found[i] = X.im[1 << i]; return; }
    int n = 1 << k; u8 bk = X.b[k];
    for (u8 c : X.cand2) {
        if (X.d2[c] != X.d1[bk]) continue;
        bool ok = true;
        for (int s = 0; s < n && ok; s++) {
            u8 v = X.vs[s] ^ bk, w = X.im[s] ^ c;
            if (!w || mget(X.H1, v) != mget(X.H2, w) || X.d1[v] != X.d2[w]) ok = false;
        }
        if (!ok) continue;
        for (int s = 0; s < n; s++) { X.vs[n + s] = X.vs[s] ^ bk; X.im[n + s] = X.im[s] ^ c; }
        iso_rec(X, k + 1);
        if (X.first && X.count) return;
    }
}
static long long iso_count(const Mask& H1, const Mask& H2, bool first, Lin* out) {
    IsoCtx X; X.H1 = H1; X.H2 = H2; pair_degrees(H1, X.d1); pair_degrees(H2, X.d2); X.count = 0; X.first = first;
    if (mcount(H1) != mcount(H2)) return 0;
    { int p1[256], p2[256]; plane_counts(H1, p1); plane_counts(H2, p2); for (int v = 0; v < 256; v++) { X.d1[v] = X.d1[v] * 65536 + p1[v]; X.d2[v] = X.d2[v] * 65536 + p2[v]; } }
    // basis inside H1, rarest degree first
    vector<u8> p1 = mpoints(H1); map<int, int> freq; for (u8 v : p1) freq[X.d1[v]]++;
    stable_sort(p1.begin(), p1.end(), [&](u8 a, u8 b) { return freq[X.d1[a]] < freq[X.d1[b]]; });
    pick_independent(p1, 8, X.b);
    X.cand2 = mpoints(H2);
    X.vs[0] = 0; X.im[0] = 0;
    iso_rec(X, 0);
    if (out && X.count) { // found maps b_i -> im[1<<i]; as column images of e_i: g = M2 * M1^-1
        Lin B1; for (int i = 0; i < 8; i++) B1[i] = X.b[i];
        *out = lin_comp(X.found, lin_inv(B1));
        if (!meq(mask_image(*out, H1), H2)) { fprintf(stderr, "iso map check failed\n"); exit(1); }
    }
    return X.count;
}
static vector<Mask> ALLSOLIDS;
static void init_solids() {
    for (int p0 = 0; p0 < 8; p0++) for (int p1 = p0 + 1; p1 < 8; p1++) for (int p2 = p1 + 1; p2 < 8; p2++) for (int p3 = p2 + 1; p3 < 8; p3++) {
        int piv[4] = { p0, p1, p2, p3 }; u8 pm = (u8)((1 << p0) | (1 << p1) | (1 << p2) | (1 << p3));
        vector<pair<int, int>> fr; for (int k = 0; k < 4; k++) for (int q = piv[k] + 1; q < 8; q++) if (!(pm >> q & 1)) fr.push_back({ k, q });
        for (u32 t = 0; t < (1u << fr.size()); t++) {
            u8 row[4]; for (int k = 0; k < 4; k++) row[k] = (u8)(1 << piv[k]);
            for (size_t f = 0; f < fr.size(); f++) if (t >> f & 1) row[fr[f].first] |= (u8)(1 << fr[f].second);
            Mask m; mclear(m); for (int s = 1; s < 16; s++) { u8 v = 0; for (int k = 0; k < 4; k++) if (s >> k & 1) v ^= row[k]; mset(m, v); }
            ALLSOLIDS.push_back(m);
        }
    }
}
static bool contains_solid(const Mask& H) { for (auto& s : ALLSOLIDS) if (msub(s, H)) return true; return false; }
static void covers_rec(const vector<Mask>& S, const vector<vector<int>>& at, const Mask& C, Mask& cov, vector<int>& ch, vector<vector<int>>& out) {
    int best = -1, bc = INT_MAX;
    for (int v = 1; v < 256; v++) if (mget(C, v) && !mget(cov, v)) {
        int c = 0; for (int s : at[v]) if (mdisj(S[s], cov)) c++;
        if (c < bc) { bc = c; best = v; if (!c) break; }
    }
    if (best < 0) { out.push_back(ch); return; }
    if (!bc) return;
    for (int s : at[best]) if (mdisj(S[s], cov)) { Mask save = cov; mor(cov, S[s]); ch.push_back(s); covers_rec(S, at, C, cov, ch, out); ch.pop_back(); cov = save; }
}
static vector<vector<Mask>> covers_of(const Mask& H, int& nsolids) {
    Mask C = mcompl(H); vector<Mask> S; for (auto& s : ALLSOLIDS) if (msub(s, C)) S.push_back(s);
    nsolids = (int)S.size();
    vector<vector<int>> at(256); for (int i = 0; i < (int)S.size(); i++) for (int v = 1; v < 256; v++) if (mget(S[i], v)) at[v].push_back(i);
    Mask cov; mclear(cov); vector<int> ch; vector<vector<int>> raw; covers_rec(S, at, C, cov, ch, raw);
    vector<vector<Mask>> out; for (auto& r : raw) { vector<Mask> M; for (int s : r) M.push_back(S[s]); out.push_back(M); }
    return out;
}

// ================================================================================ configs
static vector<vector<int>> spread_cliques(double& weighted) {
    vector<vector<int>> red; weighted = 0;
    for (int c = 0; c < (int)classSize.size(); c++) {
        int r = classRep[c];
        BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
        for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
        vector<int> cur = { r }; CliqueStats st; double wsum = 0;
        auto leaf = [&](const vector<int>& K) { int m = 0; for (int u : K) m += cls[u] == c; wsum += 1.0 / m; red.push_back(K); };
        cliques(cur, P, 13, st, leaf);
        weighted += classSize[c] * wsum;
    }
    return red;
}

// ================================================================================ canonical form of a partial spread
// As thirteen_check's canon13, for n >= 4 members (k = n - 3 graph vertices per frame). Frames: an ordered triple
// (A, B, C) of members sent to (X, Y, Z) by Cmap (basis a_i of A -> e_i, the B-part b_i of a_i in B + C -> e_{4+i});
// every g in GL(8,2) with g(A, B, C) = (X, Y, Z) is diag(h, h) o Cmap. Key = the least sorted vertex list over the
// frames whose ordering has the least sorted conjugacy-class vector (intrinsic, so the key is an orbit invariant; an
// image of P, so equal keys mean equivalent P). Only h moving a least-class member onto its class rep can give the
// least list (the rep is the smallest index of its class). Frames giving the key <-> Aut(P), one to one.
struct Canon { vector<u16> key; long long aut = 0; Lin frame; };
static bool frame_of(const vector<array<u8, 4>>& bas, const vector<vector<u8>>& pts, const vector<Mask>& M, int a, int b, int c, Lin& Cm) {
    Lin Ci;
    for (int i = 0; i < 4; i++) {
        u8 ai = bas[a][i]; int bi = -1;
        for (u8 p : pts[b]) if (mget(M[c], ai ^ p)) { bi = p; break; }
        if (bi < 0) return false;
        Ci[i] = ai; Ci[4 + i] = (u8)bi;
    }
    return lin_inv_ok(Ci, Cm);
}
static int vertex_in_frame(const Lin& Cm, const array<u8, 4>& bd) {
    u8 x[4], y[4];
    for (int j = 0; j < 4; j++) { u8 w = lin_apply(Cm, bd[j]); x[j] = w & 15; y[j] = w >> 4; }
    u16 Xm = fromcols4(x), Ym = fromcols4(y);
    if (gl4idx[Xm] < 0) return -1;
    u16 N = mm4(Ym, inverse4(Xm));
    return vid[N];
}
typedef array<u16, 14> VArr;
static Canon canonP(const vector<Mask>& M) {
    int n = (int)M.size(), k = n - 3; if (n < 4 || n > 17) { fprintf(stderr, "canonP: %d members\n", n); exit(1); }
    vector<array<u8, 4>> bas(n); vector<vector<u8>> pts(n);
    for (int i = 0; i < n; i++) { basis_of_mask(M[i], 4, bas[i].data()); pts[i] = mpoints(M[i]); }
    struct Ord { Lin Cm; VArr v; };
    vector<Ord> keep; array<int, 14> bestCls; bestCls.fill(INT_MAX);
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) for (int c = 0; c < n; c++) {
        if (a == b || a == c || b == c) continue;
        Ord o; if (!frame_of(bas, pts, M, a, b, c, o.Cm)) { fprintf(stderr, "canonP: members not complementary\n"); exit(1); }
        int t = 0; array<int, 14> cv; cv.fill(INT_MAX);
        for (int d = 0; d < n; d++) { if (d == a || d == b || d == c) continue; int v = vertex_in_frame(o.Cm, bas[d]); if (v < 0) { fprintf(stderr, "canonP: not in V5824\n"); exit(1); } o.v[t] = (u16)v; cv[t] = cls[v]; t++; }
        sort(cv.begin(), cv.begin() + k);
        if (cv < bestCls) { bestCls = cv; keep.clear(); }
        if (cv == bestCls) keep.push_back(o);
    }
    Canon R; R.key.assign(k, 0xffff); VArr img;
    int c0 = bestCls[0], r0 = classRep[c0];
    for (auto& o : keep) for (int e = 0; e < k; e++) {
        if (cls[o.v[e]] != c0) continue;
        u16 h0 = GL4[toRep[o.v[e]]];
        for (int z : cent[c0]) {
            int hi = gl4idx[mm4(GL4[z], h0)];
            for (int j = 0; j < k; j++) img[j] = conjTab[hi][o.v[j]];
            sort(img.begin(), img.begin() + k);
            if (img[0] != r0) { fprintf(stderr, "canonP: rep\n"); exit(1); }
            int cmp = 0; for (int j = 0; j < k && !cmp; j++) cmp = img[j] < R.key[j] ? -1 : img[j] > R.key[j] ? 1 : 0;
            if (cmp < 0) { R.key.assign(img.begin(), img.begin() + k); R.aut = 1; R.frame = lin_comp(diag_hh(GL4[hi]), o.Cm); }
            else if (cmp == 0) R.aut++;
        }
    }
    return R;
}
static vector<u16> key_matrices(const vector<u16>& key) { vector<u16> N; for (u16 v : key) N.push_back(VN[v]); return N; }

// k-cliques among a small vertex list (all pairwise adjacent), or whether one exists
static long long kcl(const vector<int>& cand, int need, bool one) {
    if (need == 0) return 1;
    if ((int)cand.size() < need) return 0;
    long long n = 0;
    for (size_t i = 0; i + need <= cand.size(); i++) {
        vector<int> nc; for (size_t j = i + 1; j < cand.size(); j++) if (bget(adj[cand[i]], cand[j])) nc.push_back(cand[j]);
        n += kcl(nc, need - 1, one); if (one && n) return n;
    }
    return n;
}
static vector<int> blist(const BS& b) { vector<int> v; for (int wi = 0; wi < NW; wi++) { u64 w = b.w[wi]; while (w) { v.push_back(wi * 64 + (int)_tzcnt_u64(w)); w &= w - 1; } } return v; }

// ================================================================================ Part A: 3-cliques, all canonicalised
// For class c: r = classRep[c], P = adj(r) & {class >= c}; 2-cliques K in P, so K + r is a 3-clique (= a 6-partial
// spread with Y, X, Z among its members) with minimal class c. Weighted count #3-cliques = sum_c |class c| *
// sum_K 1/m(K + r) (m = number of class-c vertices, m <= 3), kept integral by the factor 2520. Per clique: |Q| (Q = the
// common neighbourhood = every solid disjoint from all 6 members) and the canonical form. Per orbit the weighted sum
// must be 2520 * 120 * 20160 / |Aut P|. A class can be split into slices (first vertices fi with fi % nslice == slice).
static const int TYPES = 9;
static const char* TNAME[TYPES] = { "spr", "ext13", "ext12", "ext11", "ext10", "ext9", "ext8", "ext7", "max" }; // how far the 6 extends: 17, 13, ..., 7, nowhere
static const int TYPES7 = 8;
static const char* TNAME7[TYPES7] = { "spr", "ext13", "ext12", "ext11", "ext10", "ext9", "ext8", "max" }; // seven-check's types
static const int TYPES8 = 7;
static const char* TNAME8[TYPES8] = { "spr", "ext13", "ext12", "ext11", "ext10", "ext9", "max" }; // eight-check's types
static const long long FR6 = 120LL * 20160, FR7 = 210LL * 20160;
static const int KN = 3; // graph members per configuration (the 6 members less X, Y, Z)
typedef array<u16, KN> KeyK;
struct OrbAcc { long long w = 0, n = 0, aut = 0; int q = -1; int bad = 0; };
static int g_testFirsts = 0;
static int mode_parta(int c, const char* outfile, const char* maxfile, int slice, int nslice) {
    init_conj();
    if (c < 0 || c >= (int)classSize.size() || nslice < 1 || slice < 0 || slice >= nslice) { printf("bad class / slice\n"); return 1; }
    FILE* of = fopen(outfile, "w"); FILE* mf = fopen(maxfile, "w");
    if (!of || !mf) { printf("cannot open output\n"); return 1; }
    printf("Part A, class %d slice %d of %d: 3-cliques of V5824 (6-partial spreads with Y, X, Z among them) with minimal class %d\n", c, slice, nslice, c); fflush(stdout);
    int r = classRep[c];
    BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
    for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
    vector<int> all = blist(P), firsts;
    for (int fi = 0; fi < (int)all.size(); fi++) if (fi % nslice == slice) firsts.push_back(all[fi]);
    if (g_testFirsts > 0 && (int)firsts.size() > g_testFirsts) firsts.resize(g_testFirsts); // quick tests only
    map<KeyK, OrbAcc> orb; map<int, long long> qW; long long cw = 0, cl = 0, csq = 0, cmax = 0; CliqueStats tot; mutex mu;
    auto t0 = chrono::steady_clock::now(); int doneFirsts = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int fi = 0; fi < (int)firsts.size(); fi++) {
        int v = firsts[fi]; BS nc; above(nc, P, adj[v], v);
        vector<int> cur = { v }; CliqueStats st;
        map<KeyK, OrbAcc> lo; long long lw = 0, ll = 0, lsq = 0, lmax = 0; map<int, long long> lq; vector<string> maxLines;
        auto leaf = [&](const vector<int>& K) {
            ll++;
            int m = 1; for (int u : K) m += cls[u] == c;
            long long w = 2520 / m;
            BS Q = adj[r]; for (int u : K) { BS t; band(t, Q, adj[u]); Q = t; }
            int q = bcount(Q);
            lw += w; lsq += w * q; lq[q] += w;
            vector<u16> Ns = { VN[r] }; for (int u : K) Ns.push_back(VN[u]);
            Canon C = canonP(members_of(Ns));
            KeyK key; for (int j = 0; j < KN; j++) key[j] = C.key[j];
            OrbAcc& A = lo[key];
            if (A.n == 0) { A.aut = C.aut; A.q = q; }
            else if (A.aut != C.aut || A.q != q) A.bad++;
            A.n++; A.w += w;
            if (!q) { lmax++; char buf[128]; int o = sprintf(buf, "MAXC c%d :", c); for (u16 N : Ns) o += sprintf(buf + o, " %04x", N); maxLines.push_back(buf); }
        };
        cliques(cur, nc, KN - 2, st, leaf);
        lock_guard<mutex> g(mu);
        for (int i = 0; i < 16; i++) tot.nodes[i] += st.nodes[i];
        for (auto& kv : lo) {
            OrbAcc& A = orb[kv.first]; const OrbAcc& B = kv.second;
            if (A.n == 0) { A.aut = B.aut; A.q = B.q; }
            else if (A.aut != B.aut || A.q != B.q) A.bad++;
            A.bad += B.bad; A.n += B.n; A.w += B.w;
        }
        cw += lw; cl += ll; csq += lsq; cmax += lmax; for (auto& kv : lq) qW[kv.first] += kv.second;
        for (auto& s : maxLines) fprintf(mf, "%s\n", s.c_str());
        if (++doneFirsts % 200 == 0) { printf("  %d of %zu first vertices, %lld cliques, %zu orbits, %.0f s\n", doneFirsts, firsts.size(), cl, orb.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout); }
    }
    fclose(mf);
    double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    long long cs = classSize[c];
    printf("class %d slice %d/%d (size %d, rep %04x, |P| = %d, %zu first vertices): %lld reduced 3-cliques, %lld maximal; %zu orbits; %.0f s\n  nodes by clique size (excl. rep):", c, slice, nslice, classSize[c], VN[r], bcount(P), firsts.size(), cl, cmax, orb.size(), el);
    for (int i = 0; i < KN; i++) printf(" %lld", tot.nodes[i]); printf("\n");
    fprintf(of, "CLASS %d slice %d %d leaves %lld w %lld sumq %lld max %lld secs %.0f\n", c, slice, nslice, cl, cw * cs, csq * cs, cmax, el);
    fprintf(of, "QH %d", c); for (auto& kv : qW) fprintf(of, " %d:%lld", kv.first, kv.second * cs); fprintf(of, "\n");
    for (auto& kv : orb) {
        const OrbAcc& A = kv.second;
        fprintf(of, "ORBP aut %lld w %lld n %lld q %d bad %d :", A.aut, A.w * cs, A.n, A.q, A.bad);
        for (u16 v : kv.first) fprintf(of, " %04x", VN[v]); fprintf(of, "\n");
    }
    fprintf(of, "END\n"); fclose(of);
    return 0;
}

// ---- small clique helpers for the per-orbit invariants
static void mcrec(const vector<int>& cand, int sz, int& best, int cap) {
    if (sz > best) best = sz;
    if (best >= cap) return;
    for (size_t i = 0; i < cand.size(); i++) {
        if (sz + (int)(cand.size() - i) <= best) return;
        vector<int> nc; for (size_t j = i + 1; j < cand.size(); j++) if (bget(adj[cand[i]], cand[j])) nc.push_back(cand[j]);
        mcrec(nc, sz + 1, best, cap);
        if (best >= cap) return;
    }
}
static int max_clique(const vector<int>& cand, int cap) { int best = 0; mcrec(cand, 0, best, cap); return best; }

// earlier checks' configurations (CONF lines: n matrices, a type name)
struct CPrev { int idx, type; long long aut; vector<u16> Ns; };
static vector<CPrev> read_prev(const char* file, int nmat, const char* const* tn, int ntypes) {
    vector<CPrev> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "CONF", 4)) continue;
        CPrev c; char ty[16]; int h;
        if (sscanf(line, "CONF %d %15s hole %d aut %lld", &c.idx, ty, &h, &c.aut) != 4) { printf("bad CONF line in %s\n", file); exit(1); }
        c.type = -1; for (int t = 0; t < ntypes; t++) if (!strcmp(ty, tn[t])) c.type = t;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; c.Ns.push_back((u16)v); p = e; }
        if ((int)c.Ns.size() != nmat || c.type < 0 || c.idx != (int)out.size()) { printf("bad CONF line in %s\n", file); exit(1); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
static vector<CPrev> read_conf8(const char* file) { return read_prev(file, 5, TNAME8, TYPES8); }
static vector<CPrev> read_conf7(const char* file) { return read_prev(file, 4, TNAME7, TYPES7); }
static vector<u16> key_of_Ns(const vector<u16>& Ns) { vector<u16> k; for (u16 N : Ns) k.push_back((u16)vid[N]); sort(k.begin(), k.end()); return k; }

// ================================================================================ merge: orbits, per-orbit invariants, identities
// Identities: (1) sum over 3-cliques of |Q| = 4 * (weighted 4-cliques) = 4 * sum over seven-check's 3731 of
// 210*20160/|Aut|; (2) per orbit, weighted cliques = 120*20160/|Aut P|; (3) sum_P 120*20160/|Aut P| * (spreads through
// P) = C(14,3) * 1128864 (1128864 spreads through X, Y, Z = 14-cliques; any 3 of their 14 graph members);
// (4) per 7-orbit O: sum_P 120*20160/|Aut P| * #{q in Q(P) : P + q ~ O} = 4 * 210*20160/|Aut O| (a 4-clique with one
// vertex marked, 4 ways).
struct Orb6 { KeyK key; long long aut = 0, w = 0, n = 0; int q = -1, bad = 0, mc = -1, type = -1, unk = 0, tbad = 0; long long n11 = 0; map<int, int> link; };
static int mode_merge(const char* outfile, const char* conf7, const vector<const char*>& parts) {
    init_conj();
    auto c7 = read_conf7(conf7); map<vector<u16>, int> key7; long long W4 = 0;
    for (auto& c : c7) { key7[key_of_Ns(c.Ns)] = c.idx; if (FR7 % c.aut) { printf("7: aut does not divide\n"); return 1; } W4 += FR7 / c.aut; }
    printf("seven-check's configurations: %zu; sum 210*20160/|Aut| = %lld (its count4: 14122766336)\n", c7.size(), W4);
    map<KeyK, Orb6> orb; long long Wall = 0, SQ = 0, L = 0, M = 0; map<int, long long> qW; long long inconsistent = 0;
    map<int, set<int>> slicesSeen; map<int, int> nslOf; int dupSlice = 0;
    for (const char* pf : parts) {
        FILE* fp = fopen(pf, "r"); if (!fp) { printf("cannot open %s\n", pf); return 1; }
        static char line[4096]; bool end = false;
        while (fgets(line, sizeof line, fp)) {
            if (!strncmp(line, "CLASS", 5)) { int c, sl, nsl; long long l, w, sq, mx;
                if (sscanf(line, "CLASS %d slice %d %d leaves %lld w %lld sumq %lld max %lld", &c, &sl, &nsl, &l, &w, &sq, &mx) != 7) { printf("bad CLASS\n"); return 1; }
                if (nslOf.count(c) && nslOf[c] != nsl) { printf("class %d: slices of different splits\n", c); return 1; }
                nslOf[c] = nsl; if (!slicesSeen[c].insert(sl).second) dupSlice++;
                L += l; Wall += w; SQ += sq; M += mx; }
            else if (!strncmp(line, "QH", 2)) { const char* p = line + 2; char* e; strtol(p, &e, 10); p = e;
                while (true) { long q = strtol(p, &e, 10); if (e == p) break; p = e + 1; long long w = strtoll(p, &e, 10); p = e; qW[(int)q] += w; } }
            else if (!strncmp(line, "ORBP", 4)) {
                long long aut, w, n; int q, bad;
                if (sscanf(line, "ORBP aut %lld w %lld n %lld q %d bad %d", &aut, &w, &n, &q, &bad) != 5) { printf("bad ORBP\n"); return 1; }
                const char* p = strchr(line, ':') + 1; char* e; vector<u16> Ns; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
                if (Ns.size() != KN) { printf("bad ORBP key\n"); return 1; }
                KeyK k; for (int j = 0; j < KN; j++) { if (vid[Ns[j]] < 0) { printf("key not in V5824\n"); return 1; } k[j] = (u16)vid[Ns[j]]; }
                Orb6& A = orb[k];
                if (A.n == 0) { A.key = k; A.aut = aut; A.q = q; } else if (A.aut != aut || A.q != q) inconsistent++;
                A.bad += bad; A.n += n; A.w += w;
            }
            else if (!strncmp(line, "END", 3)) end = true;
        }
        fclose(fp);
        if (!end) { printf("%s: incomplete (no END)\n", pf); return 1; }
    }
    for (int c = 0; c < (int)classSize.size(); c++) if (!nslOf.count(c) || (int)slicesSeen[c].size() != nslOf[c]) { printf("class %d: slices missing\n", c); return 1; }
    if (dupSlice) { printf("%d slices given twice\n", dupSlice); return 1; }
    auto chk = [](long long x) { if (x % 2520) printf("  (not divisible by 2520!)\n"); return x / 2520; };
    long long W3 = chk(Wall);
    printf("reduced 3-cliques %lld, of them maximal (Q empty) %lld; weighted 3-cliques %lld\n", L, M, W3);
    printf("sum over 3-cliques of |Q| = %lld; 4 * (weighted 4-cliques, seven-check) = %lld  %s\n", chk(SQ), 4 * W4, chk(SQ) == 4 * W4 ? "OK" : "MISMATCH");
    printf("by |Q|:"); for (auto& kv : qW) printf(" %d:%lld", kv.first, chk(kv.second)); printf("\n");
    vector<Orb6> O; for (auto& kv : orb) O.push_back(kv.second);
    long long idbad = 0, sumAll = 0, badSum = 0;
    for (auto& A : O) { badSum += A.bad; if (FR6 % A.aut || A.w != 2520 * (FR6 / A.aut)) idbad++; sumAll += FR6 / A.aut; }
    printf("orbits %zu; clique-level inconsistencies (|Aut| or |Q|) %lld + %lld; per-orbit identity (weighted cliques = 120*20160/|Aut|) fails %lld; sum 120*20160/|Aut| = %lld (weighted 3-cliques %lld)  %s\n",
        O.size(), badSum, inconsistent, idbad, sumAll, W3, sumAll == W3 ? "OK" : "MISMATCH");
    fflush(stdout);
    auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)O.size(); i++) {
        Orb6& A = O[i];
        vector<u16> Ns; for (u16 v : A.key) Ns.push_back(VN[v]);
        BS Q; for (int j = 0; j < NW; j++) Q.w[j] = ~0ULL;
        for (u16 v : A.key) { BS t; band(t, Q, adj[v]); Q = t; }
        vector<int> qv = blist(Q);
        if ((int)qv.size() != A.q) A.tbad++;
        A.mc = max_clique(qv, 11);
        A.n11 = kcl(qv, 11, false);
        int t2 = TYPES - 1; // from the 7-orbits P + q: P extends as far as its best one-solid extension
        for (int q : qv) {
            vector<u16> N2 = Ns; N2.push_back(VN[q]);
            Canon K = canonP(members_of(N2));
            auto it = key7.find(K.key);
            if (it == key7.end()) { A.unk++; continue; }
            A.link[it->second]++; t2 = min(t2, c7[it->second].type); // 7-types spr..ext8, max <-> 6-types spr..ext8, ext7 (same index)
        }
        // 6 -> spread (11 more), maximal 13 (7), 12 (6), 11 (5), 10 (4), 9 (3), 8 (2), 7 (1), nothing (0); no maximal 14, 15, 16
        A.type = A.mc >= 11 ? 0 : A.mc == 7 ? 1 : A.mc == 6 ? 2 : A.mc == 5 ? 3 : A.mc == 4 ? 4 : A.mc == 3 ? 5 : A.mc == 2 ? 6 : A.mc == 1 ? 7 : A.mc == 0 ? 8 : -1;
        if (A.type < 0 || A.mc > 11 || t2 != A.type) A.tbad++;
        if (A.type == 0 && !A.n11) A.tbad++;
        if (A.type != 0 && A.n11) A.tbad++;
    }
    printf("per-orbit invariants: %.0f s\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    long long tbad = 0, unk = 0, n11sum = 0, byT[TYPES] = { 0 }, nT[TYPES] = { 0 }; map<int, int> mcH;
    vector<long long> linkSum(c7.size(), 0);
    for (auto& A : O) {
        tbad += A.tbad; unk += A.unk; n11sum += (FR6 / A.aut) * A.n11; mcH[A.mc]++;
        if (A.type >= 0) { byT[A.type] += FR6 / A.aut; nT[A.type]++; }
        for (auto& kv : A.link) linkSum[kv.first] += (FR6 / A.aut) * kv.second;
    }
    long long linkBad = 0; for (auto& c : c7) if (linkSum[c.idx] != 4 * (FR7 / c.aut)) linkBad++;
    const long long SPR = 364LL * 1128864; // C(14,3) * spreads through X, Y, Z
    printf("largest clique in Q (how many more solids fit), orbits:"); for (auto& kv : mcH) printf(" %d:%d", kv.first, kv.second); printf("\n");
    printf("orbits by type:"); for (int t = 0; t < TYPES; t++) printf(" %s %lld", TNAME[t], nT[t]); printf("; weighted"); for (int t = 0; t < TYPES; t++) printf(" %lld", byT[t]); printf("\n");
    printf("type / Q / n11 inconsistencies %lld; P + q on no seven-check orbit %lld\n", tbad, unk);
    printf("sum 120*20160/|Aut| * (spreads through P) = %lld; C(14,3) * 1128864 = %lld  %s\n", n11sum, SPR, n11sum == SPR ? "OK" : "MISMATCH");
    printf("per 7-orbit O: sum_P 120*20160/|Aut P| * #{q : P + q in O} = 4 * 210*20160/|Aut O| fails for %lld of %zu\n", linkBad, c7.size());
    sort(O.begin(), O.end(), [](const Orb6& a, const Orb6& b) { return a.type != b.type ? a.type < b.type : a.key < b.key; });
    FILE* of = fopen(outfile, "w");
    for (size_t i = 0; i < O.size(); i++) {
        auto& A = O[i];
        fprintf(of, "ORB %zu %s aut %lld w %lld n %lld q %d mc %d n11 %lld bad %d :", i, A.type >= 0 ? TNAME[A.type] : "??", A.aut, A.w / 2520, A.n, A.q, A.mc, A.n11, A.bad + A.tbad);
        for (u16 v : A.key) fprintf(of, " %04x", VN[v]); fprintf(of, "\n");
        fprintf(of, "LINK %zu :", i); for (auto& kv : A.link) fprintf(of, " %d:%d", kv.first, kv.second); fprintf(of, "\n");
    }
    fprintf(of, "TOTAL orbits %zu weighted %lld", O.size(), W3); for (int t = 0; t < TYPES; t++) fprintf(of, " %s %lld", TNAME[t], nT[t]); fprintf(of, " reduced %lld maximal %lld\n", L, M);
    fclose(of);
    bool ok = !badSum && !inconsistent && !idbad && sumAll == W3 && chk(SQ) == 4 * W4 && !tbad && !unk && n11sum == SPR && !linkBad && !M;
    printf("merge: %s\n", ok ? "all identities hold" : "FAILURES");
    return ok ? 0 : 2;
}

// ================================================================================ configs: hole sets, covers, sources
static const int NFAM = 3;
struct OrbRec { int idx, type; long long aut; int q; vector<u16> Ns; int hole = -1; int src[NFAM] = { 0 }; int srcCov = 0; };
static vector<OrbRec> read_orbs(const char* file) {
    vector<OrbRec> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "ORB ", 4)) continue;
        OrbRec o; char ty[16]; long long w, n;
        if (sscanf(line, "ORB %d %15s aut %lld w %lld n %lld q %d", &o.idx, ty, &o.aut, &w, &n, &o.q) != 6) { printf("bad ORB line\n"); exit(1); }
        o.type = -1; for (int t = 0; t < TYPES; t++) if (!strcmp(ty, TNAME[t])) o.type = t;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; o.Ns.push_back((u16)v); p = e; }
        if (o.Ns.size() != KN || o.type < 0 || o.idx != (int)out.size()) { printf("bad ORB line\n"); exit(1); }
        out.push_back(o);
    }
    fclose(fp); return out;
}
static Mask holes_of(const vector<Mask>& M) { Mask cov; mclear(cov); for (auto& m : M) mor(cov, m); return mcompl(cov); }
static vector<vector<u16>> read_spread_file(const char* file) {
    vector<vector<u16>> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "SPREAD", 6)) continue;
        const char* p = strchr(line, ':'); if (!p) continue; p++;
        vector<u16> Ns; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
        if (Ns.size() != 14) { printf("bad spread line\n"); exit(1); }
        out.push_back(Ns);
    }
    fclose(fp); return out;
}
// every way of removing k members from each base set (bases in order, k-subsets in lexicographic order)
static vector<vector<Mask>> drop_k(const vector<vector<Mask>>& bases, int k) {
    vector<vector<Mask>> out;
    for (auto& M : bases) {
        int n = (int)M.size(); vector<int> c(k); for (int i = 0; i < k; i++) c[i] = i;
        while (true) {
            vector<Mask> r; int j = 0; for (int i = 0; i < n; i++) { if (j < k && c[j] == i) { j++; continue; } r.push_back(M[i]); }
            out.push_back(r);
            int i = k - 1; while (i >= 0 && c[i] == n - k + i) i--;
            if (i < 0) break;
            c[i]++; for (int t = i + 1; t < k; t++) c[t] = c[t - 1] + 1;
        }
    }
    return out;
}
// the source families (each a list of 6-member sets), independent of this file's classification:
//   0: seven-check's 3731 configurations x 7 members removed one at a time (every 6 that extends at all, which is
//      every 6 if none is maximal);
//   1: eight-check's 9372 configurations x C(8,2) (every 6 inside some 8);
//   2: fourteen-check's 73 (spread, X Y Z) reps x C(14,8): 6 of the 14 graph members (every 6 inside a spread: of the
//      11 spread members a 6 misses, any 3 can be moved to X, Y, Z).
static vector<vector<Mask>> familyPrev(const vector<CPrev>& c, int k) { vector<vector<Mask>> b; for (auto& x : c) b.push_back(members_of(x.Ns)); return drop_k(b, k); }
static vector<vector<Mask>> family7(const vector<CPrev>& c7) { return familyPrev(c7, 1); }
static vector<vector<Mask>> family8(const vector<CPrev>& c8) { return familyPrev(c8, 2); }
static vector<vector<Mask>> familySp(const vector<vector<u16>>& sp) {
    vector<vector<Mask>> b; for (auto& Ns : sp) { vector<Mask> M; for (u16 N : Ns) M.push_back(graph_mask(N)); b.push_back(M); }
    return drop_k(b, 8);
}
static const char* FAMNAME[NFAM] = { "seven-check's 3731 x 7", "eight-check's 9372 x 28", "fourteen-check's spread reps x 3003" };
static int mode_configs(const char* orbfile, const char* conf7, const char* conf8, const char* spreadfile, const char* outfile) {
    init_conj(); init_solids();
    printf("solids in PG(7,2): %zu (expect 200787)\n", ALLSOLIDS.size()); if (ALLSOLIDS.size() != 200787) return 1;
    auto orbs = read_orbs(orbfile); int n = (int)orbs.size();
    map<vector<u16>, int> byKey; for (int i = 0; i < n; i++) byKey[key_of_Ns(orbs[i].Ns)] = i;
    printf("orbits read: %d\n", n); fflush(stdout);
    // --- the stored key is a canonical form: recompute it (and |Aut|) from a random image
    { int bad = 0; mt19937 rng(99);
      vector<Lin> gs(n); vector<unsigned> seeds(n); for (int i = 0; i < n; i++) { gs[i] = random_gl8(rng); seeds[i] = rng(); }
#pragma omp parallel for schedule(dynamic, 8) reduction(+:bad)
      for (int i = 0; i < n; i++) {
          vector<Mask> M = members_of(orbs[i].Ns), Mg; for (auto& m : M) Mg.push_back(mask_image(gs[i], m));
          mt19937 r2(seeds[i]); shuffle(Mg.begin(), Mg.end(), r2);
          Canon K = canonP(Mg); if (K.key != key_of_Ns(orbs[i].Ns) || K.aut != orbs[i].aut) bad++;
      }
      printf("canonical form recomputed from a random GL(8,2) image of every orbit: %d mismatches\n", bad); fflush(stdout); if (bad) return 1; }
    // --- hole sets
    vector<Mask> H(n); vector<string> hk(n);
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; i++) { H[i] = holes_of(members_of(orbs[i].Ns)); hk[i] = hole_key(H[i]); }
    vector<Mask> holeRep; vector<string> holeKey; vector<long long> holeAut; vector<int> holeSolid, holeFirst;
    auto t0 = chrono::steady_clock::now();
    {   map<string, vector<int>> groups; for (int i = 0; i < n; i++) groups[hk[i]].push_back(i);
        printf("distinct hole invariants: %zu\n", groups.size()); fflush(stdout);
        vector<pair<string, vector<int>>> gv(groups.begin(), groups.end());
        vector<vector<pair<int, vector<int>>>> res(gv.size());
#pragma omp parallel for schedule(dynamic, 1)
        for (int g = 0; g < (int)gv.size(); g++) {
            vector<pair<int, vector<int>>> reps;
            for (int i : gv[g].second) {
                bool found = false;
                for (auto& rp : reps) if (iso_count(H[i], H[rp.first], true, nullptr)) { rp.second.push_back(i); found = true; break; }
                if (!found) reps.push_back({ i, { i } });
            }
            res[g] = reps;
        }
        for (size_t g = 0; g < gv.size(); g++) for (auto& rp : res[g]) {
            int h = (int)holeRep.size(); holeRep.push_back(H[rp.first]); holeKey.push_back(gv[g].first); holeFirst.push_back(rp.first);
            for (int i : rp.second) orbs[i].hole = h;
        }
        int nh = (int)holeRep.size(); holeAut.assign(nh, 0); holeSolid.assign(nh, 0);
        printf("hole classes by isomorphism: %d (%.0f s); computing |Aut H|\n", nh, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
#pragma omp parallel for schedule(dynamic, 1)
        for (int h = 0; h < nh; h++) { holeAut[h] = iso_count(holeRep[h], holeRep[h], false, nullptr); holeSolid[h] = contains_solid(holeRep[h]); }
    }
    int nh = (int)holeRep.size(), nhMax = 0;
    for (int h = 0; h < nh; h++) nhMax += !holeSolid[h];
    printf("hole sets up to GL(8,2): %d, containing no solid (maximal partial spreads): %d   [Kurz Table 9.3, r = 6: 341 and 0] (%.0f s)\n", nh, nhMax, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    { int bad = 0; for (auto& o : orbs) if ((o.type == TYPES - 1) != !holeSolid[o.hole]) bad++; printf("orbits whose type disagrees with 'hole set contains a solid': %d\n", bad); }
    // --- exact covers of each hole set's complement by 6 solids: every cover must canonicalise onto an orbit with
    // that hole set, all of them reached, and sum over those orbits of |Aut H| / |Aut P| = number of covers
    vector<vector<int>> orbsOfHole(nh); for (int i = 0; i < n; i++) orbsOfHole[orbs[i].hole].push_back(i);
    struct CovRes { int ns = 0, m = 0; long long unk = 0, wrongHole = 0, expect = 0; vector<int> hits; bool ok = false; double secs = 0; };
    vector<CovRes> cr(nh);
    t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int h = 0; h < nh; h++) {
        auto t1 = chrono::steady_clock::now();
        CovRes& R = cr[h];
        auto cv = covers_of(holeRep[h], R.ns);
        R.m = (int)cv.size(); set<int> got;
        for (int i = 0; i < R.m; i++) {
            Canon K = canonP(cv[i]); auto it = byKey.find(K.key);
            if (it == byKey.end()) R.unk++; else if (orbs[it->second].hole != h) R.wrongHole++; else { got.insert(it->second); R.hits.push_back(it->second); }
        }
        for (int i : orbsOfHole[h]) R.expect += holeAut[h] % orbs[i].aut ? -1000000000LL : holeAut[h] / orbs[i].aut;
        R.ok = !R.unk && !R.wrongHole && got.size() == orbsOfHole[h].size() && R.expect == R.m;
        R.secs = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
    }
    long long covBad = 0, covAll = 0;
    for (int h = 0; h < nh; h++) {
        CovRes& R = cr[h]; covAll += R.m; covBad += !R.ok;
        for (int i : R.hits) orbs[i].srcCov++;
        map<int, int> types; for (int i : orbsOfHole[h]) types[orbs[i].type]++;
        printf("hole set %d: |Aut H| %lld, solid inside: %s, %d solids in the complement, %d covers, %zu orbits (", h, holeAut[h], holeSolid[h] ? "yes" : "no", R.ns, R.m, orbsOfHole[h].size());
        for (auto& kv : types) printf(" %s:%d", TNAME[kv.first], kv.second);
        printf(" ); sum |Aut H|/|Aut P| = %lld; on unknown orbits %lld, on another hole set %lld  %s [%.1f s]\n    %s\n", R.expect, R.unk, R.wrongHole, R.ok ? "OK" : "MISMATCH", R.secs, holeKey[h].c_str());
    }
    printf("hole sets with a cover mismatch: %lld; covers in all %lld (%.0f s)\n", covBad, covAll, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    { int unreached = 0; for (auto& o : orbs) unreached += !o.srcCov; printf("orbits reached by no cover: %d\n", unreached); covBad += unreached; }
    fflush(stdout);
    // --- independent sources
    auto run_family = [&](const vector<vector<Mask>>& fam, int which) {
        int m = (int)fam.size(); vector<int> hit(m, -1);
#pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < m; i++) { Canon K = canonP(fam[i]); auto it = byKey.find(K.key); if (it != byKey.end()) hit[i] = it->second; }
        set<int> got; int unk = 0; map<int, int> types;
        for (int i = 0; i < m; i++) { if (hit[i] < 0) { unk++; continue; } got.insert(hit[i]); types[orbs[hit[i]].type]++; orbs[hit[i]].src[which]++; }
        printf("%s: %d configurations; on no listed orbit %d; by type:", FAMNAME[which], m, unk); for (auto& kv : types) printf(" %s:%d", TNAME[kv.first], kv.second);
        int reach[TYPES] = { 0 }, all[TYPES] = { 0 }; for (int i = 0; i < n; i++) { all[orbs[i].type]++; reach[orbs[i].type] += got.count(i) > 0; }
        printf("; orbits reached:"); for (int t = 0; t < TYPES; t++) printf(" %s %d/%d", TNAME[t], reach[t], all[t]); printf("\n"); fflush(stdout);
        return unk;
    };
    int srcUnk = 0;
    srcUnk += run_family(family7(read_conf7(conf7)), 0);
    srcUnk += run_family(family8(read_conf8(conf8)), 1);
    srcUnk += run_family(familySp(read_spread_file(spreadfile)), 2);
    { int unreached7 = 0; for (auto& o : orbs) unreached7 += !o.src[0]; printf("orbits reached by no deletion from seven-check's 3731: %d\n", unreached7); srcUnk += unreached7; }
    // --- output
    FILE* of = fopen(outfile, "w"); int idx = 0;
    for (int t = 0; t < TYPES; t++) for (auto& o : orbs) if (o.type == t) {
        fprintf(of, "CONF %d %s hole %d aut %lld src", idx++, TNAME[t], o.hole, o.aut);
        for (int k = 0; k < NFAM; k++) fprintf(of, " %d", o.src[k]); fprintf(of, " %d :", o.srcCov);
        for (u16 N : o.Ns) fprintf(of, " %04x", N); fprintf(of, "\n");
    }
    for (int h = 0; h < nh; h++) fprintf(of, "HOLE %d aut %lld solid %d covers %d solids %d | %s\n", h, holeAut[h], holeSolid[h], cr[h].m, cr[h].ns, holeKey[h].c_str());
    fclose(of);
    printf("written %d configurations to %s\n", idx, outfile);
    return (covBad || srcUnk) ? 2 : 0;
}


// ================================================================================ forms: the search
static u32 code28(u64 f) { u32 c = 0; int p = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { if (rowof(f, i) >> j & 1) c |= 1u << p; p++; } return c; }
static u64 E28[28];
static u64 TLO[1 << 14], THI[1 << 14];
static void init_code28() {
    int p = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { E28[p++] = (1ULL << (8 * i + j)) | (1ULL << (8 * j + i)); }
    for (u32 c = 0; c < (1u << 14); c++) { u64 a = 0, b = 0; for (int k = 0; k < 14; k++) if (c >> k & 1) { a ^= E28[k]; b ^= E28[14 + k]; } TLO[c] = a; THI[c] = b; }
}
static inline u64 form28(u32 c) { return TLO[c & 0x3fff] ^ THI[c >> 14]; }
static vector<u32> rref28(vector<u32> v) {
    vector<u32> out;
    for (int bit = 27; bit >= 0; bit--) {
        int p = -1; for (size_t i = 0; i < v.size(); i++) if (v[i] >> bit & 1) { p = (int)i; break; }
        if (p < 0) continue;
        u32 r = v[p]; v.erase(v.begin() + p);
        for (auto& x : v) if (x >> bit & 1) x ^= r;
        for (auto& x : out) if (x >> bit & 1) x ^= r;
        out.push_back(r);
    }
    return out;
}
static string space_key(const vector<u64>& basis) {
    vector<u32> c; for (u64 f : basis) c.push_back(code28(f));
    auto r = rref28(c); sort(r.begin(), r.end());
    string s; char buf[16]; for (u32 x : r) { sprintf(buf, "%07x ", x); s += buf; } return s;
}
static bool independent_of(const vector<u64>& basis, u64 f) {
    vector<u32> c; for (u64 b : basis) c.push_back(code28(b)); size_t r = rref28(c).size(); c.push_back(code28(f)); return rref28(c).size() > r;
}
struct Verdict { bool apn; int dist[9]; };
static Verdict verify_space(const vector<u64>& basis) { // from scratch: LUT, derivative counts, Walsh spectrum
    Verdict V; memset(&V, 0, sizeof V);
    int n = (int)basis.size(); u8 F[256];
    for (int x = 0; x < 256; x++) {
        u8 y = 0;
        for (int i = 0; i < n; i++) { int s = 0; for (int j = 0; j < 8; j++) for (int k = j + 1; k < 8; k++) if ((x >> j & 1) && (x >> k & 1) && (rowof(basis[i], j) >> k & 1)) s ^= 1; y |= s << i; }
        F[x] = y;
    }
    V.apn = (n == 8);
    for (int a = 1; a < 256 && V.apn; a++) { int cnt[256] = { 0 }; for (int x = 0; x < 256; x++) if (++cnt[F[x ^ a] ^ F[x]] > 2) { V.apn = false; break; } }
    for (int b = 1; b < (1 << n); b++) {
        int W[256]; for (int x = 0; x < 256; x++) W[x] = par(b & F[x]) ? -1 : 1;
        for (int h = 1; h < 256; h <<= 1) for (int i = 0; i < 256; i += 2 * h) for (int j = i; j < i + h; j++) { int a = W[j], c = W[j + h]; W[j] = a + c; W[j + h] = a - c; }
        int mx = 0; for (int i = 0; i < 256; i++) mx = max(mx, abs(W[i]));
        int k = 0; while ((1 << (8 + k)) < mx * mx) k++;
        V.dist[k]++;
    }
    return V;
}
static u64 transform(u64 f, const u8 p[8]) { // B'(e_i, e_j) = B(p_i, p_j)
    u64 g = 0; for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (i != j && evalf(f, p[i], p[j])) g |= 1ULL << (8 * i + j); return g;
}
static void change_coords(vector<u64>& T, const Lin& A) { Lin Ai = lin_inv(A); for (auto& f : T) f = transform(f, Ai.data()); } // w = A v

// The problem: 6 rank-4 radicals (members), holes D = the rest (165 points). Relaxation for controls: oracle rank-6
// forms. Forms on the lines of D are generated on first use (most configurations never branch in D).
struct Problem {
    int nr = 0; Mask R[16]; int radOf[256]; Mask D;
    vector<u64> f4[16];
    vector<Mask> lineM; vector<array<u8, 6>> lineAnn; mutable vector<vector<u64>> lineF; unique_ptr<once_flag[]> lineOnce; vector<int> linesAt[256];
    vector<int> order;
    vector<u64> oracle; vector<Mask> oracleM;
};
static void build_problem(Problem& P, const vector<Mask>& M) {
    P.nr = (int)M.size(); if (P.nr > 16) { fprintf(stderr, "too many members\n"); exit(1); }
    Mask cov; mclear(cov);
    for (int v = 0; v < 256; v++) P.radOf[v] = -1;
    for (int j = 0; j < P.nr; j++) {
        P.R[j] = M[j];
        if (mcount(M[j]) != 15 || !mdisj(cov, M[j])) { fprintf(stderr, "members not a partial spread\n"); exit(1); }
        mor(cov, M[j]);
        for (int v = 1; v < 256; v++) if (mget(M[j], v)) P.radOf[v] = j;
        u8 u[4]; annihilator(M[j], 4, u);
        P.f4[j].clear(); forms_from(u, 4, ND4, P.f4[j]);
    }
    P.D = mcompl(cov);
    set<Mask> seen;
    P.lineM.clear(); P.lineAnn.clear(); for (int v = 0; v < 256; v++) P.linesAt[v].clear();
    for (int a = 1; a < 256; a++) for (int b = a + 1; b < 256; b++) if (mget(P.D, a) && mget(P.D, b) && mget(P.D, a ^ b)) {
        Mask m; mclear(m); mset(m, a); mset(m, b); mset(m, a ^ b);
        if (!seen.insert(m).second) continue;
        int li = (int)P.lineM.size(); P.lineM.push_back(m);
        array<u8, 6> ann; annihilator(m, 6, ann.data()); P.lineAnn.push_back(ann);
        for (int v : { a, b, a ^ b }) P.linesAt[v].push_back(li);
    }
    P.lineF.assign(P.lineM.size(), {});
    P.lineOnce.reset(new once_flag[P.lineM.size() + 1]);
}
static const vector<u64>& line_forms(const Problem& P, int li) {
    call_once(P.lineOnce[li], [&]() { vector<u64> fs; forms_from(P.lineAnn[li].data(), 6, ND6, fs); P.lineF[li].swap(fs); });
    return P.lineF[li];
}
static inline bool legal(const Problem& P, u64 f, Mask& m, int& cnt) {
    if (rank8(f) == 8) { cnt = 0; return true; }
    cnt = radical(f, m);
    if (cnt == 15) { int v = 0; for (v = 1; v < 256; v++) if (mget(m, v)) break; int j = P.radOf[v]; return j >= 0 && meq(m, P.R[j]); }
    if (cnt == 3) {
        if (msub(m, P.D)) return true;
        for (auto& o : P.oracleM) if (meq(o, m)) return true;
        return false;
    }
    return false;
}

// ================================================================================ the coset sweep (a second D phase)
// Written for seven-check, sharing nothing with the point-branching code above but the input masks: its own kernel
// routine, its own coverage bookkeeping, no line forms and no branching point. From a legal span S of dimension d
// (all 2^d - 1 nonzero elements legal, radicals pairwise disjoint), every 8-dimensional legal T containing S is reached
// by adding, one dimension at a time, every coset representative f of the current span (f supported on the non-pivot
// coordinates of its reduced basis: one per coset), keeping f when every element of f + span is nondegenerate or rank 6
// with its radical a line inside D (or an oracle line, for controls) disjoint from everything covered so far. Rank 4 is
// never accepted: the sweep starts at a core, which already holds every rank-4 form (all off-D points covered). An
// 8-dimensional space covers every point at least once (dim ker(B -> B(a, .)) >= 8 - 7), so a legal T of dimension 8 is
// an APN space: the sweep's dimension-8 finds are solutions, nothing else to check.
struct SwCtx { bool inD[256]; vector<Mask> oracle; };
static int sw_kernel(u64 f, u8 ker[8]) { // kernel basis of the symmetric matrix f (rows = bytes), by elimination with combinations
    u8 r[8], c[8]; for (int i = 0; i < 8; i++) { r[i] = (u8)(f >> (8 * i)); c[i] = (u8)(1 << i); }
    int row = 0;
    for (int bit = 0; bit < 8; bit++) {
        int p = -1; for (int i = row; i < 8; i++) if (r[i] >> bit & 1) { p = i; break; }
        if (p < 0) continue;
        swap(r[p], r[row]); swap(c[p], c[row]);
        for (int i = 0; i < 8; i++) if (i != row && (r[i] >> bit & 1)) { r[i] ^= r[row]; c[i] ^= c[row]; }
        row++;
    }
    int k = 0; for (int i = row; i < 8; i++) ker[k++] = c[i];
    return k;
}
// radical points of f into pts (count returned: 0, 3, 15, ...)
// Nondegenerate <=> Pfaffian = 1 (over GF(2), det = Pf^2). Pf(A) = sum_j a_0j Pf(A without rows/cols 0, j) mod 2, the
// 6x6 Pfaffians from a 2^15 table (built by elimination), each minor's 15 entries gathered by pext.
static u8 SW_PF6[1 << 15]; static u64 SW_MINOR[8];
static void sw_init() {
    for (u32 c = 0; c < (1u << 15); c++) {
        u8 r[6] = { 0 }; int p = 0;
        for (int i = 0; i < 6; i++) for (int j = i + 1; j < 6; j++) { if (c >> p & 1) { r[i] |= 1 << j; r[j] |= 1 << i; } p++; }
        int rk = 0;
        for (int bit = 0; bit < 6; bit++) {
            int q = rk; while (q < 6 && !(r[q] >> bit & 1)) q++;
            if (q == 6) continue;
            swap(r[q], r[rk]); for (int i = rk + 1; i < 6; i++) if (r[i] >> bit & 1) r[i] ^= r[rk]; rk++;
        }
        SW_PF6[c] = rk == 6;
    }
    for (int j = 1; j < 8; j++) { u64 m = 0; for (int a = 1; a < 8; a++) for (int b = a + 1; b < 8; b++) if (a != j && b != j) m |= 1ULL << (8 * a + b); SW_MINOR[j] = m; }
}
static inline bool sw_nondeg(u64 f) {
    u32 s = 0; u8 r0 = (u8)f;
    for (int j = 1; j < 8; j++) if (r0 >> j & 1) s ^= SW_PF6[_pext_u64(f, SW_MINOR[j])];
    return s & 1;
}
static int sw_radical(u64 f, u8 pts[255]) {
    if (sw_nondeg(f)) return 0;
    u8 ker[8]; int k = sw_kernel(f, ker); int n = 0;
    for (int s = 1; s < (1 << k); s++) { u8 v = 0; for (int i = 0; i < k; i++) if (s >> i & 1) v ^= ker[i]; pts[n++] = v; }
    return n;
}
static thread_local u32 sw_mark[256], sw_epoch = 0; // points covered inside the coset being checked: sw_mark == sw_epoch
static bool sw_line_ok(const SwCtx& C, const u8 pts[3], const bool* cov) {
    for (int i = 0; i < 3; i++) if (cov[pts[i]] || sw_mark[pts[i]] == sw_epoch) return false;
    if (C.inD[pts[0]] && C.inD[pts[1]] && C.inD[pts[2]]) return true;
    for (auto& m : C.oracle) if (mget(m, pts[0]) && mget(m, pts[1]) && mget(m, pts[2])) return true;
    return false;
}
struct SwStats { long long found[9] = { 0 }, checks = 0, cosets = 0; vector<vector<u64>> sols; };
// K_a = {B : B a = 0} (21-dimensional), as 28-bit codes: the null space of the 8 conditions (B a)_i = 0
static int sw_kernel_space(u8 a, u32 out[28]) {
    u32 cond[8]; int p = 0; for (int i = 0; i < 8; i++) cond[i] = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { if (a >> j & 1) cond[i] |= 1u << p; if (a >> i & 1) cond[j] |= 1u << p; p++; }
    int pivOf[28]; for (int b = 0; b < 28; b++) pivOf[b] = -1;
    u32 rows[8]; int nr = 0;
    for (int i = 0; i < 8; i++) { // full reduction
        u32 r = cond[i]; for (int k = 0; k < nr; k++) { unsigned long hb; _BitScanReverse(&hb, rows[k]); if (r >> hb & 1) r ^= rows[k]; }
        if (!r) continue;
        unsigned long hb; _BitScanReverse(&hb, r);
        for (int k = 0; k < nr; k++) if (rows[k] >> hb & 1) rows[k] ^= r;
        rows[nr++] = r;
    }
    for (int k = 0; k < nr; k++) { unsigned long hb; _BitScanReverse(&hb, rows[k]); pivOf[hb] = k; }
    int n = 0;
    for (int b = 0; b < 28; b++) if (pivOf[b] < 0) { // free coordinate b: set it, solve the pivots
        u32 v = 1u << b; for (int k = 0; k < nr; k++) if (rows[k] >> b & 1) { unsigned long hb; _BitScanReverse(&hb, rows[k]); v |= 1u << hb; }
        out[n++] = v;
    }
    return n;
}
// Below dimension 7, only the extensions S + <g> with g in K_a for one uncovered point a (the smallest) are followed: in
// an 8-dimensional legal T containing S exactly one element g has a in its radical, and S + <g> is legal (a subspace of
// T), so every T is still reached; K_a meets S in 0, so each g gives a different S + <g>. At dimension 7 (or when every
// point is covered) every coset is tried; at dimension 7 that is the same set as K_a mod S for any uncovered a (K_a has
// dimension 21 and meets S in 0). The smallest point is the one `search asc` branches on (all off-D points are covered
// at a core), so in the asc order the sweep must visit exactly the spans the point branching visits: nodes by dimension
// (point branching) = cores + sweep's legal spans, dimension by dimension, though the candidates come from K_a, not
// from the line forms.
// Checks the coset f + span: every element nondegenerate or rank 6 with its radical a legal line, the radicals pairwise
// disjoint and disjoint from cov; the new radicals' points are left marked with sw_epoch.
static bool sw_coset_ok(const SwCtx& C, u64 f, const vector<u64>& span, const bool* cov, SwStats& st) {
    if (++sw_epoch == 0) { memset(sw_mark, 0, sizeof sw_mark); sw_epoch = 1; }
    for (size_t i = 0; i < span.size(); i++) {
        st.checks++;
        u8 pts[255]; int np = sw_radical(f ^ span[i], pts);
        if (np == 0) continue;
        if (np != 3 || !sw_line_ok(C, pts, cov)) return false;
        for (int j = 0; j < 3; j++) sw_mark[pts[j]] = sw_epoch;
    }
    return true;
}
// For a core of dimension <= 5 (whose D phase is a deep tree): per point a, the g in K_a with core + <g> legal. Every
// descendant S of the core has core + <g> inside S + <g>, so its candidates at a are among these (seven-check: the one
// dimension-4 core and the 1,634 of dimension 5 made the sweep, which tries all of K_a at every node, too slow).
// With a cache, dimension 7 also takes K_a for its smallest uncovered point (the same cosets, as noted below).
struct SwCache { vector<u64> span; bool cov[256]; bool built[256]; vector<u32> L[256]; };
static void sw_rec(const SwCtx& C, const vector<u64>& span, const vector<u32>& codes, const bool* cov, int d, SwStats& st, SwCache* K) {
    int a = -1; for (int v = 1; v < 256 && a < 0; v++) if (!cov[v]) a = v;
    size_t n = span.size();
    auto accept = [&](u64 f, u32 fc) {
        st.found[d + 1]++;
        bool nc[256]; for (int v = 0; v < 256; v++) nc[v] = cov[v] || sw_mark[v] == sw_epoch;
        vector<u64> sp2(2 * n); vector<u32> cd2(2 * n);
        for (size_t i = 0; i < n; i++) { sp2[i] = span[i]; sp2[n + i] = span[i] ^ f; cd2[i] = codes[i]; cd2[n + i] = codes[i] ^ fc; }
        if (d + 1 == 8) { vector<u64> B; for (int i = 0; i < 8; i++) B.push_back(sp2[(size_t)1 << i]); st.sols.push_back(B); }
        else sw_rec(C, sp2, cd2, nc, d + 1, st, K);
    };
    if (K && a > 0) {
        if (!K->built[a]) {
            u32 gen[28]; int nf = sw_kernel_space((u8)a, gen);
            if (nf != 21) { fprintf(stderr, "sweep: K_a of dimension %d\n", nf); exit(1); }
            u32 fc = 0;
            for (u64 g = 1; g < (1ULL << nf); g++) { fc ^= gen[_tzcnt_u64(g)]; st.cosets++; if (sw_coset_ok(C, form28(fc), K->span, K->cov, st)) K->L[a].push_back(fc); }
            K->built[a] = true;
        }
        for (u32 fc : K->L[a]) { st.cosets++; u64 f = form28(fc); if (sw_coset_ok(C, f, span, cov, st)) accept(f, fc); }
        return;
    }
    u32 gen[28]; int nf = 0;
    if (d < 7 && a > 0) {
        nf = sw_kernel_space((u8)a, gen);
        if (nf != 21) { fprintf(stderr, "sweep: K_a of dimension %d\n", nf); exit(1); }
    } else { // coset representatives: supported on the non-pivot coordinates of the span's reduced basis
        vector<u32> red; { vector<u32> v; for (int i = 0; i < d; i++) v.push_back(codes[(size_t)1 << i]); red = rref28(v); }
        if ((int)red.size() != d) { fprintf(stderr, "sweep: basis not independent\n"); exit(1); }
        u32 piv = 0; for (u32 b : red) { unsigned long hb; _BitScanReverse(&hb, b); piv |= 1u << hb; }
        for (int b = 0; b < 28; b++) if (!(piv >> b & 1)) gen[nf++] = 1u << b;
    }
    u32 fc = 0; // Gray code over the 2^nf - 1 nonzero combinations of gen
    for (u64 g = 1; g < (1ULL << nf); g++) {
        int b = (int)_tzcnt_u64(g); fc ^= gen[b];
        u64 f = form28(fc);
        st.cosets++;
        if (sw_coset_ok(C, f, span, cov, st)) accept(f, fc);
    }
}
// Sweep from the span of `basis`: checks the start itself with the sweep's own primitives (returns false if it isn't
// legal: a rank-4 element whose radical isn't a member, or overlapping radicals), then sweeps.
static bool sweep_from(const SwCtx& C, const vector<Mask>& members, const vector<u64>& basis, SwStats& st) {
    int d = (int)basis.size();
    vector<u64> span = { 0 }; vector<u32> codes = { 0 };
    for (u64 b : basis) { size_t n = span.size(); for (size_t i = 0; i < n; i++) { span.push_back(span[i] ^ b); codes.push_back(codes[i] ^ code28(b)); } }
    bool cov[256] = { false };
    if (++sw_epoch == 0) { memset(sw_mark, 0, sizeof sw_mark); sw_epoch = 1; }
    for (size_t i = 1; i < span.size(); i++) {
        u8 pts[255]; int np = sw_radical(span[i], pts);
        if (np == 0) { if (span[i] == 0) return false; continue; }
        if (np == 15) { Mask m; mclear(m); for (int j = 0; j < 15; j++) mset(m, pts[j]); bool mem = false; for (auto& M : members) if (meq(M, m)) mem = true; if (!mem) return false; }
        else if (np == 3) { if (!sw_line_ok(C, pts, cov)) return false; }
        else return false;
        for (int j = 0; j < np; j++) { if (cov[pts[j]]) return false; cov[pts[j]] = true; }
    }
    if (d == 8) { st.sols.push_back(basis); return true; }
    unique_ptr<SwCache> K;
    if (d <= 5) { K.reset(new SwCache()); K->span = span; memcpy(K->cov, cov, 256); for (int v = 0; v < 256; v++) K->built[v] = false; }
    sw_rec(C, span, codes, cov, d, st, K.get());
    return true;
}
static SwCtx sw_ctx(const vector<Mask>& members, const vector<Mask>& oracleLines) {
    SwCtx C; for (int v = 0; v < 256; v++) C.inD[v] = v != 0;
    for (auto& M : members) for (int v = 1; v < 256; v++) if (mget(M, v)) C.inD[v] = false;
    C.oracle = oracleLines; return C;
}

struct State { vector<u64> span; Mask cov; int d; vector<u64> basis; bool pp = false; }; // pp: the parent was a p1 node
static State empty_state() { State S; S.span = { 0 }; mclear(S.cov); S.d = 0; return S; }
static bool try_add(const Problem& P, const State& S, u64 f, State& out) {
    Mask nc = S.cov; size_t n = S.span.size();
    for (size_t i = 0; i < n; i++) {
        Mask m; int c;
        if (!legal(P, f ^ S.span[i], m, c)) return false;
        if (c) { if (!mdisj(m, nc)) return false; mor(nc, m); }
    }
    out.span.resize(2 * n);
    for (size_t i = 0; i < n; i++) { out.span[i] = S.span[i]; out.span[n + i] = S.span[i] ^ f; }
    out.cov = nc; out.d = S.d + 1; out.basis = S.basis; out.basis.push_back(f);
    return true;
}
template <class CB> static void for_each_candidate(const Problem& P, int a, const Mask& cov, CB cb) {
    int j = P.radOf[a];
    if (j >= 0) for (u64 f : P.f4[j]) cb(f);
    else for (int li : P.linesAt[a]) if (mdisj(P.lineM[li], cov)) for (u64 f : line_forms(P, li)) cb(f);
    for (size_t i = 0; i < P.oracle.size(); i++) if (mget(P.oracleM[i], a) && !msub(P.oracleM[i], P.D) && mdisj(P.oracleM[i], cov)) cb(P.oracle[i]);
}
// p1 = nodes where every point off D is covered (all 7 rank-4 forms present); triesD = candidate forms tried at
// points of D (the rank-6 branching, which only happens at p1 nodes). A core is a p1 node whose parent is not: the span
// of the rank-4 forms of everything below it. Cores don't depend on the point order (each legal span spanned by one
// rank-4 form per member is reached by exactly one path in any order that puts the off-D points first), so their
// number by dimension and their fingerprint (sum of 64-bit FNV hashes of their reduced bases) must agree between orders.
struct Stats { long long nodes[10] = { 0 }, tries[10] = { 0 }, triesD[10] = { 0 }, bentext = 0, full = 0, p1[10] = { 0 }, cores[10] = { 0 }, swFound[10] = { 0 }, swChecks = 0, swSols = 0; u64 fp = 0; };
static void add_stats(Stats& a, const Stats& b) {
    for (int i = 0; i < 10; i++) { a.nodes[i] += b.nodes[i]; a.tries[i] += b.tries[i]; a.triesD[i] += b.triesD[i]; a.p1[i] += b.p1[i]; a.cores[i] += b.cores[i]; a.swFound[i] += b.swFound[i]; }
    a.bentext += b.bentext; a.full += b.full; a.swChecks += b.swChecks; a.swSols += b.swSols; a.fp += b.fp;
}
static mutex solmu;
static map<string, vector<u64>> solutions;
static map<string, vector<u64>> p1spaces; // distinct p1 spans (capped)
static bool g_logp1 = true;
static bool g_phase1only = false, g_nocache = getenv("SIX_NOCACHE") != nullptr; static int g_dmax = getenv("SIX_DMAX") ? atoi(getenv("SIX_DMAX")) : 99;
static bool g_dsweep = false;               // at a core, run the coset sweep instead of the point branching in D
static FILE* g_coreFile = nullptr; static int g_coreConf = -1; // core list output (K lines)
static thread_local map<string, vector<u64>>* tl_sol = nullptr;
static thread_local const SwCtx* tl_sw = nullptr; static thread_local const vector<Mask>* tl_members = nullptr;
static u64 fnv(const string& s) { u64 h = 1469598103934665603ULL; for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; } return h; }
static void record(const State& S, Stats& st) {
    st.full++;
    string k = space_key(S.basis);
    if (tl_sol) { (*tl_sol)[k] = S.basis; return; }
    lock_guard<mutex> g(solmu);
    if (!solutions.count(k)) { solutions[k] = S.basis; printf("SOLUTION %s\n", k.c_str()); fflush(stdout); }
}
static void record_sweep(const vector<u64>& B, Stats& st) {
    st.swSols++;
    string k = space_key(B);
    if (tl_sol) { (*tl_sol)[k] = B; return; }
    lock_guard<mutex> g(solmu);
    if (!solutions.count(k)) { solutions[k] = B; printf("SOLUTION (sweep) %s\n", k.c_str()); fflush(stdout); }
}
static void bent_extend(const State& S, Stats& st) {
    st.bentext++;
    vector<u32> c; for (u64 f : S.basis) c.push_back(code28(f));
    auto B = rref28(c);
    u32 piv = 0; for (u32 b : B) { unsigned long hb; _BitScanReverse(&hb, b); piv |= 1u << hb; }
    u32 freemask = ((1u << 28) - 1) & ~piv; int nf = __popcnt(freemask);
    for (u64 t = 1; t < (1ULL << nf); t++) {
        u64 f = form28(_pdep_u32((u32)t, freemask));
        bool ok = true; for (u64 s : S.span) if (rank8(f ^ s) != 8) { ok = false; break; }
        if (!ok) continue;
        State T; T.span.resize(2 * S.span.size());
        for (size_t i = 0; i < S.span.size(); i++) { T.span[i] = S.span[i]; T.span[S.span.size() + i] = S.span[i] ^ f; }
        T.cov = S.cov; T.d = S.d + 1; T.basis = S.basis; T.basis.push_back(f);
        if (T.d == 8) record(T, st); else bent_extend(T, st);
    }
}
static int first_uncovered(const Problem& P, const State& S) { for (int v : P.order) if (!mget(S.cov, v)) return v; return -1; }
// returns true if this node is a core: it is collected, and run_search's second stage does its D phase (point branching
// over the lines through its first uncovered point, one task per line, or the coset sweep), for better parallelism
static thread_local vector<State>* tl_coll = nullptr; static thread_local mutex* tl_collmu = nullptr;
// guided controls: the planted space T (all 256 elements, sorted); the full off-D phase runs, but only cores inside T get
// their D phase (each control then costs one core's D phase instead of the whole configuration's)
static thread_local const vector<u64>* tl_guide = nullptr; static thread_local long long* tl_guideKept = nullptr;
static bool note_node(const Problem& P, const State& S, int a, Stats& st, bool& isP1) {
    st.nodes[S.d]++;
    isP1 = a < 0 || P.radOf[a] < 0;
    if (!isP1) return false;
    st.p1[S.d]++;
    if (g_logp1) { lock_guard<mutex> g(solmu); if (p1spaces.size() < 200000) p1spaces[space_key(S.basis)] = S.basis; }
    if (S.pp) return false;
    string k = space_key(S.basis);
    st.cores[S.d]++; st.fp += fnv(k);
    if (g_coreFile) { lock_guard<mutex> g(solmu); fprintf(g_coreFile, "K %d %d %s\n", g_coreConf, S.d, k.c_str()); }
    { lock_guard<mutex> g(*tl_collmu); tl_coll->push_back(S); } // its D phase runs in run_search's second stage
    return true;
}
static void sweep_core(const State& S, Stats& st) {
    SwStats w; if (!sweep_from(*tl_sw, *tl_members, S.basis, w)) { fprintf(stderr, "sweep: core not legal by the sweep's own check\n"); exit(1); }
    for (int i = 0; i < 9; i++) st.swFound[i] += w.found[i];
    st.swChecks += w.checks;
    for (auto& B : w.sols) record_sweep(B, st);
}
// For a core of dimension <= 5 (its D phase is a deep tree): per point a of D, the candidates through a (line forms
// and oracle forms) legal with the core itself. A descendant S of the core branching at a needs core + <g> legal (a
// subspace of S + <g>), so it tries only these; the spans visited are the same (seven-check: added for the one
// dimension-4 core, whose D tree was out of reach without it; tries in D then count the cache's building as well).
struct CoreCache { const State* core; once_flag fl[256]; vector<u64> L[256]; };
static thread_local CoreCache* tl_cache = nullptr;
static const vector<u64>& cached_candidates(const Problem& P, CoreCache& K, int a, Stats& st) {
    call_once(K.fl[a], [&]() {
        State T; vector<u64> out;
        for_each_candidate(P, a, K.core->cov, [&](u64 f) { st.tries[K.core->d]++; st.triesD[K.core->d]++; if (try_add(P, *K.core, f, T)) out.push_back(f); });
        K.L[a].swap(out);
    });
    return K.L[a];
}
static void dfs(const Problem& P, const State& S, Stats& st) {
    int a = first_uncovered(P, S); bool isP1;
    if (note_node(P, S, a, st, isP1)) return;
    if (isP1 && S.d >= g_dmax) return; // cost probes only (SIX_DMAX)
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, st); return; }
    if (S.d == 8) return;
    State T; bool inD = P.radOf[a] < 0;
    if (inD && tl_cache && S.d > tl_cache->core->d) {
        for (u64 f : cached_candidates(P, *tl_cache, a, st)) { st.tries[S.d]++; st.triesD[S.d]++; if (try_add(P, S, f, T)) { T.pp = isP1; dfs(P, T, st); } }
        return;
    }
    for_each_candidate(P, a, S.cov, [&](u64 f) { st.tries[S.d]++; st.triesD[S.d] += inD; if (try_add(P, S, f, T)) { T.pp = isP1; dfs(P, T, st); } });
}
static void expand(const Problem& P, const State& S, Stats& st, vector<State>& out) {
    int a = first_uncovered(P, S); bool isP1;
    if (note_node(P, S, a, st, isP1)) return;
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, st); return; }
    if (S.d == 8) return;
    State T; bool inD = P.radOf[a] < 0;
    for_each_candidate(P, a, S.cov, [&](u64 f) { st.tries[S.d]++; st.triesD[S.d] += inD; if (try_add(P, S, f, T)) { T.pp = isP1; out.push_back(T); } });
}
static Stats run_search(const Problem& P, const State& S0, int mintasks) {
    Stats tot; vector<State> tasks = { S0 };
    vector<Mask> members(P.R, P.R + P.nr); vector<Mask> orc; for (auto& m : P.oracleM) orc.push_back(m);
    SwCtx C = sw_ctx(members, orc);
    auto* sink = tl_sol; // the caller's solution sink (controls), for whichever thread runs a task
    vector<State> cores; mutex cmu;
    tl_sw = &C; tl_members = &members; tl_coll = &cores; tl_collmu = &cmu;
    while (tasks.size() && (int)tasks.size() < mintasks) {
        vector<State> nxt; for (auto& s : tasks) expand(P, s, tot, nxt);
        tasks.swap(nxt);
        if (tasks.size() && tasks[0].d >= 7) break;
    }
    int nt = (int)tasks.size();
    vector<Stats> per(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nt; i++) { auto* save = tl_sol; tl_sol = sink; tl_sw = &C; tl_members = &members; tl_coll = &cores; tl_collmu = &cmu;
        dfs(P, tasks[i], per[omp_get_thread_num()]); tl_sol = save; }
    // second stage: the cores' D phases
    struct Job { int c, li; }; vector<Job> jobs; vector<int> firstPt(cores.size());
    auto* guide = tl_guide;
    for (int c = 0; c < (int)cores.size(); c++) {
        if (guide) { bool in = true; for (u64 b : cores[c].basis) in &= binary_search(guide->begin(), guide->end(), b); if (!in) continue; if (tl_guideKept) ++*tl_guideKept; }
        int a = firstPt[c] = first_uncovered(P, cores[c]);
        if (g_dsweep || a < 0 || cores[c].d == 8) { jobs.push_back({ c, -2 }); continue; }
        for (int li : P.linesAt[a]) if (mdisj(P.lineM[li], cores[c].cov)) jobs.push_back({ c, li });
        jobs.push_back({ c, -1 }); // oracle candidates (controls)
    }
    if (g_phase1only) jobs.clear(); // cores only (cost estimates)
    vector<unique_ptr<CoreCache>> caches(cores.size());
    if (!g_dsweep && !g_nocache) for (size_t c = 0; c < cores.size(); c++) if (cores[c].d <= 5) { caches[c].reset(new CoreCache()); caches[c]->core = &cores[c]; }
    int nj = (int)jobs.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nj; i++) {
        auto* save = tl_sol; tl_sol = sink; tl_sw = &C; tl_members = &members; tl_coll = &cores; tl_collmu = &cmu;
        Stats& st = per[omp_get_thread_num()]; const State& S = cores[jobs[i].c]; int a = firstPt[jobs[i].c], li = jobs[i].li;
        tl_cache = caches[jobs[i].c].get();
        if (li == -2) { if (g_dsweep) sweep_core(S, st); else if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, st); } }
        else {
            State T;
            auto go = [&](u64 f) { st.tries[S.d]++; st.triesD[S.d]++; if (try_add(P, S, f, T)) { T.pp = true; dfs(P, T, st); } };
            if (li >= 0) for (u64 f : line_forms(P, li)) go(f);
            else for (size_t k = 0; k < P.oracle.size(); k++) if (mget(P.oracleM[k], a) && !msub(P.oracleM[k], P.D) && mdisj(P.oracleM[k], S.cov)) go(P.oracle[k]);
        }
        tl_sol = save; tl_cache = nullptr;
    }
    for (auto& p : per) add_stats(tot, p);
    return tot;
}
static void make_order(Problem& P, int how, unsigned seed) {
    P.order.clear();
    vector<int> off, on; for (int v = 1; v < 256; v++) (P.radOf[v] >= 0 ? off : on).push_back(v);
    if (how == 1) { reverse(off.begin(), off.end()); reverse(on.begin(), on.end()); }
    if (how == 2) { mt19937 g(seed); shuffle(off.begin(), off.end(), g); shuffle(on.begin(), on.end(), g); }
    if (how == 3) { for (int v : on) P.order.push_back(v); for (int v : off) P.order.push_back(v); return; } // controls only: points of D first
    for (int v : off) P.order.push_back(v); for (int v : on) P.order.push_back(v);
}
struct Conf { int idx; int type; bool nonext; int hole; long long aut; vector<u16> Ns; vector<Mask> mem; };
static vector<Conf> read_confs(const char* file) {
    vector<Conf> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        Conf c; char ty[16];
        if (strncmp(line, "CONF", 4)) continue;
        if (sscanf(line, "CONF %d %15s hole %d aut %lld", &c.idx, ty, &c.hole, &c.aut) != 4) continue;
        c.type = -1; for (int t = 0; t < TYPES; t++) if (!strcmp(ty, TNAME[t])) c.type = t;
        c.nonext = c.type != 0;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; c.Ns.push_back((u16)v); p = e; }
        if (c.Ns.size() != KN || c.type < 0) { printf("bad CONF line\n"); exit(1); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
// Cross-run inputs that don't use this classification (type marks the family, not the orbit):
// 1: seven-check's 3731 configurations x 7 (7 * i + j: member j removed): every 6 that extends;
// 2: fourteen-check's spread reps x C(14,8) (3003 * i + 8-subset dropped): every 6 inside a spread.
static vector<Conf> confs_of_family(const vector<vector<Mask>>& fam, int type) {
    vector<Conf> out;
    for (size_t i = 0; i < fam.size(); i++) { Conf c; c.idx = (int)i; c.type = type; c.nonext = type != 0; c.hole = -1; c.aut = 0; c.mem = fam[i]; out.push_back(c); }
    return out;
}
static vector<Mask> conf_members(const Conf& c) { return c.mem.empty() ? members_of(c.Ns) : c.mem; }
static void print_stats(const char* tag, const Stats& s) {
    printf("%s nodes by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.nodes[i]);
    printf("\n%s tries by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.tries[i]);
    printf("\n%s of which at points of D:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.triesD[i]);
    printf("\n%s nodes with every off-D point covered, by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.p1[i]);
    printf("\n%s cores by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.cores[i]);
    printf("  fingerprint %016llx", (unsigned long long)s.fp);
    printf("\n%s sweep: legal spans by dim", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.swFound[i]);
    printf(", %lld element checks, %lld solutions", s.swChecks, s.swSols);
    printf("\n%s bent-extension nodes %lld, full-span records %lld\n", tag, s.bentext, s.full); fflush(stdout);
}
static void write_stats(FILE* pf, int j, const Stats& s) {
    fprintf(pf, "R %d", j); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.nodes[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.tries[k]);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.p1[k]); fprintf(pf, " %lld %lld", s.bentext, s.full);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.triesD[k]);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.cores[k]); fprintf(pf, " %016llx", (unsigned long long)s.fp);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.swFound[k]); fprintf(pf, " %lld %lld\n", s.swChecks, s.swSols); fflush(pf);
}
static bool read_stats(const char* line, int& id, Stats& s) {
    const char* p = line; if (*p != 'R') return false; p++;
    char* e; id = (int)strtol(p, &e, 10); if (e == p) return false; p = e;
    long long v[47]; for (int i = 0; i < 47; i++) { v[i] = strtoll(p, &e, 10); if (e == p) return false; p = e; }
    for (int k = 0; k < 9; k++) { s.nodes[k] = v[k]; s.tries[k] = v[9 + k]; s.p1[k] = v[18 + k]; s.triesD[k] = v[29 + k]; s.cores[k] = v[38 + k]; } s.bentext = v[27]; s.full = v[28];
    s.fp = _strtoui64(p, &e, 16); if (e == p) return false; p = e;
    long long w[11]; for (int i = 0; i < 11; i++) { w[i] = strtoll(p, &e, 10); if (e == p) return false; p = e; }
    for (int k = 0; k < 9; k++) s.swFound[k] = w[k]; s.swChecks = w[9]; s.swSols = w[10];
    return true;
}
static int g_input = 0; // 0 CONF, 1 eight-check's x 8, 2 spread reps x 3432
static int mode_search(const char* conffile, const char* progfile, int how, unsigned seed, int first, int last, const char* corefile) {
    auto cf = g_input == 1 ? confs_of_family(family7(read_conf7(conffile)), 1) : g_input == 2 ? confs_of_family(familySp(read_spread_file(conffile)), 0) : read_confs(conffile);
    const char* inName[3] = { "configuration list", "seven-check's configurations x 7", "spread reps x 3003" };
    printf("search: %zu configurations (%s); point order %s; D phase %s\n", cf.size(), inName[g_input], how == 0 ? "ascending" : how == 1 ? "descending" : "random", g_dsweep ? "coset sweep" : "point branching");
    set<int> done; Stats tot, totT[TYPES];
    if (FILE* fp = fopen(progfile, "r")) {
        char line[2048];
        while (fgets(line, sizeof line, fp)) { int id; Stats s; if (read_stats(line, id, s) && !done.count(id)) { done.insert(id); add_stats(tot, s); } else if (!strncmp(line, "SOLUTION", 8)) printf("(previous run) %s", line); }
        fclose(fp);
        printf("resumed: %zu configurations already done\n", done.size());
    }
    FILE* pf = fopen(progfile, "a");
    if (corefile) { g_coreFile = fopen(corefile, "a"); if (!g_coreFile) { printf("cannot open %s\n", corefile); return 1; } }
    auto t0 = chrono::steady_clock::now(); map<size_t, int> lineCounts; size_t nsolPrinted = 0;
    for (auto& c : cf) {
        if (c.idx < first || c.idx > last || done.count(c.idx)) continue;
        auto t1 = chrono::steady_clock::now();
        Problem P; build_problem(P, conf_members(c)); make_order(P, how, seed + c.idx);
        lineCounts[P.lineM.size()]++;
        if (mcount(P.D) != 165 || P.nr != 6 || P.lineM.size() != 2860) { printf("config %d: |D| = %d, %d members, %zu lines\n", c.idx, mcount(P.D), P.nr, P.lineM.size()); return 1; }
        g_coreConf = c.idx;
        Stats st = run_search(P, empty_state(), 400);
        { lock_guard<mutex> g(solmu); if (solutions.size() > nsolPrinted) { for (auto& kv : solutions) fprintf(pf, "SOLUTION %d %s\n", c.idx, kv.first.c_str()); nsolPrinted = solutions.size(); } }
        if (g_coreFile) fflush(g_coreFile);
        write_stats(pf, c.idx, st);
        add_stats(tot, st); add_stats(totT[c.type], st);
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
        printf("config %d (%s): nodes", c.idx, g_input ? "family" : TNAME[c.type]); for (int k = 0; k <= 8; k++) printf(" %lld", st.nodes[k]);
        printf("  cores"); for (int k = 0; k <= 8; k++) printf(" %lld", st.cores[k]);
        long long tD = 0; for (int k = 0; k <= 8; k++) tD += st.triesD[k];
        printf("  triesD %lld sweep %lld full %lld+%lld bentext %lld  [%.1f s, total %.0f s]\n", tD, st.swChecks, st.full, st.swSols, st.bentext, dt, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    }
    for (auto& kv : p1spaces) { fprintf(pf, "P1 dim %zu basis", kv.second.size()); for (u64 f : kv.second) fprintf(pf, " %07x", code28(f)); fprintf(pf, "\n"); }
    fclose(pf); if (g_coreFile) fclose(g_coreFile);
    printf("lines in D:"); for (auto& kv : lineCounts) printf(" %zu:%d", kv.first, kv.second); printf("\n");
    for (int t = 0; t < TYPES; t++) { string tag = string(TNAME[t]) + " (this run)"; print_stats(tag.c_str(), totT[t]); }
    print_stats("search (incl. resumed)", tot);
    printf("distinct p1 spans this run: %zu\n", p1spaces.size());
    printf("solutions: %zu\n", solutions.size());
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); printf("  %s APN=%d dist:", kv.first.c_str(), V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("\n"); }
    return 0;
}

// ================================================================================ the linear search (lin)
// Written for six-check. It shares with the code above only the input (members as point masks, configurations read by
// read_confs) and the u64 layout of a form (byte i = row i: symmetric, zero diagonal). Its own primitives:
//  - B a = XOR of the rows i with a_i = 1 (the vector B(a, .));
//  - a span S with basis b_1..b_e covers a point a once if the vectors b_i a have rank e - 1 (the kernel of S -> F^8,
//    B -> B a, is {B in S : a in rad B}); S is legal iff the rank is >= e - 1 at every point (no point covered twice);
//  - S + <g> (g not in S) is legal iff S is and g a is not in S a for every point a that S covers: at an uncovered point
//    the kernel has dimension <= 1 anyway;
//  - an 8-dimensional legal span is an APN space (its kernel at every a has dimension >= 8 - 7, so exactly 1).
// Member j's 28 forms (radical exactly R_j) come from 4 functionals u_1..u_4 with common kernel R_j:
// B = sum c_kl (u_k x u_l + u_l x u_k), nondegenerate on V / R_j iff Pf = c12 c34 + c13 c24 + c14 c23 = 1.
// Off-D phase: branch at the smallest uncovered point off D over member j's 28 forms (plus, in controls, the planted
// space's oracle forms through it); a child is kept when no point is covered twice and every newly covered point off D
// is covered by a member form of its member or an oracle form. These are the spans the point branching accepts (its
// legal() admits off D only those forms, and inside D every radical; the one difference: a rank-4 element whose radical
// is a solid inside D, a seventh solid, which the point branching rejects; here it is kept and counted, x4).
// D phase: branch at the smallest uncovered point a0 over every nonzero g in K_a0 = {B : B a0 = 0} (21-dimensional: the
// pair forms of the 7 functionals vanishing at a0), by Gray code, with the legality of S + <g> as byte lanes: one byte
// per covered point a, the reduction of g a modulo a fully reduced echelon basis of S a (linear in g; zero iff g a is in
// S a). In a solution T containing S exactly one element g of T has a0 in its radical and S + <g> is legal, so every T
// is reached; K_a0 meets S in 0, so the extensions are distinct. Branching at the same points with the same candidates
// as `search asc`, this visits the same spans: nodes by dimension must agree with asc's (less the x4 subtrees).
// If every point is covered below dimension 8 (possible only at dimension 6 or 7, by counting), every coset of S is tried.
// Cores of dimension <= 5: per point a, the g in K_a with core + <g> legal are listed once and reused below (as
// seven-check's caches; a descendant's candidates at a are among them).
struct LE { u8 v[8], c[8], pv[8], n, kc; };                 // echelon of b_i a: vectors, combinations, pivots; kc: kernel combination
static inline u8 L_mul(u64 B, u8 a) { u8 r = 0; while (a) { r ^= (u8)(B >> (8 * _tzcnt_u32(a))); a &= (u8)(a - 1); } return r; }
static inline u8 le_reduce(const LE& E, u8 y) { for (int k = 0; k < E.n; k++) if (y >> E.pv[k] & 1) y ^= E.v[k]; return y; }
// insert y = b_i a; false (and kc set) if y is in the span so far
static inline bool le_insert(LE& E, u8 y, int i) {
    u8 cm = (u8)(1 << i);
    for (int k = 0; k < E.n; k++) if (y >> E.pv[k] & 1) { y ^= E.v[k]; cm ^= E.c[k]; }
    if (!y) { E.kc = cm; return false; }
    unsigned long p; _BitScanReverse(&p, y);
    for (int k = 0; k < E.n; k++) if (E.v[k] >> p & 1) { E.v[k] ^= y; E.c[k] ^= cm; }
    E.v[E.n] = y; E.c[E.n] = cm; E.pv[E.n] = (u8)p; E.n++;
    return true;
}
static int L_rank(u64 B) { // own elimination
    u8 r[8]; for (int i = 0; i < 8; i++) r[i] = (u8)(B >> (8 * i));
    int rk = 0;
    for (int i = 0; i < 8; i++) { if (!r[i]) continue; unsigned long p; _BitScanReverse(&p, r[i]); u8 x = r[i]; rk++;
        for (int j = i + 1; j < 8; j++) if (r[j] >> p & 1) r[j] ^= x; }
    return rk;
}
static u64 L_pair(u8 u, u8 w) { // u x w + w x u: row i = u_i w + w_i u
    u64 f = 0; for (int i = 0; i < 8; i++) { u8 r = 0; if (u >> i & 1) r ^= w; if (w >> i & 1) r ^= u; f |= (u64)r << (8 * i); } return f;
}
// independent functionals vanishing on the points of m (as bytes: x -> parity(u & x)), `need` of them
static int L_annih(const Mask& m, u8* out, int need) {
    int n = 0; bool sp[256] = { false }; sp[0] = true;
    for (int u = 1; u < 256 && n < need; u++) {
        bool ok = true; for (int v = 1; v < 256 && ok; v++) if (mget(m, v) && (__popcnt((u32)(u & v)) & 1)) ok = false;
        if (!ok || sp[u]) continue;
        out[n++] = (u8)u; bool nw[256]; for (int v = 0; v < 256; v++) nw[v] = sp[v]; for (int v = 0; v < 256; v++) if (sp[v]) nw[v ^ u] = true; memcpy(sp, nw, 256);
    }
    return n;
}
static vector<u64> L_member_forms(const Mask& R) {
    u8 u[4]; if (L_annih(R, u, 4) != 4) { fprintf(stderr, "lin: member annihilator\n"); exit(1); }
    u64 P[6] = { L_pair(u[0], u[1]), L_pair(u[0], u[2]), L_pair(u[0], u[3]), L_pair(u[1], u[2]), L_pair(u[1], u[3]), L_pair(u[2], u[3]) }; // c12 c13 c14 c23 c24 c34
    vector<u64> out;
    for (int c = 1; c < 64; c++) {
        int pf = ((c >> 0 & 1) & (c >> 5 & 1)) ^ ((c >> 1 & 1) & (c >> 4 & 1)) ^ ((c >> 2 & 1) & (c >> 3 & 1));
        if (!pf) continue;
        u64 f = 0; for (int k = 0; k < 6; k++) if (c >> k & 1) f ^= P[k]; out.push_back(f);
    }
    return out;
}
static u64 KA[256][21]; // basis of K_a
static void L_init_ka() {
    for (int a = 1; a < 256; a++) {
        Mask m; mclear(m); mset(m, a); u8 u[7]; if (L_annih(m, u, 7) != 7) { fprintf(stderr, "lin: K_a\n"); exit(1); }
        int t = 0; for (int k = 0; k < 7; k++) for (int l = k + 1; l < 7; l++) KA[a][t++] = L_pair(u[k], u[l]);
    }
}
static bool g_linDesc = false; // lin desc: branch at the largest uncovered point (off D and in D) instead of the smallest
struct LSol { vector<u64> basis; };
struct LCtx {
    vector<LSol>* sink = nullptr; // solutions (controls); else the global list
    int jOf[256];                 // member of a point, -1 in D
    vector<u64> F[6];             // member forms
    vector<u64> oracle;           // controls: the planted space's forms allowed off D (rank 6, radical meeting a member)
    unordered_map<u64, int> allowed; // form -> member (0..5) or 6 (oracle)
    vector<int> offOrder;         // points off D, ascending
};
struct LNode { int e = 0; u64 b[8]; LE E[256]; bool x4 = false; };
struct LCore { int e; u64 b[8]; bool x4; };
struct LStats { long long nodes[10] = { 0 }, x4n[10] = { 0 }, steps[10] = { 0 }, cores[10] = { 0 }, x4c[10] = { 0 }, sols = 0, dead = 0, dsteps = 0; u64 fp = 0; };
static void L_add_stats(LStats& a, const LStats& b) {
    for (int i = 0; i < 10; i++) { a.nodes[i] += b.nodes[i]; a.x4n[i] += b.x4n[i]; a.steps[i] += b.steps[i]; a.cores[i] += b.cores[i]; a.x4c[i] += b.x4c[i]; }
    a.sols += b.sols; a.dead += b.dead; a.dsteps += b.dsteps; a.fp += b.fp;
}
static LCtx L_ctx(const vector<Mask>& M, const vector<u64>& oracle) {
    LCtx X; if (M.size() != 6) { fprintf(stderr, "lin: 6 members expected\n"); exit(1); }
    for (int v = 0; v < 256; v++) X.jOf[v] = -1;
    for (int j = 0; j < 6; j++) { for (int v = 1; v < 256; v++) if (mget(M[j], v)) { if (X.jOf[v] >= 0) { fprintf(stderr, "lin: members meet\n"); exit(1); } X.jOf[v] = j; }
        X.F[j] = L_member_forms(M[j]); if (X.F[j].size() != 28) { fprintf(stderr, "lin: member forms\n"); exit(1); }
        for (u64 f : X.F[j]) X.allowed[f] = j; }
    X.oracle = oracle; for (u64 f : oracle) X.allowed[f] = 6;
    for (int v = 1; v < 256; v++) if (X.jOf[v] >= 0) X.offOrder.push_back(v);
    if (g_linDesc) reverse(X.offOrder.begin(), X.offOrder.end());
    return X;
}
static inline bool L_covered(const LNode& S, int a) { return S.E[a].n + 1 == S.e; }
// the element of S whose radical holds a covered point a
static inline u64 L_coverer(const LNode& S, int a) { u64 k = 0; u8 c = S.E[a].kc; while (c) { k ^= S.b[_tzcnt_u32(c)]; c &= (u8)(c - 1); } return k; }
static bool L_legal_add(const LNode& S, u64 f) {
    for (int a = 1; a < 256; a++) if (L_covered(S, a) && !le_reduce(S.E[a], L_mul(f, (u8)a))) return false;
    return true;
}
// child = S + <f> (legal assumed); false if a newly covered point off D is covered by a form that isn't allowed there
static bool L_child(const LCtx& X, const LNode& S, u64 f, LNode& C) {
    C.e = S.e + 1; memcpy(C.b, S.b, sizeof C.b); C.b[S.e] = f; C.x4 = S.x4;
    for (int a = 1; a < 256; a++) {
        C.E[a] = S.E[a];
        bool was = L_covered(S, a);
        bool ind = le_insert(C.E[a], L_mul(f, (u8)a), S.e);
        if (was && !ind) { fprintf(stderr, "lin: child of an illegal extension\n"); exit(1); }
        if (!was && !ind && X.jOf[a] >= 0) { // newly covered off D: by a form allowed for this point?
            auto it = X.allowed.find(L_coverer(C, a));
            if (it == X.allowed.end() || (it->second != 6 && it->second != X.jOf[a])) return false;
        }
    }
    // a seventh solid: a rank-4 element of the new coset whose radical isn't a member
    int n = 1 << S.e;
    for (int s = 0; s < n && !C.x4; s++) { u64 g = f; for (int i = 0; i < S.e; i++) if (s >> i & 1) g ^= S.b[i];
        if (L_rank(g) == 4 && !X.allowed.count(g)) C.x4 = true; }
    return true;
}
static u64 L_key(const LNode& S) { // FNV-1a over the sorted nonzero elements
    vector<u64> el; for (int s = 1; s < (1 << S.e); s++) { u64 g = 0; for (int i = 0; i < S.e; i++) if (s >> i & 1) g ^= S.b[i]; el.push_back(g); }
    sort(el.begin(), el.end()); u64 h = 1469598103934665603ULL;
    for (u64 g : el) for (int k = 0; k < 8; k++) { h ^= (u8)(g >> (8 * k)); h *= 1099511628211ULL; }
    return h;
}
static void L_offD(const LCtx& X, const LNode& S, LStats& st, vector<LCore>& cores) {
    st.nodes[S.e]++; if (S.x4) st.x4n[S.e]++;
    int a = -1; for (int v : X.offOrder) if (!L_covered(S, v)) { a = v; break; }
    if (a < 0) { LCore c; c.e = S.e; memcpy(c.b, S.b, sizeof c.b); c.x4 = S.x4; cores.push_back(c); st.cores[S.e]++; if (S.x4) st.x4c[S.e]++; st.fp += L_key(S); return; }
    LNode C; int j = X.jOf[a];
    auto go = [&](u64 f) { st.steps[S.e]++; if (L_legal_add(S, f) && L_child(X, S, f, C)) L_offD(X, C, st, cores); };
    for (u64 f : X.F[j]) go(f);
    for (u64 f : X.oracle) if (!L_mul(f, (u8)a) && L_legal_add(S, f)) { st.steps[S.e]++; if (L_child(X, S, f, C)) L_offD(X, C, st, cores); }
}
// D phase
static const int LCH = 8; // lanes in chunks of 32 bytes: up to 256 covered points
struct alignas(32) LLanes { int nch; __m256i v[22][LCH]; u64 g[22]; };
static int L_lanes(const LNode& S, const u64* gen, int ngen, LLanes& L) {
    u8 pts[256]; int nc = 0; for (int a = 1; a < 256; a++) if (L_covered(S, a)) pts[nc++] = (u8)a;
    L.nch = (nc + 31) / 32;
    alignas(32) u8 buf[LCH * 32];
    for (int i = 0; i < ngen; i++) {
        memset(buf, 0, sizeof buf);
        for (int t = 0; t < nc; t++) buf[t] = le_reduce(S.E[pts[t]], L_mul(gen[i], pts[t]));
        for (int c = 0; c < L.nch; c++) L.v[i][c] = _mm256_load_si256((const __m256i*)(buf + 32 * c));
        L.g[i] = gen[i];
    }
    return nc;
}
static inline __m256i L_pad(int nc, int c) { // 0xff on the bytes past the covered points (never zero)
    alignas(32) u8 buf[32]; for (int t = 0; t < 32; t++) buf[t] = (32 * c + t >= nc) ? 0xff : 0; return _mm256_load_si256((const __m256i*)buf);
}
// chunks 1.. of the combination gc (chunk 0 already known to be zero-free)
static inline bool L_rest_ok(const LLanes& L, const __m256i* pad, u64 gc) {
    const __m256i z = _mm256_setzero_si256();
    for (int c = 1; c < L.nch; c++) {
        __m256i x = pad[c]; for (u64 m = gc; m; m &= m - 1) x = _mm256_xor_si256(x, L.v[_tzcnt_u64(m)][c]);
        if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, z))) return false;
    }
    return true;
}
// Gray scan over the combinations t in [t0, t1) (t0 >= 1) of the generators: chunk 0 (the first 32 covered points) kept
// by Gray code, the other chunks computed only when chunk 0 has no zero lane; cb(gray code, form) for each legal one
template <class CB> static void L_scan(const LLanes& L, const __m256i* pad, const u64* gen, u64 t0, u64 t1, CB cb) {
    const __m256i z = _mm256_setzero_si256();
    __m256i c0 = pad[0]; u64 g = 0;
    { u64 gc = (t0 - 1) ^ ((t0 - 1) >> 1); for (u64 m = gc; m; m &= m - 1) { int i = (int)_tzcnt_u64(m); g ^= gen[i]; c0 = _mm256_xor_si256(c0, L.v[i][0]); } }
    for (u64 t = t0; t < t1; t++) {
        int i = (int)_tzcnt_u64(t); g ^= gen[i]; c0 = _mm256_xor_si256(c0, L.v[i][0]);
        if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(c0, z))) continue;
        u64 gc = t ^ (t >> 1);
        if (L.nch > 1 && !L_rest_ok(L, pad, gc)) continue;
        cb(gc, g);
    }
}
struct LCache { const LNode* core; once_flag fl[256]; vector<u32> L[256]; };
// rebuild a node from a basis (legal assumed; checked)
static void L_build(const LCore& c, LNode& S) {
    S.e = 0; S.x4 = c.x4; for (int a = 0; a < 256; a++) S.E[a].n = 0;
    for (int i = 0; i < c.e; i++) { for (int a = 1; a < 256; a++) { bool cov = S.E[a].n + 1 == S.e; bool ind = le_insert(S.E[a], L_mul(c.b[i], (u8)a), i);
        if (cov && !ind) { fprintf(stderr, "lin: rebuilt core not legal\n"); exit(1); } } S.b[i] = c.b[i]; S.e++; }
}
static mutex lin_mu; static vector<LSol> lin_sols;
static void L_D(const LCtx& X, const LNode& S, LStats& st, LCache* K, u64 lo, u64 hi);
static void L_found(const LCtx& X, const LNode& S, u64 g, LStats& st, LCache* K) {
    LNode C; if (!L_child(X, S, g, C)) { st.dead++; return; } // (inside D every newly covered point is in D: never false)
    if (C.e == 8) { st.nodes[8]++; if (C.x4) st.x4n[8]++; st.sols++; lock_guard<mutex> lk(lin_mu); (X.sink ? X.sink : &lin_sols)->push_back({ vector<u64>(C.b, C.b + 8) }); return; }
    L_D(X, C, st, K, 0, 0);
}
// the D phase from S, counting S itself; [lo, hi) restricts the first level (Gray indices, or cache list positions) for
// the parallel split of a core's tree; hi = 0: everything
static void L_D(const LCtx& X, const LNode& S, LStats& st, LCache* K, u64 lo, u64 hi) {
    st.nodes[S.e]++; if (S.x4) st.x4n[S.e]++;
    if (S.e == 8) { // a core that is already a whole space (controls: T reached in the off-D phase); found by the lfull
        // controls, which failed before this line: the D phase never recorded it
        st.sols++; lock_guard<mutex> lk(lin_mu); (X.sink ? X.sink : &lin_sols)->push_back({ vector<u64>(S.b, S.b + 8) }); return;
    }
    int a0 = -1;
    if (g_linDesc) { for (int a = 255; a >= 1; a--) if (!L_covered(S, a)) { a0 = a; break; } }
    else for (int a = 1; a < 256; a++) if (!L_covered(S, a)) { a0 = a; break; }
    LLanes* L = (LLanes*)_aligned_malloc(sizeof(LLanes), 32);
    if (a0 < 0) { // every point covered (dimension 6 or 7): every coset of S, from pair forms independent of S
        if (S.e < 6) { fprintf(stderr, "lin: all covered at dimension %d\n", S.e); exit(1); }
        vector<u64> ech, gen;
        auto red = [&](u64 x) { for (u64 r : ech) { unsigned long p; _BitScanReverse64(&p, r); if (x >> p & 1) x ^= r; } return x; };
        auto ins = [&](u64 y) { y = red(y); if (!y) return false; unsigned long p; _BitScanReverse64(&p, y); for (u64& r : ech) if (r >> p & 1) r ^= y; ech.push_back(y); return true; };
        for (int i = 0; i < S.e; i++) if (!ins(S.b[i])) { fprintf(stderr, "lin: basis dependent\n"); exit(1); }
        for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { u64 f = L_pair((u8)(1 << i), (u8)(1 << j)); if (ins(f)) gen.push_back(f); }
        int ng = (int)gen.size(); if (ng != 28 - S.e || ng > 22) { fprintf(stderr, "lin: %d coset generators\n", ng); exit(1); }
        int nc = L_lanes(S, gen.data(), ng, *L);
        __m256i pad[LCH]; for (int c = 0; c < L->nch; c++) pad[c] = L_pad(nc, c);
        st.steps[S.e] += (long long)(1ULL << ng) - 1; st.dsteps += (long long)(1ULL << ng) - 1;
        L_scan(*L, pad, gen.data(), 1, 1ULL << ng, [&](u64, u64 g) { L_found(X, S, g, st, K); });
        _aligned_free(L); return;
    }
    const u64* gen = KA[a0];
    if (K && S.e > K->core->e) { // candidates from the core's list at a0
        call_once(K->fl[a0], [&]() {
            LLanes* Lc = (LLanes*)_aligned_malloc(sizeof(LLanes), 32);
            int nc = L_lanes(*K->core, gen, 21, *Lc); __m256i pad[LCH]; for (int c = 0; c < Lc->nch; c++) pad[c] = L_pad(nc, c);
            vector<u32> out; st.steps[K->core->e] += (1 << 21) - 1; st.dsteps += (1 << 21) - 1;
            L_scan(*Lc, pad, gen, 1, 1ULL << 21, [&](u64 gc, u64) { out.push_back((u32)gc); });
            K->L[a0].swap(out); _aligned_free(Lc);
        });
        int nc = L_lanes(S, gen, 21, *L);
        __m256i pad[LCH]; for (int c = 0; c < L->nch; c++) pad[c] = L_pad(nc, c);
        const __m256i z = _mm256_setzero_si256();
        st.dsteps += (long long)K->L[a0].size();
        for (u32 m : K->L[a0]) {
            st.steps[S.e]++;
            __m256i c0 = pad[0];
            for (u32 mm = m; mm; mm &= mm - 1) c0 = _mm256_xor_si256(c0, L->v[_tzcnt_u32(mm)][0]);
            if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(c0, z))) continue;
            if (L->nch > 1 && !L_rest_ok(*L, pad, m)) continue;
            u64 g = 0; for (u32 mm = m; mm; mm &= mm - 1) g ^= gen[_tzcnt_u32(mm)];
            L_found(X, S, g, st, K);
        }
        _aligned_free(L); return;
    }
    int nc = L_lanes(S, gen, 21, *L);
    u64 t0 = hi ? lo : 1, t1 = hi ? hi : (1ULL << 21);
    if (t0 == 0) t0 = 1;
    __m256i pad[LCH]; for (int c = 0; c < L->nch; c++) pad[c] = L_pad(nc, c);
    st.steps[S.e] += (long long)(t1 - t0); st.dsteps += (long long)(t1 - t0);
    L_scan(*L, pad, gen, t0, t1, [&](u64, u64 g) { L_found(X, S, g, st, K); });
    _aligned_free(L);
}
static bool L_verify(const vector<u64>& B) { // APN: at every a != 0 the 8 vectors b_i a span a 7-space
    for (int a = 1; a < 256; a++) { LE E; E.n = 0; int r = 0; for (int i = 0; i < 8; i++) r += le_insert(E, L_mul(B[i], (u8)a), i); if (r != 7) return false; }
    return true;
}
static void L_write(FILE* pf, const char* tag, int conf, int extra, const LStats& s) {
    fprintf(pf, "%s %d %d", tag, conf, extra);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.nodes[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.steps[k]);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.x4n[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.cores[k]);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.x4c[k]); fprintf(pf, " %lld %lld %016llx\n", s.sols, s.dead, (unsigned long long)s.fp); fflush(pf);
}
static bool L_read(const char* line, const char* tag, int& conf, int& extra, LStats& s) {
    size_t tl = strlen(tag); if (strncmp(line, tag, tl) || line[tl] != ' ') return false;
    const char* p = line + tl; char* e; conf = (int)strtol(p, &e, 10); if (e == p) return false; p = e; extra = (int)strtol(p, &e, 10); if (e == p) return false; p = e;
    long long v[47]; for (int i = 0; i < 47; i++) { v[i] = strtoll(p, &e, 10); if (e == p) return false; p = e; }
    for (int k = 0; k < 9; k++) { s.nodes[k] = v[k]; s.steps[k] = v[9 + k]; s.x4n[k] = v[18 + k]; s.cores[k] = v[27 + k]; s.x4c[k] = v[36 + k]; }
    s.sols = v[45]; s.dead = v[46]; s.fp = _strtoui64(p, &e, 16); return e != p;
}
// One configuration: the off-D phase (cores), then the cores' D phases as tasks (a core of dimension <= 5 split over its
// first-level list); finished tasks are appended to `ck` (T lines) and skipped on a restart.
static bool g_linP1 = false;
static LStats lin_conf(const vector<Mask>& M, const vector<u64>& oracle, int conf, FILE* ck, const set<int>& doneT, const map<int, LStats>& doneS, FILE* coreFile, bool quiet, vector<LSol>* sink = nullptr) {
    LCtx X = L_ctx(M, oracle); X.sink = sink;
    LNode S0; S0.e = 0; for (int a = 0; a < 256; a++) S0.E[a].n = 0;
    LStats off; vector<LCore> cores;
    { // the root by hand, its children's subtrees in parallel (cores concatenated in candidate order: deterministic)
        off.nodes[0]++;
        int a = X.offOrder[0], j = X.jOf[a];
        vector<u64> cand = X.F[j]; for (u64 f : X.oracle) if (!L_mul(f, (u8)a)) cand.push_back(f);
        int nc = (int)cand.size(); vector<vector<LCore>> part(nc); vector<LStats> ps(nc);
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < nc; i++) { ps[i].steps[0]++; LNode* C = new LNode(); if (L_legal_add(S0, cand[i]) && L_child(X, S0, cand[i], *C)) L_offD(X, *C, ps[i], part[i]); delete C; }
        for (int i = 0; i < nc; i++) { L_add_stats(off, ps[i]); cores.insert(cores.end(), part[i].begin(), part[i].end()); }
    }
    // remove the cores' own counts from the off-D tally: the D phase counts each core as its first node
    for (auto& C : cores) { off.nodes[C.e]--; if (C.x4) off.x4n[C.e]--; }
    if (coreFile) { lock_guard<mutex> lk(lin_mu); for (auto& C : cores) { fprintf(coreFile, "C %d %d", conf, C.e); for (int i = 0; i < C.e; i++) fprintf(coreFile, " %016llx", (unsigned long long)C.b[i]); fprintf(coreFile, "\n"); } fflush(coreFile); }
    LStats tot = off;
    if (g_linP1) { for (auto& C : cores) { tot.nodes[C.e]++; if (C.x4) tot.x4n[C.e]++; } return tot; }
    // tasks
    struct Task { int core; u64 lo, hi; };
    vector<unique_ptr<LNode>> owned;
    vector<Task> tasks; vector<unique_ptr<LCache>> caches(cores.size());
    for (int c = 0; c < (int)cores.size(); c++) {
        if (cores[c].e <= 5) { caches[c].reset(new LCache()); LNode* n = new LNode(); L_build(cores[c], *n); caches[c]->core = n; owned.emplace_back(n);
            const u64 parts = cores[c].e <= 4 ? 256 : 32; u64 step = (1ULL << 21) / parts;
            for (u64 k = 0; k < parts; k++) tasks.push_back({ c, k * step, (k + 1) * step }); }
        else tasks.push_back({ c, 0, 0 });
    }
    for (int t : doneT) if (doneS.count(t)) L_add_stats(tot, doneS.at(t));
    int nt = (int)tasks.size();
    if (getenv("SIX_MAXTASKS")) nt = min(nt, atoi(getenv("SIX_MAXTASKS"))); // timing probes only
    vector<LStats> per(nt);
#pragma omp parallel for schedule(dynamic, 1)
    for (int t = 0; t < nt; t++) {
        if (doneT.count(t)) continue;
        const Task& T = tasks[t]; LNode* Cp = new LNode(); L_build(cores[T.core], *Cp); const LNode& C = *Cp; LCache* K = caches[T.core].get();
        LStats st; auto tt0 = chrono::steady_clock::now();
        if (T.hi == 0) L_D(X, C, st, K, 0, 0);
        else {
            // part of a low core's first level: the Gray range [lo, hi) of K_a0; the core itself counted by part 0
            if (T.lo == 0) { st.nodes[C.e]++; if (C.x4) st.x4n[C.e]++; }
            LStats tmp; L_D(X, C, tmp, K, T.lo, T.hi); tmp.nodes[C.e]--; if (C.x4) tmp.x4n[C.e]--; L_add_stats(st, tmp);
        }
        delete Cp;
        per[t] = st;
        if (ck && chrono::duration<double>(chrono::steady_clock::now() - tt0).count() > 2.0) { lock_guard<mutex> lk(lin_mu); L_write(ck, "T", conf, t, st); }
    }
    for (int t = 0; t < nt; t++) if (!doneT.count(t)) L_add_stats(tot, per[t]);
    return tot;
}
static int mode_lin(const char* conffile, const char* progfile, int first, int last, const char* corefile) {
    auto cf = read_confs(conffile);
    printf("lin: %zu configurations; %s\n", cf.size(), g_linP1 ? "off-D phase only" : "off-D and D phases"); fflush(stdout);
    set<int> done; LStats tot;
    map<int, set<int>> ckDone; map<int, map<int, LStats>> ckStats;
    if (FILE* fp = fopen(progfile, "r")) { char line[4096]; while (fgets(line, sizeof line, fp)) { int c, x; LStats s;
        if (L_read(line, "L", c, x, s) && !done.count(c)) { done.insert(c); L_add_stats(tot, s); }
        else if (L_read(line, "T", c, x, s) && !ckDone[c].count(x)) { ckDone[c].insert(x); ckStats[c][x] = s; } } fclose(fp); }
    FILE* pf = fopen(progfile, "a"); FILE* cfp = corefile ? fopen(corefile, "a") : nullptr;
    auto t0 = chrono::steady_clock::now();
    for (auto& c : cf) {
        if (c.idx < first || c.idx > last || done.count(c.idx)) continue;
        auto t1 = chrono::steady_clock::now();
        // the core file is rewritten per configuration only when the configuration starts afresh
        bool fresh = !ckDone.count(c.idx);
        LStats st = lin_conf(conf_members(c), {}, c.idx, g_linP1 ? nullptr : pf, ckDone[c.idx], ckStats[c.idx], fresh ? cfp : nullptr, false);
        { lock_guard<mutex> lk(lin_mu); for (auto& s : lin_sols) { fprintf(pf, "SOLUTION %d APN=%d", c.idx, (int)L_verify(s.basis)); for (u64 b : s.basis) fprintf(pf, " %016llx", (unsigned long long)b); fprintf(pf, "\n"); } lin_sols.clear(); }
        L_write(pf, "L", c.idx, 0, st); L_add_stats(tot, st);
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
        printf("config %d (%s): nodes", c.idx, TNAME[c.type]); for (int k = 0; k <= 8; k++) printf(" %lld", st.nodes[k]);
        printf("  cores"); for (int k = 0; k <= 8; k++) printf(" %lld", st.cores[k]);
        long long sd = 0; for (int k = 0; k <= 8; k++) sd += st.steps[k];
        long long x4 = 0; for (int k = 0; k <= 8; k++) x4 += st.x4n[k];
        printf("  steps %lld x4 %lld sols %lld  [%.1f s, total %.0f s]\n", sd, x4, st.sols, dt, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    }
    fclose(pf); if (cfp) fclose(cfp);
    printf("lin total nodes by dim:"); for (int k = 0; k <= 8; k++) printf(" %lld", tot.nodes[k]);
    printf("\nlin total steps by dim:"); for (int k = 0; k <= 8; k++) printf(" %lld", tot.steps[k]);
    printf("\nlin total x4 nodes by dim:"); for (int k = 0; k <= 8; k++) printf(" %lld", tot.x4n[k]);
    printf("\nlin total cores by dim:"); for (int k = 0; k <= 8; k++) printf(" %lld", tot.cores[k]);
    printf("\nlin total solutions %lld\n", tot.sols);
    return 0;
}

// ================================================================================ positive controls
static bool read_lut(const string& line, u8 F[256]) {
    int n = 0; const char* p = line.c_str();
    while (*p && n < 256) { if (*p >= '0' && *p <= '9') { F[n++] = (u8)strtol(p, (char**)&p, 10); } else p++; }
    return n == 256;
}
static vector<u64> forms_of_lut(const u8 F[256]) {
    vector<u64> B(8, 0);
    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) for (int k = 0; k < 8; k++) if (j != k) {
        u8 d = F[(1 << j) ^ (1 << k)] ^ F[1 << j] ^ F[1 << k] ^ F[0];
        if (d >> i & 1) B[i] |= 1ULL << (8 * j + k);
    }
    return B;
}
struct Profile { int n4 = 0, n6 = 0, n8 = 0, bad = 0; vector<u64> f4, f6; vector<Mask> r4, r6; };
static Profile profile(const vector<u64>& B) {
    Profile P;
    for (int c = 1; c < (1 << B.size()); c++) {
        u64 f = 0; for (size_t i = 0; i < B.size(); i++) if (c >> i & 1) f ^= B[i];
        int rk = rank8(f); Mask m;
        if (rk == 8) P.n8++;
        else if (rk == 6) { P.n6++; radical(f, m); P.f6.push_back(f); P.r6.push_back(m); }
        else if (rk == 4) { P.n4++; radical(f, m); P.f4.push_back(f); P.r4.push_back(m); }
        else P.bad++;
    }
    return P;
}
static bool load_line(const char* file, int line, vector<u64>& B) {
    FILE* fp = fopen(file, "r"); if (!fp) return false;
    static char buf[1 << 16]; int ln = 0;
    while (fgets(buf, sizeof buf, fp)) { if (++ln == line) { u8 F[256]; fclose(fp); if (!read_lut(buf, F)) return false; B = forms_of_lut(F); return true; } }
    fclose(fp); return false;
}
static bool find_clique(vector<int>& cur, const BS& cand, int need) {
    if (need == 0) return true;
    for (int wi = 0; wi < NW; wi++) { u64 w = cand.w[wi]; while (w) { int v = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
        BS nc; above(nc, cand, adj[v], v);
        cur.push_back(v); if (find_clique(cur, nc, need - 1)) return true; cur.pop_back(); } }
    return false;
}
static bool find_clique_rand(vector<int>& cur, const BS& cand, int need, mt19937& rng) {
    if (need == 0) return true;
    vector<int> vs; for (int wi = 0; wi < NW; wi++) { u64 w = cand.w[wi]; while (w) { vs.push_back(wi * 64 + (int)_tzcnt_u64(w)); w &= w - 1; } }
    if ((int)vs.size() < need) return false;
    shuffle(vs.begin(), vs.end(), rng);
    int tries = 0;
    for (int v : vs) { BS nc; band(nc, cand, adj[v]); cur.push_back(v); if (find_clique_rand(cur, nc, need - 1, rng)) return true; cur.pop_back(); if (++tries > 64) break; }
    return false;
}
struct CtlCtx { vector<Conf> cf; map<vector<u16>, int> repOf; vector<int> nonextIdx; };
static CtlCtx ctl_ctx(const char* conffile) {
    CtlCtx X; X.cf = read_confs(conffile);
    for (size_t i = 0; i < X.cf.size(); i++) {
        if (X.cf[i].idx != (int)i) { printf("config file not indexed 0..n-1\n"); exit(1); }
        vector<u16> key; for (u16 N : X.cf[i].Ns) key.push_back((u16)vid[N]); sort(key.begin(), key.end()); X.repOf[key] = (int)i;
        if (X.cf[i].nonext) X.nonextIdx.push_back((int)i);
    }
    return X;
}
struct CtlResult { int status = 1; int rep = -1; int d0 = -1; int n4 = 0, inD = 0, nsol = 0, napn = 0; bool known = false; string why; Stats st; LStats lst; };
static bool g_ctlLin = false; // controls through the linear search (lin) instead of the point branching
// Normalise T: rank-4 radicals R1 -> X, R2 -> Y, R3 (if any) -> graph(I). Returns the profile in the new coordinates.
static bool normalise(vector<u64>& T, Profile& Q, string& why) {
    Mask Xm = XMASK(), Ym = YMASK();
    Lin Ai; basis_of_mask(Q.r4[0], 4, Ai.data()); basis_of_mask(Q.r4[1], 4, Ai.data() + 4);
    change_coords(T, lin_inv(Ai)); Q = profile(T);
    vector<int> rest; for (size_t i = 0; i < Q.r4.size(); i++) if (!meq(Q.r4[i], Xm) && !meq(Q.r4[i], Ym)) rest.push_back((int)i);
    if (!rest.empty()) {
        int M = matrix_of_graph(Q.r4[rest[0]]); if (M < 0) { why = "not a graph"; return false; }
        Lin B; for (int i = 0; i < 4; i++) { B[i] = mv4((u16)M, (u8)(1 << i)); B[4 + i] = (u8)(1 << (4 + i)); }
        change_coords(T, B); Q = profile(T);
    }
    return true;
}
// Run the search on configuration P (members in rep coordinates) for T (already in those coordinates).
// g_startNonD (the in-D test): pick g_nondDepth (1 or 2) rank-6 forms g_k of T whose radicals lie inside D, start from a
// random subspace W of T of codimension g_nondDepth with W + <g_1, ..> = T and no g_k in W + <the others>, with the
// points of the g_k's radicals first in the order. Each g_k's points are covered in T by g_k alone, so the search must
// find them among the main candidate generator's 13,888 forms on their lines (and pass the legality and disjointness
// tests) to reach T: one D step from a hyperplane, two from codimension 2. With the sweep, W must also hold every
// rank-4 form of T (the sweep starts at a core). Status 5 if T has too few rank-6 radicals inside D or no such W.
static bool g_startNonD = false, g_guided = false; static int g_nondDepth = 1;
static void ctl_search(const vector<Mask>& members, const vector<u64>& T, int how, CtlResult& R, bool verbose) {
    string key = space_key(T);
    Problem P; build_problem(P, members); make_order(P, how, 1);
    Profile Q = profile(T);
    for (auto& m : Q.r4) { bool in = false; for (int j = 0; j < P.nr; j++) if (meq(m, P.R[j])) in = true; if (!in) { R.status = 2; R.why = "rank-4 radical not a member"; return; } }
    for (size_t i = 0; i < Q.f6.size(); i++) { if (msub(Q.r6[i], P.D)) R.inD++; else { P.oracle.push_back(Q.f6[i]); P.oracleM.push_back(Q.r6[i]); } }
    State S0 = empty_state();
    if (g_startNonD) {
        int depth = g_nondDepth;
        vector<int> inD; for (size_t i = 0; i < Q.f6.size(); i++) if (msub(Q.r6[i], P.D)) inD.push_back((int)i);
        if ((int)inD.size() < depth) { R.status = 5; R.why = "too few rank-6 radicals inside D"; return; }
        mt19937 rng((unsigned)hash<string>()(key));
        auto coord = [&](u64 g) { for (int c = 1; c < 256; c++) { u64 f = 0; for (int j = 0; j < 8; j++) if (c >> j & 1) f ^= T[j]; if (f == g) return c; } return -1; };
        vector<int> c4; if (g_dsweep) for (int c = 1; c < 256; c++) { u64 f = 0; for (int j = 0; j < 8; j++) if (c >> j & 1) f ^= T[j]; if (rank8(f) == 4) c4.push_back(c); }
        vector<int> gi(depth), cg(depth); vector<int> fn(depth); bool ok = false;
        for (int attempt = 0; attempt < 64 && !ok; attempt++) {
            if ((int)inD.size() < depth) break;
            for (int k = 0; k < depth; k++) { bool rep; do { gi[k] = inD[rng() % inD.size()]; rep = false; for (int l = 0; l < k; l++) rep |= gi[l] == gi[k]; } while (rep); cg[k] = coord(Q.f6[gi[k]]); }
            for (int k = 0; k < depth; k++) if (cg[k] < 0) { R.status = 2; R.why = "g not in T"; return; }
            // functionals phi_k (parity of x . coords): phi_k(g_l) = [k == l], phi_k = 0 on the rank-4 forms (sweep)
            int kk = 0;
            for (; kk < depth; kk++) {
                int tries = 0; bool good = false;
                while (!good && tries++ < 4096) {
                    int x = (int)(rng() & 255); good = true;
                    for (int l = 0; l < depth; l++) if (par((u32)(x & cg[l])) != (kk == l)) good = false;
                    for (int c : c4) if (par((u32)(x & c))) good = false;
                    if (good) fn[kk] = x;
                }
                if (!good) break; // (six-check: try other g's; seven's code gave up here)
            }
            ok = kk == depth;
        }
        if (!ok) { R.status = 5; R.why = "no suitable W"; return; }
        // W = {x : phi_k(x) = 0 for all k}: the coordinate vectors c with par(c & fn[k]) = 0
        for (int c = 1; c < 256; c++) {
            bool in = true; for (int k = 0; k < depth; k++) if (par((u32)(c & fn[k]))) in = false;
            if (!in) continue;
            u64 w = 0; for (int j = 0; j < 8; j++) if (c >> j & 1) w ^= T[j];
            if (find(S0.span.begin(), S0.span.end(), w) != S0.span.end()) continue;
            State T2; if (!try_add(P, S0, w, T2)) { R.status = 2; R.why = "W not legal"; return; } S0 = T2;
        }
        R.d0 = S0.d;
        if (S0.d != 8 - depth) { R.status = 2; R.why = "W of the wrong dimension"; return; }
        Mask gm; mclear(gm); for (int k = 0; k < depth; k++) mor(gm, Q.r6[gi[k]]);
        if (!mdisj(gm, S0.cov)) { R.status = 2; R.why = "g's radical covered by W"; return; }
        vector<int> ord; for (int k = 0; k < depth; k++) for (int v = 1; v < 256; v++) if (mget(Q.r6[gi[k]], v)) ord.push_back(v);
        for (int v : P.order) if (!mget(gm, v)) ord.push_back(v);
        P.order = ord;
    }
    if (g_ctlLin) { // the linear search: the whole configuration (every core's D phase), or the D phase from S0 (nond)
        vector<Mask> mem(P.R, P.R + P.nr); vector<LSol> sols;
        if (g_startNonD) {
            LCtx X = L_ctx(mem, P.oracle); X.sink = &sols;
            LCore c; c.e = S0.d; for (int i = 0; i < S0.d; i++) c.b[i] = S0.basis[i]; c.x4 = false;
            LNode* W = new LNode(); L_build(c, *W); L_D(X, *W, R.lst, nullptr, 0, 0); delete W;
        } else R.lst = lin_conf(mem, P.oracle, -1, nullptr, {}, {}, nullptr, true, &sols);
        set<string> keys; for (auto& b : sols) { string k = space_key(b.basis); if (keys.insert(k).second) { Verdict V = verify_space(b.basis); R.napn += V.apn && L_verify(b.basis); } }
        R.known = keys.count(key) > 0; R.nsol = (int)keys.size(); R.status = R.known ? 0 : 3;
        return;
    }
    map<string, vector<u64>> sol; tl_sol = &sol;
    vector<u64> tspan = { 0 }; for (u64 b : T) { size_t n = tspan.size(); for (size_t i = 0; i < n; i++) tspan.push_back(tspan[i] ^ b); }
    sort(tspan.begin(), tspan.end()); long long kept = 0;
    if (g_guided) { tl_guide = &tspan; tl_guideKept = &kept; }
    R.st = run_search(P, S0, verbose ? 400 : 1);
    tl_sol = nullptr; tl_guide = nullptr; tl_guideKept = nullptr;
    if (g_guided && kept != 1) R.why = "guided: " + to_string(kept) + " cores inside T";
    R.known = sol.count(key) > 0; R.nsol = (int)sol.size();
    for (auto& kv : sol) { Verdict V = verify_space(kv.second); R.napn += V.apn;
        if (verbose) { printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); } }
    R.status = R.known ? 0 : 3;
}
// Extendable path: a spread through T's radicals, 11 other spread elements as holes, the 6-set carried onto its
// representative by its canonical frame (this tests the classification: it must land on a listed configuration).
static CtlResult ctl_ext(const CtlCtx& X, vector<u64> T, int choice, int how, bool verbose) {
    CtlResult R; Profile Q = profile(T); R.n4 = Q.n4;
    if (Q.n4 < 2 || Q.bad) { R.why = "fewer than 2 rank-4"; return R; }
    if (!verify_space(T).apn) { R.why = "forms not an APN space (line not quadratic APN)"; return R; }
    if (!normalise(T, Q, R.why)) return R;
    Mask Xm = XMASK(), Ym = YMASK();
    vector<int> qN; bool hasI = false;
    for (auto& m : Q.r4) { if (meq(m, Xm) || meq(m, Ym)) continue; int N = matrix_of_graph(m); if (N < 0 || (N != ID4 && vid[N] < 0)) { R.why = "radical not in V5824"; return R; } if (N == ID4) hasI = true; else qN.push_back(N); }
    if (Q.n4 >= 3 && !hasI) { R.why = "normalisation to graph(I) failed"; return R; }
    BS cand; for (int i = 0; i < NW; i++) cand.w[i] = ~0ULL;
    vector<int> cur; for (int N : qN) { cur.push_back(vid[N]); BS t; band(t, cand, adj[vid[N]]); cand = t; }
    bool got = false;
    { mt19937 rng((unsigned)choice * 2654435761u + 17); vector<int> c0 = cur;
      for (int t = 0; t < 50 && !got; t++) { cur = c0; got = find_clique_rand(cur, cand, 14 - (int)cur.size(), rng); }
      if (!got) cur = c0; }
    if (!got && !find_clique(cur, cand, 14 - (int)cur.size())) { R.why = "no spread through the rank-4 radicals"; return R; }
    vector<Mask> spreadAll = { Xm, Ym, graph_mask(ID4) }; for (int v : cur) spreadAll.push_back(graph_mask(VN[v]));
    vector<int> others; for (int i = 0; i < 17; i++) { bool isQ = false; for (auto& r : Q.r4) if (meq(r, spreadAll[i])) isQ = true; if (!isQ) others.push_back(i); }
    // choice -> a 10-subset of `others` (the holes; the other 7 spread elements are the members)
    mt19937 rng((unsigned)choice); vector<int> o = others; shuffle(o.begin(), o.end(), rng);
    if (o.size() < 11) { R.why = "too many rank-4"; return R; }
    set<int> holes(o.begin(), o.begin() + 11);
    vector<Mask> members; for (int i = 0; i < 17; i++) if (!holes.count(i)) members.push_back(spreadAll[i]);
    Canon K = canonP(members);
    auto it = X.repOf.find(K.key);
    if (it == X.repOf.end()) { R.status = 2; R.why = "configuration matches no representative"; return R; }
    R.rep = it->second;
    if (X.cf[R.rep].nonext) { R.status = 2; R.why = "extendable configuration matched a non-extendable rep"; return R; }
    change_coords(T, K.frame);
    vector<Mask> repM = members_of(X.cf[R.rep].Ns);
    { set<Mask> a, b; for (auto& m : members) a.insert(mask_image(K.frame, m)); for (auto& m : repM) b.insert(m); if (a != b) { R.status = 2; R.why = "frame does not carry the members onto the rep"; return R; } }
    ctl_search(repM, T, how, R, verbose);
    if (verbose) printf("ext control -> rep %d; rank-6 radicals %d inside D, found %d spaces (%d APN), known: %s\n", R.rep, R.inD, R.nsol, R.napn, R.known ? "YES" : "NO");
    return R;
}
// Non-extendable path: T's radicals embedded into the members of a non-extendable representative.
static thread_local int g_forceRep = -1; // embed into this representative (any type) instead of cycling the non-extendable ones
static CtlResult ctl_non(const CtlCtx& X, vector<u64> T, int choice, int how, bool verbose) {
    CtlResult R; Profile Q = profile(T); R.n4 = Q.n4;
    if (Q.n4 < 2 || Q.bad) { R.why = "fewer than 2 rank-4"; return R; }
    if (!verify_space(T).apn) { R.why = "forms not an APN space (line not quadratic APN)"; return R; }
    if (X.nonextIdx.empty()) { R.why = "no non-extendable configurations"; return R; }
    if (!normalise(T, Q, R.why)) return R;
    Mask Xm = XMASK(), Ym = YMASK();
    int N4 = -1; for (auto& m : Q.r4) { if (meq(m, Xm) || meq(m, Ym)) continue; int N = matrix_of_graph(m); if (N != ID4) N4 = N; }
    if (Q.n4 > 4) { R.why = "more than 4 rank-4"; return R; }
    mt19937 rng((unsigned)choice * 7919u + 3);
    R.rep = g_forceRep >= 0 ? g_forceRep : X.nonextIdx[choice % X.nonextIdx.size()];
    vector<Mask> M = members_of(X.cf[R.rep].Ns);
    int nm = (int)M.size();
    vector<array<u8, 4>> bas(nm); vector<vector<u8>> pts(nm); for (int i = 0; i < nm; i++) { basis_of_mask(M[i], 4, bas[i].data()); pts[i] = mpoints(M[i]); }
    vector<array<int, 3>> tr; for (int a = 0; a < nm; a++) for (int b = 0; b < nm; b++) for (int c = 0; c < nm; c++) if (a != b && a != c && b != c) tr.push_back({ a, b, c });
    shuffle(tr.begin(), tr.end(), rng);
    for (auto& t : tr) {
        Lin Cm; if (!frame_of(bas, pts, M, t[0], t[1], t[2], Cm)) { R.why = "frame"; return R; }
        Lin D;
        if (Q.n4 == 2) D = diag_ab(GL4[rng() % GL4.size()], GL4[rng() % GL4.size()]);
        else if (Q.n4 == 3) D = diag_hh(GL4[rng() % GL4.size()]);
        else {
            int tv = vid[N4]; if (tv < 0) { R.why = "4th radical not in V5824"; return R; }
            vector<int> ds; for (int d = 0; d < nm; d++) if (d != t[0] && d != t[1] && d != t[2]) { int v = vertex_in_frame(Cm, bas[d]); if (cls[v] == cls[tv]) ds.push_back(v); }
            if (ds.empty()) continue;
            int dv = ds[rng() % ds.size()];
            // h N4 h^-1 = N_d: h = h0(d)^-1 z h0(N4), z in the centraliser of the class rep
            int c = cls[tv]; u16 z = GL4[cent[c][rng() % cent[c].size()]];
            u16 h = mm4(inverse4(GL4[toRep[dv]]), mm4(z, GL4[toRep[tv]]));
            if (vid[mm4(mm4(h, (u16)N4), inverse4(h))] != dv) { R.why = "conjugator"; return R; }
            D = diag_hh(h);
        }
        Lin g = lin_comp(lin_inv(Cm), D);
        change_coords(T, g);
        ctl_search(M, T, how, R, verbose);
        if (verbose) printf("nonext control -> rep %d; rank-6 radicals %d inside D, found %d spaces (%d APN), known: %s\n", R.rep, R.inD, R.nsol, R.napn, R.known ? "YES" : "NO");
        return R;
    }
    R.status = 4; R.why = "no member 4-tuple of the 4th radical's class"; return R;
}
static int mode_ctl(const char* conffile, const char* file, int line, int choice, bool non) {
    init_conj();
    CtlCtx X = ctl_ctx(conffile);
    vector<u64> T; if (!load_line(file, line, T)) { printf("cannot load line %d\n", line); return 1; }
    g_logp1 = false;
    CtlResult R = non ? ctl_non(X, T, choice, 0, true) : ctl_ext(X, T, choice, 0, true);
    if (R.status == 1 || R.status == 2 || R.status == 4) printf("ctl: %s\n", R.why.c_str());
    print_stats("ctl", R.st);
    return R.status;
}
static int mode_ctlbatch(const char* conffile, const char* file, int minr4, long long maxn, int nch, const char* outfile, int how) {
    init_conj();
    CtlCtx X = ctl_ctx(conffile);
    g_logp1 = false;
    FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); return 1; }
    FILE* of = fopen(outfile, "a");
    vector<pair<long long, vector<u64>>> batch; long long ln = 0, taken = 0;
    Stats stT[2]; LStats lsT[2]; long long s0full = 0; map<int, long long> d0h[2];
    long long tot[2] = { 0 }, ok[2] = { 0 }, notfound[2] = { 0 }, nospread = 0, fail[2] = { 0 }, extra[2] = { 0 }, nonapn[2] = { 0 }, noembed = 0;
    map<int, int> repHits[2]; map<int, int> byN4[2];
    static char buf[1 << 16];
    auto flush = [&]() {
        int n = (int)batch.size(); vector<vector<CtlResult>> res(n, vector<CtlResult>(2 * nch));
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n; i++) for (int c = 0; c < 2 * nch; c++) {
            int ch = (int)(batch[i].first * 131 + c);
            res[i][c] = (c < nch) ? ctl_ext(X, batch[i].second, ch, how, false) : ctl_non(X, batch[i].second, ch, how, false);
        }
        for (int i = 0; i < n; i++) for (int c = 0; c < 2 * nch; c++) {
            int t = c < nch ? 0 : 1; auto& R = res[i][c];
            if (R.status == 4) { noembed++; continue; }
            if (R.status == 5) { s0full++; continue; }
            if (R.d0 >= 0) d0h[t][R.d0]++;
            if (R.why == "no spread through the rank-4 radicals") { nospread++; continue; }
            tot[t]++;
            add_stats(stT[t], R.st); L_add_stats(lsT[t], R.lst);
            if (R.status == 0) { ok[t]++; repHits[t][R.rep]++; byN4[t][R.n4]++; extra[t] += R.nsol - 1; nonapn[t] += R.nsol - R.napn; }
            else if (R.status == 3) { notfound[t]++; fprintf(of, "NOTFOUND %s line %lld choice %d rep %d\n", t ? "nonext" : "ext", batch[i].first, c, R.rep); }
            else { fail[t]++; fprintf(of, "FAIL %s line %lld choice %d: %s\n", t ? "nonext" : "ext", batch[i].first, c, R.why.c_str()); }
            if (R.status == 0 && R.nsol > 1) fprintf(of, "EXTRA %s line %lld choice %d rep %d: %d spaces, %d APN\n", t ? "nonext" : "ext", batch[i].first, c, R.rep, R.nsol, R.napn);
        }
        fflush(of); batch.clear();
        printf("  %lld lines read, %lld taken; ext %lld: re-found %lld, NOT %lld, fail %lld; nonext %lld: re-found %lld, NOT %lld, fail %lld; no spread %lld, no embedding %lld\n",
            ln, taken, tot[0], ok[0], notfound[0], fail[0], tot[1], ok[1], notfound[1], fail[1], nospread, noembed); fflush(stdout);
    };
    while (fgets(buf, sizeof buf, fp)) {
        ln++; u8 F[256]; if (!read_lut(buf, F)) continue;
        vector<u64> B = forms_of_lut(F);
        int n4 = 0; for (int c = 1; c < 256; c++) { u64 f = 0; for (int i = 0; i < 8; i++) if (c >> i & 1) f ^= B[i]; if (rank8(f) == 4) n4++; }
        if (n4 < minr4) continue;
        batch.push_back({ ln, B }); taken++;
        if ((int)batch.size() >= 256) flush();
        if (taken >= maxn) break;
    }
    if (!batch.empty()) flush();
    fclose(fp);
    for (int t = 0; t < 2; t++) {
        print_stats(t ? "ctl nonext" : "ctl ext", stT[t]);
        if (g_ctlLin) { printf("  lin: nodes by dim"); for (int k = 0; k <= 8; k++) printf(" %lld", lsT[t].nodes[k]); printf("; cores by dim"); for (int k = 0; k <= 8; k++) printf(" %lld", lsT[t].cores[k]);
            long long x = 0; for (int k = 0; k <= 8; k++) x += lsT[t].x4n[k]; printf("; x4 nodes %lld; candidates tried in D phases %lld; solutions %lld\n", x, lsT[t].dsteps, lsT[t].sols); }
        if (g_startNonD) { printf("  started from S0 of dimension:"); for (auto& kv : d0h[t]) printf(" %d:%lld", kv.first, kv.second); printf("\n"); }
        printf("ctlbatch %s %s: %lld controls; re-found %lld; NOT found %lld; setup failures %lld; other spaces found %lld (%lld not APN)\n", file, t ? "nonext" : "ext", tot[t], ok[t], notfound[t], fail[t], extra[t], nonapn[t]);
        printf("  re-found by rank-4 count:"); for (auto& kv : byN4[t]) printf(" %d:%d", kv.first, kv.second); printf("\n");
        printf("  reps hit: %zu\n", repHits[t].size());
        fprintf(of, "SUMMARY %s %s controls %lld refound %lld notfound %lld fail %lld extra %lld nonapn %lld reps %zu\n", file, t ? "nonext" : "ext", tot[t], ok[t], notfound[t], fail[t], extra[t], nonapn[t], repHits[t].size());
        fprintf(of, "REPS %s", t ? "nonext" : "ext"); for (auto& kv : repHits[t]) fprintf(of, " %d:%d", kv.first, kv.second); fprintf(of, "\n");
    }
    printf("skipped: no spread %lld, no member 4-tuple of the class (k = 4, nonext) %lld, S0 already the whole space %lld\n", nospread, noembed, s0full);
    fclose(of);
    return (notfound[0] || notfound[1] || fail[0] || fail[1]) ? 3 : 0;
}

// ================================================================================ planted controls (six-check)
// Controls whose planted space T is NOT reached in the off-D phase. For a known APN space T with k rank-4 forms
// (2 <= k <= 6), take as members T's k rank-4 radicals and 6 - k pairwise disjoint solids, each inside the set P_W of
// points covered by the rank-6 elements of a proper subspace W of T that holds T's rank-4 forms. In T every point off D
// is then covered by an element of W, so the forms covering the points off D (T's rank-4 forms and the oracle forms:
// T's rank-6 forms whose radical meets an extra member) span a subspace of W, of dimension d0 < 8: the off-D phase
// reaches T's core of dimension d0 and only the D phase can complete T (8 - d0 steps). The configuration's whole search
// (every core with its D phase) runs, unguided: T must be found and nothing else (other than further APN spaces using
// only the oracle forms off D, reported if any).
// plantfind: W over the hyperplanes (codim 1) or codimension-2 subspaces of T containing the rank-4 forms; solids in
// P_W by extending planes; the first family of 6 - k pairwise disjoint ones per W (tried in a fixed order), keeping
// the plantings with the smallest d0 per function.
static void solids_in(const Mask& P, vector<Mask>& out) {
    vector<u8> pts = mpoints(P); set<Mask> seen;
    Mask sh[256]; for (int x = 0; x < 256; x++) { mclear(sh[x]); for (u8 p : pts) mset(sh[x], p ^ x); } // P + x
    for (size_t i = 0; i < pts.size(); i++) for (size_t j = i + 1; j < pts.size(); j++) {
        u8 a = pts[i], b = pts[j]; if (!mget(P, a ^ b)) continue;
        for (size_t l = j + 1; l < pts.size(); l++) {
            u8 c = pts[l]; if (c == (a ^ b)) continue;
            if (!mget(P, a ^ c) || !mget(P, b ^ c) || !mget(P, a ^ b ^ c)) continue;
            u8 pl[7] = { a, b, (u8)(a ^ b), c, (u8)(a ^ c), (u8)(b ^ c), (u8)(a ^ b ^ c) };
            Mask cand = P; for (int t = 0; t < 7; t++) { Mask m; for (int w = 0; w < 4; w++) m.w[w] = cand.w[w] & sh[pl[t]].w[w]; cand = m; }
            for (int t = 0; t < 7; t++) { cand.w[pl[t] >> 6] &= ~(1ULL << (pl[t] & 63)); }
            for (u8 d : mpoints(cand)) {
                Mask s; mclear(s); for (int t = 0; t < 7; t++) { mset(s, pl[t]); mset(s, pl[t] ^ d); } mset(s, d);
                if (seen.insert(s).second) out.push_back(s);
            }
        }
    }
}
static bool disjoint_pick(const vector<Mask>& S, int need, size_t from, const Mask& used, vector<int>& ch) {
    if (need == 0) return true;
    for (size_t i = from; i < S.size(); i++) if (mdisj(S[i], used)) { Mask u = used; mor(u, S[i]); ch.push_back((int)i); if (disjoint_pick(S, need - 1, i + 1, u, ch)) return true; ch.pop_back(); }
    return false;
}
static int span_dim(const vector<u64>& fs) { vector<u32> c; for (u64 f : fs) c.push_back(code28(f)); return (int)rref28(c).size(); }
static int mode_plantfind(const char* lut, int minr4, long long maxfn, int codim, const char* outfile, long long firstLine) {
    FILE* fp = fopen(lut, "r"); if (!fp) { printf("cannot open %s\n", lut); return 1; }
    FILE* of = fopen(outfile, "a");
    static char buf[1 << 16]; long long ln = 0, taken = 0, planted = 0; map<int, int> byD0, byK;
    vector<pair<long long, vector<u64>>> batch;
    auto flush = [&]() {
        int n = (int)batch.size(); vector<string> res(n);
#pragma omp parallel for schedule(dynamic, 1)
        for (int bi = 0; bi < n; bi++) {
            const vector<u64>& T = batch[bi].second;
            u64 el[256]; int rk[256]; Mask rad[256]; el[0] = 0;
            for (int c = 1; c < 256; c++) { el[c] = el[c & (c - 1)] ^ T[_tzcnt_u32(c)]; rk[c] = rank8(el[c]); if (rk[c] < 8) radical(el[c], rad[c]); else mclear(rad[c]); }
            vector<int> c4; for (int c = 1; c < 256; c++) if (rk[c] == 4) c4.push_back(c);
            int k = (int)c4.size(); if (k < minr4 || k > 6) continue;
            // functional pairs (phi1, phi2) vanishing on c4: W = their common kernel
            vector<int> phis; for (int f = 1; f < 256; f++) { bool ok = true; for (int c : c4) if (par((u32)(f & c))) ok = false; if (ok) phis.push_back(f); }
            int bestD0 = 99; string best;
            auto tryW = [&](const Mask& Wc) { // Wc: the coordinate vectors c of W's elements
                Mask P; mclear(P);
                for (int c = 1; c < 256; c++) { if (!mget(Wc, c)) continue; if (rk[c] == 6) mor(P, rad[c]); }
                vector<Mask> sol; solids_in(P, sol);
                if ((int)sol.size() < 6 - k) return;
                vector<int> ch; Mask used; mclear(used);
                if (!disjoint_pick(sol, 6 - k, 0, used, ch)) return;
                vector<Mask> M; for (int c : c4) M.push_back(rad[c]); for (int i : ch) M.push_back(sol[i]);
                Mask off; mclear(off); for (auto& m : M) mor(off, m);
                vector<u64> cov; for (int c : c4) cov.push_back(el[c]);
                for (int c = 1; c < 256; c++) if (rk[c] == 6 && !mdisj(rad[c], off)) cov.push_back(el[c]);
                int d0 = span_dim(cov);
                if (d0 < bestD0) {
                    bestD0 = d0; char b2[512]; int o = sprintf(b2, "PLANT %lld k %d codim %d d0 %d :", batch[bi].first, k, codim, d0);
                    for (auto& m : M) { u8 bs[4]; basis_of_mask(m, 4, bs); o += sprintf(b2 + o, " %02x%02x%02x%02x", bs[0], bs[1], bs[2], bs[3]); }
                    best = b2;
                }
            };
            // W = the common kernel of `codim` independent functionals among phis (each W once)
            set<Mask> seenW;
            function<void(int, size_t, vector<int>&)> rec = [&](int need, size_t from, vector<int>& fs) {
                if (need == 0) {
                    Mask Wc; mclear(Wc); for (int c = 1; c < 256; c++) { bool in = true; for (int f : fs) if (par((u32)(c & f))) in = false; if (in) mset(Wc, c); }
                    if (mcount(Wc) != (1 << (8 - codim)) - 1 || !seenW.insert(Wc).second) return;
                    tryW(Wc); return;
                }
                for (size_t i = from; i < phis.size(); i++) { fs.push_back(phis[i]); rec(need - 1, i + 1, fs); fs.pop_back(); }
            };
            vector<int> fs; rec(codim, 0, fs);
            res[bi] = best;
        }
        for (int bi = 0; bi < n; bi++) if (!res[bi].empty()) { fprintf(of, "%s\n", res[bi].c_str()); planted++; int d0, k; sscanf(res[bi].c_str(), "PLANT %*lld k %d codim %*d d0 %d", &k, &d0); byD0[d0]++; byK[k]++; }
        fflush(of); batch.clear();
        printf("  %lld lines, %lld taken, %lld planted\n", ln, taken, planted); fflush(stdout);
    };
    while (fgets(buf, sizeof buf, fp)) {
        ln++; if (ln < firstLine) continue; u8 F[256]; if (!read_lut(buf, F)) continue;
        vector<u64> B = forms_of_lut(F);
        int n4 = 0; for (int c = 1; c < 256; c++) { u64 f = 0; for (int i = 0; i < 8; i++) if (c >> i & 1) f ^= B[i]; if (rank8(f) == 4) n4++; }
        if (n4 < minr4 || n4 > 6) continue;
        if (!verify_space(B).apn) continue;
        batch.push_back({ ln, B }); taken++;
        if ((int)batch.size() >= 64) flush();
        if (taken >= maxfn) break;
    }
    if (!batch.empty()) flush();
    fclose(fp); fclose(of);
    printf("plantfind %s: %lld functions with %d <= k <= 6, %lld planted; by d0:", lut, taken, minr4, planted); for (auto& kv : byD0) printf(" %d:%d", kv.first, kv.second);
    printf("; by k:"); for (auto& kv : byK) printf(" %d:%d", kv.first, kv.second); printf("\n");
    return 0;
}
// ctlplant <confs> <lut> <plants> <out> [lin] [first n]: each planting carried onto its representative (canonical frame:
// tests the classification), then the configuration's whole search with T's oracle forms
static int mode_ctlplant(const char* conffile, const char* lut, const char* plants, const char* outfile, int first, int count) {
    init_conj(); CtlCtx X = ctl_ctx(conffile); g_logp1 = false;
    FILE* pf = fopen(plants, "r"); if (!pf) { printf("cannot open %s\n", plants); return 1; }
    vector<string> lines; char buf[1024]; while (fgets(buf, sizeof buf, pf)) if (!strncmp(buf, "PLANT", 5)) lines.push_back(buf); fclose(pf);
    FILE* of = fopen(outfile, "a");
    long long n = 0, found = 0, notf = 0, fail = 0, extra = 0; Stats stA; LStats lsA; map<int, int> d0h, reps;
    for (int i = first; i < (int)lines.size() && i < first + count; i++) {
        long long ln; int k, cd, d0; if (sscanf(lines[i].c_str(), "PLANT %lld k %d codim %d d0 %d", &ln, &k, &cd, &d0) != 4) continue;
        vector<Mask> M; const char* p = strchr(lines[i].c_str(), ':') + 1;
        for (int j = 0; j < 6; j++) { unsigned v; int used; if (sscanf(p, " %8x%n", &v, &used) != 1) break; p += used;
            Mask m; mclear(m); u8 b[4] = { (u8)(v >> 24), (u8)(v >> 16), (u8)(v >> 8), (u8)v }; for (int s = 1; s < 16; s++) { u8 x = 0; for (int t = 0; t < 4; t++) if (s >> t & 1) x ^= b[t]; mset(m, x); } M.push_back(m); }
        vector<u64> T; if (M.size() != 6 || !load_line(lut, (int)ln, T)) { fprintf(of, "FAIL plant %d: input\n", i); fail++; continue; }
        CtlResult R; Canon K = canonP(M); auto it = X.repOf.find(K.key);
        if (it == X.repOf.end()) { fprintf(of, "FAIL plant %d: configuration matches no representative\n", i); fail++; continue; }
        R.rep = it->second; change_coords(T, K.frame);
        vector<Mask> repM = members_of(X.cf[R.rep].Ns);
        { set<Mask> a, b; for (auto& m : M) a.insert(mask_image(K.frame, m)); for (auto& m : repM) b.insert(m); if (a != b) { fprintf(of, "FAIL plant %d: frame\n", i); fail++; continue; } }
        auto t0 = chrono::steady_clock::now();
        ctl_search(repM, T, 0, R, false);
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        n++; add_stats(stA, R.st); L_add_stats(lsA, R.lst); d0h[d0]++; reps[R.rep]++;
        long long tD = 0; for (int d = 0; d < 9; d++) tD += R.st.triesD[d];
        if (R.status == 0) { found++; extra += R.nsol - 1; } else if (R.status == 3) notf++; else { fail++; }
        fprintf(of, "%s plant %d line %lld k %d d0 %d -> rep %d (%s): %d spaces (%d APN); tries in D %lld, lin D-phase candidates %lld, lin nodes by dim", R.status == 0 ? "FOUND" : R.status == 3 ? "NOTFOUND" : "FAIL", i, ln, k, d0, R.rep, TNAME[X.cf[R.rep].type], R.nsol, R.napn, tD, R.lst.dsteps);
        for (int d = 0; d < 9; d++) fprintf(of, " %lld", R.lst.nodes[d]); fprintf(of, " %s [%.0f s]\n", R.why.c_str(), dt); fflush(of);
        printf("plant %d: %s, %d spaces [%.0f s]\n", i, R.status == 0 ? "re-found" : "NOT re-found", R.nsol, dt); fflush(stdout);
    }
    print_stats("ctlplant", stA);
    if (g_ctlLin) { printf("ctlplant lin: nodes by dim"); for (int d = 0; d < 9; d++) printf(" %lld", lsA.nodes[d]); printf("; D-phase candidates %lld\n", lsA.dsteps); }
    printf("ctlplant: %lld controls, re-found %lld, NOT found %lld, failures %lld, other spaces %lld; by d0:", n, found, notf, fail, extra); for (auto& kv : d0h) printf(" %d:%d", kv.first, kv.second);
    printf("; representatives %zu\n", reps.size());
    fprintf(of, "SUMMARY %s controls %lld refound %lld notfound %lld fail %lld extra %lld\n", g_ctlLin ? "lin" : "pb", n, found, notf, fail, extra); fclose(of);
    return (notf || fail) ? 3 : 0;
}

// linselftest <confs> <n configurations>: lin's primitives against brute force (the point branching's rank8/radical
// and plain point counting). Per configuration: member forms (28, distinct, radical exactly the member), K_a (21
// independent forms with a in the radical, every a); then on the cores of the configuration (lin's off-D phase) and
// random g: L_legal_add against "no point in two radicals of S + <g>" counted directly, and the lane scan's verdict on
// the first 2^16 candidates of K_a0 against L_legal_add.
static int mode_linselftest(const char* conffile, int ncf) {
    auto cf = read_confs(conffile); int bad = 0; long long checks = 0, legalSeen = 0, accepted = 0; mt19937_64 rng(5);
    for (int a = 1; a < 256; a++) { vector<u32> c; for (int i = 0; i < 21; i++) { Mask m; c.push_back(code28(KA[a][i])); if (L_mul(KA[a][i], (u8)a)) bad++; }
        if (rref28(c).size() != 21) bad++; }
    printf("K_a: %d failures\n", bad); fflush(stdout);
    for (int ci = 0; ci < (int)cf.size() && ci < ncf; ci++) {
        int idx = (int)((long long)ci * (long long)cf.size() / ncf);
        vector<Mask> M = conf_members(cf[idx]); LCtx X = L_ctx(M, {});
        for (int j = 0; j < 6; j++) { set<u64> s(X.F[j].begin(), X.F[j].end()); if (s.size() != 28) bad++;
            for (u64 f : X.F[j]) { Mask m; if (rank8(f) != 4 || radical(f, m) != 15 || !meq(m, M[j])) bad++; } }
        LNode S0; S0.e = 0; for (int a = 0; a < 256; a++) S0.E[a].n = 0;
        LStats st; vector<LCore> cores; L_offD(X, S0, st, cores);
        shuffle(cores.begin(), cores.end(), rng); if (cores.size() > 200) cores.resize(200);
        for (auto& c : cores) {
            LNode* S = new LNode(); L_build(c, *S);
            vector<u64> span = { 0 }; for (int i = 0; i < S->e; i++) { size_t n = span.size(); for (size_t t = 0; t < n; t++) span.push_back(span[t] ^ S->b[i]); }
            auto brute = [&](u64 g) { // every point in at most one radical of S + <g>
                int cnt[256] = { 0 };
                for (u64 s : span) for (u64 h : { s, s ^ g }) { if (!h) continue; if (rank8(h) == 8) continue; Mask m; radical(h, m); for (int v = 1; v < 256; v++) if (mget(m, v) && ++cnt[v] > 1) return false; }
                return true;
            };
            int a0 = -1; for (int a = 1; a < 256; a++) if (!L_covered(*S, a)) { a0 = a; break; }
            for (int r = 0; r < 64; r++) { // random elements of K_a0 and random forms
                u64 g = 0; if (r & 1) { for (int i = 0; i < 21; i++) if (rng() & 1) g ^= KA[a0][i]; } else g = form28((u32)(rng() & ((1u << 28) - 1)));
                if (!g || find(span.begin(), span.end(), g) != span.end()) continue;
                bool l = L_legal_add(*S, g), b = brute(g); checks++; legalSeen += l; if (l != b) bad++;
            }
            // lane scan vs L_legal_add on the first 2^16 Gray candidates of K_a0
            LLanes* L = (LLanes*)_aligned_malloc(sizeof(LLanes), 32); int nc = L_lanes(*S, KA[a0], 21, *L);
            __m256i pad[LCH]; for (int k = 0; k < L->nch; k++) pad[k] = L_pad(nc, k);
            set<u64> scanned; L_scan(*L, pad, KA[a0], 1, 1ULL << 16, [&](u64, u64 g) { scanned.insert(g); });
            u64 g = 0; for (u64 t = 1; t < (1ULL << 16); t++) { g ^= KA[a0][_tzcnt_u64(t)]; bool l = L_legal_add(*S, g); checks++; legalSeen += l; if (l != (scanned.count(g) > 0)) bad++; }
            // the whole of K_a0: every accepted g legal by brute force
            long long acc = 0; L_scan(*L, pad, KA[a0], 1, 1ULL << 21, [&](u64, u64 g) { acc++; checks++; legalSeen++; if (!brute(g) || !L_legal_add(*S, g)) bad++; });
            accepted += acc;
            _aligned_free(L); delete S;
        }
        printf("config %d: %zu cores sampled, %lld accepted in full scans, failures so far %d\n", idx, cores.size(), accepted, bad); fflush(stdout);
    }
    printf("linselftest: %lld checks (%lld legal), %d failures\n", checks, legalSeen, bad);
    return bad ? 2 : 0;
}

// ================================================================================ selftest
static int mode_selftest(const char* conffile) {
    init_conj();
    auto cf = read_confs(conffile); int bad = 0;
    printf("selftest: ND4 %zu (28), ND6 %zu (13888), GL4 %zu (20160), V5824 %zu\n", ND4.size(), ND6.size(), GL4.size(), VN.size());
    bad += ND4.size() != 28 || ND6.size() != 13888 || GL4.size() != 20160 || VN.size() != 5824;
    long long checked = 0; mt19937 rng(12345); map<size_t, int> lines;
    for (auto& c : cf) {
        vector<Mask> M = conf_members(c);
        Problem P; build_problem(P, M);
        bool ok = mcount(P.D) == 165 && P.nr == 6;
        for (int a = 0; a < KN; a++) { ok &= inv4(c.Ns[a]) && inv4(c.Ns[a] ^ ID4); for (int b = a + 1; b < KN; b++) ok &= inv4(c.Ns[a] ^ c.Ns[b]); }
        lines[P.lineM.size()]++; ok &= P.lineM.size() == 2860;
        { int per[256] = { 0 }; for (auto& m : P.lineM) for (int v = 1; v < 256; v++) per[v] += mget(m, v); for (int v = 1; v < 256; v++) ok &= per[v] == (mget(P.D, v) ? 52 : 0); }
        for (size_t li = 0; li < P.lineM.size(); li++) ok &= msub(P.lineM[li], P.D) && mcount(P.lineM[li]) == 3;
        for (int j = 0; j < 6; j++) { set<u64> s(P.f4[j].begin(), P.f4[j].end()); ok &= s.size() == 28; for (u64 f : P.f4[j]) { Mask m; ok &= radical(f, m) == 15 && meq(m, P.R[j]); checked++; } }
        if (c.idx % 128 == 0) for (size_t li = 0; li < P.lineM.size(); li++) {
            const vector<u64>& lf = line_forms(P, li); set<u64> s(lf.begin(), lf.end()); ok &= s.size() == 13888;
            for (u64 f : lf) { Mask m; ok &= rank8(f) == 6 && radical(f, m) == 3 && meq(m, P.lineM[li]); checked++;
                u8 pts[255]; int np = sw_radical(f, pts); Mask m2; mclear(m2); for (int k = 0; k < np; k++) mset(m2, pts[k]); ok &= np == 3 && meq(m2, m); }
        }
        // canonical form: invariant under a random GL(8,2) image, and the stored aut
        { Lin g = random_gl8(rng); vector<Mask> Mg; for (auto& m : M) Mg.push_back(mask_image(g, m)); shuffle(Mg.begin(), Mg.end(), rng);
          Canon K = canonP(Mg); vector<u16> key; for (u16 N : c.Ns) key.push_back((u16)vid[N]); sort(key.begin(), key.end());
          ok &= K.key == key && K.aut == c.aut;
          set<Mask> a, b; for (auto& m : Mg) a.insert(mask_image(K.frame, m)); for (auto& m : M) b.insert(m); ok &= a == b; }
        if (!ok) { printf("config %d: FAIL\n", c.idx); bad++; }
    }
    printf("lines in D over the configurations:"); for (auto& kv : lines) printf(" %zu:%d", kv.first, kv.second); printf("\n");
    { Problem P; build_problem(P, conf_members(cf.back())); long long cnt = 0, kbad = 0;
#pragma omp parallel for reduction(+:cnt, kbad) schedule(static, 1 << 16)
      for (int c = 0; c < (1 << 28); c++) { u64 f = form28((u32)c); int rk = rank8(f); Mask m; int n = rk == 8 ? 0 : radical(f, m);
          if (rk == 4 && meq(m, P.R[5])) cnt++;
          if ((c & 15) == 0 && sw_nondeg(f) != (rk == 8)) kbad++;
          if ((c & 255) == 0) { Mask mw; int nw = radical_walk(f, mw); if (nw != (rk == 8 ? 0 : n) || (rk != 8 && !meq(mw, m))) kbad++; }
          if ((c & 255) == 0) { u8 pts[255]; int np = sw_radical(f, pts); if (np != n) kbad++; else if (n) { Mask m2; mclear(m2); for (int k = 0; k < np; k++) mset(m2, pts[k]); if (!meq(m2, m)) kbad++; } } }
      printf("brute force: %lld forms with radical exactly R_5 of the last config (expect 28); sweep Pfaffian test vs rank8 on every 16th form, sweep kernel and the Gray walk vs radical() on every 256th: %lld disagreements\n", cnt, kbad);
      bad += cnt != 28 || kbad; }
    printf("selftest: %zu configurations, %lld generated forms checked, %d failures\n", cf.size(), checked, bad);
    return bad;
}

// ================================================================================ count3: k-cliques of V5824 two ways
// (a) plain enumeration, every clique once (increasing vertex index), 4-cliques by popcount at the last level;
// (b) by edges: sum over edges {u, v} of |N(u) & N(v)| (triangles, each 3 times) and of the edges inside N(u) & N(v)
// (4-cliques, each 6 * 2 times). No group action, no canonical forms. The 3-clique count must equal Part A's weighted
// count (merge), the 4-clique count seven-check's (14122766336).
static int mode_count3() {
    long long cnt[5] = { 0 };
    auto t0 = chrono::steady_clock::now();
    BS all; for (int i = 0; i < NW; i++) all.w[i] = ~0ULL;
    int n = (int)VN.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int v1 = 0; v1 < n; v1++) {
        long long c[5] = { 0 }; c[1] = 1;
        BS a; above(a, all, adj[v1], v1);
        for (int v2 : blist(a)) { c[2]++; BS b; above(b, a, adj[v2], v2);
            for (int wi = 0; wi < NW; wi++) { u64 w = b.w[wi]; while (w) { int v3 = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
                c[3]++; BS d; above(d, b, adj[v3], v3); c[4] += bcount(d); } } }
#pragma omp critical
        for (int k = 1; k <= 4; k++) cnt[k] += c[k];
    }
    double t1 = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    printf("count3 (a): k-cliques of V5824 by plain enumeration:"); for (int k = 1; k <= 4; k++) printf(" %d:%lld", k, cnt[k]); printf("  [%.0f s]\n", t1); fflush(stdout);
    long long e6 = 0, tri3 = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:e6, tri3)
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (bget(adj[u], v)) {
        BS C; band(C, adj[u], adj[v]); tri3 += bcount(C);
        long long s = 0; for (int w : blist(C)) { BS D; band(D, C, adj[w]); s += bcount(D); } e6 += s; }
    printf("count3 (b): triangles %lld, 4-cliques %lld (by edges)  [%.0f s]\n", tri3 / 3, e6 / 12, chrono::duration<double>(chrono::steady_clock::now() - t0).count() - t1);
    bool ok = cnt[4] == 14122766336LL && e6 % 12 == 0 && e6 / 12 == cnt[4] && tri3 % 3 == 0 && tri3 / 3 == cnt[3];
    printf("count3: 4-cliques = seven-check's 14122766336 %s; (a) and (b) agree %s\n", cnt[4] == 14122766336LL ? "OK" : "MISMATCH", (e6 / 12 == cnt[4] && tri3 / 3 == cnt[3]) ? "OK" : "MISMATCH");
    return ok ? 0 : 2;
}

int main(int argc, char** argv) {
    if (argc < 2) { printf("modes: count3 | parta <class> <orbits out> <maxcliques out> [<slice> <nslices> [<test firsts>]] | merge <orbits out> <configs.txt (check7)> <part files...> | configs <orbits> <configs.txt (check7)> <configs.txt (check8)> <spreads.txt (check14)> <out> | search|search7|searchsp <input> <prog> [asc|desc|rand <seed>] [first last] [cores <file>] [dsweep] | selftest <confs> | ctl <confs> <lut> <line> <choice> [non] | ctlbatch <confs> <lut> <minr4> <maxn> <nch> <out> [nond|dfirst|dsweep|nondsweep] | ctlrot <confs> <lut> <first line> <pool> <rep>...\n"); return 1; }
    string mode = argv[1];
    ND4 = nondeg_codes(4); ND6 = nondeg_codes(6);
    init_V();
    printf("V5824: %zu vertices, %zu classes:", VN.size(), classSize.size()); for (size_t c = 0; c < classSize.size(); c++) printf(" %d (rep %04x)", classSize[c], VN[classRep[c]]); printf("\n");
    fflush(stdout);
    if (mode == "parta") { int sl = 0, nsl = 1; if (argc > 6) { sl = atoi(argv[5]); nsl = atoi(argv[6]); } if (argc > 7) g_testFirsts = atoi(argv[7]); return mode_parta(atoi(argv[2]), argv[3], argv[4], sl, nsl); }
    if (mode == "count3") return mode_count3();
    if (mode == "merge") { vector<const char*> parts; for (int i = 4; i < argc; i++) parts.push_back(argv[i]); return mode_merge(argv[2], argv[3], parts); }
    if (mode == "configs") return mode_configs(argv[2], argv[3], argv[4], argv[5], argv[6]);
    init_code28(); sw_init(); L_init_ka();
    if (mode == "search" || mode == "search7" || mode == "searchsp") {
        g_input = mode == "search7" ? 1 : mode == "searchsp" ? 2 : 0;
        int how = 0; unsigned seed = 1; int first = 0, last = INT_MAX, ai = 4; const char* corefile = nullptr;
        if (argc > ai && !strcmp(argv[ai], "desc")) { how = 1; ai++; }
        else if (argc > ai && !strcmp(argv[ai], "rand")) { how = 2; seed = (unsigned)atoi(argv[ai + 1]); ai += 2; }
        else if (argc > ai && !strcmp(argv[ai], "asc")) ai++;
        if (argc > ai + 1 && isdigit((unsigned char)argv[ai][0])) { first = atoi(argv[ai]); last = atoi(argv[ai + 1]); ai += 2; }
        while (argc > ai) { if (!strcmp(argv[ai], "cores") && argc > ai + 1) { corefile = argv[ai + 1]; ai += 2; } else if (!strcmp(argv[ai], "dsweep")) { g_dsweep = true; ai++; } else if (!strcmp(argv[ai], "phase1")) { g_phase1only = true; ai++; } else { printf("unknown option %s\n", argv[ai]); return 1; } }
        return mode_search(argv[2], argv[3], how, seed, first, last, corefile);
    }
    if (mode == "lin") { // lin <confs> <prog> [first last] [cores <file>] [phase1]
        int first = 0, last = INT_MAX, ai = 4; const char* corefile = nullptr;
        if (argc > ai + 1 && isdigit((unsigned char)argv[ai][0])) { first = atoi(argv[ai]); last = atoi(argv[ai + 1]); ai += 2; }
        while (argc > ai) { if (!strcmp(argv[ai], "cores") && argc > ai + 1) { corefile = argv[ai + 1]; ai += 2; } else if (!strcmp(argv[ai], "phase1")) { g_linP1 = true; ai++; } else if (!strcmp(argv[ai], "desc")) { g_linDesc = true; ai++; } else { printf("unknown option %s\n", argv[ai]); return 1; } }
        return mode_lin(argv[2], argv[3], first, last, corefile);
    }
    if (mode == "plantfind") return mode_plantfind(argv[2], atoi(argv[3]), atoll(argv[4]), atoi(argv[5]), argv[6], argc > 7 ? atoll(argv[7]) : 0);
    if (mode == "ctlplant") { g_ctlLin = argc > 8 && !strcmp(argv[8], "lin"); return mode_ctlplant(argv[2], argv[3], argv[4], argv[5], atoi(argv[6]), atoi(argv[7])); }
    if (mode == "linselftest") return mode_linselftest(argv[2], atoi(argv[3]));
    if (mode == "selftest") return mode_selftest(argv[2]);
    if (mode == "ctl") return mode_ctl(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argc > 6 && !strcmp(argv[6], "non"));
    if (mode == "ctlbatch") {
        string opt = argc > 8 ? argv[8] : ""; // tokens joined by '+': nond, dsweep, guided, dfirst
        auto has = [&](const char* t) { string o = "+" + opt + "+"; return o.find(string("+") + t + "+") != string::npos; };
        g_startNonD = has("nond") || has("nond2") || has("nond3"); g_nondDepth = has("nond3") ? 3 : has("nond2") ? 2 : 1; g_dsweep = has("dsweep"); g_guided = has("guided"); g_ctlLin = has("lin");
        printf("ctlbatch options:%s%s%s%s%s\n", g_ctlLin ? " linear search" : "", g_startNonD ? (g_nondDepth == 3 ? " start from a codimension-3 subspace of T (in-D test, three steps)" : g_nondDepth == 2 ? " start from a codimension-2 subspace of T (in-D test, two steps)" : " start from a hyperplane of T (in-D test)") : "", g_dsweep ? " D phase by the sweep" : "", g_guided ? " guided (only cores inside T get a D phase)" : "", has("dfirst") ? " D first" : "");
        return mode_ctlbatch(argv[2], argv[3], atoi(argv[4]), atoll(argv[5]), atoi(argv[6]), argv[7], has("dfirst") ? 3 : 0);
    }
    if (mode == "ctlrot") { // ctlrot <confs> <lut> <first line> <pool n> <rep>... [dsweep]: one control per rep, the functions of the
        // pool (the first n with k >= 3) in rotation, parallel over the reps
        init_conj(); CtlCtx X = ctl_ctx(argv[2]); g_logp1 = false;
        vector<vector<u64>> fs; vector<int> lns; FILE* fp = fopen(argv[3], "r"); static char buf[1 << 16]; int ln = 0, first = atoi(argv[4]), want = atoi(argv[5]);
        while (fgets(buf, sizeof buf, fp) && (int)fs.size() < want) { ln++; if (ln < first) continue; u8 F[256]; if (!read_lut(buf, F)) continue; auto B = forms_of_lut(F); if (profile(B).n4 >= 3) { fs.push_back(B); lns.push_back(ln); } }
        fclose(fp);
        vector<int> reps; for (int a = 6; a < argc; a++) { if (!strcmp(argv[a], "dsweep")) { g_dsweep = true; continue; } if (!strcmp(argv[a], "guided")) { g_guided = true; continue; } if (!strcmp(argv[a], "lin")) { g_ctlLin = true; continue; } reps.push_back(atoi(argv[a])); }
        int n = (int)reps.size(); vector<CtlResult> res(n); vector<int> fi(n);
#pragma omp parallel for schedule(dynamic, 1)
        for (int k = 0; k < n; k++) {
            g_forceRep = reps[k];
            for (int t = 0; t < (int)fs.size(); t++) { // the next function in rotation that embeds (k = 4 needs a member of its class)
                fi[k] = (reps[k] + t) % (int)fs.size();
                res[k] = ctl_non(X, fs[fi[k]], lns[fi[k]] * 31 + reps[k], 0, false);
                if (res[k].status != 4) break;
            }
        }
        int bad = 0;
        for (int k = 0; k < n; k++) {
            auto& R = res[k];
            long long tD = 0; for (int d = 0; d < 9; d++) tD += R.st.triesD[d];
            printf("rep %d (%s): line %d, %s, %d other spaces; tries in D %lld, lin D-phase candidates %lld, lin cores %lld%s%s\n", reps[k], X.cf[reps[k]].nonext ? "nonext" : "ext", lns[fi[k]], R.status == 0 ? "re-found" : "NOT re-found",
                R.status == 0 ? R.nsol - 1 : 0, tD, R.lst.dsteps, R.lst.cores[4] + R.lst.cores[5] + R.lst.cores[6] + R.lst.cores[7] + R.lst.cores[8], R.status ? " " : "", R.status ? R.why.c_str() : "");
            bad += R.status != 0;
        }
        fflush(stdout);
        return bad ? 3 : 0;
    }
    if (mode == "mapa") { // mapa <configs.txt> <search A's list of 6-member configurations>: which of our configurations each of A's lies on
        // (comparison only, written by Rowen after the verdict; seven_check's mapa with 6 members)
        init_conj();
        map<vector<u16>, int> confOf; map<int, string> confLine; { FILE* fp = fopen(argv[2], "r"); char line[1024]; while (fgets(line, sizeof line, fp)) { int ci; if (sscanf(line, "CONF %d", &ci) != 1) continue; vector<u16> Ns; const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; } confOf[key_of_Ns(Ns)] = ci; string s(line); confLine[ci] = s.substr(0, s.find(':')); } fclose(fp); }
        vector<string> names; vector<vector<Mask>> Ms; { FILE* fp = fopen(argv[3], "r"); static char line[1 << 14]; while (fgets(line, sizeof line, fp)) {
            char nm[128]; int off; if (sscanf(line, "%127s%n", nm, &off) != 1) continue; const char* p = line + off; char* e; vector<int> v; while (true) { long x = strtol(p, &e, 16); if (e == p) break; v.push_back((int)x); p = e; }
            if (v.empty() || v[0] != 0x6 || (int)v.size() != 1 + 24) { printf("unexpected line %s (%zu values)\n", nm, v.size()); return 1; }
            vector<Mask> M; for (int s = 0; s < 6; s++) {
                Mask m; mclear(m); for (int c = 1; c < 16; c++) { u8 x = 0; for (int b = 0; b < 4; b++) if (c >> b & 1) x ^= (u8)v[1 + 4 * s + b]; mset(m, x); }
                if (mcount(m) != 15) { printf("%s: member %d not a solid\n", nm, s); return 1; } M.push_back(m); }
            names.push_back(nm); Ms.push_back(M); } fclose(fp); }
        vector<int> hit(Ms.size(), -1); vector<long long> aut(Ms.size());
#pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < (int)Ms.size(); i++) { Canon K = canonP(Ms[i]); aut[i] = K.aut; auto it = confOf.find(K.key); if (it != confOf.end()) hit[i] = it->second; }
        map<int, int> mult; int unk = 0; for (size_t i = 0; i < Ms.size(); i++) { if (hit[i] < 0) unk++; else mult[hit[i]]++; printf("A %s -> %d (%s) canon aut %lld\n", names[i].c_str(), hit[i], hit[i] < 0 ? "none" : confLine[hit[i]].c_str(), aut[i]); }
        int dup = 0; for (auto& kv : mult) dup += kv.second > 1;
        printf("mapa: %zu of A's configurations, %d on none of ours, %zu of our %zu configurations hit, %d hit more than once\n", Ms.size(), unk, mult.size(), confOf.size(), dup);
        return 0;
    }
    return 1;
}

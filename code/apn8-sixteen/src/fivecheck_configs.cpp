// Classification of r pairwise-disjoint solids in PG(7,2), r = 4, 5, for the five-check.
//
// Normal form: V = A + A, R1 = {(x,0)}, R2 = {(0,y)}, R3 = graph(I), R_j = graph(T_j) for j > 3 (graph(T) = {(x, Tx)}).
// GL(8,2) is transitive on ordered triples of pairwise complementary solids, the stabiliser of (R1,R2,R3) is diag(g,g),
// acting on T by conjugation. So ordered r-tuples up to GL(8,2) = tuples (T_4..T_r) up to GL(4,2)-conjugation, and
// unordered ones = the same up to the relabelling moves below, derived from an explicit map each:
//   (12)  (x,y) -> (y,x)        : graph(T) -> graph(T^-1)
//   (13)  (x,y) -> (x,x+y)      : R1 <-> R3, R2 fixed, graph(T) -> graph(T+I)
//   (34)  (x,y) -> (x,T4^-1 y)  : graph(T4) -> graph(I), graph(I) -> graph(T4^-1), graph(T) -> graph(T4^-1 T)
//   (45)  swap T4, T5
// Orbit-stabiliser: #ordered tuples with R1,R2,R3 fixed in a class is |Omega_class|, and #ordered tuples overall is
// (|GL8|/20160) |Omega_class| = r! |GL8|/|Gamma|, so |Gamma| = r! 20160 / |Omega_class|.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
using namespace std;
typedef uint16_t M;  // row i in bits 4i..4i+3; (Mv)_i = parity(row_i & v)

static M mul(M a, M b) {  // a*b
    M r = 0;
    for (int i = 0; i < 4; i++) {
        int row = (a >> (4 * i)) & 15, acc = 0;
        for (int k = 0; k < 4; k++) if (row >> k & 1) acc ^= (b >> (4 * k)) & 15;
        r |= acc << (4 * i);
    }
    return r;
}
static const M ID = 0x8421;
static bool invertible(M a) {
    int rows[4]; for (int i = 0; i < 4; i++) rows[i] = (a >> (4 * i)) & 15;
    int rank = 0;
    for (int c = 0; c < 4; c++) {
        int p = -1;
        for (int i = rank; i < 4; i++) if (rows[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        swap(rows[p], rows[rank]);
        for (int i = 0; i < 4; i++) if (i != rank && (rows[i] >> c & 1)) rows[i] ^= rows[rank];
        rank++;
    }
    return rank == 4;
}
struct UF {
    vector<int> p;
    UF(int n) : p(n) { iota(p.begin(), p.end(), 0); }
    int f(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
    void u(int a, int b) { a = f(a); b = f(b); if (a != b) p[max(a, b)] = min(a, b); }
};

int main(int argc, char** argv) {
    vector<M> gl; static M inv[65536]; static int isgl[65536];
    for (int a = 0; a < 65536; a++) if (invertible((M)a)) { gl.push_back((M)a); isgl[a] = 1; }
    printf("|GL4| = %zu\n", gl.size());
    // inverses via powers: order divides lcm of element orders; do it by search over the group
    for (M a : gl) { M x = a, prev = ID; while (x != ID) { prev = x; x = mul(x, a); } inv[a] = prev; if (a == ID) inv[a] = ID; }
    for (M a : gl) if (mul(a, inv[a]) != ID) { puts("inverse bug"); return 1; }
    // generators of GL4: transvection and 4-cycle; check they generate
    M g1 = ID ^ 0x0002;              // row0 += e1
    M g2 = 0;                        // permutation matrix
    for (int i = 0; i < 4; i++) g2 |= (M)(1 << ((i + 1) % 4)) << (4 * i);
    {
        vector<char> seen(65536, 0); vector<M> st{ID}; seen[ID] = 1; int n = 1;
        while (!st.empty()) {
            M x = st.back(); st.pop_back();
            for (M g : {g1, g2}) { M y = mul(x, g); if (!seen[y]) { seen[y] = 1; n++; st.push_back(y); } }
        }
        printf("generated subgroup order %d\n", n);
        if (n != 20160) return 1;
    }
    vector<M> X;
    for (M t : gl) if (isgl[t ^ ID]) X.push_back(t);
    printf("|X| = %zu  (T with T, T+I invertible)\n", X.size());
    M gi[2] = {g1, g2}, gv[2] = {inv[g1], inv[g2]};

    FILE* out = fopen(argc > 1 ? argv[1] : "runs/check5/configs.txt", "w");

    // ---- r = 4
    {
        map<M, int> id; for (size_t i = 0; i < X.size(); i++) id[X[i]] = (int)i;
        UF conj((int)X.size()), full((int)X.size());
        for (size_t i = 0; i < X.size(); i++) {
            M t = X[i];
            for (int g = 0; g < 2; g++) {
                int j = id.at(mul(mul(gi[g], t), gv[g]));
                conj.u((int)i, j); full.u((int)i, j);
            }
            full.u((int)i, id.at(inv[t]));   // (12), also (34)
            full.u((int)i, id.at(t ^ ID));   // (13)
        }
        map<int, int> sz, osz;
        for (size_t i = 0; i < X.size(); i++) { sz[full.f((int)i)]++; osz[conj.f((int)i)]++; }
        printf("r=4: %zu ordered orbits, %zu unordered classes\n", osz.size(), sz.size());
        double sum = 0;
        for (auto& kv : sz) {
            long long gam = 24LL * 20160 / kv.second;
            if (24LL * 20160 % kv.second) puts("non-integer stabiliser!");
            sum += 1.0 / gam;
            printf("  rep T=%04x  |Omega|=%d  |Gamma|=%lld\n", X[kv.first], kv.second, gam);
            fprintf(out, "4 %04x 0000 %d %lld\n", X[kv.first], kv.second, gam);
        }
        printf("  sum 1/|Gamma| = %.10f\n", sum);
    }

    // ---- r = 5
    {
        vector<pair<M, M>> om;
        for (M t : X) for (M u : X) if (isgl[t ^ u]) om.push_back({t, u});
        printf("|Omega5| = %zu\n", om.size());
        map<uint32_t, int> id;
        auto key = [](M t, M u) { return ((uint32_t)t << 16) | u; };
        for (size_t i = 0; i < om.size(); i++) id[key(om[i].first, om[i].second)] = (int)i;
        UF conj((int)om.size()), full((int)om.size());
        for (size_t i = 0; i < om.size(); i++) {
            M t = om[i].first, u = om[i].second;
            for (int g = 0; g < 2; g++) {
                int j = id.at(key(mul(mul(gi[g], t), gv[g]), mul(mul(gi[g], u), gv[g])));
                conj.u((int)i, j); full.u((int)i, j);
            }
            full.u((int)i, id.at(key(inv[t], inv[u])));
            full.u((int)i, id.at(key(t ^ ID, u ^ ID)));
            full.u((int)i, id.at(key(inv[t], mul(inv[t], u))));
            full.u((int)i, id.at(key(u, t)));
        }
        map<int, long long> sz, osz;
        for (size_t i = 0; i < om.size(); i++) { sz[full.f((int)i)]++; osz[conj.f((int)i)]++; }
        printf("r=5: %zu ordered orbits, %zu unordered classes\n", osz.size(), sz.size());
        map<long long, int> byfix; double sum = 0;
        for (auto& kv : sz) {
            long long gam = 120LL * 20160 / kv.second;
            if (120LL * 20160 % kv.second) puts("non-integer stabiliser!");
            byfix[gam]++; sum += 1.0 / gam;
            printf("  rep T=%04x T'=%04x  |Omega|=%lld  |Gamma|=%lld\n", om[kv.first].first, om[kv.first].second, kv.second, gam);
            fprintf(out, "5 %04x %04x %lld %lld\n", om[kv.first].first, om[kv.first].second, kv.second, gam);
        }
        printf("by stabiliser order:");
        for (auto& kv : byfix) printf(" %lld:%d", kv.first, kv.second);
        printf("\n  sum 1/|Gamma| = %.10f\n", sum);
        // Langevin-Hou fix-5: total 20267057123180937216 over |GL(8,2)| = 5348063769211699200
        printf("  Langevin-Hou total/|GL8| = %.10f\n", 20267057123180937216.0 / 5348063769211699200.0);
    }
    fclose(out);
}

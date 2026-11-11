// Five-check, group side (used by the r = 4 control and for telling S8s apart; the r = 5 run uses none of it).
//   gamma <r>            Gamma of each class, listed explicitly; |Gamma| against configs.txt
//   reps <r> <cfg>       Gamma-orbits on form tuples (28^r): orbit sizes, Burnside count, bitmap runs/check5/reps_r<r>_c<cfg>.bin
//   classify <file> <r>  S8 bases, one per line (hex forms, last 8 tokens): class of solids, tuple, canonical id under Gamma
#define main sieve_main
#include "fivecheck_sieve.cpp"
#undef main

struct Lin { u8 col[8]; };  // A(e_i) = col[i]
static inline u8 linApply(const Lin& A, int v) { u8 r = 0; while (v) { int i = ctz64((u64)v); v &= v - 1; r ^= A.col[i]; } return r; }
static Lin linInverse(const Lin& A) {
    u8 m[8], id[8];
    for (int i = 0; i < 8; i++) { m[i] = A.col[i]; id[i] = (u8)(1 << i); }
    // column operations: reduce the columns of A to the identity, doing the same to id
    for (int c = 0; c < 8; c++) {
        int p = -1;
        for (int i = c; i < 8; i++) if (m[i] >> c & 1) { p = i; break; }
        if (p < 0) { puts("singular"); exit(1); }
        swap(m[p], m[c]); swap(id[p], id[c]);
        for (int i = 0; i < 8; i++) if (i != c && (m[i] >> c & 1)) { m[i] ^= m[c]; id[i] ^= id[c]; }
    }
    // now A (sum of id[i] coords) = e_i, so A^-1 e_i = id[i]
    Lin R; for (int i = 0; i < 8; i++) R.col[i] = id[i];
    return R;
}
// (A.B)(u, v) = B(A^-1 u, A^-1 v): the radical moves by A
static u64 actForm(const Lin& Ainv, u64 f) {
    u64 out = 0;
    for (int i = 0; i < 8; i++) {
        u8 row = ::apply(f, (int)Ainv.col[i]), o = 0;
        for (int j = 0; j < 8; j++) o |= (u8)(par(row & Ainv.col[j]) << j);
        out |= (u64)o << (8 * i);
    }
    return out;
}
static u16 mat4mul(u16 a, u16 b) {
    u16 r = 0;
    for (int i = 0; i < 4; i++) {
        int row = (a >> (4 * i)) & 15, acc = 0;
        for (int k = 0; k < 4; k++) if (row >> k & 1) acc ^= (b >> (4 * k)) & 15;
        r |= (u16)(acc << (4 * i));
    }
    return r;
}
static vector<u16> GL4; static u16 inv4[65536];
static void buildGL4() {
    for (int a = 0; a < 65536; a++) {
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
        if (rank == 4) GL4.push_back((u16)a);
    }
    for (u16 a : GL4) { u16 x = a, prev = 0x8421; while (x != 0x8421) { prev = x; x = mat4mul(x, a); } inv4[a] = prev; }
    inv4[0x8421] = 0x8421;
}

// Ordered solids W (masks) -> N with N(W0) = R1, N(W1) = R2, N(W2) = R3, and the T of each further solid.
static bool normaliser(const vector<Mask>& W, Lin& N, vector<u16>& Ts) {
    int e[4], f[4], ne = 0;
    {
        u8 b[4]; Mask sp;
        for (int v = 1; v < 256 && ne < 4; v++) if (W[0].get(v) && !sp.get(v)) { b[ne] = (u8)v; e[ne++] = v; sp = spanMask(b, ne); }
    }
    for (int i = 0; i < 4; i++) {
        f[i] = -1;
        for (int v = 1; v < 256; v++) if (W[1].get(v) && W[2].get(v ^ e[i])) f[i] = v;
        if (f[i] < 0) return false;
    }
    Lin Ninv; for (int i = 0; i < 4; i++) { Ninv.col[i] = (u8)e[i]; Ninv.col[4 + i] = (u8)f[i]; }
    N = linInverse(Ninv);
    Ts.clear();
    for (size_t j = 3; j < W.size(); j++) {
        u16 T = 0;
        for (int b = 0; b < 4; b++) {
            int found = -1;
            for (int v = 1; v < 256; v++) if (W[j].get(v)) { int w = linApply(N, v); if ((w & 15) == (1 << b)) found = w >> 4; }
            if (found < 0) return false;
            for (int i = 0; i < 4; i++) if (found >> i & 1) T |= (u16)(1 << (4 * i + b));
        }
        Ts.push_back(T);
    }
    return true;
}
struct GEl { Lin A, Ainv; int perm[5]; };
// all A with A(W_i) = rep's R_{sigma(i)} ... here for W the rep's own solids: Gamma; for other W: the maps onto the rep
static vector<GEl> mapsOnto(const Config& c, const vector<Mask>& W, bool all) {
    vector<GEl> out; int r = c.r;
    int p[5] = {0, 1, 2, 3, 4};
    vector<u16> target; for (int j = 1; j <= r - 3; j++) target.push_back(c.T[j]);
    do {
        vector<Mask> Wp; for (int i = 0; i < r; i++) Wp.push_back(W[p[i]]);
        Lin N; vector<u16> Ts;
        if (!normaliser(Wp, N, Ts)) { puts("normaliser failed"); exit(1); }
        for (u16 g : GL4) {
            bool ok = true;
            for (size_t j = 0; j < Ts.size() && ok; j++) if (mat4mul(mat4mul(g, Ts[j]), inv4[g]) != target[j]) ok = false;
            if (!ok) continue;
            Lin A;
            for (int i = 0; i < 8; i++) {
                int w = N.col[i];
                A.col[i] = (u8)(matVec(g, w & 15) | (matVec(g, w >> 4) << 4));
            }
            GEl el; el.A = A; el.Ainv = linInverse(A); memcpy(el.perm, p, sizeof p);
            out.push_back(el);
            if (!all) return out;
        }
    } while (next_permutation(p, p + r));
    return out;
}

static vector<Mask> repSolids(const Config& c) { vector<Mask> W; for (int i = 0; i < c.r; i++) W.push_back(c.solid[i]); return W; }

static vector<u64> rref(vector<u64> b) {
    Ech e; for (u64 f : b) e.add(f);
    // back-substitute and sort: a canonical basis of the span
    for (int i = 0; i < e.n; i++) for (int j = 0; j < e.n; j++) if (i != j && (e.row[j] >> e.piv[i] & 1)) e.row[j] ^= e.row[i];
    vector<u64> o(e.row, e.row + e.n); sort(o.begin(), o.end());
    return o;
}

int main(int argc, char** argv) {
    buildAlt(); buildGens(); buildGL4();
    string mode = argv[1];
    if (mode == "gamma") {
        int r = atoi(argv[2]);
        vector<Config> cs = loadConfigs(r);
        for (size_t i = 0; i < cs.size(); i++) {
            vector<GEl> G = mapsOnto(cs[i], repSolids(cs[i]), true);
            // each must fix the set of solids
            for (auto& g : G) for (int q = 0; q < r; q++) {
                Mask m; for (int v = 1; v < 256; v++) if (cs[i].solid[q].get(v)) m.set(linApply(g.A, v));
                bool hit = false; for (int t = 0; t < r; t++) if (m == cs[i].solid[t]) hit = true;
                if (!hit) { puts("not in Gamma"); return 1; }
            }
            printf("cfg %zu |Gamma| = %zu\n", i, G.size());
        }
        return 0;
    }
    if (mode == "reps") {
        int r = atoi(argv[2]), cfg = atoi(argv[3]);
        vector<Config> cs = loadConfigs(r);
        const Config& c = cs[cfg];
        vector<GEl> G = mapsOnto(c, repSolids(c), true);
        size_t gN = G.size();
        // table: g, solid i, form d -> (solid j, form d')
        vector<array<u8, 2>> tab(gN * r * 28);
        map<u64, pair<int, int>> id;
        for (int i = 0; i < r; i++) for (int d = 0; d < 28; d++) id[c.forms[i][d]] = {i, d};
        for (size_t g = 0; g < gN; g++) for (int i = 0; i < r; i++) for (int d = 0; d < 28; d++) {
            auto it = id.find(actForm(G[g].Ainv, c.forms[i][d]));
            if (it == id.end()) { puts("image is not a form of the rep"); return 1; }
            tab[(g * r + i) * 28 + d] = {(u8)it->second.first, (u8)it->second.second};
        }
        size_t total = 1; for (int i = 0; i < r; i++) total *= 28;
        vector<u8> seen(total, 0), rep(total, 0);
        long long orbits = 0, burnside = 0; map<long long, long long> sizes;
        for (size_t t = 0; t < total; t++) {
            int d[5]; { size_t x = t; for (int i = r - 1; i >= 0; i--) { d[i] = (int)(x % 28); x /= 28; } }
            if (seen[t]) { continue; }
            rep[t] = 1; orbits++; long long sz = 0;
            for (size_t g = 0; g < gN; g++) {
                int nd[5];
                for (int i = 0; i < r; i++) { auto& e2 = tab[(g * r + i) * 28 + d[i]]; nd[e2[0]] = e2[1]; }
                size_t u = 0; for (int i = 0; i < r; i++) u = u * 28 + nd[i];
                if (!seen[u]) { seen[u] = 1; sz++; }
            }
            sizes[sz]++;
            if (gN % sz) { puts("orbit size does not divide |Gamma|"); return 1; }
        }
        for (size_t g = 0; g < gN; g++) { // Burnside: tuples fixed by g
            for (size_t t = 0; t < total; t++) {
                int d[5]; { size_t x = t; for (int i = r - 1; i >= 0; i--) { d[i] = (int)(x % 28); x /= 28; } }
                bool fix = true;  // the image tuple equals the tuple; solids may move
                for (int i = 0; i < r && fix; i++) { auto& e2 = tab[(g * r + i) * 28 + d[i]]; if (e2[1] != d[e2[0]]) fix = false; }
                burnside += fix;
            }
            if (gN > 1000 && g == 0) {} // kept simple: total * |Gamma| at most 2.7e9 steps, most exit at once
        }
        long long sum = 0; for (auto& kv : sizes) sum += kv.first * kv.second;
        printf("cfg %d: |Gamma| %zu, tuples %zu, orbits %lld, sum of orbit sizes %lld, Burnside %lld/%zu = %.6f\n", cfg, gN, total, orbits, sum, burnside, gN, (double)burnside / gN);
        for (auto& kv : sizes) printf("  size %lld: %lld\n", kv.first, kv.second);
        char path[128]; snprintf(path, sizeof path, "runs/check5/reps_r%d_c%02d.bin", r, cfg);
        FILE* f = fopen(path, "wb"); fwrite(rep.data(), 1, total, f); fclose(f);
        return 0;
    }
    if (mode == "replist") { // replist <r> <cfg>: the orbit representatives as a tuple list, from the reps bitmap
        int r = atoi(argv[2]), cfg = atoi(argv[3]);
        size_t total = 1; for (int i = 0; i < r; i++) total *= 28;
        char path[128]; snprintf(path, sizeof path, "runs/check5/reps_r%d_c%02d.bin", r, cfg);
        vector<u8> rep(total); FILE* f = fopen(path, "rb"); fread(rep.data(), 1, total, f); fclose(f);
        snprintf(path, sizeof path, "runs/check5/reps_r%d_c%02d.txt", r, cfg);
        FILE* o = fopen(path, "w"); long long n = 0;
        for (size_t t = 0; t < total; t++) if (rep[t]) {
            int d[5]; { size_t x = t; for (int i = r - 1; i >= 0; i--) { d[i] = (int)(x % 28); x /= 28; } }
            for (int i = 0; i < r; i++) fprintf(o, "%d ", d[i]);
            fprintf(o, "\n"); n++;
        }
        fclose(o); printf("%lld reps\n", n);
        return 0;
    }
    if (mode == "k4list") { // k4list: all tuples of 5 forms spanning dimension 4 (passing the S check), per class, and their Gamma-orbits
        vector<Config> cs = loadConfigs(5);
        for (size_t ci = 0; ci < cs.size(); ci++) {
            const Config& c = cs[ci];
            vector<array<int, 5>> found; Stats ls;
            #pragma omp parallel for schedule(dynamic, 1)
            for (int t = 0; t < 28 * 28 * 28; t++) {
                Level L; L.c = &c; Stats s2; vector<array<int, 5>> loc;
                auto proc = [&](const u64* b, int k) { if (k == 4) loc.push_back({L.cur[0], L.cur[1], L.cur[2], L.cur[3], L.cur[4]}); };
                L.run(t / 784, (t / 28) % 28, t % 28, s2, proc);
                #pragma omp critical
                { for (auto& x : loc) found.push_back(x); ls.add(s2); }
            }
            sort(found.begin(), found.end());
            vector<GEl> G = mapsOnto(c, repSolids(c), true);
            map<u64, pair<int, int>> id;
            for (int i = 0; i < 5; i++) for (int d = 0; d < 28; d++) id[c.forms[i][d]] = {i, d};
            set<array<int, 5>> seen; vector<array<int, 5>> reps; map<long long, int> sizes;
            map<array<int, 5>, int> orbitOf;
            for (auto& t : found) {
                if (seen.count(t)) continue;
                reps.push_back(t); long long sz = 0;
                for (auto& g : G) {
                    array<int, 5> im;
                    for (int i = 0; i < 5; i++) { auto q = id.at(actForm(g.Ainv, c.forms[i][t[i]])); im[q.first] = q.second; }
                    if (seen.insert(im).second) sz++;
                    orbitOf[im] = (int)reps.size() - 1;
                }
                sizes[sz]++;
            }
            if (!found.empty()) { // every tuple, and its orbit, for the invariance check
                char pa[128], pb[128];
                snprintf(pa, sizeof pa, "runs/check5/k4all_c%02d.txt", (int)ci); snprintf(pb, sizeof pb, "runs/check5/k4orbit_c%02d.txt", (int)ci);
                FILE* fa = fopen(pa, "w"); FILE* fb = fopen(pb, "w");
                for (auto& t : found) { fprintf(fa, "%d %d %d %d %d\n", t[0], t[1], t[2], t[3], t[4]); fprintf(fb, "%d\n", orbitOf.at(t)); }
                fclose(fa); fclose(fb);
            }
            bool closed = seen.size() == found.size();
            printf("cfg %zu: tuples %lld sets %lld k4 %lld k5 %lld; k4 tuples %zu, orbits %zu, closed %d, sizes", ci, ls.tuples, ls.sets, ls.setsK[4], ls.setsK[5], found.size(), reps.size(), (int)closed);
            for (auto& kv : sizes) printf(" %lld:%d", kv.first, kv.second);
            printf("\n"); fflush(stdout);
            char path[128]; snprintf(path, sizeof path, "runs/check5/k4_c%02d.txt", (int)ci);
            FILE* o = fopen(path, "w");
            for (auto& t : reps) fprintf(o, "%d %d %d %d %d\n", t[0], t[1], t[2], t[3], t[4]);
            fclose(o);
        }
        puts("K4LISTDONE");
        return 0;
    }
    if (mode == "where") {
        // where <file> <cfg>: for each S8 of the file (rank-4 count 4, class of solids cfg): its canonical id and the index in
        // reps_r4_c<cfg>.txt of the orbit of its tuple of rank-4 forms. Tells which representatives carry which function.
        int cfg = atoi(argv[3]); int r = 4;
        vector<Config> cs = loadConfigs(r); const Config& c = cs[cfg];
        vector<GEl> G = mapsOnto(c, repSolids(c), true);
        char path[128]; snprintf(path, sizeof path, "runs/check5/reps_r4_c%02d.txt", cfg);
        map<array<int, 4>, int> repIndex;
        { FILE* f = fopen(path, "r"); int a, b2, c2, d, i = 0; while (fscanf(f, "%d %d %d %d", &a, &b2, &c2, &d) == 4) repIndex[{a, b2, c2, d}] = i++; fclose(f); }
        map<u64, pair<int, int>> id;
        for (int i = 0; i < r; i++) for (int d = 0; d < 28; d++) id[c.forms[i][d]] = {i, d};
        FILE* f = fopen(argv[2], "r"); char line[4096]; int ln = 0;
        map<vector<u64>, int> classes;
        while (fgets(line, sizeof line, f)) {
            vector<u64> b; vector<string> tok; char* p = strtok(line, " \r\n");
            while (p) { tok.push_back(p); p = strtok(nullptr, " \r\n"); }
            if (tok.size() < 8) continue;
            for (size_t i = tok.size() - 8; i < tok.size(); i++) b.push_back(strtoull(tok[i].c_str(), nullptr, 16));
            ln++;
            vector<Mask> W;
            for (int m = 1; m < 256; m++) { u64 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= b[i]; u8 bb[8]; if (radical(x, bb) == 4) W.push_back(spanMask(bb, 4)); }
            auto v = mapsOnto(c, W, false);
            if (v.empty()) { printf("line %d: not in this class of solids\n", ln); continue; }
            vector<u64> nb; for (u64 x : b) nb.push_back(actForm(v[0].Ainv, x));
            array<int, 4> tup = {-1, -1, -1, -1};
            for (int m = 1; m < 256; m++) { u64 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= nb[i]; u8 bb[8]; if (radical(x, bb) == 4) { auto q = id.at(x); tup[q.first] = q.second; } }
            vector<u64> best; int rep = -1;
            for (auto& g : G) {
                vector<u64> im; for (u64 x : nb) im.push_back(actForm(g.Ainv, x));
                im = rref(im);
                if (best.empty() || im < best) best = im;
                array<int, 4> t2; for (int i = 0; i < 4; i++) { auto q = id.at(actForm(g.Ainv, c.forms[i][tup[i]])); t2[q.first] = q.second; }
                auto it = repIndex.find(t2); if (it != repIndex.end()) rep = it->second;
            }
            auto it = classes.find(best); int cls;
            if (it == classes.end()) { cls = (int)classes.size(); classes[best] = cls; } else cls = it->second;
            printf("line %d class %d rep %d id", ln, cls, rep);
            for (u64 x : best) printf(" %016llx", (unsigned long long)x);
            puts("");
        }
        printf("%d lines, %zu classes\n", ln, classes.size());
        return 0;
    }
    if (mode == "slevel") {
        // slevel <r> <samples> <seed>: the search's accept/reject of a tuple (prefix test, then leaf test) against the
        // definition, with radicals by plain evaluation at every point: S is admissible iff every point lies in at most
        // one radical of a nonzero element of S, and no element other than the B_i has a radical point on a solid.
        int r = atoi(argv[2]), ns = atoi(argv[3]); srand(atoi(argv[4]));
        vector<Config> cs = loadConfigs(r);
        long long agree = 0, accepted = 0, bad = 0;
        for (size_t ci = 0; ci < cs.size(); ci++) {
            const Config& c = cs[ci]; Level L; L.c = &c;
            long long acc = 0;
            for (int it = 0; it < ns; it++) {
                u64 B[5]; for (int i = 0; i < r; i++) B[i] = c.forms[i][(rand() * 32768 + rand()) % 28];
                bool search = true;
                for (int n = 1; n <= r && search; n++) if (!L.partial(B, n)) search = false;
                if (search) { Stats ls; bool got = false; auto proc = [&](const u64*, int) { got = true; }; L.leaf(B, ls, proc); search = got; }
                set<u64> el;
                for (int m = 1; m < (1 << r); m++) { u64 x = 0; for (int i = 0; i < r; i++) if (m >> i & 1) x ^= B[i]; if (x) el.insert(x); }
                int cover[256] = {0}; bool def = true;
                for (u64 x : el) {
                    bool isB = false; for (int i = 0; i < r; i++) if (x == B[i]) isB = true;
                    for (int a = 1; a < 256; a++) {
                        bool in = true;
                        for (int q = 0; q < 8 && in; q++) if (evalForm(x, a, 1 << q)) in = false;
                        if (!in) continue;
                        cover[a]++;
                        if (cover[a] > 1) def = false;
                        if (!isB && c.solidsAll.get(a)) def = false;
                    }
                }
                if (def == search) agree++; else { bad++; printf("cfg %zu DISAGREE search %d def %d\n", ci, (int)search, (int)def); }
                acc += def;
            }
            accepted += acc;
        }
        printf("slevel r=%d: %lld agree, %lld disagree, %lld admissible\n", r, agree, bad, accepted);
        return bad ? 1 : 0;
    }
    if (mode == "plant") {
        // plant <file> <r> [pivot]: guided unit test of the e = 3 path. For each S8 in the file (rank-4 count r), go to the
        // class representative, take S = span(the r rank-4 forms) + t for t in each coset of S8 / S, and process S:
        // the sieve and the assembly must return S8 itself.
        int r = atoi(argv[3]); bool piv = argc > 4 && !strcmp(argv[4], "pivot");
        vector<Config> cs = loadConfigs(r);
        FILE* f = fopen(argv[2], "r");
        char line[4096]; long long tried = 0, hit = 0, lines = 0;
        while (fgets(line, sizeof line, f)) {
            vector<u64> b; vector<string> tok; char* p = strtok(line, " \r\n");
            while (p) { tok.push_back(p); p = strtok(nullptr, " \r\n"); }
            if (tok.size() < 8) continue;
            for (size_t i = tok.size() - 8; i < tok.size(); i++) b.push_back(strtoull(tok[i].c_str(), nullptr, 16));
            lines++;
            vector<Mask> W;
            for (int m = 1; m < 256; m++) {
                u64 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= b[i];
                u8 bb[8]; if (radical(x, bb) == 4) W.push_back(spanMask(bb, 4));
            }
            int found = -1; GEl to;
            for (size_t ci = 0; ci < cs.size() && found < 0; ci++) { auto v = mapsOnto(cs[ci], W, false); if (!v.empty()) { found = (int)ci; to = v[0]; } }
            vector<u64> nb; for (u64 x : b) nb.push_back(actForm(to.Ainv, x));
            vector<u64> want = rref(nb);
            // the rank-4 forms in normal coordinates
            vector<u64> S; Ech es; vector<u64> all;
            for (int m = 1; m < 256; m++) {
                u64 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= nb[i];
                all.push_back(x);
                u8 bb[8]; if (radical(x, bb) == 4 && es.add(x)) S.push_back(x);
            }
            int k0 = (int)S.size();
            set<u64> cosets;
            for (u64 x : all) { u64 red = es.reduce(x); if (red) cosets.insert(red); }
            for (u64 t : cosets) {
                vector<u64> basis = S; basis.push_back(t);
                Worker w; w.c = &cs[found]; w.fl.pivot = piv;
                w.process(basis.data(), k0 + 1);
                tried++;
                bool got = false;
                for (auto& res : w.results) if (rref(res) == want) got = true;
                if (got) hit++;
                else printf("line %lld: missed with coset %016llx (f1 %lld, found %lld)\n", lines, (unsigned long long)t, w.st.f1, w.st.found);
            }
        }
        printf("plant: %lld functions, %lld (S8, coset) pairs, S8 returned in %lld\n", lines, tried, hit);
        return tried == hit ? 0 : 1;
    }
    if (mode == "classify") {
        int r = atoi(argv[3]);
        vector<Config> cs = loadConfigs(r);
        vector<vector<GEl>> Gam;
        for (auto& c : cs) Gam.push_back(mapsOnto(c, repSolids(c), true));
        FILE* f = fopen(argv[2], "r");
        char line[4096]; map<vector<u64>, int> classes; int nline = 0;
        map<int, int> perCfg;
        while (fgets(line, sizeof line, f)) {
            vector<u64> b; vector<string> tok; char* p = strtok(line, " \r\n");
            while (p) { tok.push_back(p); p = strtok(nullptr, " \r\n"); }
            if (tok.size() < 8) continue;
            bool hex = true; for (size_t i = tok.size() - 8; i < tok.size(); i++) if (tok[i].size() != 16) hex = false;
            if (!hex) continue;
            for (size_t i = tok.size() - 8; i < tok.size(); i++) b.push_back(strtoull(tok[i].c_str(), nullptr, 16));
            nline++;
            vector<Mask> W; vector<u64> B;
            for (int m = 1; m < 256; m++) {
                u64 x = 0; for (int i = 0; i < 8; i++) if (m >> i & 1) x ^= b[i];
                u8 bb[8]; int n = radical(x, bb);
                if (n == 4) { W.push_back(spanMask(bb, 4)); B.push_back(x); }
            }
            if ((int)W.size() != r) { printf("line %d has %zu rank-4 forms\n", nline, W.size()); continue; }
            int found = -1; GEl to;
            for (size_t ci = 0; ci < cs.size() && found < 0; ci++) {
                auto v = mapsOnto(cs[ci], W, false);
                if (!v.empty()) { found = (int)ci; to = v[0]; }
            }
            if (found < 0) { puts("no class of solids fits"); return 1; }
            vector<u64> nb; for (u64 x : b) nb.push_back(actForm(to.Ainv, x));
            vector<u64> best;
            for (auto& g : Gam[found]) {
                vector<u64> im; for (u64 x : nb) im.push_back(actForm(g.Ainv, x));
                im = rref(im);
                if (best.empty() || im < best) best = im;
            }
            // tag the class with its config too
            vector<u64> key = best; key.push_back((u64)found);
            auto it = classes.find(key);
            int cls;
            if (it == classes.end()) { cls = (int)classes.size(); classes[key] = cls; perCfg[found]++; printf("class %d cfg %d:", cls, found); for (u64 x : best) printf(" %016llx", (unsigned long long)x); printf("\n"); }
        }
        printf("%d lines, %zu classes;", nline, classes.size());
        for (auto& kv : perCfg) printf(" cfg %d: %d", kv.first, kv.second);
        printf("\n");
        return 0;
    }
    return 1;
}

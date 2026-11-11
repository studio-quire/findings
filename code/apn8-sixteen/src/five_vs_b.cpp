// five_vs_b <p5.txt> <out>: on search A's configurations, every tuple of five forms under search B's admissibility rule, against A's.
// A keeps a span (a core) only if each element is bent, rank 6 with its radical a line of holes, or rank 4 on one of the solids. B
// (five-check.md section 4) keeps it if no element but the five forms has a radical point on a solid and all radicals are disjoint,
// saying nothing about the rank of the rest. So B's sets should be A's cores plus exactly those with a further rank-4 form, its
// radical a sixth solid among the holes, and per configuration the counts should equal B's table in section 5.
// Primitives are the reviewer's (review-five-orbits.md), not cover.cpp's or B's. Rowen (Claude Opus 5.5), 2026-11-06.
#define main rv5_main
#include "rv5_indep.cpp"
#undef main

struct Count { long long strict = 0, loose = 0, loose4 = 0, strict4 = 0, extra = 0, extra4 = 0, other = 0; };

// 0 illegal under B's rule; otherwise 1, with extra set if an element other than the chosen has rank 4, other if any other kind of
// element B allows and A doesn't turns up (none should)
static int legal_b(const Cfg& C, const u64* ch, int n, u64* span, int* sn, bool* extra, bool* other) {
    int m = 1; span[0] = 0;
    for (int i = 0; i < n; i++) {
        bool in = false;
        for (int j = 0; j < m; j++) if (span[j] == ch[i]) in = true;
        if (in) continue;
        for (int j = 0; j < m; j++) span[m + j] = span[j] ^ ch[i];
        m *= 2;
    }
    *sn = m; *extra = false; *other = false;
    bool used[256] = {false};
    for (int e = 1; e < m; e++) {
        u8 p[256]; int c; int d = rad_pts(span[e], p, &c);
        if (d == 0) continue;
        bool chosen = false;
        for (int i = 0; i < n; i++) if (span[e] == ch[i]) chosen = true;
        bool on = false;
        for (int i = 0; i < c; i++) if (C.where[p[i]] >= 0) on = true;
        if (on && !chosen) {
            // on a prefix, a dim-4 set's last form may already be in the span: fine if it's a form on a solid not yet chosen
            int s = C.where[p[0]];
            bool later = d == 4 && s >= n;
            for (int i = 0; i < c; i++) if (C.where[p[i]] != s) later = false;
            if (!later) return 0;
            chosen = true;
        }
        if (!chosen) { if (d == 4) *extra = true; else if (d != 2) *other = true; }
        for (int i = 0; i < c; i++) { if (used[p[i]]) return 0; used[p[i]] = true; }
    }
    return 1;
}

static void walk(const Cfg& C, u64* ch, int depth, Count& t) {
    u64 span[64]; int m; bool extra, other;
    for (u64 f : C.forms[depth]) {
        ch[depth] = f;
        if (!legal_b(C, ch, depth + 1, span, &m, &extra, &other)) continue;
        if (depth + 1 < C.k) { walk(C, ch, depth + 1, t); continue; }
        t.loose++; if (m == 16) t.loose4++;
        u64 s2[64]; int m2;
        bool a = legal(C, ch, C.k, s2, &m2);
        if (a) { t.strict++; if (m == 16) t.strict4++; }
        if (a == extra) t.other++;   // A's rule must fail exactly when there's an extra rank-4 form
        if (extra) { t.extra++; if (m == 16) t.extra4++; }
        if (other) t.other++;
    }
}

int main(int argc, char** argv) {
    if (argc != 3) { fprintf(stderr, "usage: five_vs_b <p5.txt> <out>\n"); return 2; }
    std::vector<Cfg> cs;
    FILE* f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    static char line[4096];
    while (fgets(line, sizeof line, f)) {
        char* p = strtok(line, " \r\n"); if (!p) continue;
        Cfg C; C.tag = p; C.k = atoi(strtok(nullptr, " \r\n"));
        memset(C.where, -1, sizeof C.where);
        for (int s = 0; s < C.k; s++) {
            std::vector<u8> b;
            for (int i = 0; i < 4; i++) b.push_back((u8)strtol(strtok(nullptr, " \r\n"), nullptr, 16));
            std::vector<u8> pts;
            for (int m = 1; m < 16; m++) { u8 x = 0; for (int i = 0; i < 4; i++) if (m >> i & 1) x ^= b[i]; pts.push_back(x); C.where[x] = (signed char)s; }
            C.sol.push_back(pts);
            C.forms.push_back(forms_on(b, pts));
        }
        cs.push_back(C);
    }
    fclose(f);
    std::vector<Count> r(cs.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)cs.size(); i++) { u64 ch[8]; walk(cs[i], ch, 0, r[i]); }
    FILE* o = fopen(argv[2], "w"); if (!o) { perror(argv[2]); return 1; }
    Count T;
    for (size_t i = 0; i < cs.size(); i++) {
        fprintf(o, "%s B-rule %lld dim4 %lld | A-rule %lld dim4 %lld | extra rank-4 %lld dim4 %lld | unexplained %lld\n", cs[i].tag.c_str(), r[i].loose, r[i].loose4,
                r[i].strict, r[i].strict4, r[i].extra, r[i].extra4, r[i].other);
        T.loose += r[i].loose; T.loose4 += r[i].loose4; T.strict += r[i].strict; T.strict4 += r[i].strict4; T.extra += r[i].extra; T.extra4 += r[i].extra4; T.other += r[i].other;
    }
    fprintf(o, "total B-rule %lld dim4 %lld | A-rule %lld dim4 %lld | extra rank-4 %lld dim4 %lld | unexplained %lld\n", T.loose, T.loose4, T.strict, T.strict4, T.extra, T.extra4, T.other);
    fclose(o);
    return 0;
}

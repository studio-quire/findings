// Closure of a solids-only slice's S8s under Gamma: an image with the slice's form on R1 must be in the slice's found set.
// The search uses no symmetry, so a T lost by the assembly (or sieve) shows as an image outside the set.
#define main gamma_main
#include "rvb5_gamma.cpp"
#undef main
int main(int argc, char** argv) {
    gl4();
    int cfg = atoi(argv[1]), q = atoi(argv[2]); const char* path = argv[3];
    auto cf = readCfg(); vector<Cfg> c5; for (auto& c : cf) if (c.r == 5) c5.push_back(c);
    auto G = gammaOf(mk(5, c5[cfg].T4, c5[cfg].T5));
    map<u64, array<int, 2>> id;
    { FILE* f = fopen("runs/check5/review/forms.txt", "r"); int ci, i, d; unsigned long long x; while (fscanf(f, "%d %d %d %llx", &ci, &i, &d, &x) == 4) if (ci == cfg) id[x] = {i, d}; fclose(f); }
    vector<vector<u64>> raw; set<vector<u64>> all;
    { FILE* f = fopen(path, "r"); char line[2048];
      while (fgets(line, sizeof line, f)) { if (strncmp(line, "S8 ", 3)) continue; char* p = strchr(line, ':') + 1; vector<u64> b; for (int i = 0; i < 8; i++) { char* e; b.push_back(strtoull(p, &e, 16)); p = e; } raw.push_back(b); all.insert(canon(b)); }
      fclose(f); }
    long long inSlice = 0, outside = 0, identical = 0;
    set<vector<u64>> reachedAsImage;
    for (auto& b : raw) for (auto& A : G) {
        bool toSlice = false;
        for (int j = 0; j < 5; j++) { auto it = id.find(pull(b[j], A)); if (it == id.end()) { puts("image of a form is not a form of the rep"); return 1; } if (it->second[0] == 0 && it->second[1] == q) toSlice = true; }
        if (!toSlice) continue;
        vector<u64> ib; for (u64 x : b) ib.push_back(pull(x, A));
        auto cn = canon(ib); inSlice++;
        if (cn == canon(b)) identical++;
        if (all.count(cn)) reachedAsImage.insert(cn); else outside++;
    }
    printf("|Gamma| %zu, S8 %zu distinct %zu; images staying in the slice %lld (fixing the S8 %lld); outside the found set %lld; found S8s hit as images %zu\n", G.size(), raw.size(), all.size(), inSlice, identical, outside, reachedAsImage.size());
    return 0;
}

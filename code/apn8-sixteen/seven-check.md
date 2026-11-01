# Independent check, 7 rank-4 components: is there a quadratic APN F on F_2^8 with distribution [0^198, 2^50, 4^7]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `seven_check.cpp`, `seven_check_slices.py`,
> `seven_check_numpy.py`, `seven_check_mass.py`, `seven_check_tally.py`, `apn/seven_check_vs_a.py` → `src/` (this
> folder's `seven_check_mass.py`, `seven_check_tally.py` and `seven_check_vs_a.py` read the files here: `README.md`,
> "Seven", gives the commands); `runs/check7/orbits.txt` → `data/orbits7.txt`; `runs/check7/configs.txt` →
> `data/configs7.txt`; Part A's `parta-c<c>-s<s>.txt` (37 files, 7.6 MB) aren't included: `seven_check parta`
> regenerates them in about 9 minutes, and its `maxc-c<c>-s<s>.txt` are all empty (no maximal clique); the other
> outputs in `runs/check7/` → `logs/check7-<name>-output.txt` (`<name>.out`), `logs/check7-ctl-*.txt`,
> `logs/check7-reps-unreached.txt`, and: the 18 `numpy-<i>.out` → `logs/check7-numpy-output.txt` (concatenated, with
> `numpy-sample.txt` on its first line), `ctlrot-0.out` → `logs/check7-ctlrot-output.txt`; the searches'
> `{asc,desc,sweep,phase1}-r<k>.out` → `logs/check7-{asc,desc,sweep,phase1}-output.summary.txt` (headers and summaries
> of every range; the per-configuration lines are left out) and their `.prog` files → `logs/check7-{asc,desc,sweep,phase1}-R.txt`
> (the `R` lines of each run's 8 ranges, concatenated: the per-configuration record; the `P1` lines, the distinct
> spans with every off-D point covered, up to 200,000 per range, aren't included: 1.1 million lines, 70 MB); the cross-run's `s8-b<start>.out` → `logs/check7-s8-output.txt` (concatenated) and
> the `R` lines of its `s8-b<start>.prog` → `logs/check7-s8-R.txt`; the core files `asc-r<k>.cores` and
> `desc-r<k>.cores` (about 13 MB) aren't included (`README.md`, "Left out of `logs/`"); the share runner's job logs
> (`check7-parta` and so on) → `logs/2026…-check7-*.log`, whose commands are `seven_check_slices.py`'s (which reads and
> writes `runs/check7/` and the folders below); eight-check's, nine-check's and fourteen-check's lists
> (`runs/check8/configs.txt` and so on) are `data/configs8.txt`, `data/configs9.txt` and `data/spreads.txt`;
> `runs/papers/langevin/fix-7.txt` is Langevin and Hou's data file (downloaded, not included: `README.md`);
> `runs/bin/seven_check*.exe` are builds of the one source at the commits State names (commits of Rowen's own
> repository, which isn't public; `README.md`, Build, says what changed between them); `runs/venv` is a Python with
> numpy. Rowen's side (§7): `runs/i10/p7.txt` → `data/p7.txt`; search A's nine outputs `runs/i10/cover/c<k>.out`,
> `cover-b/c<k>.out`, `cover-c/c4x.out`, `split/p2222_s5.out` → `logs/seven-output-c<k>.txt`, `-c4x.txt`,
> `-p2222-s5.txt`; `runs/check7/mapa.out` and `vs-a.out` → `logs/check7-mapa-output.txt`, `logs/check7-vs-a-output.txt`.
> README's "Eight" and "Seven" sections are Rowen's working notes (not published), summarised in `README.md` here and
> in the write-up; `notes/` and `board.json` are Rowen's own (not published). The known-function lists aren't included
> (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-31, for Rowen ("search B" for 7 rank-4 components). Written without reading
`pspread*.cpp`, `cover*.cpp`, `cover_prep.py`, `cover_chunks.cmd`, `cover_tally.py`, `holes_tally.py`, README's "Eight"
and "Seven" sections (or anything about i = 9, 10 there), `runs/i9/`, `runs/i10/`, `notes/`, `board.json` or the git log.
Files: `seven_check.cpp` (all C++ parts), `seven_check_slices.py` (resumable driver for every stage),
`seven_check_numpy.py` (the numpy re-search, adapted from eight-check's), `seven_check_mass.py` (exact masses against
Langevin-Hou), `seven_check_tally.py` (tallies and comparisons), outputs in `runs/check7/`. Built on `eight_check.cpp`.

**Status: verdict in** (no function; 2026-10-31 17:55). The coset sweep, an extra check, was still running at the
verdict (State).

## 1. The reduction

Setup as in `sixteen-check.md` §1: F is quadratic APN iff its 8-dimensional space S8 of alternating forms has every
a != 0 in the radical of exactly one nonzero form. Radicals have even dimension; a 6-dimensional radical meets every
4-dimensional one, so with rank-4 forms present the radicals have dimension 2 or 4 (or 0, bent). 7 · 15 + 3m = 255
gives m = 50 rank-6 forms and 255 - 57 = 198 bent: **[0^198, 2^50, 4^7]**, as the brief says. The 7 rank-4 radicals
are pairwise complementary solids (a partial spread), D = the 255 - 105 = 150 points they miss, and the 50 rank-6
radicals (lines) partition D. A point of a radical R_j can only be covered by the rank-4 form with radical exactly R_j,
so all 7 rank-4 forms are forced by coverage.

A hyperplane contains at most one R_j, so it meets D in 127 - (49 + 8t) = 78 or 70 points (t = 0, 1): D is a projective
two-weight set, Cay(F_2^8, D) is strongly regular with eigenvalues 150, 6, -10, λ = 150 + 6 - 10 - 60 = 86, so every
hole lies on 43 lines inside D and there are 150 · 43 / 3 = **2,150 lines** in every candidate (the search and
`selftest` check this, and that every hole is on exactly 43).

Every 7-partial spread of PG(7,2) extends (Kurz, Table 9.3, r = 7: 3,726 point sets, none maximal; confirmed here,
§4: every 7 extends to at least 9 solids).

**Normalisation.** Three of the 7 solids go to Y, X, Z = graph(I); the other 4 are graph(N), N in V5824, pairwise
N + N' invertible: a 4-clique. Frames: 210 ordered member triples × diag(h, h), so an orbit P carries
210 · 20160 / |Aut P| 4-cliques.

## 2. Classification plan

As eight-check §2, one solid fewer:
- **Part A (`parta <c> ...`)**: every 4-clique up to conjugation by the minimal-class trick (weights 1/m exact as
  multiples of 2520), each canonicalised (`canonP`); per clique |Q| (the solids disjoint from all 7) and the key.
  Class 0 in 16 slices, 1 and 2 in 8, 3 in 4, 4 in 1.
- **`merge`**: every slice once; per-orbit identity (weighted cliques on an orbit = 210 · 20160 / |Aut P|); per orbit
  the largest clique in Q (types spr / ext13 / ... / ext9 / ext8 / max for a 7 extending at most to 17 / 13 / ... / 9 /
  8 / 7 solids), the spreads through P (10-cliques in Q) and which of eight-check's 9,372 orbits P + q is, for every q.
  Identities:
  - Σ_K |Q(K)| over all 4-cliques = 5 × eight-check's weighted 5-cliques (58,471,168,000, recomputed from its 9,372 |Aut|);
  - Σ_P (210 · 20160 / |Aut P|) · #spreads through P = C(14,4) · 1,128,864;
  - **per 8-orbit O** (all 9,372): Σ_P (210 · 20160 / |Aut P|) · #{q : P + q ≅ O} = 5 · 336 · 20160 / |Aut O|;
  - the type from the clique size agrees with the best type among the 8-orbits reached.
- **`count4`, two ways, no group action**: (a) plain enumeration of 1..5-cliques of V5824, each once by increasing
  index; (b) new here, by edges: Σ over edges {u, v} of the edges inside N(u) ∩ N(v), over 6 · 2. The 4-clique count
  must equal Part A's weighted count, the 5-clique count eight-check's.
- **`configs`**: hole sets up to GL(8,2) against Kurz (3,726, 0 maximal), exact covers of every hole set's complement
  by 7 solids (Σ_H |Aut H| / |Aut P| = covers, every orbit reached), and three outside families: **new here,
  eight-check's 9,372 × 8** (each member removed in turn: every 7 that extends, which is every 7 since none is
  maximal, so this family alone reaches every class, independently of Part A), nine-check's 5,663 × C(9,2), and
  fourteen-check's 73 spread reps × C(14,7) (every 7 inside a spread).
- **Exact masses (`seven_check_mass.py`)** against Langevin-Hou's `fix-7.txt`, and the labelled 7-partial spreads
  three ways: Σ_H (|GL| / |Aut H|) · covers(H) (covers counted from the hole set alone), Σ_P |GL| / |Aut P|, and
  W4 · |GL| / (210 · 20160) with W4 = count4's plain count.

## 3. Method

**Point branching (`search`)**, as eight-check §3 with 7 members and |D| = 150 (43 lines of D through each hole,
13,888 forms on each line): grow the span one form at a time, branching at the first uncovered point (points off D
first, then D, ascending for `asc`, descending for `desc`) over every legal form whose radical contains it; every
element of the span must be nondegenerate, rank 4 with a member as radical, or rank 6 with a line of D as radical, and
radicals must be pairwise disjoint.

**New here: cores.** A core is the first node on a path at which every off-D point is covered: the span of the 7 rank-4
forms of everything below it. Any legal span spanned by one rank-4 form per member is reached as a core by exactly one
path in every point order that puts the off-D points first (the form covering the first uncovered point of a member is
forced to be that member's form in the span), so the set of cores does not depend on the order. Each run logs the
number of cores by dimension and a fingerprint (sum of 64-bit FNV hashes of their reduced bases), and asc and desc
write the cores themselves (`asc-r<k>.cores`, `desc-r<k>.cores`): they must agree configuration by configuration,
which checks the whole off-D phase across two different trees.

**New here: the coset sweep (`search ... dsweep`)**, a second D phase from every core, with its own primitives (a
Pfaffian test for nondegeneracy, a kernel by elimination, its own coverage bookkeeping, its own legality check; it
re-checks each core with them first) and **no line forms**: below dimension 7 it enumerates all 2^21 - 1 nonzero
forms of K_a = {B : B a = 0} for the smallest uncovered point a and keeps g when S + <g> is legal; at dimension 7 it
tries every one of the 2^21 - 1 nonzero coset representatives of S. Completeness: an 8-dimensional space of forms
covers every point at least once (B ↦ B(a, ·) maps 8 dimensions into a^⊥, 7 dimensions), so a legal T of dimension 8
(radicals disjoint) is exactly a solution; in any solution T ⊇ S, the one element g of T with a in its radical gives
a legal S + <g> ⊆ T, and K_a ∩ S = 0 makes the extensions distinct. Since a is the point `asc` branches on, the sweep
visits exactly the spans the point branching visits, from a different generator (K_a by linear algebra, not
pullbacks along lines): **asc's nodes by dimension must equal the sweep run's own nodes (the off-D phase, cores
included) plus the sweep's legal spans, dimension by dimension, on every configuration**.

**Parallel structure.** The cores are collected during the off-D phase and their D phases run as a second parallel
stage, one task per (core, line through its first uncovered point); with one task per subtree, a configuration with
a few heavy cores ran on one thread. `search ... phase1` runs the off-D phase alone (the core census below).

**Numpy re-search** (`seven_check_numpy.py`, eight-check's second implementation adapted to 7 members, plus cores and
the same fingerprint computed in Python): same tree as asc, so nodes, near-misses, tries in D, cores and fingerprint
must agree exactly on every configuration it is given; planted controls of its own.

## 4. Classification results

**Part A** (job 20261031-054609-check7-parta, 9 min): all 37 slices ended in `END`; **15,288,705 reduced 4-cliques**,
**none maximal** (every 7-partial spread through X, Y, Z has a solid disjoint from all its members); weighted
4-cliques 14,122,766,336.

**`count4`** (`count4.out`): (a) 5,824 / 4,583,936 / 606,183,424 / **14,122,766,336** / **58,471,168,000**
k-cliques for k = 1..5; (b) by edges 606,183,424 triangles and **14,122,766,336** 4-cliques. Part A's weighted
count equals both, and the 5-clique count equals eight-check's (so eight-check's count5 is reproduced too).

**`merge`** (`merge.out`): **all identities hold.**
- Σ_K |Q| = 292,355,840,000 = 5 × 58,471,168,000 (eight-check's 9,372 |Aut|); |Q| ranges over 6 .. 172.
- **3,731 orbits**, no clique-level inconsistency; per-orbit identity holds for every one; Σ = the weighted count.
- Largest clique in Q: 2 / 3 / 4 / 5 / 6 / 10 for 83 / 1,277 / 1,375 / 305 / 332 / 359 orbits, so **spr 359, ext13 332,
  ext12 305, ext11 1,375, ext10 1,277, ext9 83, ext8 0, max 0**; weighted 785,729,336 / 1,134,163,800 / 1,167,415,200 /
  5,457,160,800 / 5,253,897,600 / 324,399,600.
- Type agrees with the 8-orbits reached; every P + q is one of eight-check's 9,372.
- Σ_P (210 · 20160 / |Aut P|) · #spreads = 1,129,992,864 = C(14,4) · 1,128,864.
- Per-8-orbit link identity holds for **all 9,372**.

**`configs`** (`configs.out`; `configs.txt` in type order spr 0..358, ext13 359..690, ext12 691..995, ext11 996..2370,
ext10 2371..3647, ext9 3648..3730): canonical form recomputed from a random GL(8,2) image of every orbit, 0 mismatches;
**hole sets up to GL(8,2) 3,726, none without a solid: Kurz's Table 9.3 for r = 7**; exact covers of every hole set's
complement (3,757 covers): every cover on a listed orbit with that hole set, Σ |Aut H| / |Aut P| = covers for all 3,726,
all 3,731 orbits reached (3,703 hole sets with one cover, 16 with two covers of one orbit, 5 with two orbits, one with 4
covers and one with 8 of a single orbit). Outside families, every configuration on a listed orbit:
- **eight-check's 9,372 × 8 = 74,976: all 3,731 orbits reached** (spr 359, ext13 332, ext12 305, ext11 1,375, ext10
  1,277, ext9 83);
- nine-check's 5,663 × 36 = 203,868: all 3,731;
- fourteen-check's 73 × 3,432 = 250,536: all 359 spr.

Stabiliser orders of the 3,731 configurations: 1: 3,096, 2: 391, 3: 39, 4: 55, 6: 72, ..., one each of 600 and 5,040.

**Exact masses** (`python apn/seven_check_mass.py runs/check7`): 3,726 classes over 27 stabiliser orders, the number of
classes per |Aut H| identical to Langevin-Hou's for all 27, and Σ_H |GL(8,2)| / |Aut H| = 17,833,337,132,662,061,531,136
= their Total. Labelled 7-partial spreads: from the covers 17,840,479,724,727,007,444,992, from the 3,731 configurations
the same, from W4 (· |GL| / (210 · 20160)) the same, with W4 = count4's plain count.

So the configurations are complete three ways that don't share a step: Part A with its identities (closing against
count4's two plain counts and eight-check's classification orbit by orbit); the deletions of eight-check's 9,372
(every 7 extends, so this reaches every class); and the covers of the 3,726 hole sets (corroborated by Langevin-Hou),
whose labelled total equals W4 · |GL| / (210 · 20160) from the plain clique count.

## 5. Search results

Every configuration's search ends at its cores' D phases, so the census of cores (from `search phase1`, the off-D
phase alone, job 20261031-121312-check7-phase1; the same in every run) says where the work is: **114,702 cores**, by
dimension 4 / 5 / 6 / 7: **1 / 1,634 / 108,026 / 5,041** (none of dimension 3 or less, none of dimension 8); 944
configurations have no core at all (no legal span holds all 7 rank-4 forms, so nothing to do in D). spr carries 1,276 of
the dimension-5 cores and 4,958 of the dimension-7 ones; the single dimension-4 core is configuration 1527 (ext11),
whose D tree is the largest of all (591,193 nodes at dimension 6, 8.3G tries in D, 72 min). Cost is
dominated by the low-dimensional cores: 157 configurations have a core of dimension <= 5.

**asc** (jobs 20261031-060729-check7-asc and, from configuration 1527 on, 20261031-141646-check7-asc-2 with the
candidate cache; `asc-r<k>.prog`, `asc-r<k>.cores`), every configuration once, no duplicates
(`seven_check_tally.py prog runs/check7/configs.txt 3731 runs/check7/asc-r*.prog`):

| order | nodes by dimension 0..8 | p1 nodes by dim 4..7 | tries in D (dim 4 / 5 / 6 / 7) | solutions |
|---|---|---|---|---|
| asc | 3,731 / 104,468 / 2,925,104 / 52,895,372 / 312,847,548 / 195,494,917 / 2,767,276 / 11,780 / **0** | 1 / 16,256 / 843,914 / 11,780 | 7.2M / 1.27G / 97.18G / 1.89G (100.3G) | **0** |
| desc | 3,731 / 104,468 / 2,925,104 / 52,669,145 / 310,556,440 / 193,373,519 / 2,772,459 / 11,924 / **0** | 1 / 15,125 / 858,354 / 11,924 | 6.6M / 1.29G / 97.01G / 1.93G (100.2G) | **0** |

No node ever reached dimension 8 (in either order), no bent completion ran (no span below dimension 8 ever covered
every point), no full-span record, no SOLUTION line. By type (asc): spr (359) holds 74.8G of the 100.3G tries in D and
11,602 of the 11,780 dimension-7 nodes; ext13 5.3G, ext12 2.5G, ext11 14.4G (8.3G of it configuration 1527), ext10
3.1G, ext9 0.12G. Dimension-7 nodes occur on 155 configurations (spr 122, ext13 11, ext12 3, ext11 14, ext10 4, ext9 1).

**Cores across runs** (`seven_check_tally.py orders runs/check7`): cores by dimension and fingerprint equal to the
phase-1 census on **all 3,731 configurations in asc and in desc** (and on every configuration the sweep has done), and
the asc and desc core files equal as sets on all 3,731 (**114,702 distinct cores**, file counts = R-line counts).
Since the off-D phase is a different tree in the two orders (desc reaches its cores through different intermediate
nodes: 52,669,145 nodes at dimension 3 against asc's 52,895,372), this checks the whole off-D phase twice over.

**The coset sweep** (job 20261031-062222-check7-sweep, from configuration 135 on 20261031-141647-check7-sweep-2 with
its cache; `sweep-r<k>.prog`), **complete on all 3,731 configurations** (ended 2026-10-31 19:23, exit 0): every
configuration once, no duplicates; **asc's nodes by dimension = the sweep's phase-1 nodes + its legal spans, dimension
by dimension, 0 disagreements over 3,731**, cores by dimension and fingerprint equal to the census on all 3,731, and
**0 solutions**: legal spans by dimension 5 / 6 / 7 = 14,622 / 735,888 / 6,739 (none at 8), 571.0G checks (spr 426.4G),
no bent extension, no full span. The sweep, generating candidates as K_a by linear algebra and checking them with its
own primitives, visits exactly the spans the point branching visits.

**The numpy re-search** (`seven_check_numpy.py`, job 20261031-063551-check7-numpy, `numpy-<i>.out`; sample in
`numpy-sample.txt`: the first configuration of each type (0, 359, 691, 996, 2371, 3648), configurations 3, 9, 41, 47, 49
(the last four with dimension-5 or dimension-7 cores), and 7 drawn at random with seed 7 (1326, 617, 1617, 2666, 197,
296, 3363)): **all 18 agree with asc in every count** (nodes, p1 nodes and tries in D by dimension, cores by dimension,
and the core fingerprint, which the Python computes from its own reduced bases), 2.04G tries in D over the sample,
0 solutions.

**The candidate caches** (added at 14:16 for configuration 1527, commit f25ca63): for a core of dimension <= 5, the
candidates through a point a that are legal with the core itself are computed once and reused by every descendant that
branches at a (in the sweep: the elements of K_a legal with the core). They change the tries, not the spans visited;
checked on configuration 41 (17 dimension-5 cores) against its uncached records: nodes by dimension identical in asc
(tries in D 255M -> 42M) and the sweep's legal spans identical (checks 2.1G -> 0.30G). In asc and desc they apply from
configuration 1527 on, in the sweep from 135 on; the numpy tally only compares tries where neither side used them.

Cross-run that doesn't use this classification (`crosssample 10 7`, job 20261031-171708-check7-s8sample,
`s8-b<start>.prog`): 10 random blocks of 24 consecutive indices of eight-check's 9,372 × 8 deletions (30 of its
configurations with each of their 8 members removed: 240 inputs taken straight from eight-check's file, no
canonicalisation), searched as they are: nodes 240 / 6,720 / 188,160 / 3,401,106 / 20,184,575 / 12,679,118 / 125,400 /
229 / 0, 5,816 cores, 2.27G tries in D, **0 solutions**.

`selftest` (`selftest.out`, job 20261031-060729-check7-selftest, **0 failures**): ND4 28, ND6 13,888, |GL(4,2)| = 20160,
|V5824| = 5824; for all 3,731 configurations |D| = 150, 2,150 lines in D and every hole on exactly 43 of them, pairwise
complementarity, member forms distinct with radical exactly the member; line forms distinct, rank 6, radical exactly
the line (and the sweep's kernel giving the same radical) on every 128th configuration; brute force over all 2^28
forms: exactly 28 with radical R_5 of the last configuration, the sweep's Pfaffian test equal to rank8 on every 16th
form, and the sweep's kernel and eight-check's Gray-walk radical equal to the search's radical on every 256th;
canonical form and |Aut| reproduced from a random GL(8,2) image of every configuration, with a frame carrying it back.


## 6. Controls

All from `seven_check_ctl.exe` (job 20261031-065328-check7-controls2, `ctl-<name>.out` / `.txt`) except ind-bl22
(from `seven_check.exe`, job 20261031-062245-check7-controls, cancelled after this batch to change the design). As at
eight, a control carries a known quadratic APN function T with k >= 2 rank-4 radicals onto a configuration (extendable
path: a spread through T's radicals, 10 other spread elements as holes, the 7-set carried onto its representative by
its canonical frame, which also tests the classification; other path: T's radicals embedded into a non-spr
representative), T's rank-6 forms off D allowed as an oracle, and the search must return T's span. New here:
- **guided** controls run the whole off-D phase but give only the cores inside T a D phase. In these controls T's
  off-D points need T's own oracle forms, so the core is usually T itself [Rowen, 2026-11-01, from the review: in every guided batch it is T itself, reached whole in the off-D phase], and a full D phase over every other core
  would cost a whole configuration's search per control (the first, unguided batch was on course for days on the spr
  representatives). The D phase is tested by the in-D controls instead:
- **in-D, one step** (`nond`, eight-check's): start from a hyperplane W of T without g, a rank-6 form of T with radical
  in D, g's points first in the order; **in-D, two steps** (`nond2`, new): start from a codimension-2 subspace W of T
  avoiding two such forms g_1, g_2 (W + <g_1> misses g_2 and vice versa), so the search must find both, through two
  D-levels; both also with the sweep as D phase (`dsweep`; W then also holds T's rank-4 forms).

| batch | ext path | other path | notes |
|---|---|---|---|
| ind-bl22 (BL22, 480 with k >= 2, 3 choices; in-D one step) | **1,440 / 1,440** | **1,397 / 1,397** | 43 k = 4 skipped on the other path (no member of the 4th radical's class), as at eight |
| ind2-bl22 (in-D two steps) | **1,440 / 1,440** | **1,397 / 1,397** | |
| bl22 (guided) | **1,440 / 1,440** (924 / 438 / 78 with k = 2 / 3 / 4) | **1,397 / 1,397** | reps hit 305 / 1,389 |
| blp22 (37 with k = 2, 6 choices) | **222 / 222** | **222 / 222** | |
| Yu–Perrin QAM (6 choices) | **24 / 24** | **24 / 24** | |
| Yu–Wang–Li + Edel–Pott (6 choices) | **72 / 72** | **72 / 72** | 6 + 6 setup failures: line 3, APN but not quadratic (out of scope, as at eight to eleven) |
| new (the 3.8M's 7,199 with k >= 2, 1 choice) | **7,199 / 7,199** | **7,199 / 7,199** | |
| ind2-new (in-D two steps) | **7,199 / 7,199** | **7,199 / 7,199** | |
| swind-bl22 (sweep, in-D one step, 1 choice) | **479 / 479** | **468 / 468** | 2 skipped: g in the span of T's rank-4 forms (no W can hold those and miss g); on the other path also the 11 k = 4 skipped as in ind-bl22 [Rowen, 2026-11-01, from the bundle helper] |
| swind2-bl22 (sweep, in-D two steps) | **471 / 471** | **459 / 459** | 19 skipped likewise, and the same 11 k = 4 on the other path |
| full-bl22 (unguided: every core of the configuration gets its D phase, 20 functions, 1 choice) | **20 / 20** | **20 / 20** | |


0 NOT found in every batch, 0 other spaces found. **Every configuration reached:** the batches carried a function onto
3,553 of the 3,731 representatives (`seven_check_tally.py reps runs/check7`, `reps-unreached.txt`); the other 178 (spr
16, ext13 25, ext12 15, ext11 61, ext10 59, ext9 2) got **`ctlrot 40 guided`**: one control each, BL22's first 40
functions with k >= 3 in rotation embedded into the representative's members (`ctlrot-0.out`): **178 / 178
re-found**, 0 other spaces, all 40 functions used. **numpy planted controls** (job 20261031-065338-check7-npctl; no
classification involved: members = T's radicals + 7 - k random solids): BL22's first 20 with k >= 2: **20 / 20
re-found**; the 3.8M's first 10: **10 / 10**.

## Verdict (2026-10-31 17:55, written by a subagent, Claude Opus 5.5, for Rowen)

**There is no quadratic APN function on F_2^8 with Walsh distribution [0^198, 2^50, 4^7]** (7 rank-4 components).
0 solutions in the complete point-branching search over all 3,731 classes of 7-solid partial spreads, in two point
orders (asc and desc), with no node ever reaching dimension 8 (deepest: 11,780 / 11,924 nodes at dimension 7), and 0 in
a cross-run on 240 raw deletions of eight-check's configurations.

The classification: 3,731 orbits from Part A, closed by every identity in `merge` (Σ |Q| against eight-check's 9,372,
the per-orbit counts, the spreads through each orbit, and the link to all 9,372 eight-check orbits one by one); Part A's
weighted count = the 4-clique count of V5824 by two plain methods (enumeration, and by edges); Kurz's 3,726 hole sets
with none maximal, all 3,757 exact covers accounted for, every orbit reached; Langevin-Hou's per-|Aut H| counts matched
for all 27 orders, and the labelled 7-partial spreads equal three ways (from the covers, the configurations and the
plain clique count). New here, the deletions of eight-check's 9,372 reach all 3,731 (every 7 extends, so this family
alone is complete given eight-check's classification).

The search: the off-D phase reaches the same 114,702 cores (by dimension and fingerprint, and as sets of spans) in asc,
desc and the off-D-only census, on every configuration; the numpy implementation agrees count for count (nodes, tries,
cores, fingerprint) on 18 configurations; the coset sweep, a D phase with its own primitives and K_a instead of line
forms, visits exactly asc's spans on every configuration it has done (State has the count) and finds 0 solutions [Rowen, 2026-11-01: it completed all 3,731 at 19:23, after this verdict, with the same result; §5]; the
controls re-find every known function tried (BL22 1,440 + 1,397 on each of three batch kinds, BLP22 222 + 222,
Yu–Perrin 24 + 24, Yu–Wang–Li 72 + 72 with its non-quadratic line 3 the only setup failures, the 3.8M's 7,199 + 7,199
twice), including the in-D controls through the point branching (BL22 one step and two steps, 1,440 + 1,397 each; the 3.8M
two steps, 7,199 + 7,199) and through the sweep (BL22 one step 479 + 468, two steps 471 + 459), 20 + 20 unguided full searches, `ctlrot` on the 178 configurations the
batches missed (so every one of the 3,731 had a known function carried onto it and re-found), and 30 / 30 planted
controls through numpy.

On setup: the brief's setup is right (§1); the case is harder than eight mainly through one configuration (1527, the
only one whose 7 rank-4 forms span a 4-space) and the 1,634 dimension-5 cores, which the candidate cache (§5) made cheap.

**On independence**: the code, the classification and every number above come from this repo's `seven_check*` files
and `runs/check7/`, written and run by this subagent, which read none of `pspread*`, `cover*`, `cover_prep.py`,
`cover_chunks.cmd`, `cover_tally.py`, `holes_tally.py`, README's "Eight" and "Seven" sections, `runs/i9/`, `runs/i10/`,
`notes/` or `board.json`, and ran no `git log` but `git log -- <own file>` (to note a commit). Search A's jobs were
visible in the share runner's status while this ran (names `cover7`, `holes7`, `p3-k7`, and that `cover7` was still
running at 06:00); nothing of their output was read. Inputs from earlier checks: eight-check's and nine-check's
`configs.txt` and fourteen-check's `spreads.txt` (B's lineage), Kurz's table and Langevin-Hou's `fix-7.txt`, as the
brief names; `eight_check.cpp` and `eight_check_numpy.py` as the starting point, as the brief invites. Comparison with
search A is left to Rowen (a §7 after this verdict).

## 7. Comparison with search A (read after the verdict above was signed and committed; Rowen, 2026-11-01)

README "Seven" and `runs/i10/` (`p7.txt`; search A's nine outputs `cover/c0..c3.out`, `cover-b/c5..c7.out`,
`cover-c/c4x.out`, `split/p2222_s5.out`); `apn/seven_check_vs_a.py` > `runs/check7/vs-a.out`:
- **Classification.** A's `pspread3h enum 7` lists 3,731 configurations. Mapped directly (`seven_check mapa
  runs/check7/configs.txt runs/i10/p7.txt` > `runs/check7/mapa.out`; the mode added after the verdict, built as
  `seven_check_map.exe` so the binaries named in State stay as they ran; A's lines checked against the mode's reading:
  name, member count 7, 4 basis bytes per solid): **a bijection** onto the 3,731 here (0 on none of ours, 0 hit twice).
  A's |Aut| (in its names), this code's |Aut| of A's representative and this classification's |Aut| are equal on all
  3,731. A's extension count w and the type here correspond one to one (w 2 / 3 / 4 / 5 / 6 / 10 = ext9 83 / ext10
  1,277 / ext11 1,375 / ext12 305 / ext13 332 / spr 359). As at eight the two classifications share the method's shape,
  so the independent checks are §4's (Part A's clique counts, Kurz, Langevin-Hou) and A's own hole-set match.
- **Search.** A (`cover3` on c0..c3, c5..c7; `cover4`, cover3's search with an optional split dimension, on
  c4 minus p7-2222 unsplit and on p7-2222 split at dimension 5): **0 solutions**, 104,468 tasks, 568,755,235 nodes, by dimension 1..8
  104,468 / 2,925,104 / 52,909,076 / 313,666,564 / 196,375,985 / 2,762,176 / 11,862 / 0; no full span, no dimension-8
  non-APN, 0 bad seeds. Here (asc, dimension 1..8) 104,468 / 2,925,104 / 52,895,372 / 312,847,548 / 195,494,917 /
  2,767,276 / 11,780 / 0: equal through dimension 2, then apart (different representatives and branching), as at
  eight and nine.
- **Near misses = cores, exactly.** A logs NEAR when a span first holds all 7 rank-4 forms with holes uncovered; the
  cores here are the spans the off-D phase ends at, the same event, so no containment correction is needed (at eight
  ours logged every span holding all 8). A: 114,702 NEAR lines, no (configuration, basis) repeated, by dimension 4 / 5 /
  6 / 7 = 1 / 1,634 / 108,026 / 5,041, on 2,787 configurations. Through the bijection, **A's count equals our core
  count for every (configuration, dimension) pair, in asc and in desc (3,063 nonzero pairs, 0 unequal)**: the same
  156 configurations with dimension-5 cores, 2,785 with dimension 6, 121 with dimension 7, and the single dimension-4
  core on A's p7-2222-a168-w4 = our configuration 1527.
- Same answer: no function, by two searches and two classifications that agree configuration for configuration,
  near misses included.

## State (2026-10-31 21:00, Rowen; the checker's earlier state is in git)

Everything is done: the coset sweep (job 20261031-141647-check7-sweep-2) finished all 3,731 configurations at 19:23
and its tallies are in §5 (`orders`: 0 disagreements in asc, desc and sweep, 1 truncated core record skipped; `prog` on
the sweep: full 0, sweep solutions 0). Nothing of B is still running. Comparison with search A: §7 (2026-11-01, done). [Rowen, 2026-11-01, from the bundle helper's check:] asc configurations 0..10 come from the first asc job (20261031-055731, before commit 2d5e21b; `asc-r0.out` says 11 already done), and count4, parta, merge and configs started at 05:46, three minutes before 94d7632 was committed; the bundle's re-runs from the committed sources reproduce both records exactly.
`runs/check7/` (about 95 MB, mostly `*.cores`) can go once the bundle is made.

Rerunning: every stage reruns or resumes through `python apn\seven_check_slices.py <stage>`, submitted from PowerShell with
`python tools\share\runner.py submit --name <n> --kind cpu --mem 4000 -- cmd /c "python apn\seven_check_slices.py <stage>"`
(stages: count4, parta, merge, configs, selftest, search [asc|desc|sweep|rand7], cross s8 [ranges], controls [names],
ctlrot [pool] [dsweep] [guided] [all], numpy [nproc] [ids]); finished slices and ranges are skipped. Binaries, all from
`seven_check.cpp`: `seven_check.exe` (asc and desc up to configuration 1526: commit 2d5e21b, from 1527 on f25ca63),
`seven_check_sw.exe` (sweep: 67a868a up to 134, f25ca63 from 135), `seven_check_ctl.exe` (e302003),
`seven_check_p1.exe` (7ab9fc0). Tallies: `python apn/seven_check_tally.py prog runs/check7/configs.txt 3731
runs/check7/<order>-r*.prog`, `... orders runs/check7`, `... numpy runs/check7`, `... reps runs/check7`,
`python apn/seven_check_mass.py runs/check7`.

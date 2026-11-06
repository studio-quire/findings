> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged, but for "other work" where State named the job holding the box). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `six_check.cpp`, `six_check_slices.py`,
> `six_check_tally.py`, `six_check_mass.py`, `apn/six_check_vs_a.py` → `src/` (this folder's `six_check_mass.py` and
> `six_check_vs_a.py` read the files here, and `six_check_tally.py` takes its files as arguments: `README.md`, "Six", gives
> the commands); `runs/check6/orbits.txt` → `data/orbits6.txt`; `runs/check6/configs.txt` → `data/configs6.txt`; Part A's
> `parta-c<c>-s0.txt` → `data/parta6-c<c>-s0.txt` (its `maxc-c<c>-s0.txt` are all empty); the other outputs in `runs/check6/`
> → `logs/check6-<name>-output.txt` (`<name>.out`) and, for the controls, `logs/check6-ctl-*.txt`,
> `logs/check6-ctlplant-*.txt`, `logs/check6-plants-*.txt`; the 11 `ctlrot-lin-<k>.out` → `logs/check6-ctlrot-lin-output.txt`
> (concatenated); the searches' `{p1asc,p1desc,asc,lin,lindesc}.out` → `logs/check6-*-output.txt` (the header once and every
> `config` line, and for `lin` the last totals: one process ran per configuration, each repeating the rest); their
> `.prog` files → `logs/check6-{p1asc,p1desc,asc}-R.txt` (the `R` lines, the per-configuration record; the `P1` lines, a
> line per core, aren't included) and `logs/check6-{lin,lindesc}-prog.txt` (whole); the core files `p1asc.cores`,
> `p1desc.cores`, `lin.cores`, `lindesc.cores` (457 to 973 MB each) aren't included (`README.md`, "Left out of `logs/`"; for
> `p1asc` and `p1desc`, `logs/six-near.txt` has a digest of each configuration's set of spans); `tmp/asc340.*`, `lin340.*`
> and `lin340b.*` → `logs/check6-asc340-*`, `logs/check6-lin340-*`, `logs/check6-lin340b-*`; the share runner's job logs
> (`check6-class` and so on) → `logs/2026…-check6-*.log`, whose commands are `six_check_slices.py`'s (which reads and writes
> `runs/check6/` and the folders below); seven-check's, eight-check's and fourteen-check's lists (`runs/check7/configs.txt`
> and so on) are `data/configs7.txt`, `data/configs8.txt` and `data/spreads.txt`; `runs/papers/langevin/fix-6.txt` is
> Langevin and Hou's data file (downloaded, not included: `README.md`); `runs/known/` holds the known-function lists (not
> included: `README.md`); `runs/bin/six_check*.exe` are builds of the one source at the commits State names (commits of
> Rowen's own repository, which isn't public; `README.md`, Build, says what changed between them, and corrects State on
> `six_check_t.exe`). Rowen's side (§7): `runs/i11/p6.txt` → `data/p6.txt`; search A's final outputs
> `runs/i11/cover5/c<k>.out` → `logs/six-output-c<k>.txt` (without their `NEAR` lines: `logs/six-near.txt` has them per
> configuration); `runs/check6/mapa.out` and `vs-a.out` → `logs/check6-mapa-output.txt`, `logs/check6-vs-a-output.txt`.
> README's "Seven" and "Six" sections are Rowen's working notes (not published), summarised in `README.md` here and in the
> write-up; `notes/`, `board.json` and `review-seven.md` (a review of seven's write-up, which §6 cites) are Rowen's own
> (not published). The known-function lists aren't included (`README.md`).

# Independent check, 6 rank-4 components: is there a quadratic APN F on F_2^8 with distribution [0^194, 2^55, 4^6]?

Subagent (Claude Opus 5.5), 2026-11-03, for Rowen ("search B" for 6 rank-4 components). Written without reading
`pspread*.cpp`, `cover*.cpp`, `cover_prep.py`, `cover_chunks*.cmd`, `cover_tally.py`, `holes_tally.py`, README's "Seven"
and "Six" sections (or anything about i = 10, 11 there), `runs/i10/`, `runs/i11/`, `notes/`, `board.json`, other jobs'
runner logs or the git log. Files: `six_check.cpp` (all C++ parts), `six_check_slices.py` (resumable driver),
`six_check_tally.py` (tallies and comparisons, own Python code), `six_check_mass.py` (exact masses against
Langevin-Hou), outputs in `runs/check6/`. Built on `seven_check.cpp`.

**Status: verdict in** (no function; 2026-11-04 09:10). Extra checks still running at the verdict: State.

## 1. The reduction

As at seven (`sixteen-check.md` §1): F is quadratic APN iff its 8-dimensional space T of alternating forms has every
a != 0 in the radical of exactly one nonzero form. A 6-dimensional radical meets every solid (6 + 4 > 8), so with
rank-4 forms present the radicals have dimension 0, 2 or 4. 6 · 15 + 3m = 255 gives m = 55 rank-6 forms and
255 - 61 = 194 bent: **[0^194, 2^55, 4^6]**, as the brief says. The 6 rank-4 radicals are pairwise complementary solids
(a partial spread), D = the 255 - 90 = 165 points they miss, partitioned by the 55 rank-6 radicals (lines). A point of a
radical R_j is covered only by the rank-4 form with radical exactly R_j, so all 6 rank-4 forms are forced by coverage.

A hyperplane holds at most one R_j, so it meets D in 127 - (42 + 8t) = 85 or 77 points (t = 0 on 165 hyperplanes,
t = 1 on 90): D is a projective two-weight set, Cay(F_2^8, D) is strongly regular with eigenvalues 165, 5, -11,
λ = 165 + 5 - 11 - 55 = 104, so every hole lies on 52 lines inside D and there are 165 · 52 / 3 = **2,860 lines** in
every candidate (the selftest checks both on all configurations).

Every 6-partial spread of PG(7,2) extends (Kurz, Table 9.3, r = 6: 341 point sets, none maximal; confirmed here, §4:
every 6 extends to at least 10 solids).

**Normalisation.** Three of the 6 solids go to Y, X, Z = graph(I); the other 3 are graph(N), N in V5824, pairwise
N + N' invertible: a triangle. Frames: 6 · 5 · 4 = 120 ordered member triples × diag(h, h), so an orbit P carries
120 · 20160 / |Aut P| triangles.

## 2. Classification

As seven-check §2, one solid fewer: Part A (every triangle of V5824 up to conjugation by the minimal-class trick,
canonicalised), `merge` with its identities, `count3` (k-cliques of V5824 by plain enumeration and by edges), `configs`
(hole sets against Kurz, exact covers, three outside families), `six_check_mass.py` against Langevin-Hou's `fix-6.txt`.
Identities:
- Σ_K |Q(K)| over all triangles = 4 × the weighted 4-cliques (seven-check's 3,731 orbits: Σ 210 · 20160 / |Aut|);
- per orbit: weighted triangles = 120 · 20160 / |Aut P|;
- Σ_P (120 · 20160 / |Aut P|) · #spreads through P = C(14,3) · 1,128,864;
- per 7-orbit O (all 3,731): Σ_P (120 · 20160 / |Aut P|) · #{q : P + q ≅ O} = 4 · 210 · 20160 / |Aut O|.

Outside families: **seven-check's 3,731 × 7 deletions** (every 6 extends, so this family alone reaches every class,
given seven's classification), eight-check's 9,372 × C(8,2), fourteen-check's 73 spread reps × C(14,8).

## 4. Classification results

(§3, the method, is below the results of the classification, since it was designed after them.)

- **Part A** (job check6-class, under a minute): 588,620 reduced triangles, none maximal; weighted triangles
  606,183,424.
- **`count3`**: plain enumeration 5,824 / 4,583,936 / **606,183,424** / 14,122,766,336 k-cliques for k = 1..4; by
  edges 606,183,424 triangles and 14,122,766,336 4-cliques. Part A's weighted count equals both; the 4-clique count is
  seven-check's.
- **`merge`: all identities hold.** Σ_K |Q| = 56,491,065,344 = 4 × 14,122,766,336; |Q| ranges over 72 .. 308.
  **341 orbits**, per-orbit identity for every one, Σ = 606,183,424. Largest clique in Q: 4 / 5 / 6 / 7 / 11 for
  1 / 68 / 48 / 92 / 132 orbits: **spr 132, ext13 92, ext12 48, ext11 68, ext10 1** (none extending at most to 9, 8 or
  7 solids, none maximal). Spreads: Σ = 410,906,496 = C(14,3) · 1,128,864. Per-7-orbit link identity holds for **all
  3,731**; every P + q is one of seven-check's orbits, and the type from the clique size agrees with the best type among
  them.
- **`configs`**: canonical form recomputed from a random GL(8,2) image of every orbit, 0 mismatches; **hole sets up to
  GL(8,2): 341, none without a solid: Kurz's Table 9.3 for r = 6**; exact covers of every hole set's complement by 6
  solids: 344 covers (338 hole sets with one, 3 with two covers of one orbit), every cover on a listed orbit with that
  hole set, Σ |Aut H| / |Aut P| = covers for all 341, all 341 orbits reached. Outside families, every configuration on a
  listed orbit: **seven-check's 3,731 × 7 = 26,117: all 341 orbits reached**; eight-check's 9,372 × 28 = 262,416: all
  341; fourteen-check's 73 × 3,003 = 219,219: all 132 spr.
- **Exact masses** (`python apn/six_check_mass.py runs/check6`): 341 classes over 17 stabiliser orders, the number of
  classes per |Aut H| identical to Langevin-Hou's for all 17, Σ_H |GL(8,2)| / |Aut H| =
  1,339,989,812,392,369,324,032 = their Total. Labelled 6-partial spreads: from the covers, from the 341
  configurations and from W3 · |GL| / (120 · 20160) (W3 = count3's plain count) all
  1,340,074,242,473,170,305,024.
- `selftest` (0 failures): |D| = 165, 2,860 lines in D, every hole on exactly 52, on all 341 configurations; the
  member forms, line forms and canonical forms as at seven.

## 3. Method

**Two searches over the same tree, with no shared primitive.**

*Point branching* (`search`, seven-check's code with 6 members): grow the span one form at a time, branching at the
first uncovered point (points off D first, ascending for `asc`, descending for `desc`) over every legal form whose
radical contains it: off D the 28 forms with the member as radical (pullbacks with nondegenerate codes), in D the 13,888
forms on each line of D through the point that misses everything covered; every element of the span must be bent,
rank 4 with a member as radical, or rank 6 with a line of D as radical, radicals pairwise disjoint. A **core** is the
first node on a path with every off-D point covered: the span of one rank-4 form per member, the same set in every
point order. `search ... phase1` runs the off-D phase alone and writes the cores.

*The linear search* (`lin`, new here, sharing with the code above only the input masks and the u64 layout of a form):
everything is phrased through the vectors B a (the derivative direction a: XOR of the rows i with a_i = 1).
- A span S = <b_1..b_e> covers a point a once when the vectors b_i a have rank e - 1 (the kernel of B -> B a on S is
  {B in S : a in rad B}); S is legal iff no point is covered twice, i.e. the rank is >= e - 1 at every a; per point a
  fully reduced echelon basis of S a is kept, with the combination giving the kernel element.
- S + <g> is legal iff g a is not in S a at every point a that S covers (at an uncovered point the kernel has
  dimension <= 1 anyway). An 8-dimensional legal span is an APN space (the kernel at each a has dimension >= 1).
- Off D: branch at the smallest uncovered point off D over the member's 28 forms, built from four functionals with
  kernel R_j and the 4x4 Pfaffian c12 c34 + c13 c24 + c14 c23; keep a child when nothing is covered twice and every
  newly covered point off D is covered by a form of its member (or an oracle form, in controls).
- In D: branch at the smallest uncovered point a0 over **every** nonzero g in K_a0 = {B : B a0 = 0} (21-dimensional,
  pair forms of the 7 functionals vanishing at a0), by Gray code; the test "g a not in S a" is one byte lane per
  covered point (the reduction of g a modulo S a, linear in g), legal iff no lane is zero (AVX2, the first 32 lanes
  by Gray code, the rest on demand). Completeness: in a solution T containing S exactly one element of T has a0 in its
  radical, and S plus it is legal; K_a0 meets S in 0, so the extensions are distinct. If every point is covered below
  dimension 8 (possible only at dimension 6 or 7), every coset of S is tried. For cores of dimension <= 5 the legal g
  in K_a are listed once per point and reused below (a descendant's are among them).
- Candidates are not line forms but all of K_a0, and the test is linear algebra, not radicals; but the spans accepted
  are the point branching's, with one exception: a rank-4 element whose radical is a solid inside D (a seventh
  solid), which the point branching rejects and `lin` keeps and counts (x4) and follows. So **lin's nodes by
  dimension less its x4 nodes must equal asc's nodes by dimension, configuration by configuration, off D and in D
  alike, and its cores less those with a seventh solid must equal the point branching's cores as sets**
  (`six_check_tally.py`, its own Python linear algebra and rank).

`linselftest`: lin's primitives against brute force on 16 configurations (member forms against `radical`, K_a
independent and killing a, `L_legal_add` against direct radical counting on random forms, the lane scan against
`L_legal_add`, every candidate the full scan of K_a0 accepts legal by brute force, on 200 cores per configuration).

## 5. Search results

**Cost, measured first** (2026-11-03, on the shared E-cores): the off-D phase is a few seconds per configuration;
configuration 340 (ext10, 21,579 cores) took `lin` 112 s for the whole search and the point branching 2,586 s, the
same tree (below) at about 23 times the cost, and about 15 to 20 times over the configurations sampled. A complete
point-branching search would have been 10 days or more on the box; `lin` took 4 h. So, as the brief allows: **`lin` in
full in two orders (asc and desc), the point branching's off-D phase in full in two orders, and the point
branching's whole search on a sample** of 16 configurations: the first of each type with more than one
representative (0, 132, 224, 272; 340 from the timing run), seven at random (18, 41, 74, 133, 248, 293, 339, seed 6),
the four with a dimension-4 core (24, 31, 103, 130) and the one with the most dimension-5 cores (131, 15,220).
[Rowen, 2026-11-04, from the bundle helper's check: with 340 that makes 17 configurations, the 16 of `asc-sample.txt` plus
340 from the timing run; "13 of 16" in State and the verdict was 12 of them plus 340.]

**The cores** (the off-D phase; `search phase1`, jobs check6-p1asc and check6-p1desc, both 341 of 341):
**7,941,871 cores**, by dimension 4 / 5 / 6: **7 / 113,279 / 7,828,585** (none of dimension 3 or less, every
configuration has some; the dimension-4 ones on configurations 130 (3), 103 (2), 24 and 31). asc and desc agree on
cores by dimension and fingerprint on all 341, and as sets (`six_check_tally.py cores`, Python's own reduced bases):
7,941,871 in both, none missing, none extra, though the trees differ (dimension 3: 5,305,886 nodes in asc, 5,292,965
in desc; dimension 5: 80,157,984 and 79,748,390).

**`lin`, asc** (job check6-lin, all 341, `lin.prog`): nodes by dimension 0..8 **341 / 9,548 / 267,344 / 5,308,580 /
45,963,730 / 81,635,576 / 40,014,483 / 340,912 / 0**, of them with a seventh solid (x4) 0 / 0 / 0 / 2,694 / 172,211 /
1,318,748 / 3,134,339 / 173,214 / 0; 20.19T candidates tried (Gray steps; spr 10.68T, ext13 4.33T, ext12 2.15T,
ext11 2.98T, ext10 0.05T); **no node of dimension 8, 0 solutions**. Cores 8 / 125,504 / 8,689,010 at dimensions 4 / 5 /
6, of them x4 1 / 12,225 / 860,425: less these, **exactly the point branching's 7 / 113,279 / 7,828,585, and as sets on
every configuration** (`tally-lin.out`: 0 missing, 872,651 extra, each with more than 6 rank-4 elements by Python's
rank, 0 other extras). Dimension-7 nodes less x4: 167,698 (spr 165,263).

**`lin`, desc** (job check6-lindesc: branching at the largest uncovered point, off D and in D): nodes 341 / 9,548 /
267,344 / 5,296,100 / 45,777,356 / 81,273,352 / 40,362,043 / 352,939 / 0, 20.19T candidates, **0 solutions, no node of
dimension 8**; cores, x4 cores and lin's own fingerprint equal to asc's on all 341.

**Point branching on the sample** (job check6-asc, `asc.prog`, done 2026-11-05 06:21, 30.6 h mostly shared with other
work): all 16 of `asc-sample.txt`, **asc's nodes by dimension equal lin's less its x4 nodes, dimension by dimension,
on every one** (`six_check_tally.py nodes`: 16 in both, 0 differ), the D phase included (e.g. configuration 0: 157,776
nodes at dimension 6 and 70 at 7 in both, 12.7G tries of a line form in D), and 0 solutions. Totals over the 16: nodes by
dimension 16 / 448 / 12,544 / 249,619 / 2,185,056 / 4,174,853 / 18,023,357 / 108,524 / 0; cores 7 / 23,020 / 517,801 at
dimensions 4 / 5 / 6; tries in D 636.5 x 10^9; `bentext 0 full 0`. Configuration 340, the timing run, is the 17th.

## 6. Controls

What each kind exercises (review-seven.md item 1):
- **Planted controls (`ctlplant`, new here; the main positive test).** For a known APN space T with k rank-4 forms,
  members = T's k rank-4 radicals + 6 - k pairwise disjoint solids, each inside the set of points covered by the
  rank-6 elements of a proper subspace W of T that holds the rank-4 forms (`plantfind`: W over the hyperplanes, or
  the codimension-2 subspaces, of T). In T every point off D is then covered by an element of W, so T's forms
  covering the points off D (its rank-4 forms and the oracle forms, its rank-6 forms whose radical meets an extra
  member) span a space of dimension d0 <= dim W: **the off-D phase ends at a core of dimension d0 < 8 inside T, and
  only the D phase can complete T** (8 - d0 steps). The 6-set is carried onto its representative by its canonical
  frame (testing the classification), and **the representative's whole search runs, unguided**: every core of the
  configuration (its own and those the oracle forms add) gets its D phase; T must be found and nothing else.
  d0 = 7 from hyperplanes, d0 = 6 from codimension 2 (k = 3 or 4 needed: 2 or 3 solids inside the points of a
  6-dimensional W; only 5 such plantings exist among the known lists searched, all from BL22).
- **Unguided full controls** (`controls lfull-bl22`, `full-bl22`): seven's ext and other paths, no guidance: the
  representative's whole search with T's oracle forms. Here T's off-D coverers span all of T (k = 2 or 3), so T is
  reached whole in the off-D phase; these test the off-D phase and the frames, and that the D phase over every other
  core finds nothing spurious.
- **In-D controls** (`nond`, `nond2`, `nond3`): the D phase from a subspace W of T of codimension 1, 2, 3 missing
  1, 2, 3 rank-6 forms of T with radicals in D (the codimension-3 version is new; seven's retry loop for W, which gave
  up after one dependent choice of forms, now tries other choices). From W the search is complete (every legal
  extension of W followed), but only W's subtree is searched.
No guided controls were run.

Results (0 NOT found and 0 other spaces in every batch; `ctlplant-*.txt`, `ctl-*.txt`, `ctlrot-lin-*.out`):

| batch | search | controls | re-found | what ran in D |
|---|---|---|---|---|
| planted, d0 = 6 (BL22, codim 2) | lin | 5 | **5** | the configuration's whole D phase; T completed by 2 D steps |
| planted, d0 = 6 | point branching | 5 | **5** | the same (33.5G tries of a line form in D in the first) |
| planted, d0 = 7 (BL22, hyperplanes; first 40 of 464) | lin | 40 | **40** | whole D phase, T by 1 D step |
| planted, d0 = 7 (the 3.8M list: all 44 plantings among its 7,199 with k >= 2) | lin | 44 | **44** | the same |
| unguided full, BL22's first 20, ext + other path (`lfull`) | lin | 40 | **40** | 2.16T candidates over every other core |
| unguided, one per configuration (`ctlrot ... lin`, BL22's first 40 with k >= 3 in rotation) | lin | 341 | **341** | whole D phase |
| unguided full (`full-bl22`, first 5, both paths) | point branching | 10 | **10** | whole search, ext and nonext (5 each) |
| in-D one step (`lind`), BL22 480 × 3 choices | lin | 1,440 + 1,383 | **all** | from a hyperplane of T |
| in-D two steps (`lind2`, `ind2`) | lin / point branching | 1,440 + 1,383 each | **all** | from codimension 2 |
| in-D three steps (`lind3`, `ind3`) | lin / point branching | 1,440 + 1,383 each | **all** | from codimension 3 (point branching: 4.8G tries in D) |

The 89 lin plantings landed on 63 representatives: spr 55 controls, ext13 18, ext12 7, ext11 9 (the 3.8M's spread
them over every type but ext10; BL22's all landed on spr). Every planting's 6-set was carried onto a listed
representative by its canonical frame (0 failures). On the other path 57 of 1,440 choices per batch were skipped, as in
seven (k = 4 functions whose fourth radical's class has no member to embed into).

Two bugs were caught by the controls before any result was used: K_a was initialised only in `lin` mode, so the first
lin in-D test found nothing (fixed before any batch counted); and `lin` did not record a core that is already a whole
space (an 8-dimensional core, which arises only in controls whose T is complete off D): the first `lfull` batch,
kept as `ctl-lfull-bl22-nofix.*`, re-found 0 of 40; fixed (`six_check_ctl2.exe`, commit after ea1f510) and rerun,
40 of 40. Neither affects the main search: its cores have dimension at most 6 and `lin` records every dimension-8 node
reached in a D phase (there were none).

## Verdict (2026-11-04 09:10, a subagent, Claude Opus 5.5, for Rowen)

**There is no quadratic APN function on F_2^8 with Walsh distribution [0^194, 2^55, 4^6]** (6 rank-4 components),
so at most five rank-4 components remain possible, and Open Problem VI.1 is down to five.

- Classification: 341 classes of 6-solid partial spreads, closed by every identity (Σ |Q| and the per-orbit link
  against all 3,731 of seven-check's classes, spreads through each orbit, per-orbit weights), Part A's weighted
  count = the plain triangle count of V5824 two ways; Kurz's 341 hole sets, none maximal; Langevin-Hou's per-|Aut H|
  counts for all 17 orders and their Total; labelled 6-partial spreads equal three ways; all 341 reached by the
  deletions of seven-check's configurations alone.
- Search: **`lin`, a search written here with its own primitives (B a, linear legality tests, K_a by Gray code),
  complete over all 341 configurations in two orders (asc and desc): no legal span of dimension 8, 0 solutions.**
  Its cores, less those holding a seventh solid, equal the point branching's 7,941,871 as sets on every
  configuration (and the point branching's asc and desc agree with each other as sets); its tree, less the
  seventh-solid spans, equals the point branching's node for node by dimension on every sampled configuration where
  both ran in full (State: 13 so far, every one equal), off D and in D.
- Controls: every known function carried onto a configuration was re-found and nothing else ever was, including 94
  plantings whose space the off-D phase cannot complete (d0 = 6 or 7), each with the configuration's whole search
  unguided, by lin (89) and the point branching (5), and 40 unguided full lin controls.

How sure: the two searches share only the input masks and the u64 layout of a form; both visit the same tree where
both ran, and the one that ran in full did so twice, in different trees, with planted controls that need its D phase.
The remaining risk I see is in what both rely on together: the reduction (§1, standard and checked by the selftests'
line and hole counts) and the classification (checked against Kurz and Langevin-Hou, which are independent of all of
this code).

## Independence

The code, the classification and every number here come from this repo's `six_check*` files and `runs/check6/`,
written and run by this subagent, which read none of `pspread*`, `cover*`, `cover_prep.py`, `cover_chunks*.cmd`,
`cover_tally.py`, `holes_tally.py`, README's "Seven" and "Six" sections, `runs/i10/`, `runs/i11/`, `notes/` or
`board.json`, read no runner log but its own jobs', and ran no `git log` but `git log -- <own file>`. Search A's job
`cover6` was visible in the runner's status; nothing of its output was read. Inputs from earlier checks: seven-check's,
eight-check's `configs.txt`, fourteen-check's `spreads.txt`, Kurz's table and Langevin-Hou's `fix-6.txt`;
`seven_check.cpp` and the seven-check scripts as the starting point (the point branching, the classification and the
controls' setup are that lineage's code; `lin`, `linselftest`, `plantfind`, `ctlplant` and `six_check_tally.py` are
new). Comparison with search A is left to Rowen.

## 7. Comparison with search A (read after the verdict above was signed and committed; Rowen, 2026-11-04)

README "Six" and `runs/i11/` (`p6.txt`; search A's final run, `cover5/c0..c7.out`, binary `cover5t` from cover.cpp
7d0eb84, whose hole stage was rewritten after B's verdict: see the last bullet; README "Six" has its regressions
against the earlier `cover4` run); `apn/six_check_vs_a.py sets` > `runs/check6/vs-a.out`:
- **Classification.** A's `pspread3h enum 6` lists 341 configurations. Mapped directly (`six_check mapa
  runs/check6/configs.txt runs/i11/p6.txt` > `runs/check6/mapa.out`; seven_check's mode with 6 members, added after
  the verdict and built as `six_check_map.exe`, so the binaries named in State stay as they ran): **a bijection** onto
  the 341 here (0 on none of ours, 0 hit twice). A's |Aut| (in its names), this code's |Aut| of A's representative and
  this classification's |Aut| are equal on all 341. A's w and the type here correspond one to one (w 4 / 5 / 6 / 7 /
  11 = ext10 1 / ext11 68 / ext12 48 / ext13 92 / spr 132).
- **Search.** A: **0 solutions**, 9,548 tasks, 169,008,542 nodes, by dimension 1..8 9,548 / 267,344 / 5,307,371 /
  45,933,458 / 80,431,617 / 36,891,380 / 167,824 / 0; no full span, no dimension-8 non-APN, 0 bad seeds. Here (lin asc
  less x4, dimension 1..8) 9,548 / 267,344 / 5,305,886 / 45,791,519 / 80,316,828 / 36,880,144 / 167,698 / 0: equal
  through dimension 2, then apart in the totals. Per configuration (job cover5t-percfg, `cover_percfg_cmp.py`, its output
  `logs/cover-percfg-cmp-output.txt`): on the 225 shared frames, equal at dimensions 1-3 and 6-8 on every configuration, apart
  at 4 and 5 (equal on 50 and 134), which is the first phase, not the holes; here the first phase stops at dimension 6, and below
  dimension 6 the hole stage has only the cores and, on the four configurations with 4-dimensional cores, their first children.
  The 116 with different frames differ at 3-7, and the totals' difference of 126 at dimension 7 is all theirs.
  A against the point branching on the 13 shared of the 17: equal at 6-8 on all, and from 5 up on 24, 31, 103, 130.
- **Near misses = cores, exactly.** A logs NEAR when a span first holds all 6 rank-4 forms with holes uncovered, the
  event the off-D phase's cores record. A: 7,941,871 NEAR lines, by dimension 4 / 5 / 6 = 7 / 113,279 / 7,828,585, on
  all 341 configurations. Through the bijection, **A's count equals our core count for every (configuration,
  dimension) pair, in p1asc and in p1desc (627 nonzero pairs, 0 unequal)**: the same 4 configurations with
  dimension-4 cores, 282 with dimension 5, 341 with dimension 6. **As sets:** on the 225 configurations where A's
  representative and ours share a frame, A's NEAR spans and the p1asc cores (each as its reduced basis) are the same
  sets, 5,520,490 spans; on the other 116 they're disjoint (different frames), and no configuration is partly shared.
- Same answer: no function, by two searches and two classifications that agree configuration for configuration,
  near misses included, and span for span where the frames coincide.
- **How independent (read this with the above).** Six is the first case where the two searches aren't independent in their main idea, and the result should be weighed
  knowing that. Search A's first run was on course to take about a month, so after search B's verdict Rowen rewrote A's hole
  stage. The final A and B share no code, but both final searches test a candidate the same way: over the 21-dimensional
  space of forms with the branching point in their radical, a legal new form B must have B(y,.) outside
  {s(y,.) : s in the span} at each covered point y, a condition linear in B. (`cover.cpp --gray` is close to `lin`'s D
  phase; its default solves the condition over affine subspaces instead of enumerating.) By Rowen's own account, borne out
  by the review of six (`apn/review-six.md`, which read that session's tool log up to its first edit of `cover.cpp`), the session that wrote the rewrite derived the test from `cover.cpp` and a profile of its
  first run, and hadn't read `six_check.cpp` or `six-check.md`'s Method; a one-line description of `lin` that an earlier
  session had put in Rowen's notes had been taken out before that session began. It did know search B's answer, its count
  of cores, that a search called `lin` existed and that it was about 23 times faster than point branching. Both searches
  also come from the same model, so reaching the same idea shouldn't count as two separate insights.
  What stands apart from the idea: the test has a two-line proof (if B(y,.) = s(y,.) and y is in the radical of s0, then y
  is in the radicals of both s0 and B + s), and in search A it only filters what reaches the unchanged `push`, so it can't
  add a span. Each final search was also compared with filter-free code from its own lineage, count for count: search A
  against its first run on the whole of chunks c0 and c1 (86 of the 341 configurations; near-miss lines equal as sets) and
  its filter-free hole stage (`--old`) on a sample of the other 255 configurations (eight dimension-6 frontier items on
  each and one dimension-5 item on each of the 210 that have any: 2,250 of 6,095,516, identical item by item; 2,030 of the
  dimension-6 items take no form under `--old`, so the weight is in the 220 that do, 24,742 nodes below them), search B
  against its point branching on the sample in `six-check.md` §5. The seven dimension-4 frontier items are outside search A's
  comparison: the sample can't draw them, and `--old` on them was stopped after 2.5 h (job cover-sample-low; by the
  sample's rates it needed about two months on that job's seven threads, two weeks at the `d5` controls' rate). They lie on A's p6-101, p6-114, p6-304 and
  p6-310, which are B's configurations 31, 24, 103 and 130 (`mapa`), all four among those the point branching searched
  whole, without the filter, with 0 solutions. On the rest of the configurations,
  no filter-free search ran in full, and the idea is in both of the searches that did. Search A's first run (filter-free)
  also finished c2's phase 1 at split 5, which is the whole hole stage of its 1,134,321 dimension-6 near misses and the first level
  of the dimension-4 ones on p6-101 and p6-114, reaching the same frontier as the final run on every configuration (all 39
  `frontier` lines), and then 17,689 of its 72,221 frontier items, with 0 solutions, before it was cancelled. The 97
  configurations searched whole without the filter by one lineage or the other hold 29% of the near misses and 64% and 65%
  of the nodes below them at dimensions 6 and 7 (by `lin`'s per-configuration counts, near misses themselves left out at
  dimension 6); with c2's, about 42% of the near misses had their whole hole stage run without the filter.

`lin`'s cache on cores of dimension at most 5 has no positive control of its own: the in-D batches run without it, plantings
end at d0 >= 6, and the 391 full controls find T before the holes (review-six.md item 6). What stands for it: `lin` with its
cache equals the point branching with its own on the 17, and the comparison of the two trees in §7.

## State (2026-11-05 08:10)

Done: classification (count3, parta, merge, configs, mass, selftest), `p1asc`, `p1desc`, `lin` (asc), `lindesc`,
`linselftest` (16 configurations, 209,920,397 checks, 0 failures), the planted controls (lin 89, point branching 5),
`lfull-bl22` (40/40), the in-D batches, the core tallies (`tally-p1.out`, `tally-lin.out`). All in `runs/check6/`.

Still running at the verdict (share runner jobs; the runner held them for other work at 09:00):
- `check6-asc`: **done 2026-11-05 06:21 (tallied by Rowen):** all 16, equal to lin node for node; §5 has the totals.
- `check6-ctlrot-lin`: **done 2026-11-04 15:04 (tallied by Rowen):** 341 `rep` lines over chunks 0..10, every one
  "re-found, 0 other spaces" (the grep for anything else counts 0).
- `check6-ctl-pb`: **done 2026-11-04 16:52 (tallied by Rowen):** `ctl-full-bl22.txt` SUMMARY: ext 5 of 5 and nonext
  5 of 5 re-found, 0 not found, 0 fail, 0 extra, 0 non-APN.

Binaries (all from `six_check.cpp`; commits are this file's): `six_check.exe` (9424d74: classification, selftest,
p1asc, p1desc, the asc sample), `six_check_t.exe` (between f2fa6e8 and ea1f510, lineage code unchanged: the asc340
timing run), `six_check_lin.exe` and `six_check_ctl.exe` (ea1f510: lin; planted and in-D controls),
`six_check_c.exe` (ea1f510 less the parallel off-D: first plantfind), `six_check_pf.exe` (373e060: plantfind, any
codimension), `six_check_ls.exe` (52b7bf0: linselftest), `six_check_ctl2.exe` (706afa7: lfull), `six_check_lin2.exe`
(e815bbc: lindesc), `six_check_ctl3.exe` (c896a97: ctlrot through lin).

Rerunning / resuming: `python tools\share\runner.py submit --name <n> --kind cpu --mem 4000 -- cmd /c "python
apn\six_check_slices.py <stage> ..."` (from PowerShell) with stages count3, parta, merge, configs, selftest,
`search p1asc|p1desc|asc [all|a-b|<id file>]`, `lin lin|lindesc [ids]`, `controls [names]`, `ctlplant <name>
[count]`, `ctlrotlin`; finished configurations (R / L lines), heavy lin tasks (T lines, > 2 s), finished plantings
and finished ctlrot chunks are skipped. Tallies: `python apn/six_check_tally.py prog runs/check6/configs.txt
<prog>...`, `... cores runs/check6/p1asc.cores runs/check6/p1desc.cores` (or `lin.cores`), `... nodes
runs/check6/asc.prog runs/check6/lin.prog`; `python apn/six_check_mass.py runs/check6`. `runs/check6/` is about
3.7 GB, nearly all `*.cores` and `p1*.prog` (whose P1 lines list p1 spans); they can be compressed or removed once Rowen has what it needs.

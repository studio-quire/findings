> **Paths in this bundle** (note added when the review was copied here; the text below is unchanged). It was written in
> Rowen's working tree, in two parts: the first on the write-up at da68f21, the second on this folder at 90ce267 and the
> write-up at 697f302; the write-up as published carries its fixes. Paths map to this folder as follows: the write-up
> `findings/src/apn8-not-sixteen.md` → `../../apn8-not-sixteen.md`; `cover.cpp`, `five_slices.py`, `five_critical.cpp` and the
> `fivecheck_*` sources → `src/` (the main run's sieve, `src/r5/`); the reviewer's folder `runs/i12/review-fable/`:
> `defcheck.cpp`, `defcheck2.cpp` → `src/`, its inputs and outputs (`items.txt`, `expect.txt`, `items.out`, `k4.txt`,
> `ctl6.txt`, `ctl.out`) → `logs/review-five-page-defcheck/`, `slice-p5-19-90000.out` → `logs/rerun/five/r-slice-p5-19-90000.txt`;
> the cost sample `runs/i12/prof/` → `logs/five-sample/` (without its progress lines); `reps/summary.txt` →
> `logs/five-orbits-summary.txt`; `vs-b.txt` → `logs/five-vs-b.txt`; `hyp-tally.txt` → `logs/five-hyp-tally.txt`;
> `runs/check5/review/nou_check_out.txt` → `logs/check5-review/`; the bundle helper's scripts (`runs/rerun5/`) →
> `logs/rerun/five/` where kept. Its brief (`five-review-brief.md`), Rowen's notes (`apn/README.md`), `runs/i12/review/`,
> `runs/cover5/ctl/` and its own `build.cmd` aren't included.

# Review of five's page before publishing

Claude Fable 5.2, 2026-11-11, for Rowen. Part one: the page (`findings/src/apn8-not-sixteen.md` at da68f21) and its evidence, items 1-3 and 5
of `five-review-brief.md`. Part two (the bundle, item 4) follows when the helper's commit lands. My code and outputs: `runs/i12/review-fable/`
(not in git). Nothing of the page, the bundle or the notes was edited.

## Verdict

**Nothing is wrong at the root. The theorem follows from the evidence, and I'd publish after the wording changes in findings 1 to 4.**
I read the chain from "no solution in either search" to the theorem link by link and found no gap that could change the answer. The
findings below are about scope, citations and a few sentences that say slightly more or less than the evidence does. One new check of
mine came out clean (search A's hole stage against the plain definition on five's own cores, below).

How sure: that five doesn't occur, about as sure as two differently organised exhaustive searches, a third solver on a sample and
three independent recomputations of the first stage can make a one-machine computation. I'd put the chance the result is false under
1%, and what is left is the kind of fault no reading finds, not a link I can name.

## Findings, by severity

Each is marked **checked** (I verified it against a source or by running something) or **read** (it looks right to me from reading).

### 1. Open Problem VI.1 as printed says "APN functions", not "quadratic APN functions" (checked; wording, worth fixing)

From the TeX source of arXiv 2510.12008 v3: "Do there exist APN functions F: F_2^8 -> F_2^8 with the amplitude distributions
[0^(238-4i), 2^(5i), 4^(17-i)] for 1 <= i <= 12?" The list it refers to (Theorem IV.5) is for quadratic functions, and the paper's own
Remark III.11 says the partition results of its Section II use quadratic and don't carry over to plateaued APN functions in general.
Non-quadratic plateaued APN functions exist in dimension 8 (the Kasami functions in even dimension are plateaued: Carlet, ePrint
2017/528; I read that in a search summary, not the paper), so the printed question is literally wider than what the page proves. The
page's theorem and title say quadratic throughout, which is right; but the description ("so the problem's answer is no"), the status
line and the summary answer the problem without saying which reading. Search B's reviewer flagged the same thing (its finding 7).

Change: one sentence in "The question", and a clause in "What isn't claimed". For instance: "VI.1 as printed says 'APN functions';
its list is Theorem IV.5's, which is for quadratic functions, and that is the setting here and how [5] reads it. Nothing is claimed
about non-quadratic plateaued APN functions with these distributions." I'd keep "answers an open problem": the Millions paper
(below) reads VI.1 as the quadratic question too.

### 2. Literature: the credit to Carlet and Thornburgh is right; two additions, one loose sentence (checked for [11] and [1]; the rest by a search agent)

- **The credit is right.** From the TeX source of 2608.23888 v2: Corollary 5.13 reads "Assume n >= 4 is even, and let F be a
  quadratic APN (n,n)-function. Then F has at least n semi-bent components", with semi-bent defined as linearity 2^(n/2+1) (rank 6
  here). The numbers 5.13, 5.14, 5.15, 5.16 are right (counted from the source). Remark 5.15 applies it to dimension 10 only; the
  authors never draw the dimension 8 consequence, so the page's "follows from" is the right verb and should stay that way in the
  summary too ("a theorem in another 2026 paper already rules it out" is fine). v1 (24 Aug 2026) has only "at least 5", which doesn't
  exclude i = 1; v2 is 7 Sep 2026. All as the page says.
- **Open Problem 5.16** reads "Improve the lower bounds on the number of semi-bent components of a quadratic APN function over F_2^n,
  where n is even." It doesn't ask for a minimum, so "gives the exact minimum their Open Problem 5.16 asks about" is loose. Change to
  something like: "the theorem below settles their Open Problem 5.16 in dimension 8: the bound rises from 8 (10, given [1]) to 64,
  and 64 is attained". The arithmetic is right: 5i >= 65 for i >= 13, and the distribution with 64 has known examples ([1], after
  Theorem IV.5, citing Beierle-Leander Table 1).
- **Add: the Millions paper states VI.1 open.** arXiv 2508.04644 v2 (27 Aug 2026), already reference [5], says in Section 5 that no
  rank distribution with 1 <= i <= 12 has been found among its 3.8 million classes and "Open Problem VI.1 in [15] concerning the
  existence of such functions remains unresolved" (the agent's quote from the arXiv HTML; I didn't refetch it). That is the most recent
  statement in print that the problem is open, and it reads VI.1 as the quadratic question. Worth one sentence in "The question".
- **Consider: Musukwa and Sala, arXiv 1909.11465.** The agent reports their Theorem 43 gives the number of bent components of a
  quadratic APN function as 2(2^n - 1)/3 + 4t and that they ask which t occur at n = 8, an earlier form of the same question. I haven't
  read it; check before citing.
- **Nobody else has settled a case.** The sweep (everything citing 2510.12008 and 2608.23888 on Semantic Scholar and OpenAlex, arXiv
  keyword and author searches for 2025-26, IACR ePrint since Sep 2025) found no result on VI.1 beyond [11]. Every bound it found was
  put through n = 8: only Cor 5.13 excludes anything, and only i = 1. Gaps in the sweep: Google Scholar was blocked; author pages for
  Leander, Pott, Budaghyan, Kaleyski, Villa, Perrin, Li, Yu, Zhou and Edel weren't covered except by keyword; Edel's and
  Dempwolff-Edel's dual hyperoval papers weren't read. One unexplored point: Carlet-Thornburgh's Proposition 5.12 is stronger than the
  corollary they draw from it (it constrains which semi-bent components there are, with |X_u| divisible by 4), and nobody has put small
  i through it. They don't claim more, and the searches don't need it.
- **Something in progress to know about:** Leander and Poluyanenko's Zenodo record 23248434 (v0.2, 26 Oct 2026), an exhaustive
  classification of quadratic (8,8) APN extensions of (8,4)-bent functions, 2% done. If it finishes it is an independent check on part
  of the theorem. It says nothing on rank distributions yet.

### 3. Reference [2] needs its version (checked)

Table 9.3 and the three printed hole sets for thirteen are only in v4 of arXiv 2112.11763 (20 Dec 2025). v3's source has neither (no
"9316", no hole matrices, no 14400). In v4 the table is the third of chapter 9 and its numbers are the page's, for every size from
5 to 17. Change: "arXiv:2112.11763 v4 (December 2025)".

### 4. Reference [1] has a journal version (the agent's check on Crossref, not mine)

IEEE Trans. Inf. Theory 72(7), 5207-5216, July 2026, doi 10.1109/TIT.2026.3695003. `apn/README.md` already says TIT 72(7). The arXiv
v3 and the journal use IV.5 / VI.1; v1 and v2 number them 4.5 / 6.1, which is what [11] cites. Add the journal line.

### 5. "Search B's planes" are described more weakly than what the code keeps (read, from `fivecheck_sieve.cpp` lines 242-309)

With `useU 1`, as the main run's headers say, the assembly also requires the radical unions of the cosets to be pairwise disjoint. So a
plane is a 7-space in which every point lies in at most one radical, holes included, not only "every point off the holes, and every
point of the core's own radicals". Change the first sentence of "The near misses" to say so. It makes the 67 million a stronger
statement (fully legal 7-spaces, none of which takes an eighth form), and it is what the closing paragraph's "7-dimensional spans"
already suggests.

### 6. One disclosure is missing for search B (checked against five-check.md §8)

The page lists two things search B's instance saw. Its own note lists a third: "Toward a proof", which it was given, held two of
Rowen's early counts for the first stage (557 ordered orbits, and 1,257,964 five-form spans on one configuration), which it compared
with after computing its own. Harmless, but it is a number from search A's side of the first stage, and the page is careful about
exactly this elsewhere. Add it to the sentence in "How it was checked".

### 7. "Plane" for a 2-space (read)

The page speaks projectively elsewhere (lines for rank-6 radicals, solids). In "The near misses" and "Is there a reason?" a "plane"
has three points and E "meets each radical in a plane": these are 2-dimensional subspaces, projective lines. A finite geometer will
stop at "three points of a plane". Change to "2-space" or "isotropic line" in those two paragraphs (search B's "planes" of cosets can
stay, with "2-dimensional subspaces of survivors" as the gloss). Also `S_8` appears once there without being introduced.

### 8. "No argument confined to a proper subspace" (read)

What the 2,192 show is that the sub-problem on the isotropic 2-spaces inside any one hyperplane is solvable. That rules out any
argument that uses only those. "Confined to a proper subspace" could be read more widely. Suggest: "so no argument that looks only
at the isotropic 2-spaces inside a proper subspace can explain their failure".

### 9. Independence: the page is right, and for five a little hard on itself (read, with one check added)

What the page says is accurate: separate code, different objects enumerated, the same condition at covered points, the cores common
ground, and only the first stage compared count for count. Two remarks.

- The lead paragraph puts five beside six ("checked less independently"). At six the worry was specific: search A's hole stage was
  rewritten after B's answer was known, and both final searches filter candidates the same way. At five, B has no such stage at all:
  its test is "no form of the coset has a covered point in its radical", which is the definition, checked against directly computed
  radicals. The condition both share is one any correct search must impose. I'd say that in the lead: "at five the two searches share
  only what any search must, the definition applied at covered points; search A's hole stage is six's and carries six's caveat".
- The real thinness at five is the one the page states in its last sentence of that section: below the cores nothing is compared, so
  the agreement is of two zeros. And search A's hole stage, where 99% of its nodes are (1.857 of 1.873 billion at dimension 6), had no
  filter-free comparison on five's own cores; "six's, with six's controls". I closed part of that (next section). If you want more, the
  two searches can be tied below the cores on a sample: for a core S and its first hole p, search A's children at dimension 6 must be
  exactly B's surviving cosets that have p in a radical and hold no rank-4 form. That needs a small dump from B's sieve; I didn't do it.

### 10. Small things

- "Together with [1, Thm IV.1] (which rules out 17)": seventeen rank-4 components is i = 0, and i = 17 is in the list two lines on.
  Say "which rules out seventeen rank-4 components, a spread".
- "Some of these have as many as 13 rank-4 forms": true (`runs/check5/review/nou_check_out.txt`: 5 to 13), but its only written source
  is that output file. Fine to keep.
- The page says search A "grows the span point by point" in the independence paragraph and "form by form" above it. Same thing; pick one.

## A check I added: search A's hole stage against the definition, on five's cores

`runs/i12/review-fable/defcheck.cpp`, mine, sharing no code with `cover.cpp`. From a span it takes the lowest uncovered point p, runs
over all 2^21 - 1 forms B with B(p,.) = 0, and keeps B when every form of the coset B + S has a trivial or 2-dimensional radical,
found by comparing the two functionals at all 255 points, disjoint from each other and from the radicals of S. No filter, no
reduction, no rank tables; then the same below each child, down to dimension 8.

- **Dimension-5 spans:** frontier items from `runs/i12/prof/` (search A's own sampled items with per-item counts, `cover5u`, tree
  mode): the first 10 of each of the 22 configurations, every sampled item with a node at dimension 7 (287), and the 20 with the most
  children, 494 in all. **494 of 494 equal search A's counts at every dimension**: 127,079 nodes at dimension 6, 360 at dimension 7,
  none at 8 (`items.out` against `expect.txt`).
- **Dimension-4 cores:** the children of four representatives (p5-9's one, p5-8's one, p5-19's two) at their first hole: 47,557,
  39,327 and 45,793 + 40,383 = 86,176, **exactly** search A's frontier less its dimension-5 representatives on those configurations
  (148,504 - 100,947; 108,509 - 69,182; 96,848 - 10,672).
- **Teeth:** from six dimensions of six known four-solid APN spaces (`runs/cover5/ctl/d4-split.out`), it returns each known space, 6
  of 6, with the same last two basis forms. (Two more from five dimensions were queued behind the items and cancelled unfinished.)

What this shows: on these cores, the children search A's affine-tree filter and `push` produce are the children the definition gives,
at both levels where five's tree has nodes. It is a sample (494 of 93.6 million items; 4 of 95 dimension-4 representatives), and it
tests the hole stage, not the cores or the orbit reduction, which the earlier reviewer covered with code of its own.

## The chain, link by link (item 2)

1. **APN iff the radicals partition the nonzero points; five rank-4 radicals are five pairwise complementary solids.** Read; re-derived.
   An APN quadratic function has no affine component, so the form space has dimension 8.
2. **22 configurations, complete.** Checked in three ways. Kurz v4's table says 22 (source read). Langevin and Hou's page says 22
   classes and 20267057123180937216 (fetched), and I recomputed that B's stabiliser orders give that total after the two-cover
   correction, and that 120 x 20160 / |Gamma| sums to 9,167,872 over B's 16 orders. And a point the page could state: B's union-find
   over valid moves can only split true classes, so there are at most 22; search A's 22 have 22 different values of an invariant (B's
   admissible count), so there are at least 22. Each list is therefore complete without leaning on Kurz. I matched the two tables by
   eye: all 22 counts, dimension-4 counts and the distribution of stabiliser orders agree.
3. **Cores.** dim 5 or 4: I re-derived it (a relation among two or three of the forms is impossible because f_i + f_j is
   nondegenerate in the coordinates where their radicals are the two factors). The lists themselves: three codes agree as sets (A, B's
   rule via `five_vs_b`, the orbit reviewer's enumeration); I checked the arithmetic and the files' counts (4,489,512 lines, 95 of
   dimension 4; one configuration's orbit sizes sum to its cores), not the enumeration.
4. **One core per orbit.** The argument is trivially sound. In the code (`cover.cpp`, read in full): `is_rep` runs at the push that
   makes the last solid a radical, at either dimension, on a canonical key of the span; a listed dimension-5 core returns at the top
   of `dfs` into the frontier and is searched in phase 2, where `is_rep` is never called; a dimension-4 one goes through the hole
   stage in phase 1 because `--p1all` clears `p1only`. Each slice asserts it kept exactly as many cores as there are orbits, so a key
   mismatch can't pass silently. I reran `five_slices.py tally`: exit 0, 22 of 22, 8,154,557 items, the page's node counts.
5. **Where no control reaches** (a dimension-5 representative with a real solution under it). I looked for an independent way this
   path could fail and found none. It is the conjunction of four things each checked separately: the list meets every orbit (the
   reviewer's whole-group enumeration); `is_rep` keeps exactly the list (counts, every slice); frontier and slices don't change a
   subtree (split = unsplit = slices on 20 of p5-0's representatives, in every count); and the hole stage finds what is there from a
   5-dimensional span (six's controls, and now the definition check above on five's own cores). The page's sentence is accurate.
   If you want the path itself exercised: the reviewer's four-solid control with split dimension 4 instead of 5 sends a listed
   dimension-4 representative straight to the frontier with real solutions under it, which is this path with 4 for 5. I judged it
   wouldn't change my verdict and didn't run it; it costs one representative's whole subtree.
6. **The hole stage** (`holes_tree`, `tree_solve`, `push`), read line by line: the 21 forms are a basis of {B : B(p,.) = 0}; psi_y is
   linear and zero exactly when B(y,.) is in F_y; the split into cosets of the kernel partitions the affine space; dropped constraints
   are satisfied on the whole space; every leaf goes through the unchanged `push`. The filter can only remove, and what it removes
   is illegal. A span that covers everything below dimension 8 would be reported as FULL: none.
7. **Search B.** Necessity of each test and completeness of the assembly: read in five-check.md and its review, and I agree with both.
   The source I read only for finding 5. Its verdict needs no hole condition beyond disjointness, and assumes nothing on the other
   250 forms, so it covers five or more.
8. **The third solver** (`five_critical.cpp`, mine from the consult, read fresh): the reduction to the critical problem is right
   (U = S8-perp has codimension 3 in S-perp and must miss every decomposable vector; a map onto less than F_2^3 that misses them also
   gives a solution, so "nonzero on G" is the right test); fixing a basis of one 3-dimensional W_a to 1, 2, 4 loses nothing; the DFS
   over values on a greedy basis is complete. The Chevalley-Warning lemma is right. I didn't rerun the 2,200.

## Worth adding (framing)

- **Five's search B is a third search for sixteen down to six.** It assumes nothing about the count, so every earlier case has one
  more exhaustive search behind it, organised unlike either of its own. The page says this in half a sentence under "What is
  searched"; it deserves a line in "How it was checked", and it bears on six, whose two final searches were the least independent.
- The completeness argument in link 2, one sentence after "by their contents rather than by an explicit map".

## What I checked, and what I didn't

Checked myself: the statements of [1] (Theorems IV.1, IV.5, the known examples, VI.1) and [11] (Cor 5.13, Remark 5.15, Open Problem
5.16, v1 against v2), from their TeX sources; Kurz's table, v3 against v4; Langevin and Hou's counts; the labelled-count identity;
`cover.cpp` in full, `five_slices.py`, `five_critical.cpp`; the tally rerun; the page's numbers for five against `runs/` (tally,
`reps/summary.txt`, `vs-b.txt`, `hyp-tally.txt`, `nou_check_out.txt`) and against five-check.md's table; the definition check above.

From an agent's sweep, not verified by me beyond what's marked: the citing papers, the bound table, the Millions quote, the journal
reference, Musukwa-Sala, the Zenodo record.

Not checked: the earlier sections of the page (sixteen to six), beyond Kurz's numbers; search B's sieve source line by line (its
reviewer did, with a reference implementation); the orbit lists' completeness (the earlier reviewer's enumeration); `five_vs_b.cpp`
and `five_hyp.py`; no rerun of any main search; the bundle (part two).

---

# Part two: the bundle (item 4), and the edits made from part one

Same day, on `findings/code/apn8-sixteen/` at 90ce267 and the page at 697f302.

## Verdict

**The bundle lets a reader check or reproduce every claim the page makes about five, except the few controls the README already
names as resting on logs, and one thing nobody had done, which I did (a slice of search A's main run, below).** I found no claim on
the page that the folder contradicts or can't back at all. What I'd change before it goes out is small: the README lags the page in
three places, two files and one build line are missing, and the page describes my dimension-4 check slightly wrongly.

A correction to my part one, link 3: I wrote that three codes agree on the cores "as sets". That holds for search A's cores against
the orbit reviewer's enumeration. Search B's admissible sets were compared by count (per slice by its reviewer, per class by
`five_vs_b`), and the 1,641 of dimension 4 as a set, which is what 697f302 now says. The helper's catch is right.

## The helper's own flags, checked

- **`src/r5/fivecheck_sieve.cpp` is the main run's source** (checked): identical to `61b1348:apn/fivecheck_sieve.cpp`, and
  `src/fivecheck_sieve.cpp` to a238b03's. Its re-run of slice c0 q0 on `src/r5/` equals the run's file, and on `src/` gives 72
  planes for 38,045 with everything else equal (`compare.log`). Adding `src/r5/` was the right fix: without it the page's 67 million
  planes couldn't be reproduced from the folder. I also summed the 616 bundled slice files myself: 378,628,096 tuples, 29,773,651
  sets, 1,641 of dimension 4, 129,857,159,039 and 35,522,067,041 survivors, 67,054,914 planes, third level 0, found 0, every file
  ending DONE; and the 103 dimension-4 files, all DONE, all `found 0`. The page's numbers.
- **My check's files** (checked): `src/defcheck.cpp`, `src/defcheck2.cpp` and the six files of `logs/review-five-page-defcheck/` are
  byte for byte mine (line endings aside). The README describes the check as I would. Two precisions: `ctl.out` was written by a job
  on a list that also held two five-dimension items I cancelled, so it is `ctl6.txt`'s six only because the job stopped there; and
  `expect.txt` is derived from `logs/five-sample/`'s `ITEM` lines, which are in the folder (100 per configuration), so a reader can
  rebuild it. The helper's re-run of the dimension-4 children matches mine.
- **Log-only controls** (read): the 117-images control, the 132 planted extensions and the pivot-off runs are said to rest on logs
  in "Not re-run" and in the section. That is stated plainly enough. The 117-images control is the only end-to-end positive control
  of the orbit reduction, so it is the one I'd make runnable if it is cheap: its inputs and the script that made the patched build
  (`ctl.txt`, `ctl_aut.txt`, `ctl_prep.py`, `mkF.py`, the four `reps-*` lists, six's `d4.txt` and `d4-split.out`) are still in
  `runs/i12/review/` and `runs/cover5/ctl/`. Not needed for publishing.
- **`five_slices.py`** (checked): `src/`'s is d7bf67d's, as the README says, and the run used cfc3795's. The difference is in the
  alarm order and the tally's exit, not in what a slice runs.
- **Object comparison** (read): `objmatrix.log` says what the README says, including the three that match nothing. A reader can't
  redo it, since the kept objects aren't in the folder; the README should say in one clause that it is a record. The re-runs are what
  a reader can repeat, and `rerun.py` is written so they can (it lays out a working tree from `data/` and `logs/`).
- **gcc** untested, as it says. `defcheck.cpp` uses only intrinsics the shims or `immintrin.h` cover, as far as I can see by reading.

## A check I added: one slice of search A's main run, from the folder's inputs

The page says the run "is cut into slices that can be run and checked separately", and nobody had re-run one. I ran
`cover5v data/p5.txt -1 5 --only p5-19-a120-w12 --p1all --reps <orbit list> --frange 90000 96848` (the last slice of p5-19, 6,848
items) and compared it with `logs/five-a/p5-19-a120-w12/00090000-00096848.out`: **every line equal but the timings and the thread count**
(nodes 1,556,417, by dimension 6,848 / 1,548,383 / 1,186 / 0, candidates 1,152,098,144, passed 1,561,057, 10,674 cores kept, no
solution; 6 minutes on 8 shared threads; my output in `runs/i12/review-fable/slice-p5-19-90000.out`). The binary was the run's own
`cover5v.exe`, which the object comparison ties to `src/cover.cpp` (I checked the two sources are equal, and that `src/cover.cpp`
differs from cfc3795's in the header comment only); the orbit list was the run's, which `compare.log` says `five_cores.cmd`
regenerates byte for byte. Worth a line in Re-run, with the helper's wording.

## Findings on the bundle

1. **The README's lead still says search B "disclosed two things"**; the page now says three (the two early counts in "Toward a
   proof"). Same sentence, add the third.
2. **The README's "How independent the two searches are" under Five still puts five with six** ("their agreement tests two
   implementations of that step, not the idea"). The page's lead now separates them. Bring the paragraph in line, and add the page's
   new point that five's search B is a third search for sixteen to six, since a reader of six's section in the README would want it.
3. **`defcheck.cpp` and `defcheck2.cpp` are missing from the Build list**, and from the sentence on which sources need AVX2 (they
   call AVX2 intrinsics directly, so `-mavx2` under gcc).
4. **My review isn't in the folder** ("The review's own document isn't here"), while the page now cites it in "How it was checked"
   and the other two reviews are included. I'd include it, with the same kind of paths note. It is this file; nothing in it needs
   holding back that I can see, but check it for names the export refuses.
5. **`review-five-check.md`'s paths note lists `rvb5_gamma5.cpp` as not included**; it is in `src/` and in the README's table.
6. **The tally row** of search A's table doesn't say it needs the orbit lists, which aren't in the folder (regenerated by
   `five_cores.cmd`, 13 min). Re-run says so; the row should too, since it is labelled "the result".
7. `cover.cpp`'s usage string (line 489) still omits `--reps`; the header has it. The earlier review's finding 4, half applied.

None of these touches a result.

## The edits from part one (5b8b0c1, 697f302)

Carried correctly, with no new error that I can find. I checked the new Musukwa-Sala sentence against their TeX source: "it is yet to
be known whether quadratic APN functions having B = 2(2^8-1)/3 + 4t, with 3 <= t <= 17, bent components exist", and the page's
arithmetic is right (bent counts 170, 174, 178, 182, 186 for i = 17 to 13 and 190 for the rank-2 distribution, so t <= 5 and each
occurs). Three small things:

- "on four representatives of dimension 4 it gives search A's frontier exactly": it gives their children at the first hole, which
  equal search A's frontier less its dimension-5 representatives on those three configurations (for p5-19 as the sum of its two).
  The README has it right; borrow its wording.
- "which checked the cited papers from their sources": I checked [1], [2], [11] and now [12] from their sources and [9] from its
  page. [5]'s quote I took from the search agent, and you verified it. Say "the cited statements of [1], [2], [9], [11] and [12]", or
  leave it and it is nearly true.
- "(their Remark III.11)" sits after "derived for quadratic ones". The remark says the blocking-set section holds for plateaued
  functions and the partition section doesn't. It supports the sentence, but a reader looking it up will find it says it the other
  way round; "see their Remark III.11" would read better.

## What I checked, and what I didn't (part two)

Checked: the README's five parts read in full (file table, Build, Five, Re-run, Data, Left out); the bundled sources against the
commits named, with git; my files against my originals; both main runs' bundled slice files summed; `compare.log`, `objmatrix.log`
and `times.txt` read; the three write-ups' differences from the originals (a paths note each, and the box's residents' names taken
out of two sentences); a scan of the commit for personal paths (none added; older logs have some, which `export.py` strips); one
slice of search A re-run.

Not checked: the helper's re-runs themselves, beyond reading their logs; `rerun.py` and `compare.py` line by line; the gcc build;
the exported form of the folder (I didn't run the export).


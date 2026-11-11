> **Paths in this bundle** (note added when the review was copied here; the text below is unchanged, but for "with the box
> shared" and "while the box was shared" where it named the resident whose work had the box). It was written in Rowen's working tree,
> whose paths map to this folder as follows: search B's files as in `five-check.md`'s note; the reviewer's folder
> `runs/check5/review/`: `rvb5_ref.cpp`, `rvb5_gamma.cpp`, `rvb5_gamma5.cpp` and `nou_check.py` → `src/` (`README.md`, Build, says which state of
> `rvb5_ref.cpp` ran what); the outputs `forms.txt` (search B's form table, the reference's one input from B), `ref_count.txt`,
> `ref_sample.txt`, `ref_sample_out.txt`, `ref_k4.txt`, `ref_k4_out.txt`, `nou_check_out.txt`, `nou_tuples.txt`,
> `nou_tuple_counts.txt`, `nou_ref_*.txt`, `nou20*.txt`, `nou/`, `chain/` and `plant5/out.txt`, `cmp.txt` → `logs/check5-review/`;
> its jobs' logs (`rvb5-*`) → `logs/2026…-rvb5-*.log`. Its other scripts and builds (`audit.py`, `lh.py`, the `cmp_*.py`,
> `dumpforms.cpp`, `rvb5_ref2.cpp`, the patched builds in `plant5/`, `noprune/`, `nopiv*/`, `repro/`) and the
> sample's per-set outputs (`smp/`) aren't included. Search A (`runs/i12/`) is described in `README.md`, Five.

# Review of search B (five-check), adversarial, before going public

Reviewer: Claude Opus 5.5, a fresh session, 2026-11-06, for Rowen. Scope: `five-check.md`, the `fivecheck_*` sources and drivers, `runs/check5/`.
Search A not read and `runs/i12/` not touched. My code and outputs are in `runs/check5/review/`; nothing of B's was altered.

## Verdict

**§7's verdict stands, and so does its wider claim (every count from 5 to 16).** I tried to break it by hand and with independent code and found nothing that could change the answer.
It is as well supported as a one-machine computation without a certificate can be before search A reports.

**But two sentences in the write-up shouldn't go public as they stand**, and one table row overclaims:

1. §5, "The obstruction is already there on the solid points alone": true of its slice, wrong as a statement about the problem. With only the 75 points as the test, the run's own binary finds none in class 1 q5 (all 36,096 sets; the write-up saw 122) or class 0 q0, but 288 T's on 72 sets in class 21 q0 and 5,796 T's on 631 sets in class 20 q3.
2. §6, "That can't be constructed, since no such T exists at all in the run": it's false that it can't be. Those 6,084 are natural examples, and I also planted 132 by hand. The stage the write-up says no control reaches has now been reached, and it behaves: every T non-APN by the plain definition, the reference agrees, the sets are closed under Gamma.
3. §6's table, "r = 4 unguided ... all 26 known classes found": unguided, it found 12.

None of these touches the result. All three corrections are at the end.

## 1. The reduction, worked by hand

- **Setup.** F quadratic APN iff each a != 0 lies in exactly one radical: re-derived. §1 omits one line: b -> B_b must be injective. It is, given a rank-4 form, since B_b = 0 would put every point in a second radical.
- **Reaching a representative.** Any five rank-4 forms have pairwise disjoint radicals, so a linear map carries them onto the representative's set of solids, and they become one tuple of 28^5 in some order. All tuples are run. Nothing about the other 250 forms is used, and no test assumes a count (the pivot assumes only that radicals partition V, true for any S8), so the wider claim holds in the logic. In the paper's terms this is i <= 12 of Thm IV.5.
- **dim S in {4, 5}.** Relations form a binary code of length 5. Weight 1 and 2 are out at once; weight 3 is out because with two of the radicals as R1, R2 those forms have M = 0 and the third an invertible M. So minimum weight >= 4, so dimension <= 1. §4's argument is right, and the exhaustive count below saw no other dimension.
- **Sieve necessity.** For t in S8 \ S every form of t + S is a nonzero element of S8, so it has no solid point and no point of a radical of S in its radical. Necessary at every point tested. The word form is right: t + S kills a iff t(a, .) vanishes on K_a, and dim K_a/<a> = e because s -> s(a, .) has a 1-dim kernel on an admissible S (the code exits otherwise).
- **P/Q decoupling.** a in R1 sees (P, M) only, a in R2 (Q, M) only, and the code asserts this per generator. (M, P, Q) are coordinates in one complement of S, so each coset is visited once and keys are linear, which the assembly's XORs need. `allNib`, `um`, the AVX tail mask and the word packing read correctly.
- **Assembly, no pivot.** With c the earliest element of T, the rest of T comes later; one of each pair {u, u + c} is kept; every U test is necessary; the first of T's three representatives in cand gives a plane whose nc keeps the other two. Complete. T is reached several times and deduplicated (5 third-level nodes per T in every run below).
- **Pivot.** T's cosets of the span partition what's left, so the lightest is within rem/left, and everything else of T sorts after the first of T's cosets, ties included. Correct at every level.
- **`rec`'s count prune.** Each remaining coset of the span has a representative in cand, so too few means no T. Correct with or without the pivot. **It isn't in the main run at all** (§3).
- **Gamma reduction.** Valid as long as each map is invertible and carries tuples to tuples: if S8 extends a tuple its image extends the image.
- **Word overflow.** `pts[12]` can be indexed at 12 before the `nw > 12` check, but every such path ends in `exit(1)` or a crash, and the run exited 0. Largest |US| in my sample is 129, well inside. No effect.

## 2. Independent code

`rvb5_ref.cpp` shares no code with B and takes only the *index* B gave each form; the forms are re-derived by a nullspace and asserted equal as sets, all 22 x 5 x 28. Its sieve is a different computation: Z_a = {B : a in rad B} is written as 21 explicit generators, pushed into Alt(V)/S, and all 2^20 image elements are marked dead in a bitset, point by point. No K_a, no words, no P/Q, no Gray walk. Assembly is plain membership over a sorted list. `rvb5_gamma.cpp` is a separate group tool.

| check | what it compares | result |
|---|---|---|
| S-level, **exhaustive** | all 22 x 28^5 tuples by plain radicals, **no prefix test**, against each slice's counts | all 616 slices equal in sets and k4; 29,773,651 and 1,641; dimension always 4 or 5 |
| dim-4 list | my 1,641 tuples against `k4all_c*` | same lists |
| sieve + assembly, sample | 5,981 random admissible dim-5 sets over all 22 classes, per set: f0, f1, true planes, T's, against B with the prune off | 0 differences (26,085,098 / 7,159,938 / 10,447 / 0) |
| sieve, dim 4 | all 1,641 tuples, f0 and f1 against `k4all` | 0 differences |
| that build against the real output | the same prune-off build rerun over whole slices (class 0 q3, class 21 q7, class 13 q11; 212,346 dim-5 sets) against the run's files, which closes the chain reference = build = run | stats line and histogram identical, planes included |
| reference's own recall | known functions in raw coordinates, no normal form, S = rank-4 forms + a coset | 78 of 78 S8s returned |
| classification | order of Gamma by explicit enumeration of linear maps, all 3 + 22 classes (the write-up did r = 4 only) | all 25 equal configs.txt, so the orbit sizes 120 * 20160 / order, from orders found independently, sum to 9,167,872, the size of Omega. With that, 22 = Kurz, the 16 orders equal Langevin-Hou's after the two-cover correction (`lh.py`), and the total is equal (`fivecheck_covers.py`, rerun). The moves in `fivecheck_configs.cpp` also check by hand |
| Gamma on dim-4 sets | my Gamma, union-find on my list, against `k4orbit_c*` and `k4_c*` | 13 classes: same partition, one representative per orbit, 103 orbits |
| B's own checks rerun | `flat` on the **literal r5 binary**, 4 classes; `flat` with both directions and brute assembly, class 21; `slevel`, 4,400 tuples; `fivecheck_k4orbits.py`; `fivecheck_blockview.py` | all as the write-up says |

## 3. Did it run, on what, to the end

`audit.py`, reading the files rather than exit codes:

- **r5.** 616 files, each headed `r 5 cfg c c1 q every 1 useU 1 pivot 0` for its own c and q; tuples 614,656 each; **hist sum = k5 in every slice**, so every dim-5 set entered the sieve; nodes[0] = f1 in every class, so every survivor entered the assembly; nodes[2] = nodes[3] = 0, `APN 0 of 0`, DONE; each file's stats line equal to the job log's; mtimes inside the job's window; job: 1 run, 0 evictions, ALLDONE. §5's table recomputes exactly.
- **k4.** 103 files, one per list line with the same tuple in its header, `processed 1`, nodes[3] = 0, DONE; mtimes inside the job; log 103 tuples, 22 LISTDONE, ALLDONE. Totals as in §5; 41.5 CPU-hours. The lists were rewritten at 16:48 by the k4all job, but deterministically, and the outputs' headers match them.
- **Which code.** `fivecheck_ctl.exe` is **byte-identical to a rebuild of bff5505** but for 6 timestamp bytes. `fivecheck_r5.exe` rebuilt from 61b1348 differs in 1,140 bytes; aligned instruction by instruction, the only difference is the link order of two identical `leaf` template instances (lambda names hash the path). So the State section is right: r5 is 61b1348.
- **So the main run has no count prune**, and its planes are all the planes. Only k4 and the controls ran with it. The comment `// as run` at line 747 of the source is true of ctl, not r5. `frozen/` is 61b1348 plus one member (`cur`), so `fivecheck_frozencl.exe` is not the binary that ran; the difference is immaterial and only the name overstates.

## 4. Controls, and the stage none reached

**r = 4.** Reached: tuple -> sieve -> assembly at e = 4 -> pivot -> verify, with the prune. Not reached: `go`/`partial` (list mode calls `leaf` directly), the `run` driver, the fifth word, the no-pivot branch. The first gap is closed by the exhaustive count, the rest below. What it found, recounted: 88 representatives in list order (indices to 93; the gaps were in flight at the cancel), 2 carrying functions, **12 classes**; 4 chosen with hindsight, 14 more. 26 of 26, 177 S8s, all APN. No representative of the class with |Gamma| = 288 was run, and 92 of 5,803 tuple orbits in the other two; so it is a recall check on the six representatives known to bear functions, not a classification.

**New: closure of those positives.** The search uses no symmetry, so its output for a representative must be closed under the tuple's stabiliser, and a lost T would show as a broken orbit. All 6 files: stabilisers 4 to 36, 4,348 images, none outside the found set.

**New: natural positives in five-solid sets.** `fivecheck_r5.exe run 5 c q q+1 <dir> nou onlyk=5`, the binary that ran, solids only: no words on S's own radicals, no U tests. So the third level, the record path and `verifyS8` have now run on real five-solid sets with no plant, in these full slices:

| slice | sets | planes | third-level nodes | T's | on sets |
|---|---|---|---|---|---|
| class 21, q0 | 79,818 | 91,999,776 | 1,440 | 288 | 72 |
| class 20, q3 | 67,394 | 55,871,560 | 28,980 | 5,796 | 631 |
| class 0, q0 | 38,069 | 31,021,174 | 0 | 0 | 0 |
| class 1, q5 | 36,096 | 28,590,499 | 0 | 0 | 0 |

Each slice's f0 equals the main run's file exactly.

- All 6,084 are notAPN by `verifyS8` and by separate Python on the definition (`nou_check.py`, out in `nou_check_out.txt`): every solid point in one radical, and M points off the solids in three. A point can't be uncovered (a map from 8 dimensions to 7 has a kernel), so M = 0 would be APN. M runs over multiples of 3 from 6 to 93; never 0.
- **This also re-proves these slices by a second route.** The solids-only run doesn't use the S-radical words or the U tests, so an APN S8 in these slices would be among its T's. None is. The main run's zero, on these 221,377 sets, doesn't depend on those two stages.
- The reference, on all T-carrying sets and random others of each slice: exactly the same 288 on the 72 sets of class 21 and the same 5,796 on the 631 of class 20; none on 34 random other sets; and with the full test, none on any of the 737 sets, as the run found.
- Closure under Gamma, images whose form on R1 stays the slice's: class 21, 622,080 images; class 20, 374,328; **none outside the found set**, every found S8 hit. This is a recall test of the exact main-run path on the binary that ran, the five-solid analogue of the r = 4 closure above.

Worth Rowen's look, not a finding: in class 20 some of these have rank counts 186/64/5, five rank-4 forms and only 6 bad points, a quadratic function APN at all but 6 nonzero a. I looked at four slices, so none of this is known of the classes at large: T's in two slices of the two most structured classes, none in one slice each of classes 0 and 1. Why M stops at 6 is open.

**New: planted T in five-solid sets** (`plant5/`, the exact 61b1348 source with one patch letting a word hold fewer than 15 points). Real admissible S, a real plane of survivors, the best of 400,000 random third generators, and only the solid points where T fails dropped: 59 to 66 of 75 points stay, in all five words. **132 of 132 returned**, with 665 third-level nodes, across all 22 classes, 6 sets each. The reference on the same sets and point sets: f1, planes (3,215,465 in all) and T's (133, one more than planted) equal on all 132 lines.

**Pivot and prune on dim-4 sets.** The pivot has no independent check on a five-solid set in the write-up, and it is cheap to remove: without it six of the cheapest representatives took 78 to 254 s, against 72 to 120 s with it. Nine representatives (classes 2, 3, 12, 13) rerun with `fivecheck_ctl.exe` and no `pivot` (six), or a rebuild with the prune off as well, so a plain subspace search (four, one overlapping): no T, up to 110,460 three-spaces. Those last took up to about an hour, on capped E-cores while the box was shared, so timing there says little.

## 5. Findings, ranked by whether they could change the answer

**None could.** The ranked list is of what's wrong or unproven.

1. **§5's "obstruction on the solid points alone" doesn't generalise.** The basis was 122 sets, 1 in 300 of one slice in class 1. That slice in full does have none, but the first slices tried in classes 20 and 21 have thousands of T's. Nothing rests on it, but it reads as a claim about the mathematics in a text about to be published. Fix below.
2. **§6's "can't be constructed" is false**, for the same reason, and the "what the controls don't reach" paragraph should be replaced by §4 above.
3. **§6's table overclaims the r = 4 control** (unguided, "all stages", "all 26"). The text under the table is honest; the table isn't.
4. **§7's "on the code that ran"** is true only through the source diff: the two-direction flat check and the brute assembly exist only in later builds, and the r5 binary has only the first direction. I confirmed the diff touches nothing on the no-pivot path (debug print, `FC_NOASM`, `cur`, list mode). Say that.
5. **The checks' outputs weren't kept.** flat, brute assembly, slevel, gamma and reps left no files; only `plant`, the job logs and the runs did. I reran spot checks, so the claims stand, but the evidence should be files.
6. **A and B appear to share their first stages.** I didn't read A, but README's "five-form spans" for spread 0 is B's count, and `five_vs_b` applies B's rule to A's configurations. If A also reduces to classes of five solids and form tuples, its agreement says little about those stages, though §7 suggests comparing there first; §1 and §2 here carry them. Independence would begin at T.
7. **Small.** §1's missing injectivity line. "Class 0 and 1 of four solids weren't run beyond the one representative": class 0 wasn't run at all. §4's "cuts the first levels": the pivot acts at every level. VI.1 as printed says "APN functions" and takes "quadratic" from Thm IV.5; keep "quadratic" in any public statement, as §7 does.

**Suggested, in order of value:** run the remaining 94 dim-4 representatives with the pivot off (a rough 100 to 200 CPU-hours; the ratio is unmeasured on the costly ones; then the dim-4 leg rests on nothing but a plain subspace search; `runs\check5\review\noprune\bin_noprune.exe list 5 <c> <tuples> <dir>`, resumable); keep `nou` plus the closure test as a regression for the main path, and run it on more slices, since each is a second route to that slice's zero (about 3 CPU-hours a slice; all of it, up to about 1,800 CPU-hours, isn't warranted); save the check logs.

## 6. Corrected wording

§5, replacing "The obstruction is already there on the solid points alone in a sample: with only the 75 points as the test, 122 sets of class 1 gave planes and no T.":

> The solid points alone don't always exclude T. With only the 75 points as the test, full slices of classes 0 and 1 have none, but one of class 21 has 288 T's on 72 sets and one of class 20 has 5,796 on 631. All are non-APN, with from 6 to 93 points off the solids in three radicals. It's the test on all points that leaves none.

§6, table row 1: label it "r = 4, in part guided (below)"; reaches: "tuple list onward (`leaf`, not `partial`), sieve, assembly at e = 4, pivot, prune, verification; not the fifth word, the no-pivot branch or the `run` driver"; result: "12 classes in the 88 representatives taken in list order; the other 14 from 4 representatives chosen with hindsight; 26 of 26 in all; every S8 APN (177 of 177). Not a classification."

§6, replacing the paragraph "What the controls don't reach...":

> The stage of a T in a five-solid set is reached in two ways, neither in the main run itself, where no T survives the full test. With only the solid points as the test the run's binary finds 6,084 T's in two full slices of classes 20 and 21, all non-APN by definition, confirmed by independent code and closed under Gamma. And T's planted in real five-solid sets, on 59 to 66 of the 75 points, are returned.

§7, replacing "and each stage is tested against a definition or a known answer (§6), on the code that ran":

> and each stage is tested against a definition or a known answer (§6), on the run's source or on later builds that differ from it, on that path, only by checks and debug output.

§7, adding: "The classification and the admissible sets are common ground with search A, so A's agreement doesn't test them; they are checked separately in `review-five-check.md`."

## 7. Not checked

- The reference covered 5,981 sets of dim 5 at random and 703 T-carrying ones, not all 29.8M; the rest rests on §1 and on one code path having matched wherever compared. A whole slice in the reference itself I started and cancelled after 4.6 CPU-hours on capped E-cores, with the box shared; the three-slice chain above stands in for it. About half an hour with the box open: `bin_rvb5_ref.exe slice 0 3 4` from the repo root, against `r5_c00_q03.txt`.
- Dim-4 assembly on five-solid sets is B's code throughout; the reference checked only its sieve there. The pivot is checked by hand, by the r = 4 control and closure, and by nine of 103 representatives without it; the other 94 aren't.
- The reference is independent in code, not in idea: it imposes the same necessary conditions. Their necessity rests on §1 here.
- One machine, one compiler, one run: no rerun of the main run, so no guard against a hardware fault. Search A is the guard.
- Burnside counts 5,477 and 326, Kurz's table itself, and anything of search A's. The block view isn't load-bearing and I only reran it.

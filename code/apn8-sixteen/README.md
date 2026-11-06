# Quadratic APN functions on F₂⁸ have at most five rank-4 components: code and logs

The searches behind the finding "No quadratic APN function in dimension 8 has more than five rank-4 components"
([write-up](../../apn8-not-sixteen.md); Rowen's findings, Studio Quire). Each case was settled by two
searches, independent except as "Six" below says for `i = 11`: sixteen rank-4 components (`i = 1`), fifteen (`i = 2`), fourteen (`i = 3`), thirteen (`i = 4`), twelve
(`i = 5`), eleven (`i = 6`), ten (`i = 7`), nine (`i = 8`), eight (`i = 9`), seven (`i = 10`) and six (`i = 11`). Search A is Rowen's
(Claude Opus 5.5); search B is by separate Claude Opus 5.5 subagents (one for `i = 1, 2`, a fresh one each for `i = 3`,
`i = 4`, `i = 5`, `i = 6`, `i = 7`, `i = 8`, `i = 9`, `i = 10` and `i = 11`), each written without reading search A, with its own
write-up here (for `i = 7`, Rowen tallied search B's results and wrote its verdict, as that write-up says; for `i = 9`,
search B's author saw search A's verdict and configuration count in the repository's commit messages before starting,
as its write-up discloses; for `i = 10`, search B's author read none of search A's code or output and no commit history
beyond its own files, but saw search A's job names in the machine's shared job queue, and Rowen tallied search B's
coset sweep, which finished after the verdict; for `i = 11`, search B's author likewise read nothing of search A but saw
its job name in the queue, and search A's final run came after search B's verdict, on a hole stage Rowen rewrote then:
Six, below, says what that rewrite knew). The write-up has the reductions; here is how to re-run them. Eight's part of
this folder (its files, tables and re-runs) was put together by a separate Claude Opus 5.5 subagent for Rowen, and
seven's and six's by others.

| File | What |
|---|---|
| `src/sixteen.cpp` | Search A for `i = 1, 2, 3`: every spread set containing 0 and I (14-cliques), then a DFS assigning the 28 possible forms per radical. Modes: none (`i = 1`), `i2`, `i3`. |
| `src/extend6.cpp` (+ `src/apn8.cpp`) | Every APN 8-space containing a given 6-space, no caps; run on search A's `SHORT` lines for `i = 2`. |
| `src/cover.cpp` | Search A for `i = 4` to `i = 11`: every APN space whose rank-4 radicals are a given set of solids, rank-6 radicals in the holes; point branching, no symmetry used. Prints nodes by span dimension and `NEAR` lines (spans holding a form on every solid). Usage `cover <input> [maxlines\|-1] [splitdim]`: the split dimension D (added for `i = 10`) stops the tasks at every node of span dimension D that holds a form on every solid and leaves holes uncovered, then searches those nodes' subtrees in parallel: the same tree, the same counts. For `i = 11` the hole stage (below all solids) was rewritten: only forms passing a linear necessary condition reach the unchanged `push`, found by default as an affine-subspace tree (`holes_tree`); options `--gray` (the same filter by Gray code), `--old` (every line × 13,888 forms, as before), `--only TAG`, `--frange A B` (profiling). The same children, so the same counts; the file's header has the argument, and Build, below, the checks. |
| `src/cover_chunks.cmd`, `src/cover_tally.py` | Search A for `i = 9` and `i = 10`: the `cover3` build of `cover.cpp` run over `data/p8.txt` (and `data/p7.txt`) in 8 chunks, skipping finished ones on a restart (its paths are Rowen's working tree: `runs\bin\cover3.exe`, chunks `c0.txt` … `c7.txt` in the folder given), and one summary over the chunks' outputs. |
| `src/cover_chunks_split_exe.cmd` | Search A for `i = 11`: `cover_chunks.cmd` with the build given and a split dimension (`cover_chunks_split_exe.cmd DIR 5 cover5t` runs `runs\bin\cover5t.exe DIR\cK.txt -1 5`; the first run, on `cover4`, used a copy with the build fixed). |
| `src/cover_prep.py` | cover.cpp's inputs: `ext` (extendable case, 25,850 covers), `max` (Kurz's three hole sets, 11,004 covers), `ctl` (positive controls from a known function). |
| `src/i4_holes.py` | Checks Kurz's three hole sets for `i = 4` (`data/kurz13.txt`): distinct points, hyperplane counts, solids in the complement and exact covers by 13 of them, lines in the holes. |
| `src/spreads_to_cover.py` | Written for this folder: converts search B's 73 `i = 3` representatives (`data/spreads.txt`) to cover.cpp input; reproduces `data/i3-reps.txt` exactly. |
| `src/sixteen_check.cpp` | Search B for `i = 1, 2`: grows the span of forms point by point; no spread sets. Modes `selftest`, `main`, `i2`, `ctlB`, `i2ctl`. |
| `src/sixteen_check_level2.py`, `src/sixteen_check_level2_i2.py` | Numpy recounts of search B's first two levels for `i = 1` (5316 / 5280 / 5480) and `i = 2` (representatives 0, 1, 50, 120, 199: 1616 / 1659 / 1650 / 1493 / 1724). |
| `src/sixteen_check_short.py` | Checks search A's `i = 2` `SHORT` spans coordinate-free (rank counts, disjoint radicals, holes = two 4-spaces). |
| `src/fourteen_check.cpp` | Search B for `i = 3`: every 14-partial spread extends (`holes`), the 73 configurations up to GL(8,2) (`spreads`), point-branching search (`search`), `selftest`, known-function controls (`ctl`, `ctlbatch`). |
| `src/fourteen_check_levels.py` | Recount of the first levels of `fourteen_check search` for one configuration, sharing no code with it. |
| `src/thirteen_check.cpp` | Search B for `i = 4`: the 13-partial spreads classified (`holes`, `configs`: 190 configurations), Kurz's covers (`covers`), search (`search`, `searchsp`, `searchmem`), `selftest`, controls (`ctlbatch`, `ctlreps`), sizing for `i = 5` (`holesk`). |
| `src/thirteen_check_levels.py` | Recount of the first levels of `thirteen_check search` for one configuration, sharing no code with it. |
| `src/pspread.cpp` | Search A's classification for `i = 5` and `i = 6`: partial spreads of K solids up to GL(8,2) (`enum K`: canonical form over frames, extension class w, orbit-sum check; `selftest K file`). Run for K = 12 and 11, and for K = 17, 14, 13 as checks against the cases above. Output is cover.cpp input. |
| `src/holes_iso12.py` | The hole sets of search A's 12 maximal configurations up to GL(8,2), by explicit isomorphism (5 classes, as Kurz's table says). |
| `src/twelve_check.cpp` | Search B for `i = 5`: every 9-clique canonicalised (`parta`: 376 configurations, extension counts), hole sets and outside sources (`configs`), search (`search`, `searchsp`, `search13`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and the comparison with search A's list (`maprowen`). |
| `src/twelve_check_levels.py` | Recount of the first levels of `twelve_check search` for one configuration, sharing no code with it. |
| `src/levels_all12.py` | Runs `twelve_check_levels.py` to the bottom of the tree for all 376 configurations and compares every count with `twelve_check`'s. |
| `src/eleven_check.cpp` | Search B for `i = 6`: every 8-clique canonicalised, one minimal class per run (`parta`), the classes merged with per-orbit identities and the link to `twelve_check`'s 376 (`merge`: 817 configurations), hole sets, exact covers and outside sources (`configs`), search (`search`, `search12`, `search13`, `searchsp`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and, for comparison only, `orbitcov` (which configuration each listed clique lies on) and `mapa` (search A's list mapped onto search B's). |
| `src/eleven_check_levels.py` | Recount of the first levels of `eleven_check search` for one configuration, sharing no code with it. |
| `src/pspread3.cpp` | Search A's classification for `i = 7` to `i = 11`: `pspread.cpp`'s enumeration unchanged, plus per-task checkpoints (`--ckpt`, so an interrupted run resumes), the extension class w computed once per orbit (an orbit invariant) and spot-checked on every N-th clique (`--wcheck N`), and a time split. Run for K = 10, 9, 8, 7 and 6, and for K = 12 and 14 as checks against `pspread.cpp`'s outputs. Mode `holes K configs out` (added for `i = 9`): for each listed configuration, every exact cover of its covered points by K solids, canonicalised and looked up in the list. |
| `src/holes_tally.py` | Search A's hole sets from `pspread3 holes`: configurations grouped by shared hole set, `\|Aut H\|` from the covers, the classes and maximal ones counted (Kurz), and the `\|Aut H\|` distribution and number of point sets against Langevin and Hou's data file (downloaded, not included). |
| `src/ten_check.cpp` | Search B for `i = 7`, built on `eleven_check.cpp`: every 7-clique canonicalised, a minimal class per run, split into slices (`parta`), the slices merged with per-orbit identities and the link to `eleven_check`'s 817 (`merge`: 1,777 configurations), hole sets, exact covers and outside sources (`configs`), search (`search`, `search11`, `search12`, `search13`, `searchsp`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and, for comparison only, `orbitcov` and `mapa`. |
| `src/ten_check_levels.py` | Recount of the first levels of `ten_check search` for one configuration, sharing no code with it. |
| `src/ten_check_tally.py` | Tallies over `ten_check`'s outputs: configurations no control reached (`reps`), and the `orbitcov` counts over all of Part A's maximal cliques against Part A's own (`maxall`). |
| `src/near_misses10.py` | Written for this folder: search A's and search B's near misses for `i = 7` compared configuration by configuration, from the logs here. |
| `src/nine_check.cpp` | Search B for `i = 8`, built on `ten_check.cpp`: every 6-clique canonicalised, a minimal class per run, split into slices (`parta`), the slices merged with per-orbit identities and the link to `ten_check`'s 1,777 (`merge`: 5,663 configurations), hole sets, exact covers and outside sources (`configs`), search (`search`, `search10` to `search13`, `searchsp`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and, for comparison only, `orbitcov` and `mapa`. |
| `src/nine_check_levels.py` | Recount of the first levels of `nine_check search` for one configuration, sharing no code with it. |
| `src/nine_check_slices.py` | The driver that ran search B's jobs for `i = 8` (`parta`, `merge`, `configs`, `search`, `maxsample`, `cross`, `orbitcov`, `controls`, `ctlreps`): resumable, skipping finished slices and ranges. Its paths are those of Rowen's working tree (`nine-check.md`'s note maps them); the tables below give the commands it ran. |
| `src/nine_check_tally.py` | Tallies over `nine_check`'s `.prog` files (`prog`): each configuration's `R` line counted once, totals by type, distinct near-miss spans; used for `eight_check`'s `.prog` files too. `ten_check_tally.py` does the controls' and `orbitcov`'s tallies for nine and eight too. |
| `src/nine_check_compare.py` | Search A's and search B's near misses for `i = 8` compared configuration by configuration (`nine-check.md` §7), reading the logs here. |
| `src/nine_check_percfg.py` | Written for this folder by the reviewer: the same comparison for every (configuration, dimension) pair, assigning search B's `P1` spans to configurations by their radicals. |
| `src/nine_check_langevin.py` | Written for this folder: search B's hole sets for `i = 8` against Langevin and Hou's class count and stabiliser orders for 9-solid partial spreads (their data file is downloaded, not included). |
| `src/eight_check.cpp` | Search B for `i = 9`, built on `nine_check.cpp`: every 5-clique canonicalised, a minimal class per run, split into slices (`parta`), the slices merged with per-orbit identities and the link to `nine_check`'s 5,663 (`merge`: 9,372 configurations), the 5- and 6-cliques counted plainly with no symmetry (`count5`), hole sets, exact covers and outside sources (`configs`), search (`search`, `search9` to `search13`, `searchsp`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`; `ctlrot`: one control per configuration, functions in rotation), and, for comparison only, `orbitcov` and `mapa`. |
| `src/eight_check_slices.py` | The driver that ran search B's jobs for `i = 9` (`count5`, `parta`, `merge`, `configs`, `selftest`, `search`, `cross`, `maxsample`, `orbitcov`, `controls`, `ctlreps`, `ctlrot`, `numpy`): resumable, skipping finished slices and ranges. Its paths are those of Rowen's working tree (`eight-check.md`'s note maps them); the tables below give the commands it ran. |
| `src/eight_check_numpy.py` | Search B's second implementation of its whole search for `i = 9`, in Python with numpy, written from the mathematics rather than from the C++ (forms with a given radical as a null space, radicals of whole batches at once): node, near-miss and try counts must equal `eight_check search ... asc`'s on every configuration given; and planted controls with no classification involved (`control`). |
| `src/eight_check_mass.py` | Search B's exact masses for `i = 9`: hole classes per `\|Aut H\|` and the number of point sets against Langevin and Hou's `fix-8.txt` (downloaded, not included), and the labelled 8-partial spreads counted three ways. This folder's copy reads the files here. |
| `src/eight_check_compare.py` | Search B's own cross-checks for `i = 9`: near misses across its three point orders (`orders`), and the numpy re-search against the C++ (`numpy`). This folder's copy reads the files here. |
| `src/eight_check_vs_a.py` | Search A's and search B's near misses for `i = 9` compared (`eight-check.md` §7), by configuration and by containment. This folder's copy reads the files here. |
| `src/eight_check_percfg.py` | Written for this folder: `nine_check_percfg.py` carried over to `i = 9`, comparing search B's near-miss spans with search A's `NEAR` lines for every (configuration, dimension) pair, in all three orders. |
| `src/seven_check.cpp` | Search B for `i = 10`, built on `eight_check.cpp`: every 4-clique canonicalised, a minimal class per run, split into slices (`parta`), the slices merged with per-orbit identities and the link to `eight_check`'s 9,372 (`merge`: 3,731 configurations), the 1- to 5-cliques counted plainly and the 4-cliques again by edges, with no symmetry used (`count4`), hole sets, exact covers and outside sources (`configs`), search (`search`, with `phase1` for the part off the holes alone, `cores` to write the cores, and `dsweep` for the coset sweep in the holes; `search8` for eight_check's configurations with a member removed, `searchsp` for spread representatives, unused here), `selftest`, controls (`ctl`; `ctlbatch` with `nond`, `nond2`, `guided` and `dsweep`; `ctlrot`), and, for comparison only, `mapa`. The environment variable `SEVEN_NOCACHE` turns off the point branching's candidate caches (below). |
| `src/seven_check_slices.py` | The driver that ran search B's jobs for `i = 10` (`count4`, `parta`, `merge`, `configs`, `selftest`, `search`, `crosssample`, `controls`, `ctlrot`, `numpy`): resumable, skipping finished slices and ranges. Its paths are those of Rowen's working tree (`seven-check.md`'s note maps them); the tables below give the commands it ran. |
| `src/seven_check_numpy.py` | Search B's second implementation for `i = 10`, in Python with numpy (`eight_check_numpy.py` carried over to seven members, with the cores and their fingerprint): node, near-miss, try and core counts and the fingerprint must equal `seven_check search ... asc`'s on every configuration given; and planted controls with no classification involved (`control`). |
| `src/seven_check_mass.py` | Search B's exact masses for `i = 10`: hole classes per `\|Aut H\|` and the number of point sets against Langevin and Hou's `fix-7.txt` (downloaded, not included), and the labelled 7-partial spreads counted three ways. This folder's copy reads the files here. |
| `src/seven_check_tally.py` | Tallies over `seven_check`'s outputs: every configuration once, totals by type (`prog`); cores by dimension and fingerprint across runs, and the coset sweep against `asc` (`orders`); the numpy re-search against `asc` (`numpy`); configurations no control reached (`reps`). This folder's copy reads the files here. |
| `src/seven_check_vs_a.py` | Search A's and search B's classifications and near misses for `i = 10` compared configuration by configuration (`seven-check.md` §7). This folder's copy reads the files here. |
| `src/six_check.cpp` | Search B for `i = 11`, built on `seven_check.cpp`: every 3-clique canonicalised, one slice per minimal class (`parta`), the slices merged with per-orbit identities and the link to `seven_check`'s 3,731 (`merge`: 341 configurations), the 1- to 4-cliques counted plainly and by edges with no symmetry used (`count3`), hole sets, exact covers and outside sources (`configs`), the point branching (`search`, with `phase1` and `cores` as in `seven_check`), a second search written new, sharing only the input masks and the layout of a form (`lin`, with `desc`, `phase1`, `cores`), `selftest` and `linselftest`, controls (`ctlbatch` with `nond`, `nond2`, `nond3` and `lin`; `ctlrot` with `lin`; `plantfind` and `ctlplant`, the planted controls), and, for comparison only, `mapa`. |
| `src/six_check_slices.py` | The driver that ran search B's jobs for `i = 11` (`count3`, `parta`, `merge`, `configs`, `selftest`, `search`, `lin`, `controls`, `ctlplant`, `ctlrotlin`): resumable, one process per configuration in the searches. Its paths are those of Rowen's working tree (`six-check.md`'s note maps them); the tables below give the commands it ran. |
| `src/six_check_mass.py` | Search B's exact masses for `i = 11`: hole classes per `\|Aut H\|` and the number of point sets against Langevin and Hou's `fix-6.txt` (downloaded, not included), and the labelled 6-partial spreads counted three ways. This folder's copy reads the files here. |
| `src/six_check_tally.py` | Tallies over `six_check`'s outputs: every configuration once, totals by type (`prog`); core files compared as sets per configuration, with `lin`'s seventh-solid cores set aside by Python's own rank (`cores`); `lin`'s nodes less those against the point branching's (`nodes`). |
| `src/six_check_vs_a.py` | Search A's and search B's classifications and near misses for `i = 11` compared configuration by configuration (`six-check.md` §7). This folder's copy reads `logs/six-near.txt` for the spans, which are too large to include. |
| `src/six_near.py` | Written for this folder: reduces a run's near-miss spans (search A's `NEAR` lines, search B's `P1` lines) to a count by dimension and a digest of the set, per configuration; made `logs/six-near.txt`. |
| `sixteen-check.md`, `fourteen-check.md`, `thirteen-check.md`, `twelve-check.md`, `eleven-check.md`, `ten-check.md`, `nine-check.md`, `eight-check.md`, `seven-check.md`, `six-check.md` | Search B's write-ups for `i = 1, 2`, `i = 3`, `i = 4`, `i = 5`, `i = 6`, `i = 7`, `i = 8`, `i = 9`, `i = 10` and `i = 11`: the reduction re-derived, method, results, controls, and (at the end) the comparison with search A. |
| `data/` | Small inputs (below). |
| `logs/` | The runs' logs and outputs (below); `logs/rerun/` holds the re-runs made while preparing this folder. |
| `compat/` | Shims mapping the MSVC intrinsics onto gcc/clang builtins (from `../apn8-symmetry`, with `__popcnt64`, `_tzcnt_u64`, `_BitScanReverse`, `_pdep_u32` and, for `seven_check.cpp`, `_pext_u64`, for `six_check.cpp`, `_BitScanReverse64` and `_aligned_malloc` and, for `cover.cpp`, `_lzcnt_u32` added for the sources here). |

## Build

Windows (MSVC, from a Developer Command Prompt): `cl /O2 /openmp /EHsc /std:c++17 /arch:AVX2 src\sixteen.cpp /Fe:sixteen.exe`,
the same for `extend6.cpp` (which includes `apn8.cpp`), `cover.cpp`, `sixteen_check.cpp`, `fourteen_check.cpp`,
`thirteen_check.cpp`, `pspread.cpp`, `twelve_check.cpp`, `eleven_check.cpp`, `pspread3.cpp`, `ten_check.cpp`,
`nine_check.cpp`, `eight_check.cpp`, `seven_check.cpp` and `six_check.cpp`. The runs' binaries were built with
`/arch:AVX2`; the sources build without it too, though `six_check.cpp` calls AVX2 intrinsics directly, so its binary needs
a CPU that has them.

Linux / macOS on x86-64 (g++ or clang++ with OpenMP; untested, see `../apn8-symmetry/README.md` on the shims; the
shims added here were checked against the MSVC intrinsics on 20 million random inputs, `_pext_u64`'s on 2026-11-01 and
`_BitScanReverse64`'s and `_lzcnt_u32`'s on 2026-11-04, under MSVC, with `__lzcnt64` and `__lzcnt` for the `__builtin_clz`
family, as this machine has no gcc):

    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/sixteen.cpp -o sixteen

and the same for the other fourteen, with `-mavx2` for `six_check.cpp`. The Python scripts need Python 3 (`sixteen_check_level2*.py`,
`eight_check_numpy.py` and `seven_check_numpy.py` also numpy).

`src/sixteen.cpp` is the current source. The runs used builds with 32-bit `long` counters (MSVC), which is all that
changed since: the `i = 3` run's node count wrapped around and is printed negative in its log (a statistic only; the
spread-set, solution and SHORT counts are far below 2³¹).

`src/cover.cpp` is the current source; it has had two edits since the build that ran most of these, so what's said below
about it is as of the date given. The build that ran twelve, eleven, ten, nine and eight and most of seven (`cover3`) is
this source without those edits, the optional split dimension (Seven, below) and the rewritten hole stage (Six, below);
`cover3` adds the per-dimension node
counts and the `NEAR` lines to the one that ran thirteen, and those only report: on search B's 73 `i = 3`
representatives it visits the same 2,284,642 nodes as before. `src/pspread.cpp` is the current source; the K = 12, 13
and 14 runs used a build from before one later edit, which changed the cap on w (the largest clique in the common
neighbourhood) from a fixed 5 to 17 − K, for K below 12. A partial spread of K solids takes at most 17 − K more, so for
K ≥ 12 both caps are at least the largest possible w and the output is the same. The K = 11 run used the later build:
there w goes up to 6 (an 11 that extends to a spread). `src/pspread3.cpp` ran K = 10, 9 and 8; on K = 12 and 14 its
output is byte for byte `data/p12.txt` and `data/p14.txt` (Ten, below). Its `holes` mode was added afterwards, for
eight (below).

Which source built which binary was checked by rebuilding with the Build line's flags (`/arch:AVX2` included) and
comparing disassembly (2026-10-26, and for eleven 2026-10-27): `src/cover.cpp` gives the build that ran twelve and
eleven, instruction for instruction; `src/pspread.cpp` the later build (the one that ran K = 11), and `pspread.cpp` as
it stood before that edit the build that ran K = 12, 13 and 14; `src/twelve_check.cpp` the last of search B's three
twelve builds (the earlier two differ only by the added modes, as its write-up says); `src/eleven_check.cpp` the last
build of search B for eleven, which added the comparison-only modes `orbitcov` and `mapa` after the classification,
searches and controls had run (the source's history shows no other change since its first commit, made while Part A
was running). MSVC's code generation depends on the source file's name, so the comparison only matches when the file
is compiled under its own name. For `eleven_check.cpp` it also needed the original's compile line, from the source's
own folder with the object file named `eleven_check_` (`/Fo<dir>\eleven_check_`); compiled as `src\eleven_check.cpp`,
or to `eleven_check.obj`, the functions come out in a different order (the same instructions, counted by opcode).

For ten (2026-10-28) the binaries' object files were compared, function by function with the order ignored and the
compiler's hashed lambda names set aside (they vary from build to build): `src/pspread3.cpp` gives the build that ran
K = 10, 12 and 14, every function identical. `ten_check.cpp` changed once after its first commit (made as Part A was
submitted): the comparison-only `mapa` mode, which came from `eleven_check` reading 11 members, was changed to read
10 for the comparison with search A. The first commit's source gives the build that ran everything else (Part A to the
controls), and `src/ten_check.cpp` the rebuild that ran `mapa`, every function identical in both cases; the
executable that ran the rest was overwritten by that rebuild, its object file kept. In each comparison the only
difference is the size of a debug section, which records the build folder's path.

Nine (2026-10-29) ran on the same `cover3` and `pspread3` binaries as ten. `nine_check.cpp` has changed since the build
that ran search B for nine only in a comment (on the line where `mapa` reads search A's list); compiled under its own
name with that build's object name, `src/nine_check.cpp` gives its object file, every function identical, with the
same debug-section difference.

Eight (checked 2026-10-30) ran on the same `cover3` binary and, for `pspread3 enum 8`, the same `pspread3` binary as
ten and nine. `src/pspread3.cpp` has changed once since that build: the `holes` mode, added for eight, which inserts
its block into `main` before the enumeration and two `#include`s (`<set>`, `<array>`), and changes no other line. The
`holes` runs used a build of the current source (`pspread3h`). Compiled under their own names with the original
compile line, `src/pspread3.cpp` gives `pspread3h` and `src/cover.cpp` gives `cover3` with their code sections byte
for byte (the only differences are the hashed lambda names in the data section and the debug directory's timestamp).
`eight_check.cpp` changed once after its first commit: the `ctlrot` mode was added, with `g_forceRep` made
`thread_local` for its parallel loop. The first commit's source gives the build that ran everything else (Part A to the
controls, the cross-runs and the sample), and `src/eight_check.cpp` the second build (`eight_check_rot`, the same source
under another name, built while the first was busy with the `search9` cross-run), which ran `ctlrot` and `mapa`. As the
first build's object file was overwritten by a later compile, these were compared by the executables' disassembly,
with branch and call targets set aside: the same instructions in both cases, with whole functions placed in a
different order (the hashed lambda names again).

Seven (checked 2026-11-01) ran `pspread3 enum 7` and `holes 7` on `pspread3h` (`src/pspread3.cpp` is unchanged since
eight), seven of its eight cover chunks on `cover3`, and the rest on `cover4`, a build of `src/cover.cpp` as it then was.
The edit between the two adds the split dimension (the `splitdim` argument, each task's frontier list, a second
parallel loop over the frontier, the progress line's lock now taken for every task to append the frontier) and leaves
the search alone. It was checked on eleven's configurations before seven's runs, `cover4 data/p11.txt -1`, `-1 5` and
`-1 6`: 30,505,655 nodes, identical at every dimension, and the same 58 `NEAR` lines each time (at dimension 5 no node
has all eleven rank-4 forms, so the split at 5 finds no frontier there; at 6 it finds 22;
`logs/cover4-regress-{nosplit,s5,s6}-output.txt`), and again on one of seven's chunks for this folder, unsplit and split
at 5 (Re-run, below). Compiled under their own names with the Build line, `src/cover.cpp` gives `cover4`'s object file
and `src/pspread3.cpp` `pspread3h`'s, every function identical (compared as for ten, function by function, the order
and the hashed lambda names set aside).

`seven_check.cpp` was committed when search B started and changed six times as it ran, each change adding something,
so its runs came from several builds of the one source (`seven-check.md`'s State names them by commit; the commits
are in Rowen's own repository, which isn't public). In order: the first build (its source committed three minutes after its first jobs started) ran
`count4`, `parta`, `merge` and `configs` and configurations 0 to 10 of `search ... asc` (a first `asc` job, cancelled after ten minutes; the next
resumed after those 11); the second (a faster radical, a Pfaffian test in the sweep, the cores' D phases run as a
second parallel stage) ran `selftest`, the `ind-bl22` control batch and `asc` and `desc` up to configuration 1526; the
third (the sweep generating its candidates as K_a below dimension 7) the coset sweep up to configuration 134; the
fourth (guided and two-step controls; `seven_check_ctl.exe`) every other control batch and `ctlrot`; the fifth (the
`phase1` option; `seven_check_p1.exe`) the core census; the sixth (`SEVEN_DMAX`, a cost probe whose run was cancelled
and isn't here); the seventh (the candidate caches, `SEVEN_NOCACHE` to turn them off in the point branching) `asc`
and `desc` from configuration 1527 on, the sweep from 135 on and the cross-run on eight's deletions; and the current
source, which adds the comparison-only `mapa` mode after both verdicts, ran `mapa`. The caches change the tries in D
(the sweep's, its count of element checks) and nothing else: configuration 41, run with and without them, visits the
same spans (`seven-check.md` §5), and so does every configuration re-run for this folder with the current source
(below; `SEVEN_NOCACHE` set for the point branching's ranges that originally ran without; the sweep's cache has no
switch). Compiled under its own name with the Build line, the seventh state gives the object files of
the builds that ran the late `asc` and `desc` ranges and the sweep, the fourth that of `seven_check_ctl.exe`, and
`src/seven_check.cpp` that of the build that ran `mapa`, every function identical; the fifth gives `seven_check_p1.exe`,
compared by disassembly with addresses set aside (its object file wasn't kept). The builds of the first three states
were overwritten by later ones, so for those the check is the re-runs: the current source reproduces their records
exactly (below), apart from the sweep's count of element checks on configuration 1, the one re-run sweep
configuration with a core of dimension 5, where the cache lowers it (920,973,542 then, 242,798,460 now).

Six (checked 2026-11-04) ran `pspread3 enum 6` and `holes 6` on `pspread3h` (`src/pspread3.cpp` is unchanged since eight),
the first cover run on `cover4` and the final run on `cover5t`, a build of `src/cover.cpp` as it is now. The edit between
the two (`cover.cpp`'s header has the argument) leaves `push`, which checks a candidate against the whole span, and the
part above the holes alone, and puts a filter before `push` in the hole stage: for a covered point y, a legal new form B
has B(y,.) outside F_y = {s(y,.) : s in the span}, since otherwise y lies in the radicals of two elements. That is a
necessary condition, and linear in B, so over the 21-dimensional space of forms with the branching point in their radical
the forms that pass can be found directly, as a tree of affine subspaces (`holes_tree`, the default); `--gray` enumerates
the same space by Gray code, and `--old` keeps the original lines × 13,888 forms. In six only about 1 candidate of 2
million passes. Compiled under their own names with the Build line, `src/cover.cpp` gives `cover5t`'s object file, the
source before this edit (its commit 5265bc8, for seven) `cover4`'s, and `src/pspread3.cpp` `pspread3h`'s, every function
identical (compared as for ten, function by function, the order and the hashed lambda names set aside). `cover5`, the
build that ran the Gray filter's first profile and a regression (`logs/cover5-*`), is a state of the source between the
two that wasn't committed; `src/cover.cpp --gray` gives its counts (Re-run, below).

`six_check.cpp` was committed when search B started and changed eight times, each change adding something, so
its runs came from several builds of the one source (`six-check.md`'s State names them by commit; the commits are in
Rowen's own repository, which isn't public). Built from each commit under its own name with the Build line and compared
as above with the object files kept: `six_check.exe` (classification, `selftest`, the off-D phases, the `asc` sample) is
the first commit (9424d74); `six_check_lin.exe` (`lin`), byte for byte also `six_check_ctl.exe` (the planted and in-D
controls), is the third (ea1f510); `six_check_ls.exe` (`linselftest`), byte for byte a build whose object was kept, the
fourth (52b7bf0); `six_check_pf.exe` (`plantfind` for the other lists) the fifth (373e060); `six_check_ctl2.exe` (`lfull`)
the sixth (706afa7); `six_check_lin2.exe` (`lindesc`) the seventh (e815bbc); `six_check_ctl3.exe` (`ctlrot`) the eighth
(c896a97); and `src/six_check.cpp`, which adds the comparison-only `mapa` mode after both verdicts, `six_check_map.exe`:
every function identical. Two builds match no commit. `six_check_t.exe` ran the timing runs on configuration 340, `asc340`
and `lin340b`; State places it between the second and third commits, but its object file is from between the first and
second, and differs from the first's only in `lin`'s functions, so its point branching is `six_check.exe`'s. (`lin340`,
on `six_check.exe`, and `lin340b` were timing runs of earlier `lin` builds: `logs/check6-lin340-L.txt` has the full run's
nodes on 340, with fewer candidates counted at dimension 5.)
`six_check_c.exe` wrote `plants-bl22-c1.txt`; its object file wasn't kept, its instructions are nearest the second commit's
build (about 650 of 121,040 differing either way), and the current source writes the same list (Re-run, below). The comparisons of object files, each kept object against every fresh build: `logs/rerun/six/objmatrix.log`.

Search B's write-ups mention successive builds (`fourteen_check2.exe` and so on): each is a build of the one source in
`src/`, whose later edits added control modes and left the search path alone (the write-ups say what changed). The
logs show which build each run used.

## Runs and expected output

Times are wall-clock on a shared 20-thread desktop (Core Ultra 7 265) under a job runner that gave way to other work,
often to the E-cores only, so an idle machine is faster. Logs are in `logs/`. The job logs (`2026…-<name>.log`) are
that runner's; its own lines (command, pauses, timing) are dropped from the published copies, so for runs whose
program output went to a file (Rowen's, and search B's for twelve), the job log is nearly empty and the output is the
`*-output.txt` beside it.

### Sixteen (`i = 1`)

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen` | `spread sets (cliques) 1128864, semifield 19936, field 336 (expect 336)`; `solutions 0, short 0` | 57 min | `logs/20261024-141134-sixteen.log`, `logs/sixteen-output.txt` |
| `sixteen_check selftest` | all pass (28 / 13,888 nondegenerate forms, subspace counts, 336 subfields) | seconds | in `sixteen-check.md`; `logs/rerun/s16-selftest.log` |
| `sixteen_check main main.prog` | `solutions ...: 0` (3 starting pairs) | ~8 min | `logs/20261024-150959-check16-main.log` |
| `sixteen_check main main-desc.prog desc` | `0` (reversed point order) | ~8 min | `logs/20261024-150959-check16-main-desc.log` |
| `sixteen_check main main-all.prog all` | `0` (all 28 starting forms, no symmetry reduction) | ~60 min | `logs/20261024-150959-check16-main-all.log` |
| `sixteen_check ctlB apn_8bit.txt 12631 3` | `ctlB: 1 of 1 verified APN by LUT` (a known function re-found) | minutes | `logs/20261024-151000-check16-ctlb3.log` |
| `python src/sixteen_check_level2.py` | 5316 / 5280 / 5480 | seconds | in `sixteen-check.md` |

### Fifteen (`i = 2`)

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen i2 > fifteen.txt` | `field 336 (expect 336)`; `solutions 0, short 36` (36 `SHORT dim 6` lines) | ~3 h | `logs/20261024-163557-fifteen.log`, `logs/fifteen-output.txt` |
| `extend6 fifteen.txt 13` | `cases 36, skipped 0, APN extensions 0` | seconds | `logs/fifteen-ext-output.txt` |
| `sixteen_check i2 i2.prog 0 199` | `solutions: 0` (200 orbit representatives) | ~4 min | `logs/20261024-164052-check16-i2.log` |
| `sixteen_check i2 i2-desc.prog 0 199 desc` | `0` (reversed point order) | ~10 min | `logs/20261024-164517-check16-i2-desc.log` |
| `sixteen_check i2 i2-nopr.prog 0 199 nopr` | `0` (orbit-index pruning off) | ~20 min | `logs/20261024-164517-check16-i2-nopr.log` |
| `sixteen_check i2 i2-p1.prog 0 199 p1` | `0`; logs the 6 near-miss spaces (every non-hole point covered, no rank-6 form fits) | ~4 min | `logs/20261024-170915-check16-i2-p1.log` |
| `sixteen_check i2ctl apn_8bit.txt 15 3` / `... 12631 3` | `i2ctl: 1 of 1 verified APN by LUT` (known functions re-found) | ~10 min | `logs/*-check16-i2ctl3-*.log` |

### Fourteen (`i = 3`)

Search A, `sixteen.cpp` (holes Y, X and each of the 14 other graphs in turn; no orbit reduction), and cover.cpp as a
check on cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen i3 > fourteen.txt` | `holes 3 (14 rank-4 components)`; `spread sets (cliques) 1128864, semifield 19936, field 336 (expect 336)`; `solutions 0, short 0` | 15.3 h | `logs/20261024-165047-fourteen.log`, `logs/fourteen-output.txt` |
| `cover data/i3-reps.txt` | `done: covers 73, tasks 2044, nodes 2284642, solutions 0, full (dim < 8) 0` (search B's 73 representatives; uses B's reduction, so not independent of it) | 19 s | `logs/20261024-215424-i3-cover-reps.log`, `logs/i3-reps-output.txt` |

Search B, `fourteen_check` (`fourteen-check.md`). It writes `.prog` files to resume from; give each run a fresh one,
since an existing one makes it skip the configurations listed there.

| Command | Expected | Time | Log |
|---|---|---|---|
| `fourteen_check holes 11 holes11.txt` | `total 11-cliques (weighted count): 410906496.000; reduced leaves 863186; without extension: 0` (every partial spread of 14 solids extends) | 4 s | `logs/20261024-195125-check14-holes11.log` |
| `fourteen_check spreads spreads.txt` | `73 G-orbits; sum of orbit sizes 1128864 (weighted count 1128864.000); 1 orbit reps are fields`; writes `spreads.txt` (= `data/spreads.txt`) and `spreads.txt.all` (the 2,873 conjugation-reduced cliques) | 40 s | `logs/20261024-194724-check14-spreads.log`, `logs/20261024-195500-check14-spreads2.log` |
| `fourteen_check selftest spreads.txt` | `selftest: 73 configurations, 16694216 generated forms checked, 0 failures` | 40 s | `logs/20261024-200618-check14-selftest.log` |
| `fourteen_check search spreads.txt asc.prog asc` | `solutions: 0`; `search nodes by dim: 73 2044 57232 740398 1196814 64384 845 0 0`; no node with every off-D point covered | 23 s | `logs/20261024-195115-check14-asc.log` |
| `fourteen_check search spreads.txt desc.prog desc` | `solutions: 0` (descending point order) | 25 s | `logs/20261024-195837-check14-desc.log` |
| `fourteen_check search spreads.txt rand7.prog rand 7` | `solutions: 0` (random point order) | 26 s | `logs/20261024-200618-check14-rand.log` |
| `fourteen_check search spreads.txt.all noS3.prog asc` | `solutions: 0` (all 2,873 conjugation-reduced cliques, no S3 reduction) | 14 min | `logs/20261024-200446-check14-nos3.log` |
| `python src/fourteen_check_levels.py spreads.txt 0 4` | `[1, 28, 784, 10234, 15917]`, the C++ counts for spread 0 | 12 s | `logs/20261024-200202-check14-levels0.log` |
| `fourteen_check ctl spreads.txt apn_8bit.txt 12631 0 0` | `ctl: 1 spaces found ...; known space among them: YES` | seconds | `logs/20261024-195459-check14-ctl-12631.log` |
| `fourteen_check ctlbatch spreads.txt apn_8bit.txt 2 100000000 3 ctl-bl22.txt` | `1440 controls; re-found 1440; NOT found 0; ...; other spaces found 0` | 6 min | `logs/20261024-195740-check14-ctlbatch-bl22.log`, `logs/check14-ctl-bl22.txt` |
| `... apn_8bit.txt 3 100000000 24 ctl-wide.txt randspread` | `4128 controls; re-found 4128`; `spread reps hit: 68 of 73` | 23 min | `logs/20261024-200812-check14-ctl-wide.log`, `logs/check14-ctl-wide.txt` |
| `... blp22_8bit.txt 2 100000000 6 ctl-blp22.txt` | `222 controls; re-found 222` | 1 min | `logs/20261024-200446-check14-ctl-blp22.log` |
| `... new_apns.txt 2 100000000 1 ctl-new.txt` | `7199 controls; re-found 7199` | 62 min | `logs/20261024-200619-check14-ctl-new.log` |
| `... yuperrin_qam.txt 2 100000000 6 ctl-yuperrin.txt` | `24 controls; re-found 24` | 1 min | `logs/20261024-200618-check14-ctl-qam2.log` |
| `... yuwangli_qam.txt 2 100000000 6 ctl-yuwangli.txt` | `78 controls; re-found 72; NOT found 0; ...; setup failures 6` (exit code 3): the 6 failures are the list's line 3, which isn't a quadratic APN function as read | 1 min | `logs/20261024-203006-check14-ctl-qam1b.log`, `logs/check14-ctl-yuwangli.txt` |

### Thirteen (`i = 4`)

Search A, cover.cpp: the hole sets of a non-extendable 13-partial spread are Kurz's three (*Divisible codes*,
arXiv 2112.11763, §9; search B proves this independently), and the extendable case has holes Y, X, graph(I),
graph(T), one T per conjugacy class:

| Command | Expected | Time | Log |
|---|---|---|---|
| `python src/i4_holes.py data/kurz13.txt` | for the three sets: `60 distinct`, `points per hyperplane: {28: 195, 36: 60}`, `solids inside the complement: 1264` / `1280` / `1280`, `exact covers ...: 552` / `10440` / `12`, `solids inside the holes: 0`, `lines inside the holes: 200; partitions of the holes into 20 lines: 1000000` (the count stops at that cap) | 12 min, nearly all in the capped count | `logs/rerun/i4-holes.log` |
| `python src/cover_prep.py max max.txt data/kurz13.txt` | `config 0: 552 covers`, `config 1: 10440 covers`, `config 2: 12 covers`, `lines 11004` | 1 min | `logs/rerun/prep-max.log` |
| `cover max.txt` | `done: covers 11004, tasks 308112, nodes 374898296, solutions 0, full (dim < 8) 0, dim 8 non-APN 0` | 59 min | `logs/20261024-201233-i4-max.log`, `logs/i4-max-output.txt` |
| `python src/cover_prep.py ext ext.txt` | spread sets containing 0, I, T: `15852` / `2746` / `1760` / `2746` / `2746` (minimal polynomials 111, 10011, 10101, 11001, 11111); `lines 25850`. Weighted by class size (112, 1344, 1680, 1344, 1344) that is 15,804,096 = 14 × 1,128,864 | 22 min | `logs/20261024-194956-i4-ext-prep.log` |
| `cover ext.txt` | `done: covers 25850, tasks 723800, nodes 768769492, solutions 0, full (dim < 8) 0, dim 8 non-APN 0` | 1 h 57 min | `logs/20261024-201233-i4-ext.log`, `logs/i4-ext-output.txt` |
| `cover data/ctl-12631-k4.txt` (`k3`, `k2`) | `SOLUTION ctl-12631 dim 8 ranks 8:186 6:65 4:4 ...`; `solutions 1`: the known function (Beierle–Leander line 12631) re-found as the only space, seeded with one rank-4 form and 4, 3 or 2 forms completing a basis, so cover.cpp finds the other rank-4 forms and none, one or two rank-6 levels itself. The inputs come from `python src/cover_prep.py ctl OUT apn_8bit.txt 12631 1 KEEP` | 1 s | `logs/rerun/cover-ctl-k*.log` |

Search B, `thirteen_check` (`thirteen-check.md`; fresh `.prog` files as for fourteen):

| Command | Expected | Time | Log |
|---|---|---|---|
| `thirteen_check holes holes10.txt` | `reduced leaves 2358714 (non-extendable 141165, partial 0)`; `10-cliques, weighted: extendable 1129992864, non-extendable 62558496, common nbhd nonempty but no 4-clique 0` | 8 s | `logs/20261024-220843-check13-holes.log` |
| `thirteen_check configs holes10.txt data/kurz13.txt configs.txt` | Kurz's sets with automorphism groups of order 576, 14400, 96 and `552`, `10440`, `12` covers; `configurations up to GL(8,2): 159 extendable, 31 non-extendable`; both counts `MATCH`; writes `configs.txt` (= `data/configs.txt`) | 11 min | `logs/20261024-221011-check13-configs.log` |
| `thirteen_check selftest configs.txt` | `selftest: 190 configurations, 33400360 generated forms checked, 0 failures` | 2.5 min | `logs/20261024-222133-check13-selftest.log` |
| `thirteen_check search configs.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 190 5320 148960 1972468 3468491 783269 4006 0 0`; 4 nodes with every off-D point covered, all at dim 6 (configurations 45 and 154), dead at the first hole point | 111 s | `logs/20261024-222133-check13-asc.log` |
| `thirteen_check search configs.txt desc.prog desc` | `solutions: 0`; the same 4 near-miss nodes | 2 min | `logs/20261024-222637-check13-desc.log` |
| `thirteen_check search configs.txt rand7.prog rand 7` | `solutions: 0`; the same 4 | 2 min | `logs/20261024-222638-check13-rand.log` |
| `thirteen_check searchsp data/spreads.txt sp14.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 14 fourth holes (1,022 extendable configurations, without this program's classification) | 12 min | `logs/20261024-222638-check13-sp14.log` |
| `thirteen_check covers data/kurz13.txt covers.txt`, then `thirteen_check searchmem covers.txt covers.prog asc` | `Kurz set 1: 552 covers` (`10440`, `12`); then `solutions: 0` over all 11,004 Kurz covers unreduced, nodes by dim `11004 308112 8627136 115373106 207621076 42791038 177828 0 0`: 374,909,300 in all, cover.cpp's 374,898,296 plus the 11,004 roots | 53 min | `logs/20261024-224107-check13-covers.summary.log` |
| `python src/thirteen_check_levels.py configs.txt 0 3` / `... 170 3` | `[1, 28, 784, 10573]` / `[1, 28, 784, 10288]`, the C++ counts for configurations 0 and 170 | seconds | `logs/rerun/f13-levels0.log`, `logs/rerun/f13-levels170.log` |
| `thirteen_check ctlbatch configs.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440`, `nonext: 1440 controls; re-found 1440`, `other spaces found 0` | 20 min | `logs/20261024-222757-check13-ctl-bl22.log`, `logs/check13-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again, exit code 3), on both paths | 7 min | `logs/20261024-222758-check13-ctl-small.log`, `logs/check13-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `7199 controls; re-found 7199` on both paths | 69 min | `logs/20261024-222803-check13-ctl-new.log` |
| `thirteen_check ctlreps configs.txt apn_8bit.txt 1 40 50 122 ... 158 45` | `40 of 40 re-found, 0 other spaces` for each of the 19 configurations the batches missed | 2 min | `logs/20261025-002003-check13-ctlreps.log` |
| `thirteen_check holesk 10`, `holesk 9` | sizing for `i = 5` (`thirteen-check.md` §8) | seconds | `logs/20261024-231121-check13-size10.log` |

### Twelve (`i = 5`)

Kurz's table gives only counts for twelve solids (42 hole sets, 5 maximal), so both searches classify the 12-partial
spreads themselves. Search A, `pspread.cpp` then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread enum 12 p12.txt` | `w=0: reduced 100783, weighted 49896000, orbits 12, orbit sum 49896000 MATCH`; `w=1: reduced 1040086, weighted 502840800, orbits 84, ... MATCH`; `w=5: reduced 4128284, weighted 2258825184, orbits 280, ... MATCH`; `orbits 376`; writes `p12.txt` (= `data/p12.txt`) | 72 min | `logs/20261025-082434-pspread12.log`, `logs/pspread12-output.txt` |
| `pspread enum 13 p13.txt`, `enum 14 p14.txt` | K = 13: `w=0: reduced 141165, weighted 62558496, orbits 31`, `w=4: reduced 2217549, weighted 1129992864, orbits 159` (thirteen-check's numbers); K = 14: `w=3: reduced 863186, weighted 410906496, orbits 73` (fourteen-check's); writes `data/p13.txt`, `data/p14.txt` | 60 min, 38 min | `logs/20261025-082433-pspread1{3,4}.log`, `logs/pspread1{3,4}-output.txt` |
| `pspread enum 17 p17.txt` | `orbits 8` (the 8 spreads of PG(7,2), the translation planes of order 16; the Desarguesian one with `\|Aut P\|` 244,800), 2,873 reduced, orbit sum 1,128,864 (the spread sets of `sixteen`) | 21 s | `logs/rerun/pspread-p17.log`, `logs/rerun/p17.txt` (the original run's log wasn't kept; this is the re-run's) |
| `pspread selftest 12 data/p12.txt` (`13 data/p13.txt`, `14 data/p14.txt`) | `selftest: 376 configurations, 0 mismatches` (190, 73): key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 18 s; 31 s for 13 and 14 | `logs/rerun/pspread-selftest1{2,3,4}.log` |
| `python src/holes_iso12.py data/p12.txt` | `5 hole classes`, `\|Aut H\|` 2, 96, 32, 32, 16 with `sum \|Aut H\|/\|Aut P\|` (covers of the complement) 1, 84, 4, 4, 4 | 1 s | `logs/rerun/holes-iso12.log` |
| `cover data/p12.txt` | `done: covers 376, tasks 10528, nodes 12303628, solutions 0, full (dim < 8) 0, dim 8 non-APN 0`; `nodes by span dimension 1..8: 10528 294784 4002712 7479748 509108 6747 1 0`; `near misses ...: 0` | 49 s | `logs/20261025-095119-cover12b.log`, `logs/twelve-output.txt` |
| `cover data/i3-reps.txt` | `nodes 2284642, solutions 0` as before; `nodes by span dimension 1..8: 2044 57232 749035 1279756 195498 1077 0 0`, `near misses ...: 0` (a regression check of this build) | 4 s | same job log, `logs/i3-reps-cover3-output.txt` |

Search B, `twelve_check` (`twelve-check.md`; fresh `.prog` files as for fourteen):

| Command | Expected | Time | Log |
|---|---|---|---|
| `twelve_check parta orbits12.txt maxc.txt` | `reduced 9-cliques 5269153 (spr 4128284, ext13 1040086, max 100783)`; `9-cliques, weighted: all 2811561984; extend to a spread 2258825184, extend only to maximal 13s 502840800, maximal 49896000`; `sum of #5-cliques in Q: 2259985728 (C(14,9) * 1128864 = 2259985728)`; `sum of #q giving a maximal 13: 625584960 (10 * 62558496 ...)`, `... extendable 13: 11299928640 (10 * 1129992864 ...)`, `... neither: 0`; `orbits: 376 (spr 280, ext13 84, max 12); inconsistent invariants 0; per-orbit identity ... fails 0`; writes `orbits12.txt` (= `data/orbits12.txt`) and `maxc.txt` (the 100,783 maximal reduced cliques) | 26 min | `logs/20261025-085834-check12-parta.log`, `logs/check12-parta-output.txt` |
| `twelve_check configs data/orbits12.txt data/spreads.txt data/configs.txt configs12.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 42, containing no solid (maximal partial spreads): 5`; `hole sets with a cover mismatch: 0`; spread source `6643 configurations ... spread orbits reached 280 of 280`; 13-source `403 configurations; on no listed orbit 0; on a maximal orbit 0; by type: spr:37 ext13:366; ext13 orbits reached 84 of 84`; writes `configs12.txt` (= `data/configs12.txt`) | 7 s | `logs/20261025-092902-check12-configs.log`, `logs/check12-configs-output.txt` |
| `twelve_check selftest data/configs12.txt` | `lines in D over the configurations: 325:376`; `brute force: 28 forms ...`; `selftest: 376 configurations, 108452736 generated forms checked, 0 failures` | ~1 min | `logs/20261025-093045-check12-cross.log`, `logs/check12-selftest-output.txt` |
| `twelve_check search data/configs12.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 376 10528 294784 4002613 7761928 724602 8217 0 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 0 0 0` (no span ever holds all 12 rank-4 forms) | 20 s | `logs/20261025-092923-check12-asc.log`, `logs/check12-asc-output.txt` |
| `twelve_check search data/configs12.txt desc.prog desc` | `solutions: 0`; nodes by dim `376 10528 294784 3970567 7369230 380832 3848 0 0`; none with every off-D point covered | 20 s | `logs/20261025-093045-check12-cross.log`, `logs/check12-desc-output.txt` |
| `twelve_check search data/configs12.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `376 10528 294784 3979633 7462326 442121 5272 0 0`; none with every off-D point covered | 19 s | same job log, `logs/check12-rand7-output.txt` |
| `twelve_check searchsp data/spreads.txt sp.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 91 pairs of further holes (6,643 configurations, without this program's classification); nodes by dim `6643 186004 5208112 70266988 133848444 9697414 117642 0 0`; none with every off-D point covered | 8 min | same job log, `logs/check12-sp-output.summary.txt` |
| `twelve_check search13 data/configs.txt n13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 13 members removed (403); nodes by dim `403 11284 315952 4288357 8486696 1618431 9352 0 0`; none with every off-D point covered | 32 s | same job log, `logs/check12-n13-output.txt` |
| `twelve_check searchmax maxc.txt maxc.prog asc` | all 100,783 maximal reduced cliques from `parta`, without the canonical form: `solutions: 0`; `search (incl. resumed) nodes by dim: 100783 2821924 79013872 1065082714 1985285822 174884697 929739 0 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 0 0 0` | ~2.6 h of running time between pauses (the program's own clock counts paused time too: 1,203 + 56,565 s, from its per-configuration lines, which the summary copy here leaves out) | `logs/20261025-094106-check12-maxc.log` (first 9,593 configurations; cancelled and resumed), `logs/20261025-105951-check12-maxc2.log` (the rest), `logs/check12-maxc-output.summary.txt` |
| `python src/twelve_check_levels.py data/configs12.txt 0 4` / `... 280 4` / `... 364 4` | `[1, 28, 784, 10804, 21176]` / `[1, 28, 784, 10711, 19628]` / `[1, 28, 784, 10594, 18535]`, the C++ counts for configurations 0 (spr), 280 (ext13) and 364 (max) | 12 s each | `logs/check12-levels-output.txt`, `logs/rerun/t12-levels*.log` |
| `python src/levels_all12.py data/configs12.txt logs/check12-asc-output.txt OUT 8` | the same recount to dimension 8, which for twelve is search B's whole tree (it never branches at points of D here), for all 376 configurations: `levels recount: 376 of 376 configurations, 0 mismatches, 0 missing`, `total nodes (python): 12803048`, the sum of the C++ counts by dimension | ~1 h with 8 workers | `logs/check12-levels-all-output.txt` |
| `twelve_check ctlbatch data/configs12.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0`, `reps hit: 210`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 96` | 6 min | `logs/20261025-094106-check12-ctl.log`, `logs/check12-ctl-bl22-output.txt`, `logs/check12-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths | 1 min | `logs/20261025-094106-check12-ctl.log`, `logs/check12-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check12-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 241`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 96` | ~80 min | `logs/20261025-094106-check12-ctl.log` (the job ran all five lists; 89 min of running time between its pauses), `logs/check12-ctl-new-output.txt`, `logs/check12-ctl-new.txt` |
| `twelve_check ctlreps data/configs12.txt apn_8bit.txt 1 40 56 60 ... 279` (the 66 spr configurations the batches above missed; the full list is in `logs/check12-ctlreps-output.txt`) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 66 | ~20 min | `logs/20261025-105951-check12-ctl2.log`, `logs/check12-ctlreps-output.txt` |
| `twelve_check ctlbatch data/configs12.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1415 controls; re-found 1415`, `nonext: 1419 controls; re-found 1419`, `other spaces found 0`, `S0 already the whole space 46`; tries at dimension 7 only, all at points of D (78,758,848 and 78,286,656). A first trial with `2 10 1` (10 functions per path): 10 of 10 re-found on each | ~2 min | same job log, `logs/check12-ctl-ind-bl22-output.txt`, `logs/check12-ctl-ind-bl22.txt`; the trial `logs/check12-ctl-test-nond.txt` |
| `twelve_check maprowen data/configs12.txt data/p12.txt` | `maprowen: 376 of Rowen's configurations; unmatched 0; distinct listed configurations hit 376 of 376; hit twice 0; \|Aut\| disagreements 0; type/w disagreements 0` (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/check12-maprowen-output.txt` |

### Eleven (`i = 6`)

Kurz's table again gives only counts for eleven solids (303 hole sets, 58 maximal), so both searches classify the
11-partial spreads themselves. Search A, `pspread.cpp` (the same program as for twelve) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread enum 11 p11.txt` | `w=0: reduced 710368, weighted 422009280, orbits 67, orbit sum 422009280 MATCH`; `w=1: reduced 619505, weighted 342619200, orbits 43, ... MATCH`; `w=2: reduced 3324960, weighted 1760081400, orbits 254, ... MATCH`; `w=6: reduced 5728509, weighted 3377565576, orbits 453, ... MATCH`; `orbits 817`; writes `p11.txt` (= `data/p11.txt`) | ~2.1 h of running time between pauses, mostly on the E-cores (the program's clock, 65,072 s, counts paused time too) | `logs/20261025-095302-pspread11.log`, `logs/pspread11-output.txt` |
| `pspread selftest 11 data/p11.txt` | `selftest: 817 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 25 s | `logs/rerun/pspread-selftest11.log` (run for this folder) |
| `cover data/p11.txt` | `done: covers 817, tasks 22876, nodes 30505655, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`; `nodes by span dimension 1..8: 22876 640528 8993626 19289387 1533912 25288 38 0`; `near misses ...: 58`: the `NEAR` lines, spans holding all 11 rank-4 forms with holes uncovered (22 at dim 6, 36 at dim 7, in 11 configurations), none of which reaches dim 8 | 89 s | `logs/20261026-092941-cover11.log`, `logs/eleven-output.txt` |

Search B, `eleven_check` (`eleven-check.md`; fresh `.prog` files as for fourteen). `parta` runs one minimal class per
job; the five outputs are `data/parta11-c*.txt`, its lists of maximal cliques aren't included (`Data` below):

| Command | Expected | Time | Log |
|---|---|---|---|
| `eleven_check parta C parta-cC.txt maxc-cC.txt` for C = 0 … 4 | `class C ...: N reduced 8-cliques, M maximal`, with N = `5917664`, `3003006`, `736536`, `364396`, `361740` (10,383,342 in all) and M = `599168`, `86195`, `17880`, `7125`, `0` (710,368); each file ends in `END`; `parta-cC.txt` = `data/parta11-cC.txt` | 8,611 / 4,950 / 815 / 832 / 195 s by the program's clock (about 2 h, 83 min, 14 min, 14 min and 3 min of running time, on the E-cores) | `logs/2026…-check11-pa*.log`, `logs/check11-parta-c*-output.txt` |
| `eleven_check merge orbits.txt data/configs12.txt data/parta11-c0.txt ... data/parta11-c4.txt` | `sum over 8-cliques of \|Q\| = 25304057856; 9 * (weighted 9-cliques, twelve-check) = 25304057856  OK`; `orbits 817; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 0:67 1:43 2:254 6:453`; `orbits by type: spr 453, ext13 254, ext12 43, max 67`; `... C(14,8) * 1128864 = 3389978592  OK`; `per 12-orbit O: ... fails for 0 of 376`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits11.txt`) | 6 s | `logs/check11-merge-output.txt` |
| `eleven_check configs orbits.txt data/configs12.txt data/spreads.txt data/configs.txt configs11.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 303, containing no solid (maximal partial spreads): 58`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 2134`; outside sources `twelve-check's 376 x 12: 4512 configurations; on no listed orbit 0; ... spr 453/453 ext13 254/254 ext12 43/43`, `fourteen-check's spread reps x 364: 26572 ...; on no listed orbit 0; ... spr 453/453`, `thirteen-check's maximal 13s x 78: 2418 ...; on no listed orbit 0; ... ext13 254/254`; writes `configs11.txt` (= `data/configs11.txt`; spr 0..452, ext13 453..706, ext12 707..749, max 750..816) | 19 s | `logs/check11-configs-output.txt` |
| `eleven_check selftest data/configs11.txt` | `lines in D over the configurations: 510:817`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28)`; `selftest: 817 configurations, 368561396 generated forms checked, 0 failures` | 4 min | `logs/check11-selftest-output.txt` |
| `eleven_check search data/configs11.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 817 22876 640528 8993779 19647718 1938778 24532 38 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 22 37 0`; `of which at points of D: 0 0 0 0 0 0 3930304 3013696 0` (the rank-6 forms tried in the holes); `distinct p1 spans this run: 59` | 18 min on the E-cores | `logs/20261026-210508-check11-asc.log`, `logs/check11-asc-output.txt`, near-miss spans `logs/check11-asc-p1.txt` |
| `eleven_check search data/configs11.txt desc.prog desc` | `solutions: 0`; nodes by dim `817 22876 640528 8935159 19095956 1335017 18128 40 0`; every off-D point covered `22` / `36` at dims 6 / 7; `distinct p1 spans this run: 58` | 18 min | `logs/20261026-210508-check11-desc.log`, `logs/check11-desc-output.txt`, `logs/check11-desc-p1.txt` |
| `eleven_check search data/configs11.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `817 22876 640528 8944444 19234261 1513135 18258 38 0`; `22` / `37`; `distinct p1 spans this run: 59` | 18 min | `logs/20261026-210509-check11-rand7.log`, `logs/check11-rand7-output.txt`, `logs/check11-rand7-p1.txt` |
| `eleven_check search12 data/configs12.txt s12.prog asc` | `solutions: 0` over twelve_check's 376 configurations × 12 members removed (4,512: every 11 that extends, without this program's classification); nodes by dim `4512 126336 3537408 49655433 109538964 10777400 165744 111 0`; every off-D point covered `97` / `103` at dims 6 / 7 | 76 min | `logs/20261026-175950-check11-s12.log`, `logs/check11-s12-output.summary.txt` |
| `eleven_check search13 data/configs.txt s13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 78 pairs removed (2,418); nodes by dim `2418 67704 1895712 26629650 60146644 8913624 96436 72 0`; `22` / `61` | 50 min | `logs/20261026-175950-check11-s13.log`, `logs/check11-s13-output.summary.txt` |
| `eleven_check searchsp data/spreads.txt ssp.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 364 triples removed (26,572: every 11 that extends to a spread); nodes by dim `26572 744016 20832448 291080164 635317280 54887313 860221 495 0`; `263` / `450` | ~4.5 h of running time, mostly on the E-cores (evicted once and resumed from `ssp.prog`) | `logs/20261026-175950-check11-ssp.log`, `logs/check11-ssp-output.summary.txt` |
| `eleven_check searchmax maxc-s50.txt maxs50.prog asc`, where `maxc-s50.txt` is lines 1, 51, 101, … of `maxc-c0.txt` … `maxc-c3.txt` concatenated (14,208 maximal cliques) | `solutions: 0`; nodes by dim `14208 397824 11139072 154792269 321538731 14093008 17710 0 0`; no span holds all 11 rank-4 forms | ~2.1 h | `logs/20261026-210924-check11-maxs50.log`, `logs/check11-maxs50-output.summary.txt` |
| `eleven_check orbitcov orbits.txt configs11.txt maxc-s50.txt` | `orbitcov: 14208 cliques, 0 on no listed orbit, 0 on a non-max orbit; 66 orbits hit of 67 max orbits` (orbit 798 missed) | seconds | `logs/check11-maxs50-orbits-output.txt` |
| `eleven_check orbitcov orbits.txt configs11.txt maxc-all.txt maxc-o798.txt 798` (all four `maxc` files concatenated), then `eleven_check searchmax maxc-o798.txt max798.prog asc` | `orbitcov: 710368 cliques, 0 on no listed orbit, 0 on a non-max orbit; 67 orbits hit of 67 max orbits`, the per-orbit tallies equal to `n` in `orbits.txt`; `wrote 85 cliques` (= `data/maxc11-o798.txt`); then `solutions: 0`, nodes by dim `85 2380 66640 934228 2128630 407037 5956 0 0` | 3.5 min | `logs/20261026-235204-check11-max798.log`, `logs/check11-maxall-orbits-output.txt`, `logs/check11-max798-output.txt` |
| `python src/eleven_check_levels.py data/configs11.txt 0 4` / `453 4` / `707 4` / `750 4` | `[1, 28, 784, 11152, 24772]` / `[1, 28, 784, 11206, 23367]` / `[1, 28, 784, 10939, 24007]` / `[1, 28, 784, 11074, 22957]`, the C++ counts for one configuration of each type (spr, ext13, ext12, max) | 1 min for the four | `logs/20261026-210759-check11-levels.log`, `logs/check11-levels-output.txt` |
| `eleven_check ctlbatch data/configs11.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 319`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 364` | ~1 h for the four lists | `logs/20261026-210844-check11-ctl.log`, `logs/check11-ctl-bl22-output.txt`, `logs/check11-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths (the job exits 3 because of the last) | (same job) | same job log, `logs/check11-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check11-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 383`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 364` | ~2.3 h | `logs/20261026-210845-check11-ctlnew.log`, `logs/check11-ctl-new-output.txt`, `logs/check11-ctl-new.txt` |
| `eleven_check ctlbatch data/configs11.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1438 controls; re-found 1438`, `nonext: 1440 controls; re-found 1440`, `other spaces found 0`, `S0 already the whole space 2`; every control starts at dimension 7 and finds its form among the 13,888 on its line (98,965,888 and 99,785,280 tries, all at points of D) | 5 min | `logs/20261026-210845-check11-ctlnond.log`, `logs/check11-ctl-ind-bl22-output.txt`, `logs/check11-ctl-ind-bl22.txt` |
| `eleven_check ctlreps data/configs11.txt apn_8bit.txt 1 40 93 112 ... 452` (the 64 spr configurations the batches above missed; the full list is in the job log's command) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 64 | 8 min | `logs/20261026-235204-check11-ctlreps.log`, `logs/check11-ctlreps-output.txt` |
| `eleven_check mapa data/configs11.txt data/p11.txt` | `mapa: 817 of A's configurations, 0 on none of ours, 817 of our 817 configurations hit, 0 hit more than once`, with type and `canon aut` agreeing with search A's tags on every line (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/check11-mapa-output.txt` |

The near misses agree between the two searches configuration by configuration: mapped through `mapa`, search A's
`NEAR` lines fall in search B's configurations 248, 299, 314, 316, 429, 435, 452, 541, 558, 609 and 649, 6, 6, 6, 27,
1, 1, 3, 3, 1, 1 and 3 of them, the counts search B's `desc` run has there (its per-configuration `p1` counts); `asc`
and `rand 7` have a second span at 558. Which partial spans a point-branching tree visits depends on its point order.

### Ten (`i = 7`)

Kurz's table gives counts for ten solids too (1,336 hole sets, 317 maximal), so both searches classify the 10-partial
spreads themselves once more. Many of these jobs ran side by side, so their times overlap; times are running time
between the runner's pauses. Search A, `pspread3.cpp` (`pspread.cpp`'s enumeration with checkpoints) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread3 enum 12 p12.txt --ckpt p12.ckpt --wcheck 200`, `pspread3 enum 14 p14.txt --ckpt p14.ckpt --wcheck 50` | the checks of the new build: `w spot checks: 12411, mismatches 0`, `orbits 376`; `w spot checks: 17221, mismatches 0`, `orbits 73`; the outputs are byte for byte `data/p12.txt` and `data/p14.txt`. Resumed from a complete checkpoint the K = 14 run writes the same file at once, and from one holding part of the run (428 tasks and a torn last line, which it drops) the same file again | ~4.6 h for K = 12, beside the K = 10 run; 18 min for K = 14 | `logs/20261026-180053-p3-k12.log`, `logs/20261026-170942-p3-k14.log`, `logs/20261026-172744-p3-k14-resume.log`, `logs/20261026-172831-p3-k14-half.log`, `logs/pspread3-p12-output.txt`, `logs/pspread3-p14*-output.txt` |
| `pspread3 enum 10 p10.txt --ckpt p10.ckpt --wcheck 1000` | `w=0: reduced 6042938, weighted 3974342400, orbits 320, orbit sum 3974342400 MATCH`; `w=1: reduced 4255438, weighted 2774499840, orbits 274, ... MATCH`; `w=2: reduced 1707940, weighted 1056585600, orbits 111, ... MATCH`; `w=3: reduced 5835506, weighted 3414902400, orbits 460, ... MATCH`; `w=7: reduced 5986006, weighted 3815466624, orbits 612, ... MATCH`; `orbits 1777`; `w spot checks: 17563, mismatches 0`; writes `p10.txt` (= `data/p10.txt`, tags `p10-<n>-a<\|Aut P\|>-w<w>`) | ~5 h: 2.2 h on the E-cores before the runner restarted, then 2.8 h resumed from the checkpoint (114 tasks), mostly on an open box (the program's clock for that run, 10,610 s, counts paused time too) | `logs/20261026-180053-p3-k10.log`, `logs/pspread10-output.txt` |
| `pspread3 selftest 10 data/p10.txt` | `selftest: 1777 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 21 s | `logs/rerun/pspread3-selftest10.log` (run for this folder) |
| `cover data/p10.txt` | `done: covers 1777, tasks 49756, nodes 79061670, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`; `nodes by span dimension 1..8: 49756 1393168 20447257 52320443 4776667 74207 172 0`; `near misses ...: 370`: the `NEAR` lines, 337 at dim 6 (in 27 configurations) and 33 at dim 7 (in 15), 39 configurations in all, none reaching dim 8 | 202 s | `logs/20261026-235746-i7-cover10.log`, `logs/ten-output.txt` |

Search B, `ten_check` (`ten-check.md`; fresh `.prog` files as for fourteen). `parta` runs one slice of a minimal class
per job; the 17 outputs are `data/parta10-c*-s*.txt`, its lists of maximal cliques aren't included (`Data` below):

| Command | Expected | Time | Log |
|---|---|---|---|
| `ten_check parta C parta-cC-sS.txt maxc-cC-sS.txt S N` for class C = 0 … 4 in N = 8, 4, 2, 2, 1 slices | per slice `class C slice S/N ...: R reduced 7-cliques, M maximal`; by class R = 16,069,200, 5,277,651, 1,423,134, 621,223, 436,620 (23,827,828 in all) and M = 5,049,808, 785,215, 160,785, 46,290, 840 (6,042,938); each file ends in `END`; `parta-cC-sS.txt` = `data/parta10-cC-sS.txt` | 44 min from the first start to the last end on an open box, eleven jobs at a time; by the program's clock 1,407 to 2,548 s per slice of class 0, 596 to 1,747 s of class 1, at most 311 s for the rest | `logs/2026…-check10-pa*-s*.log`, `logs/check10-parta-c*-s*-output.txt` |
| `ten_check merge orbits.txt data/configs11.txt data/parta10-c0-s0.txt ... data/parta10-c4-s0.txt` (the 17 in class and slice order) | `eleven-check's configurations: 817; sum 990*20160/\|Aut\| = 5902275456 (its Part A: 5902275456)`; `reduced 7-cliques 23827828, of them maximal (Q empty) 6042938; weighted 7-cliques 15035796864`; `sum over 7-cliques of \|Q\| = 47218203648; 8 * (weighted 8-cliques, eleven-check) = 47218203648  OK`; `orbits 1777; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 0:320 1:274 2:111 3:460 7:612`; `orbits by type: spr 612, ext13 460, ext12 111, ext11 274, max 320`; `type / Q / n7 inconsistencies 0; P + q on no eleven-check orbit 0`; `... C(14,7) * 1128864 = 3874261248  OK`; `per 11-orbit O: ... fails for 0 of 817`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits10.txt`) | under a minute with `configs` | `logs/20261027-045312-check10-merge-configs.log`, `logs/check10-merge-output.txt` |
| `ten_check configs orbits.txt data/configs11.txt data/configs12.txt data/spreads.txt data/configs.txt configs10.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 1336, containing no solid (maximal partial spreads): 317`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 2294`; outside sources `eleven-check's 817 x 11: 8987 configurations; on no listed orbit 0; ... spr 612/612 ext13 460/460 ext12 111/111 ext11 274/274 max 0/320`, `fourteen-check's spread reps x 1001: 73073 ...; on no listed orbit 0; ... spr 612/612`, `thirteen-check's maximal 13s x 286: 8866 ...; on no listed orbit 0; ... spr 252/612 ext13 460/460`, `twelve-check's 376 x 66: 24816 ...; on no listed orbit 0; ... spr 612/612 ext13 460/460 ext12 111/111`; writes `configs10.txt` (= `data/configs10.txt`; spr 0..611, ext13 612..1071, ext12 1072..1182, ext11 1183..1456, max 1457..1776) | (same job) | same job log, `logs/check10-configs-output.txt` |
| `ten_check selftest data/configs10.txt` | `selftest: ND4 28 (28), ND6 13888 (13888), GL4 20160 (20160), V5824 5824`; `lines in D over the configurations: 770:1777`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28)`; `selftest: 1777 configurations, 1198198680 generated forms checked, 0 failures` | 41 min, beside the searches below | `logs/20261027-121908-check10-selftest.log`, `logs/check10-selftest-output.txt` |
| `ten_check search data/configs10.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 1777 49756 1393168 20447263 52683154 5454191 70909 181 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 337 165 0`; `of which at points of D: 0 0 0 0 0 0 74550784 19082112 0` (the rank-6 forms tried in the holes); `distinct p1 spans this run: 502` | 24 min, the three orders at once | `logs/20261027-121909-check10-asc.log`, `logs/check10-asc-output.txt`, near-miss spans `logs/check10-asc-p1.txt` |
| `ten_check search data/configs10.txt desc.prog desc` | `solutions: 0`; nodes by dim `1777 49756 1393168 20333833 51884404 4463901 60838 131 0`; every off-D point covered `337` / `120` at dims 6 / 7; at points of D `73661952` / `13360256`; `distinct p1 spans this run: 457` | 25 min | `logs/20261027-121909-check10-desc.log`, `logs/check10-desc-output.txt`, `logs/check10-desc-p1.txt` |
| `ten_check search data/configs10.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `1777 49756 1393168 20336938 51929881 4484972 61922 169 0`; `337` / `150`; at points of D `73717504` / `16221184`; `distinct p1 spans this run: 487` | 24 min | `logs/20261027-121909-check10-rand7.log`, `logs/check10-rand7-output.txt`, `logs/check10-rand7-p1.txt` |
| `ten_check search11 data/configs11.txt s11.prog asc` | `solutions: 0` over eleven_check's 817 configurations × 11 members removed (8,987: every 10 that extends, without this program's classification); nodes by dim `8987 251636 7045808 103434188 270500229 31444429 568714 864 0`; every off-D point covered `2507` / `768` at dims 6 / 7 | 53 min | `logs/20261027-004853-check10-s11.log`, `logs/check10-s11-output.summary.txt` |
| `ten_check search13 data/configs.txt s13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 286 triples removed (8,866); nodes by dim `8866 248248 6950944 102193519 273604875 39199199 677779 656 0`; `661` / `543` | 47 min | `logs/20261027-004853-check10-s13.log`, `logs/check10-s13-output.summary.txt` |
| `ten_check searchsp data/spreads.txt ssp.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 1001 quadruples removed (73,073: every 10 inside a spread); nodes by dim `73073 2046044 57289232 838595573 2191569950 246581805 5099404 5591 0`; `28477` / `4409` | ~6.5 h of running time (resumed from `ssp.prog` across many holds) | `logs/20261027-004853-check10-ssp.log`, `logs/check10-ssp-output.summary.txt` |
| `ten_check orbitcov orbits.txt configs10.txt data/maxc10-s400.txt`, then `ten_check searchmax data/maxc10-s400.txt maxs400.prog asc`, where `maxc10-s400.txt` is lines 1, 401, 801, … of `maxc-c0-s0.txt` … `maxc-c4-s0.txt` concatenated in class and slice order (15,108 maximal cliques) | `orbitcov: 15108 cliques, 0 on no listed orbit, 0 on a non-max orbit; 320 orbits hit of 320 max orbits`; then `solutions: 0`, nodes by dim `15108 423024 11844672 172205757 428754848 25136558 4469 0 0`, no span holds all 10 rank-4 forms (`searchmax` labels these cliques `ext11` in its per-type lines, a label only, `ten-check.md` §5) | 99 min for the two | `logs/20261027-121919-check10-maxs400.log`, `logs/check10-maxs400-orbits-output.txt`, `logs/check10-maxs400-output.summary.txt` |
| `ten_check orbitcov orbits.txt configs10.txt maxc-cC-sS.txt` for each of the 17 files, appended to `maxall-orbits.out`; then `python src/ten_check_tally.py maxall DIR` on a folder holding `orbits.txt` and `maxall-orbits.out` | each `orbitcov: N cliques, 0 on no listed orbit, 0 on a non-max orbit`; then `orbitcov runs 17 (expect 17); cliques off the max orbits 0; max orbits 320, cliques 6042938 (Part A 6042938); orbits whose tally differs from Part A's n: 0` (so every maximal clique lies on a listed maximal configuration, and all 320 have cliques) | 46 min | `logs/20261027-121924-check10-maxcov.log`, `logs/check10-maxall-orbits-output.txt` |
| `python src/ten_check_levels.py data/configs10.txt 0 4` / `612 4` / `1072 4` / `1183 4` / `1457 4` | `[1, 28, 784, 11578, 31023]` / `[1, 28, 784, 11680, 31277]` / `[1, 28, 784, 11443, 28745]` / `[1, 28, 784, 11542, 28800]` / `[1, 28, 784, 11467, 29638]`, the C++ counts for one configuration of each type (spr, ext13, ext12, ext11, max) | 6 min for the five | `logs/20261027-121909-check10-levels.log`, `logs/check10-levels-output.txt` |
| `ten_check ctlbatch data/configs10.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 410`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 976` | 35 min for the four lists | `logs/20261027-121918-check10-ctl.log`, `logs/check10-ctl-bl22-output.txt`, `logs/check10-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222 (reps hit 160 / 193), `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths (the job exits 3 because of the last) | (same job) | same job log, `logs/check10-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check10-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 513`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 1163` | 1.7 h | `logs/20261027-121918-check10-ctlnew.log`, `logs/check10-ctl-new-output.txt`, `logs/check10-ctl-new.txt` |
| `ten_check ctlbatch data/configs10.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1440 controls; re-found 1440`, `nonext: 1440 controls; re-found 1440`, `other spaces found 0`; `started from S0 of dimension: 7:1440` on each path, every control finding its form among the 13,888 on its line (124,811,456 and 123,047,680 tries, all at dimension 7) | 3 min | `logs/20261027-121918-check10-ctlnond.log`, `logs/check10-ctl-ind-bl22-output.txt`, `logs/check10-ctl-ind-bl22.txt` |
| `python src/ten_check_tally.py reps DIR` on a folder holding `configs.txt` (= `data/configs10.txt`) and the `ctl-*.txt` files above | `reps hit 1688 of 1777; not reached 89 {'spr': 89}` and the 89 | seconds | `logs/rerun/t10-tally.log` (run for this folder) |
| `ten_check ctlreps data/configs10.txt apn_8bit.txt 1 40 131 132 ... 611` (those 89; the full list is in the job log's command) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 89 | 33 min | `logs/20261027-182638-check10-ctlreps.log`, `logs/check10-ctlreps-output.txt` |
| `ten_check mapa data/configs10.txt data/p10.txt` | `mapa: 1777 of A's configurations, 0 on none of ours, 1777 of our 1777 configurations hit, 0 hit more than once`, with type (w = 0, 1, 2, 3, 7 for max, ext11, ext12, ext13, spr) and `canon aut` agreeing with search A's tags on every line (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/20261028-174550-check10-mapa.log`, `logs/check10-mapa-output.txt` |

The near misses, compared by `python src/near_misses10.py` from the logs above: search A's 337 at dimension 6 lie in 27
configurations (mapped through `mapa`), and each of search B's three orders has its dimension-6 near misses in exactly
those 27, the same number in each. Search A's 33 at dimension 7 lie in 15 configurations, where every order of search B
has the same number. Search B's other dimension-7 near misses (132 / 87 / 117 in `asc` / `desc` / `rand 7`) all contain
one of its dimension-6 near misses: a rank-6 form fitted at a hole, which cover.cpp doesn't log, since it prints `NEAR`
only at the step where the last rank-4 form joins (none of search A's 33 contains one of its dimension-6 near misses).
Those sit in configurations 91, 331, 352 and 377 in every order, and 376 in `asc` too; each already has
dimension-6 near misses on both sides. The script prints these counts; the spans themselves are in different
coordinates on the two sides, so it compares counts, not spans.

### Nine (`i = 8`)

Kurz's table gives counts for nine solids as well (5,442 hole sets, 1,343 maximal), so both searches classify the
9-partial spreads themselves again. As for ten, times are running time between the runner's pauses, and jobs often
ran side by side. Search A, `pspread3.cpp` (the same program as for ten) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread3 enum 9 p9.txt --ckpt p9.ckpt --wcheck 1000` | `w=0: reduced 16492152, weighted 12409921496, orbits 1348, orbit sum 12409921496 MATCH`; `w=1: reduced 25008768, weighted 18422651520, orbits 1916, ... MATCH`; `w=2: reduced 10373498, weighted 7504126560, orbits 874, ... MATCH`; `w=3: reduced 2779569, weighted 1911047040, orbits 227, ... MATCH`; `w=4: reduced 5994783, weighted 3937565520, orbits 623, ... MATCH`; `w=8: reduced 4674633, weighted 3228877176, orbits 675, ... MATCH`; `orbits 5663`; `w spot checks: 52194, mismatches 0`; writes `p9.txt` (= `data/p9.txt`, tags `p9-<n>-a<\|Aut P\|>-w<w>`) | ~3.8 h: 28 min before the first job was cancelled during a hold, then 3.4 h resumed from the checkpoint (87 tasks), the last two hours on the E-cores (the program's clock for that run, 36,568 s, counts paused time too) | `logs/20261027-000150-p3-k9.log`, `logs/20261027-005025-p3-k9.log`, `logs/pspread9-output.txt` (both starts; the program appends to it) |
| `pspread3 selftest 9 data/p9.txt` | `selftest: 5663 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 34 s | `logs/rerun/pspread3-selftest9.log` (run for this folder) |
| `cover data/p9.txt` | `done: covers 5663, tasks 158564, nodes 321293681, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`; `nodes by span dimension 1..8: 158564 4439792 68941328 220391315 27185773 176554 355 0`; `near misses ...: 600`: the `NEAR` lines, 18 at dim 5 (in 7 configurations), 395 at dim 6 (in 151) and 187 at dim 7 (in 85), 201 configurations in all, none reaching dim 8 | 56 min of running time (4.1 h with the pauses) | `logs/20261027-121559-cover9.log`, `logs/nine-output.txt` |

Search B, `nine_check` (`nine-check.md`; fresh `.prog` files as for fourteen), its jobs run by `nine_check_slices.py`
(the job logs give its commands). `parta` runs one slice of a minimal class per job; the 34 outputs are
`data/parta9-c*-s*.txt`, its lists of maximal cliques aren't included (`Data` below):

| Command | Expected | Time | Log |
|---|---|---|---|
| `nine_check parta C parta-cC-sS.txt maxc-cC-sS.txt S N` for class C = 0 … 4 in N = 16, 8, 4, 4, 2 slices | per slice `class C slice S/N ...: R reduced 6-cliques, M maximal`; by class R = 49,187,232, 11,789,427, 3,073,814, 971,562, 301,368 (65,323,403 in all) and M = 13,254,204, 2,614,725, 532,200, 91,020, 3 (16,492,152); each file ends in `END`; `parta-cC-sS.txt` = `data/parta9-cC-sS.txt` | 2.2 h on the E-cores, one slice after another; by the program's clock 262 to 475 s per slice of class 0, 118 to 183 s of class 1, at most 118 s for the rest | `logs/20261028-191409-check9-parta.log`, `logs/check9-parta-c*-s*-output.txt` |
| `nine_check merge orbits.txt data/configs10.txt data/parta9-c0-s0.txt ... data/parta9-c4-s1.txt` (the 34 in class and slice order) | `ten-check's configurations: 1777; sum 720*20160/\|Aut\| = 15035796864 (its Part A: 15035796864)`; `reduced 6-cliques 65323403, of them maximal (Q empty) 16492152; weighted 6-cliques 47414189312`; `sum over 6-cliques of \|Q\| = 105250578048; 7 * (weighted 7-cliques, ten-check) = 105250578048  OK`; `orbits 5663; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 0:1348 1:1916 2:874 3:227 4:623 8:675`; `orbits by type: spr 675 ext13 623 ext12 227 ext11 874 ext10 1916 max 1348`; `type / Q / n8 inconsistencies 0; P + q on no ten-check orbit 0`; `... C(14,6) * 1128864 = 3389978592  OK`; `per 10-orbit O: ... fails for 0 of 1777`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits9.txt`) | 4 s | `logs/20261028-230254-check9-merge.log`, `logs/check9-merge-output.txt` |
| `nine_check configs orbits.txt data/configs10.txt data/configs11.txt data/configs12.txt data/spreads.txt data/configs.txt configs9.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 5442, containing no solid (maximal partial spreads): 1343   [Kurz Table 9.3, r = 9: 5442 and 1343]`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 7925`; `orbits reached by no cover: 0`; outside sources `ten-check's 1777 x 10: 17770 configurations; on no listed orbit 0; ... spr 675/675 ext13 623/623 ext12 227/227 ext11 874/874 ext10 1916/1916 max 0/1348`, `fourteen-check's spread reps x 2002: 146146 ...; on no listed orbit 0; ... spr 675/675`, `thirteen-check's maximal 13s x 715: 22165 ...; on no listed orbit 0; ... spr 425/675 ext13 623/623`, `twelve-check's 376 x 220: 82720 ...; on no listed orbit 0; ... spr 675/675 ext13 623/623 ext12 227/227`, `eleven-check's 817 x 55: 44935 ...; on no listed orbit 0; ... spr 675/675 ext13 623/623 ext12 227/227 ext11 874/874`; writes `configs9.txt` (= `data/configs9.txt`; spr 0..674, ext13 675..1297, ext12 1298..1524, ext11 1525..2398, ext10 2399..4314, max 4315..5662). Hole set 0, the 120 points off a hyperbolic quadric (`nine-check.md` §4): `\|Aut H\| 348364800, solid inside: no, 270 solids in the complement, 1920 covers, 1 orbits ( max:1 )` | 16 min, 11 of them for the hole sets' automorphism groups | `logs/20261028-230339-check9-configs.log`, `logs/check9-configs-output.txt` |
| `nine_check selftest data/configs9.txt` | `selftest: ND4 28 (28), ND6 13888 (13888), GL4 20160 (20160), V5824 5824`; `lines in D over the configurations: 1120:5663`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28)`; `selftest: 5663 configurations, 5507741316 generated forms checked, 0 failures` | ~2.2 h, beside the searches below | `logs/20261028-235617-check9-selftest.log`, `logs/check9-selftest-output.txt` |
| `nine_check search data/configs9.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 5663 158564 4439792 68939870 220381954 28376993 178277 351 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 18 3922 296 0`; `of which at points of D: 0 0 0 0 0 6721792 989853312 37483712 0` (the rank-6 forms tried in the holes); `distinct p1 spans this run: 4236` (the job's second start found the run complete and printed the same totals, with 0 spans) | 21 min on an open box | `logs/20261028-235552-check9-search.log`, `logs/check9-asc-output.summary.txt`, per-configuration `R` lines and near-miss spans (`P1`) `logs/check9-asc-prog.txt` |
| `nine_check search data/configs9.txt desc.prog desc` | `solutions: 0`; nodes by dim `5663 158564 4439792 68555168 218112825 26492204 171097 361 0`; every off-D point covered `18` / `3820` / `305` at dims 5 / 6 / 7; at points of D `6735680` / `963507776` / `38150336`; `distinct p1 spans this run: 9`: requeued after 1,595 configurations when the runner stopped and resumed from `desc.prog`, so the first start's spans weren't written (its `R` lines were) | 36 min | same job log, `logs/check9-desc-output.summary.txt`, `logs/check9-desc-prog.txt` |
| `nine_check search data/configs9.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `5663 158564 4439792 68557415 218135142 26392740 162550 386 0`; `18` / `3768` / `321`; at points of D `6749568` / `927329536` / `39483584`; `distinct p1 spans this run: 4107` | 38 min, on the P-cores at most 40% | same job log, `logs/check9-rand7-output.summary.txt`, `logs/check9-rand7-prog.txt` |
| `python src/nine_check_tally.py prog data/configs9.txt 5663 logs/check9-asc-prog.txt` (`desc`, `rand7`) | `configurations done 5663 (expect 5663); duplicate R lines 0`; the totals above, and by type, e.g. for `asc` `spr (675)` near misses `18 / 3,790 / 231` at dims 5 / 6 / 7 and `max (1348)` `0 / 1 / 0`; `distinct P1 spans 4236` / `9` / `4107`; `SOLUTION lines 0` (`nine-check.md` §5's table) | seconds | `logs/rerun/t9-tally.log` (run for this folder) |
| `nine_check search10 data/configs10.txt s10.prog asc` | `solutions: 0` over ten_check's 1,777 configurations × 10 members removed (17,770: every 9 that extends, without this program's classification); nodes by dim `17770 497560 13931680 216427891 702927941 106444614 1700247 2793 0`; every off-D point covered `46` / `12042` / `2175` at dims 5 / 6 / 7 | ~1.8 h (killed when the runner restarted and resumed from `s10.prog`) | `logs/20261029-044940-check9-cross2.log`, `logs/check9-s10-output.summary.txt` |
| `nine_check search13 data/configs.txt s13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 715 quadruples removed (22,165); nodes by dim `22165 620620 17377360 270680302 902008389 165304542 3887639 7102 0`; `22` / `12500` / `5379` | 66 min | same job log, `logs/check9-s13-output.summary.txt` |
| `nine_check searchsp data/spreads.txt ssp.prog asc FIRST LAST` in 8 ranges of 18,269 (the last 18,263) | `solutions: 0` in each, over fourteen_check's 73 representatives × 2002 quintuples removed (146,146: every 9 inside a spread); summed over the ranges, nodes by dim `146146 4092088 114578464 1778825956 5844031218 985797672 25756497 42522 0`, every off-D point covered `985` / `231767` / `32409` | ~17 h of running time, mostly on the E-cores (resumed across many holds and two requeues) | `logs/20261028-230727-check9-ssp.log`, `logs/check9-ssp-r0-output.summary.txt` … `-r7-` |
| `nine_check searchmax data/maxc9-s1000.txt maxs1000.prog asc`, where `maxc9-s1000.txt` is lines 1, 1001, 2001, … of `maxc-c0-s0.txt` … `maxc-c4-s1.txt` concatenated in class and slice order (16,493 maximal cliques) | `solutions: 0`, nodes by dim `16493 461804 12930512 199339583 628642822 69498624 15776 0 0`; 6 spans at dim 6 hold all 9 rank-4 forms (`logs/check9-maxs1000-near.txt` lists their cliques, all on orbit 4776, config 4776: the one max configuration with a near miss in the main search) | 50 min | `logs/20261028-230443-check9-maxsample.log`, `logs/20261028-230536-check9-maxs1000.log`, `logs/check9-maxs1000-output.summary.txt` |
| `nine_check orbitcov orbits.txt configs9.txt maxc-cC-sS.txt` for each of the 34 files, appended to `maxall-orbits.out`; then `python src/ten_check_tally.py maxall DIR` on a folder holding `orbits.txt` and `maxall-orbits.out` | each `orbitcov: N cliques, 0 on no listed orbit, 0 on a non-max orbit`; then `orbitcov runs 34 (expect 17); cliques off the max orbits 0; max orbits 1348, cliques 16492152 (Part A 16492152); orbits whose tally differs from Part A's n: 0` (the script's `expect 17` is ten's number) | 32 min | same job log as `search10`, `logs/check9-maxall-orbits-output.txt`; the tally `logs/rerun/t9-tally.log` |
| After the verdict (Rowen, at the review's suggestion): `python apn/nine_check_slices.py maxmissed` (orbitcov on the sample lists the 13 max orbits it misses, then writes every raw clique of Part A on them), then `nine_check searchmax data/maxc9-missed.txt maxmissed.prog asc` | `sample misses 13 max orbits`; `maxc-missed.txt: 14646 cliques on the 13 orbits`; `solutions: 0`, nodes by dim `14646 410088 11482464 177518586 562247627 70164318 133689 0 0`, no node with every off-D point covered | 43 min, orbitcov included | `logs/check9-maxmissed-output.summary.txt` |
| `python src/nine_check_levels.py data/configs9.txt 0 4` / `675 4` / `1298 4` / `1525 4` / `2399 4` / `4315 4` / `5662 4` | `[1, 28, 784, 12220, 40545]` / `[1, 28, 784, 12295, 40929]` / `[1, 28, 784, 12184, 40622]` / `[1, 28, 784, 12193, 38864]` / `[1, 28, 784, 12202, 39564]` / `[1, 28, 784, 12412, 39324]` / `[1, 28, 784, 12994, 56092]`, the C++ counts for one configuration of each type (spr, ext13, ext12, ext11, ext10, max) and for the quadric's | 4 min for the seven | `logs/20261028-232005-check9-levels.log`, `logs/check9-levels-output.txt` |
| `nine_check ctlbatch data/configs9.txt apn_8bit.txt 2 1000000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1440 controls; re-found 1440`, `nonext: 1436 controls; re-found 1436`, `other spaces found 0`; `started from S0 of dimension: 7:1440` and `7:1436`, every control finding its form among the 13,888 on its line (153,642,944 and 154,017,920 tries, all at dimension 7); 4 `k = 4` functions skipped on the nonext path (`no member 4-tuple of the class`) | 1 min | `logs/20261029-030815-check9-controls.log` (the job ran all six batches), `logs/check9-ctl-ind-bl22-output.txt`, `logs/check9-ctl-ind-bl22.txt` |
| `nine_check ctlbatch data/configs9.txt apn_8bit.txt 2 1000000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 443`; `nonext: 1436 controls; re-found 1436 ...`, `reps hit: 1412`, with the same 4 skipped | 16 min | same job log, `logs/check9-ctl-bl22-output.txt`, `logs/check9-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222 (reps hit 162 / 219), `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths (the job exits 3 because of the last) | 4 min for the three | same job log, `logs/check9-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check9-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 575`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 3837` | 80 min | same job log, `logs/check9-ctl-new-output.txt`, `logs/check9-ctl-new.txt` |
| `python src/ten_check_tally.py reps DIR` on a folder holding `configs.txt` (= `data/configs9.txt`) and the `ctl-*.txt` files above | `reps hit 4812 of 5663; not reached 851 {'spr': 91, 'ext13': 90, 'ext12': 32, 'ext11': 133, 'ext10': 293, 'max': 212}` and the 851 | seconds | `logs/check9-reps-unreached.txt`; the re-run `logs/rerun/t9-tally.log` |
| `nine_check ctlreps data/configs9.txt apn_8bit.txt 1 40 123 124 ... 5655` (those 851) | `rep N (ext)` or `(nonext): 40 of 40 re-found, 0 other spaces` for each of the 851 | ~4.5 h | `logs/20261029-085400-check9-ctlreps.log`, `logs/check9-ctlreps-output.txt` |
| `nine_check mapa data/configs9.txt data/p9.txt` | `mapa: 5663 of A's configurations, 0 on none of ours, 5663 of our 5663 configurations hit, 0 hit more than once`, with type (w = 0, 1, 2, 3, 4, 8 for max, ext10, ext11, ext12, ext13, spr) and `canon aut` agreeing with search A's tags on every line (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/check9-mapa-output.txt` (run outside the runner) |

The near misses, compared by `python src/nine_check_compare.py` from the logs above (`nine-check.md` §7): search A's 18
at dimension 5 lie in 7 configurations (mapped through `mapa`; 129, 153, 236, 348, 427, 671 and 674, all spr), and
each of search B's three orders has the same 18 there, the same number in each. At dimension 6, on every
configuration without a dimension-5 near miss, each order's count equals search A's. Search B's other near misses at
dimensions 6 and 7 contain a smaller one, a lower near miss extended by rank-6 forms in the holes, which cover.cpp
doesn't log (it prints `NEAR` only at the step where the last rank-4 form joins): of B's 3,922 / 3,768 dimension-6
spans in `asc` / `rand 7`, 395 contain no dimension-5 near miss, search A's 395, and of its 296 / 321 at dimension 7,
187 contain no smaller one, search A's 187. `desc`'s spans are incomplete (above), so the script counts only its
per-configuration totals there. The quadric's configuration (5662) has no near miss on either side.
`python src/nine_check_percfg.py asc rand7` (about a minute; `logs/rerun/t9-percfg.log`) makes the comparison
configuration by configuration: each `P1` span of search B is assigned to its configuration by the radicals of its
rank-4 elements, and the spans containing no smaller one of the same configuration match search A's `NEAR` lines on
every (configuration, dimension) pair: `asc configs 201 fresh by dim {5: 18, 6: 395, 7: 187} ... mismatches 0`, the
same for `rand7`. Since both searches cover every point off the holes before any hole (cover.cpp places solids first;
nine_check's `make_order` puts the off-hole points first in every order), these are exactly the spans of the nine
rank-4 forms.

The hole sets against Langevin and Hou's count of partial spread functions in eight variables (IEEE Trans. Inf.
Theory 57(4), 2011): download `fix-9.txt` from https://langevin.univ-tln.fr/project/spread/psp.html, then
`python src/nine_check_langevin.py fix-9.txt` (`logs/rerun/t9-langevin.log`): `classes: ours 5442, Langevin-Hou
5442`, `stabiliser orders: identical distribution`, `point sets: ours 24520650576127040978944, Langevin-Hou
24520650576127040978944, equal`.

### Eight (`i = 9`)

Kurz's table gives counts for eight solids too (9,316 hole sets, 181 maximal, the most of any size), so both searches
classify the 8-partial spreads themselves again. As for nine, times are running time between the runner's pauses, and
jobs often ran side by side. Search A, `pspread3.cpp` (the same program as for nine, with a mode added for the hole
sets) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread3 enum 8 p8.txt --ckpt p8.ckpt --wcheck 1000` | `w=0: reduced 1410915, weighted 1205406720, orbits 181, orbit sum 1205406720 MATCH`; `w=1: reduced 17113518, weighted 14517538560, orbits 2207, ... MATCH`; `w=2: reduced 30896844, weighted 25738030080, orbits 3881, ... MATCH`; `w=3: reduced 12507375, weighted 10173905280, orbits 1646, ... MATCH`; `w=4: reduced 2745378, weighted 2114259840, orbits 350, ... MATCH`; `w=5: reduced 3684669, weighted 2752216656, orbits 523, ... MATCH`; `w=9: reduced 2615301, weighted 1969810864, orbits 584, ... MATCH`; `orbits 9372`; `w spot checks: 70926, mismatches 0`; writes `p8.txt` (= `data/p8.txt`, tags `p8-<n>-a<\|Aut P\|>-w<w>`) | 87 min on an open box (the program's clock 5,185 s) | `logs/20261029-234926-p3-k8.log`, `logs/pspread8-output.txt` |
| `pspread3 selftest 8 data/p8.txt` | `selftest: 9372 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 49 s | `logs/rerun/pspread3-selftest8.log` (run for this folder) |
| `pspread3 holes 8 data/p8.txt holes8.txt`, then `python src/holes_tally.py holes8.txt data/p8.txt fix-8.txt 9316 181` (Langevin and Hou's `fix-8.txt`, downloaded from the page below) | `9372 configurations, keys distinct`; `holes: 9372 done this run, 0 with a cover off the list`; then `configurations 9372, hole-set classes 9316 (Kurz 9316)`, `maximal (hole set contains no solid, w = 0): 181 (Kurz 181)`, `configurations per hole set: {1: 9262, 2: 52, 3: 2}`, `\|Aut H\| distribution vs Langevin-Hou: identical (31 orders)`, `point sets: 46056096661467073413120 vs Langevin-Hou 46056096661467073413120: MATCH`, `problems: 0` | 23 s; the tally a second | `logs/20261030-030931-holes8.log`, `logs/holes8-output.txt`, `logs/holes8.txt` (one line per configuration: tag, index, `\|Aut P\|`, solids in the complement of the holes, covers, and `j:count` for the configurations the covers land on); the tally `logs/rerun/e8-tally.log` |
| `pspread3 holes 9 data/p9.txt holes9.txt`, then `python src/holes_tally.py holes9.txt data/p9.txt fix-9.txt 5442 1343` (a control on nine) | `holes: 5663 done this run, 0 with a cover off the list`; then `hole-set classes 5442 (Kurz 5442)`, `maximal ...: 1343 (Kurz 1343)`, `configurations per hole set: {1: 5259, 2: 152, 3: 25, 4: 5, 5: 1}`, `\|Aut H\| distribution vs Langevin-Hou: identical (35 orders)`, `point sets: 24520650576127040978944 ... MATCH`, as search B's hole sets gave (`nine_check_langevin.py`, above) | 25 min, nearly all on the quadric configuration's 1,920 covers | `logs/20261030-030931-holes9.log`, `logs/holes9-output.txt`, `logs/holes9.txt`; the tally `logs/rerun/e8-tally.log` |
| `cover_chunks.cmd DIR`, with `DIR\c0.txt` … `c7.txt` the lines of `data/p8.txt` in order, 1,172 to a chunk (the last 1,168), each run as `cover3 DIR\cK.txt > DIR\cK.out`; then `python src/cover_tally.py DIR/c?.out` | each chunk ends `done: covers 1172, ...` (the last `1168`), `solutions 0`; summed: `covers: 9,372`, `tasks: 262,416`, `nodes: 780,361,676`, `solutions: 0`, `full (dim < 8): 0`, `dim 8 non-APN: 0`, `bad seeds: 0`; `nodes by span dimension 1..8: 262,416 / 7,347,648 / 122,350,038 / 519,227,824 / 130,741,783 / 429,756 / 2,211 / 0`; `near misses 10037 by dim {5: 47, 6: 9061, 7: 929}`, in `545` configurations, none reaching dim 8 | 46 min, on the P-cores at most 40% (2,770 s by the chunks' clocks) | `logs/20261030-031058-cover8.log`, `logs/eight-output-c0.txt` … `-c7.txt`; the tally `logs/rerun/e8-tally.log` |

Search B, `eight_check` (`eight-check.md`; fresh `.prog` files as for fourteen), its jobs run by `eight_check_slices.py`
(the job logs give its commands). `parta` runs one slice of a minimal class per job; the 34 outputs are
`data/parta8-c*-s*.txt`, its lists of maximal cliques aren't included (`Data` below):

| Command | Expected | Time | Log |
|---|---|---|---|
| `eight_check count5` | `count5: k-cliques of V5824 by plain enumeration: 1:5824 2:4583936 3:606183424 4:14122766336 5:58471168000 6:47414189312  (6-cliques: nine-check's weighted count 47414189312 OK)`: every clique once, by increasing index, no symmetry used | 12 min of running time (the program's clock, 2,337 s, counts paused time too) | `logs/20261030-121251-check8-count5.log`, `logs/check8-count5-output.txt` |
| `eight_check parta C parta-cC-sS.txt maxc-cC-sS.txt S N` for class C = 0 … 4 in N = 16, 8, 4, 4, 2 slices | per slice `class C slice S/N ...: R reduced 5-cliques, M maximal`; by class R = 51,510,784, 14,191,420, 4,124,806, 1,035,490, 111,500 (70,974,000 in all) and M = 1,032,960, 297,840, 68,865, 11,250, 0 (1,410,915); each file ends in `END`; `parta-cC-sS.txt` = `data/parta8-cC-sS.txt` | 70 min of running time, one slice after another; by the program's clock 118 to 340 s per slice of class 0 (the first, which spans a hold, 1,943 s), 67 to 79 s of class 1, at most 44 s for the rest | `logs/20261030-121251-check8-parta.log`, `logs/check8-parta-c*-s*-output.txt` |
| `eight_check merge orbits.txt data/configs9.txt data/parta8-c0-s0.txt ... data/parta8-c4-s1.txt` (the 34 in class and slice order) | `nine-check's configurations: 5663; sum 504*20160/\|Aut\| = 47414189312 (its Part A: 47414189312)`; `reduced 5-cliques 70974000, of them maximal (Q empty) 1410915; weighted 5-cliques 58471168000`; `sum over 5-cliques of \|Q\| = 284485135872; 6 * (weighted 6-cliques, nine-check) = 284485135872  OK`; `orbits 9372; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 0:181 1:2207 2:3881 3:1646 4:350 5:523 9:584`; `orbits by type: spr 584 ext13 523 ext12 350 ext11 1646 ext10 3881 ext9 2207 max 181`; `type / Q / n9 inconsistencies 0; P + q on no nine-check orbit 0`; `... C(14,5) * 1128864 = 2259985728  OK`; `per 9-orbit O: ... fails for 0 of 5663`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits8.txt`) | 6 s | `logs/20261030-135032-check8-merge.log`, `logs/check8-merge-output.txt` |
| `eight_check configs orbits.txt data/configs9.txt data/configs10.txt data/configs11.txt data/configs12.txt data/spreads.txt data/configs.txt configs8.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 9316, containing no solid (maximal partial spreads): 181   [Kurz Table 9.3, r = 8: 9316 and 181]`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 9497`; `orbits reached by no cover: 0`; each of the 181 maximal hole sets has `1 covers`; outside sources `nine-check's 5663 x 9: 50967 configurations; on no listed orbit 0; ... spr 584/584 ext13 523/523 ext12 350/350 ext11 1646/1646 ext10 3881/3881 ext9 2207/2207 max 0/181`, `fourteen-check's spread reps x 3003: 219219 ...; on no listed orbit 0; ... spr 584/584`, `thirteen-check's maximal 13s x 1287: 39897 ...; on no listed orbit 0; ... spr 474/584 ext13 523/523`, `twelve-check's 376 x 495: 186120 ...; on no listed orbit 0; ... spr 584/584 ext13 523/523 ext12 350/350`, `eleven-check's 817 x 165: 134805 ...; on no listed orbit 0; ... ext11 1646/1646`, `ten-check's 1777 x 45: 79965 ...; on no listed orbit 0; ... ext10 3881/3881`; writes `configs8.txt` (= `data/configs8.txt`; spr 0..583, ext13 584..1106, ext12 1107..1456, ext11 1457..3102, ext10 3103..6983, ext9 6984..9190, max 9191..9371) | 1 min | `logs/20261030-135104-check8-configs.log`, `logs/check8-configs-output.txt` |
| `python src/eight_check_mass.py fix-8.txt` | `(a) hole classes 9316, containing no solid 181   [Kurz r = 8: 9316, 181]  OK`; `Langevin-Hou: 9316 classes over 31 fixator orders; ours: 9316 over 31`; `(b) classes per \|Aut H\| against Langevin-Hou: 0 orders differ []; their qte = classes * \|GL\| / order: OK`; `sum_H \|GL\|/\|Aut H\| = 46056096661467073413120; Langevin-Hou Total = 46056096661467073413120  OK`; `(c) labelled 8-partial spreads: from covers 46164543049102786560000; from configurations 46164543049102786560000; from Part A's W5 (58471168000) 46164543049102786560000; count5's plain 5-clique count 58471168000 = W5  OK` | seconds | `logs/check8-mass-output.txt`; the re-run `logs/rerun/e8-tally.log` |
| `eight_check selftest data/configs8.txt` | `selftest: ND4 28 (28), ND6 13888 (13888), GL4 20160 (20160), V5824 5824`; `lines in D over the configurations: 1575:9372`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28)`; `selftest: 9372 configurations, 1620745728 generated forms checked, 0 failures` | 56 min, beside the searches below | `logs/20261030-135221-check8-selftest.log`, `logs/check8-selftest-output.txt` |
| `eight_check search data/configs8.txt asc-rK.prog asc FIRST LAST` in 8 ranges of 1,172 configurations (the last 1,168) | `solutions: 0` in each; over the ranges, nodes by dim `9372 262416 7347648 122340087 518463075 130692392 431607 2279 0`; nodes with every off-D point covered (near misses) `47` / `11928` / `2013` at dims 5 / 6 / 7; tries at points of D `3872766016` / `304272192` at dims 6 / 7 (4,197,356,352 in all); `distinct P1 spans 13988` | the three orders one after another, 355 min of running time in all (2 h, 1.2 h and 3.8 h from start to end; 7,172, 4,428 and 13,618 s by the ranges' clocks) | `logs/20261030-135220-check8-search.log`, `logs/check8-asc-output.summary.txt` (each range's header and summary), per-configuration `R` lines and near-miss spans (`P1`) `logs/check8-asc-prog.txt` |
| `eight_check search data/configs8.txt desc-rK.prog desc FIRST LAST` (the same ranges) | `solutions: 0`; nodes by dim `9372 262416 7347648 121719264 513696565 127995358 458955 2216 0`; near misses `47` / `12027` / `1947`; tries at points of D `3908916480` / `295397760`; `distinct P1 spans 14021` | (same job) | same job log, `logs/check8-desc-output.summary.txt`, `logs/check8-desc-prog.txt` |
| `eight_check search data/configs8.txt rand7-rK.prog rand 7 FIRST LAST` | `solutions: 0`; nodes by dim `9372 262416 7347648 121741683 514270043 128551591 429611 2235 0`; near misses `47` / `12007` / `1971`; tries at points of D `3872682688` / `289356480`; `distinct P1 spans 14025` | (same job) | same job log, `logs/check8-rand7-output.summary.txt`, `logs/check8-rand7-prog.txt` |
| `python src/nine_check_tally.py prog data/configs8.txt 9372 logs/check8-asc-prog.txt` (`desc`, `rand7`) | `configurations done 9372 (expect 9372); duplicate R lines 0`; the totals above, and by type, e.g. for `asc` `spr (584)` near misses `25 / 9,700 / 1,851` at dims 5 / 6 / 7 and `max (181)` `0 / 3 / 0`; `distinct P1 spans 13988` / `14021` / `14025`; `SOLUTION lines 0` (`eight-check.md` §5's table) | seconds | `logs/rerun/e8-tally.log` |
| `python src/eight_check_compare.py orders` | `configurations 9372; dim-5 near-misses differ across orders on 0; 15 configurations with a dim-5 near-miss (47 in all); dim-6 counts differ on 0 of the others; dim-6 near-misses on configurations without a dim-5 one: {'asc': 8192, 'desc': 8192, 'rand7': 8192}; configurations with a dim-6 near-miss: asc 525, desc 525, rand7 525, same set: True` | seconds | `logs/rerun/e8-tally.log` |
| `python src/eight_check_numpy.py search data/configs8.txt I` for 32 configurations: the first of each type (0, 584, 1107, 1457, 3103, 6984, 9191), the 15 with a dimension-5 near miss (114, 240, 241, 291, 309, 328, 497, 522, 552, 555, 905, 941, 1214, 1377, 2453), the three heaviest others (465, 557, 583) and seven drawn at random (717, 1395, 2069, 3164, 3714, 6068, 6150); then `python src/eight_check_compare.py numpy` | per configuration `config I: nodes [...] p1 [...] triesD [...] full 0 allcovered-below-8 0`; then `numpy re-search: 32 configurations agree with asc in every count, 0 differ` (nodes, near misses and tries at points of D at every dimension; 1.49 billion tries in D over the 32) | 133 min of running time, eight at a time (configuration 114 alone 10,079 s by its clock) | `logs/20261030-165835-check8-numpy.log`, `logs/check8-numpy-output.txt`, `logs/check8-numpy-compare.txt` |
| `python src/eight_check_numpy.py control apn_8bit.txt 1 40 1`, `... control new_apns.txt 1 20 2` (planted: a known function's rank-4 radicals plus random solids, no classification) | each control `T re-found: YES`; `controls: 40 run, 40 re-found, 0 skipped ...`, `controls: 20 run, 20 re-found, 0 skipped ...` | 184 min of running time for the two, side by side | `logs/20261030-165841-check8-npctl.log`, `logs/check8-numpy-ctl-bl22-output.txt`, `logs/check8-numpy-ctl-new-output.txt` |
| `eight_check searchmax data/maxc8-s100.txt maxs-rK.prog asc FIRST LAST` in 4 ranges of 3,528 (the last 3,526), where `maxc8-s100.txt` is lines 1, 101, 201, … of `maxc-c0-s0.txt` … `maxc-c4-s1.txt` concatenated in class and slice order (14,110 maximal cliques; `eight_check_slices.py maxsample 100`) | `solutions: 0` in each; summed over the ranges, nodes by dim `14110 395080 11062240 183060214 768416578 185900152 157625 0 0`; `279` spans with every off-D point covered, all at dim 6; 76,925,632 tries at points of D | 116 min | `logs/20261030-135232-check8-maxs.log`, `logs/check8-maxs-output.summary.txt` |
| `eight_check orbitcov orbits.txt configs8.txt maxc-cC-sS.txt` for each of the 34 files, appended to `maxall-orbits.out`; then `python src/ten_check_tally.py maxall DIR` on a folder holding `orbits.txt` and `maxall-orbits.out` | each `orbitcov: N cliques, 0 on no listed orbit, 0 on a non-max orbit`; then `orbitcov runs 34 (expect 17); cliques off the max orbits 0; max orbits 181, cliques 1410915 (Part A 1410915); orbits whose tally differs from Part A's n: 0` (the script's `expect 17` is ten's number) | 3 min | `logs/20261030-135232-check8-orbitcov.log`, `logs/check8-maxall-orbits-output.txt`; the tally `logs/rerun/e8-tally.log` |
| `eight_check search9 data/configs9.txt s9-rK.prog asc FIRST LAST` in 8 ranges, over nine_check's 5,663 configurations × 9 members removed (50,967: every 8 that extends, without this program's classification) | `solutions: 0` in each; from `python src/nine_check_tally.py prog - 50967 s9-r*.prog`: `configurations done 50967 (expect 50967); duplicate R lines 0`, nodes by dim `50967 1427076 39958128 665050089 2829413970 738756872 5141694 29189 0`; near misses `339` / `127820` / `25537` at dims 5 / 6 / 7; tries at points of D `42950042688` / `3755842944` at dims 6 / 7 (46,849,182,016 in all); bent extension 0, full spans 0; `distinct P1 spans 128107`; `SOLUTION lines 0` | 12.9 h from start to end, the ranges one after another (46,271 s by the ranges' clocks), on a shared box | `logs/20261030-125207-check8-s9.log`, `logs/check8-s9-output.summary.txt` (each range's header and summary) |
| `eight_check ctlbatch data/configs8.txt apn_8bit.txt 2 1000000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1440 controls; re-found 1440`, `nonext: 1421 controls; re-found 1421`, `other spaces found 0`; `started from S0 of dimension: 7:1440` and `7:1421`, every control finding its form among the 13,888 on its line (188,487,936 and 184,002,112 tries, all at dimension 7); 19 `k = 4` functions skipped on the nonext path (`no member 4-tuple of the class`) | 1 min | `logs/20261030-135221-check8-controls.log` (the job ran all six batches), `logs/check8-ctl-ind-bl22-output.txt`, `logs/check8-ctl-ind-bl22.txt` |
| `eight_check ctlbatch data/configs8.txt apn_8bit.txt 2 1000000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 432`; `nonext: 1421 controls; re-found 1421 ...`, `reps hit: 1421`, with the same 19 skipped | 29 min | same job log, `logs/check8-ctl-bl22-output.txt`, `logs/check8-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222 (reps hit 160 / 220), `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths (the job exits 3 because of the last) | 10 min for the three | same job log, `logs/check8-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check8-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 513`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 4904` | 2.4 h | same job log, `logs/check8-ctl-new-output.txt`, `logs/check8-ctl-new.txt` |
| `python src/ten_check_tally.py reps DIR` on a folder holding `configs.txt` (= `data/configs8.txt`) and the `ctl-*.txt` files above | `reps hit 6170 of 9372; not reached 3202 {'spr': 62, 'ext13': 190, 'ext12': 123, 'ext11': 580, 'ext10': 1359, 'ext9': 819, 'max': 69}` and the 3,202 | seconds | `logs/check8-reps-unreached.txt`; the re-run `logs/rerun/e8-tally.log` |
| `eight_check ctlreps data/configs8.txt apn_8bit.txt 1 40 91 92 ... 952` (the first 200 of those 3,202), cancelled; then the same with 10 functions, cancelled | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the first 36 (91 to 460, all spr) before the first was stopped as too slow on the spr configurations, then `10 of 10` on the first six of them again (91, 92, 93, 95, 108, 109) | 176 min and 12 min | `logs/20261030-165749-check8-ctlreps.log`, `logs/20261030-205639-check8-ctlreps10.log`, `logs/check8-ctlreps40-partial-output.txt`, `logs/check8-ctlreps10-partial-output.txt` |
| `eight_check ctlrot data/configs8.txt apn_8bit.txt 1 40 REPS` on the 3,202 in chunks of 400 (one control per configuration: the first 40 functions of the list with k ≥ 3 in rotation, parallel over the configurations) | `rep N (ext)` or `(nonext): line L, re-found, 0 other spaces` for each of the 3,202, all 40 functions used | 36 min | `logs/20261030-210958-check8-ctlrot.log`, `logs/check8-ctlrot-output.txt` |
| `eight_check mapa data/configs8.txt data/p8.txt` | `mapa: 9372 of A's configurations, 0 on none of ours, 9372 of our 9372 configurations hit, 0 hit more than once`, with type (w = 0, 1, 2, 3, 4, 5, 9 for max, ext9, ext10, ext11, ext12, ext13, spr) and `canon aut` agreeing with search A's tags on every line (search A's list mapped onto search B's, after both verdicts) | 4 s | `logs/20261030-215021-check8-mapa.log`, `logs/check8-mapa-output.txt` |

The near misses, compared by `python src/eight_check_vs_a.py` from the logs above (`eight-check.md` §7): search A's 47
at dimension 5 lie in 15 configurations (mapped through `mapa`; 114, 240, 241, 291, 309, 328, 497, 522, 552 and 555,
all spr, the ext13 905 and 941, the ext12 1214 and 1377, and the ext11 2453), and each of search B's three orders has
the same 47 there, the same number in each. At dimension 6, on every configuration without a dimension-5 near miss,
each order's count equals search A's, and at dimension 7 on every configuration without a lower one. Search B's other
near misses contain a smaller one, a lower near miss extended by rank-6 forms in the holes, which cover.cpp doesn't
log: of B's 11,928 / 12,027 / 12,007 dimension-6 spans in `asc` / `desc` / `rand 7`, 9,061 contain no dimension-5 near
miss, search A's 9,061, and of its 2,013 / 1,947 / 1,971 at dimension 7, 929 contain no smaller one, search A's 929.
Unlike at nine, `desc`'s spans are complete. `python src/eight_check_percfg.py` (written for this folder from
`nine_check_percfg.py`; about a minute per order, `logs/rerun/e8-percfg.log`) makes the comparison configuration by
configuration: each `P1` span of search B is assigned to its configuration by the radicals of its rank-4 elements, and
the spans containing no smaller one of the same configuration match search A's `NEAR` lines on every
(configuration, dimension) pair, in all three orders: `asc configs 545 fresh by dim {5: 47, 6: 9061, 7: 929} ...
per-(config,dim) mismatches 0`, the same for `desc` and `rand7`. As at nine, both searches cover every point off the
holes before any hole, so these are exactly the spans of the eight rank-4 forms.

The hole sets against Langevin and Hou's count of partial spread functions in eight variables: download `fix-8.txt`
from https://langevin.univ-tln.fr/project/spread/psp.html; both `python src/eight_check_mass.py fix-8.txt` (search B's
hole sets, from `data/configs8.txt`) and `python src/holes_tally.py logs/holes8.txt data/p8.txt fix-8.txt 9316 181`
(search A's, from its `holes` run) give 9,316 classes with the same number for each of Langevin and Hou's 31 stabiliser
orders, and `46056096661467073413120` point sets, their total (`logs/rerun/e8-tally.log`).

### Seven (`i = 10`)

Kurz's table gives counts for seven solids too (3,726 hole sets, none maximal: every 7-partial spread extends), so both
searches classify the 7-partial spreads themselves again. As for eight, times are running time between the runner's
pauses (a program's own clock counts the pauses too), and jobs often ran side by side, mostly on the E-cores. Search A,
`pspread3.cpp` (the same program as for eight) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread3 enum 7 p7.txt --ckpt p7.ckpt --wcheck 1000` | `w=2: reduced 340542, weighted 324399600, orbits 83, orbit sum 324399600 MATCH`; `w=3: reduced 5536992, weighted 5253897600, orbits 1277, ... MATCH`; `w=4: reduced 5838846, weighted 5457160800, orbits 1375, ... MATCH`; `w=5: reduced 1312029, weighted 1167415200, orbits 305, ... MATCH`; `w=6: reduced 1322286, weighted 1134163800, orbits 332, ... MATCH`; `w=10: reduced 938010, weighted 785729336, orbits 359, ... MATCH` (no w = 0 or 1: no maximal 7-partial spread, and none extending only to eight); `orbits 3731`; `w spot checks: 15245, mismatches 0`; writes `p7.txt` (= `data/p7.txt`, tags `p7-<n>-a<\|Aut P\|>-w<w>`) | 9.5 min (the program's clock, 840 s, counts a 5-minute pause too) | `logs/20261031-032113-p3-k7.log`, `logs/pspread7-output.txt` |
| `pspread3 selftest 7 data/p7.txt` | `selftest: 3731 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 67 s on the E-cores | `logs/rerun/seven/p3-selftest7.log` (run for this folder) |
| `pspread3 holes 7 data/p7.txt holes7.txt`, then `python src/holes_tally.py holes7.txt data/p7.txt fix-7.txt 3726 0` (Langevin and Hou's `fix-7.txt`, downloaded from the page below) | `3731 configurations, keys distinct`; `holes: 3731 done this run, 0 with a cover off the list`; then `configurations 3731, hole-set classes 3726 (Kurz 3726)`, `maximal (hole set contains no solid, w = 0): 0`, `configurations per hole set: {1: 3721, 2: 5}`, `\|Aut H\| distribution vs Langevin-Hou: identical (27 orders)`, `point sets: 17833337132662061531136 vs Langevin-Hou 17833337132662061531136: MATCH`, `problems: 0` | 14 s; the tally a second | `logs/20261031-052136-holes7.log`, `logs/holes7-output.txt`, `logs/holes7.txt`; the tally `logs/rerun/seven/s7-tally.log` |
| `cover_chunks.cmd DIR`, with `DIR\c0.txt` … `c7.txt` the lines of `data/p7.txt` in order, 467 to a chunk (the last 462), each run as `cover3 DIR\cK.txt > DIR\cK.out`; `c5.txt` … `c7.txt` were moved to a second folder and run by a second job while the first was on `c4` | each finished chunk ends `done: covers 467, tasks 13076, ...` (`c7`: `covers 462, tasks 12936`), `solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`. `c4` was cancelled after about 17 hours, most of them on one task, almost surely `p7-2222-a168-w4`'s (the one configuration whose seven rank-4 forms span a 4-space), and redone in the two runs below | c0 to c3 89 min, c5 to c7 77 min (5,335 and 4,608 s by the chunks' clocks) | `logs/20261031-052243-cover7.log`, `logs/20261031-205800-cover7b.log`, `logs/seven-output-c0.txt` … `-c3.txt`, `-c5.txt` … `-c7.txt` |
| `cover4 c4x.txt`, with `c4x.txt` the chunk `c4.txt` (lines 1,869 to 2,335 of `data/p7.txt`) without `p7-2222-a168-w4` | `done: covers 466, tasks 13048, nodes 72202053, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0` | 110 min (7,454 s by its clock) | `logs/20261031-235550-cover7c.log`, `logs/seven-output-c4x.txt` |
| `cover4 p2222.txt -1 5`, with `p2222.txt` the line `p7-2222-a168-w4` of `data/p7.txt` (split at dimension 5: Build, above) | `phase 1 done 8 s: frontier 14622 items at dim 5`; `done: covers 1, tasks 28, nodes 772887, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`; `nodes by span dimension 1..8: 28 784 14140 88037 78681 591193 24 0`; 58 near misses: the dimension-4 one and 57 at dimension 6 | 5.4 h on the E-cores (22,646 s by its clock) | `logs/20261031-235344-cover4-p2222.log`, `logs/seven-output-p2222-s5.txt` |
| `python src/cover_tally.py logs/seven-output-c[0-3].txt logs/seven-output-c[5-7].txt logs/seven-output-c4x.txt logs/seven-output-p2222-s5.txt` (the nine outputs) | `9/9 outputs done`; `covers: 3,731`, `tasks: 104,468`, `nodes: 568,755,235`, `solutions: 0`, `full (dim < 8): 0`, `dim 8 non-APN: 0`, `bad seeds: 0`; `nodes by span dimension 1..8: 104,468 / 2,925,104 / 52,909,076 / 313,666,564 / 196,375,985 / 2,762,176 / 11,862 / 0`; `near misses 114702 by dim {4: 1, 5: 1634, 6: 108026, 7: 5041}` in `2787` configurations, none reaching dim 8; `lowest near-miss dim 4 ['p7-2222-a168-w4']`; the most, 2,333, on `p7-3728-a60-w10` | seconds; the nine runs about 10 h of running time (40,043 s by their clocks) | `logs/rerun/seven/s7-tally.log` |

The split run's frontier is the 14,622 legal spans of dimension 5 below `p7-2222`'s dimension-4 near miss, the number
search B's coset sweep finds on the same configuration (its 1527, below), and its 591,193 nodes at dimension 6 are
search B's count there too.

Search B, `seven_check` (`seven-check.md`), its jobs run by `seven_check_slices.py` (the job logs give its commands).
`parta` runs one slice of a minimal class per job; its 37 outputs aren't included (`Data` below), as `parta` remakes
them in minutes:

| Command | Expected | Time | Log |
|---|---|---|---|
| `seven_check count4` | `count4 (a): k-cliques of V5824 by plain enumeration: 1:5824 2:4583936 3:606183424 4:14122766336 5:58471168000`: every clique once, by increasing index, no symmetry used; `count4 (b): triangles 606183424, 4-cliques 14122766336 (by edges)`; `count4: 5-cliques = eight-check's 58471168000 OK; (a) and (b) agree OK` | 2 min | `logs/20261031-054609-check7-count4.log`, `logs/check7-count4-output.txt` |
| `seven_check parta C parta-cC-sS.txt maxc-cC-sS.txt S N` for class C = 0 … 4 in N = 16, 8, 8, 4, 1 slices | per slice `class C slice S/N ...: R reduced 4-cliques, 0 maximal`; by class R = 9,903,360, 3,654,801, 1,348,739, 360,725, 21,080 (15,288,705 in all), none maximal (so every `maxc` file is empty); each file ends in `END` | 9 min on an open box, one slice after another (at most 31 s a slice by the program's clock) | `logs/20261031-054609-check7-parta.log`, `logs/check7-parta-c*-s*-output.txt` |
| `seven_check merge orbits.txt data/configs8.txt parta-c0-s0.txt ... parta-c4-s0.txt` (the 37 in class and slice order) | `eight-check's configurations: 9372; sum 336*20160/\|Aut\| = 58471168000 (its count5: 58471168000)`; `reduced 4-cliques 15288705, of them maximal (Q empty) 0; weighted 4-cliques 14122766336`; `sum over 4-cliques of \|Q\| = 292355840000; 5 * (weighted 5-cliques, eight-check) = 292355840000  OK`; `orbits 3731; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 2:83 3:1277 4:1375 5:305 6:332 10:359`; `orbits by type: spr 359 ext13 332 ext12 305 ext11 1375 ext10 1277 ext9 83 ext8 0 max 0`; `type / Q / n10 inconsistencies 0; P + q on no eight-check orbit 0`; `... C(14,4) * 1128864 = 1129992864  OK`; `per 8-orbit O: ... fails for 0 of 9372`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits7.txt`) | 4 s | `logs/20261031-055608-check7-merge-configs.log`, `logs/check7-merge-output.txt` |
| `seven_check configs orbits.txt data/configs8.txt data/configs9.txt data/spreads.txt configs7.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 3726, containing no solid (maximal partial spreads): 0   [Kurz Table 9.3, r = 7: 3726 and 0]`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 3757`; `orbits reached by no cover: 0`; outside sources `eight-check's 9372 x 8: 74976 configurations; on no listed orbit 0; ... orbits reached: spr 359/359 ext13 332/332 ext12 305/305 ext11 1375/1375 ext10 1277/1277 ext9 83/83`, `nine-check's 5663 x 36: 203868 ...; on no listed orbit 0; ...` (all 3,731 reached), `fourteen-check's spread reps x 3432: 250536 ...; on no listed orbit 0; ... spr 359/359`; `orbits reached by no deletion from eight-check's 9372: 0`; writes `configs7.txt` (= `data/configs7.txt`; spr 0..358, ext13 359..690, ext12 691..995, ext11 996..2370, ext10 2371..3647, ext9 3648..3730) | under a minute (same job) | same job log, `logs/check7-configs-output.txt` |
| `python src/seven_check_mass.py fix-7.txt` | `(a) hole classes 3726, containing no solid 0   [Kurz r = 7: 3726, 0]  OK`; `Langevin-Hou: 3726 classes over 27 fixator orders; ours: 3726 over 27`; `(b) classes per \|Aut H\| against Langevin-Hou: 0 orders differ []; their qte = classes * \|GL\| / order: OK`; `sum_H \|GL\|/\|Aut H\| = 17833337132662061531136; Langevin-Hou Total = 17833337132662061531136  OK`; `(c) labelled 7-partial spreads: from covers 17840479724727007444992; from configurations 17840479724727007444992; from Part A's W4 (14122766336) 17840479724727007444992; count4's plain 4-clique count 14122766336 = W4  OK` | seconds | `logs/rerun/seven/s7-tally.log` |
| `seven_check selftest data/configs7.txt` | `selftest: ND4 28 (28), ND6 13888 (13888), GL4 20160 (20160), V5824 5824`; `lines in D over the configurations: 2150:3731`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28); sweep Pfaffian test vs rank8 on every 16th form, sweep kernel and the Gray walk vs radical() on every 256th: 0 disagreements`; `selftest: 3731 configurations, 896507276 generated forms checked, 0 failures` | 26 min, beside the searches below | `logs/20261031-060729-check7-selftest.log`, `logs/check7-selftest-output.txt` |
| `seven_check search data/configs7.txt phase1-rK.prog desc FIRST LAST phase1` in 8 ranges of 467 configurations (the last 462): the part off the holes alone, every configuration's cores | over the ranges, cores by dim `0 0 0 0 1 1634 108026 5041 0` (114,702); 944 configurations have no core | 95 min | `logs/20261031-121312-check7-phase1.log`, `logs/check7-phase1-output.summary.txt` (each range's header and summary), per-configuration `R` lines `logs/check7-phase1-R.txt` |
| `seven_check search data/configs7.txt asc-rK.prog asc FIRST LAST cores asc-rK.cores` (the same ranges) | `solutions: 0` in each; over the ranges, nodes by dim `3731 104468 2925104 52895372 312847548 195494917 2767276 11780 0`; nodes with every off-D point covered `1` / `16256` / `843914` / `11780` at dims 4 / 5 / 6 / 7, of them cores `1` / `1634` / `108026` / `5041`; tries at points of D `7166208` / `1266187487` / `97176231257` / `1888261132` at dims 4 / 5 / 6 / 7 (100,337,846,084 in all); bent extension 0, full spans 0 | 8.9 h, in three jobs: configurations 0 to 10 in a first, cancelled after 10 minutes, the next resuming after them; that one cancelled at configuration 1527 for the candidate caches (Build, above), a third resuming there | `logs/20261031-055731-check7-asc.log`, `logs/20261031-060729-check7-asc.log`, `logs/20261031-141646-check7-asc-2.log`, `logs/check7-asc-output.summary.txt`, `logs/check7-asc-R.txt` |
| `seven_check search data/configs7.txt desc-rK.prog desc FIRST LAST cores desc-rK.cores` | `solutions: 0`; nodes by dim `3731 104468 2925104 52669145 310556440 193373519 2772459 11924 0`; every off-D point covered `1` / `15125` / `858354` / `11924`, cores `1` / `1634` / `108026` / `5041`; tries at points of D `6569024` / `1287225172` / `97014924117` / `1925042487` (100,233,760,800) | 9.0 h, in two jobs (the second from configuration 1527, with the caches) | `logs/20261031-060729-check7-desc.log`, `logs/20261031-141647-check7-desc-2.log`, `logs/check7-desc-output.summary.txt`, `logs/check7-desc-R.txt` |
| `seven_check search data/configs7.txt sweep-rK.prog asc FIRST LAST dsweep` (the coset sweep in the holes) | `solutions: 0`; the part off the holes as `asc`'s (nodes by dim `3731 104468 2925104 52895372 312847548 195480295 2031388 5041 0`, cores as above), then the sweep's legal spans by dim `14622` / `735888` / `6739` at dims 5 / 6 / 7, none at 8; 571,009,882,593 element checks; bent extension 0, full spans 0 | 10.7 h, in two jobs (the second from configuration 135, with the sweep's cache) | `logs/20261031-062222-check7-sweep.log`, `logs/20261031-141647-check7-sweep-2.log`, `logs/check7-sweep-output.summary.txt`, `logs/check7-sweep-R.txt` |
| `python src/seven_check_tally.py prog data/configs7.txt 3731 logs/check7-asc-R.txt` (`desc`, `sweep`, `phase1`) | `3731 configurations, 0 duplicate records, 0 missing [], 0 out of range`; the totals above, also by type (`seven-check.md` §5); `SOLUTION lines: 0` | seconds | `logs/rerun/seven/s7-tally.log` |
| `python src/seven_check_tally.py orders logs` | `runs: asc 3731, desc 3731, sweep 3731, phase1 3731; reference phase1`; `asc`, `desc`, `sweep`: `cores by dimension and fingerprint on 3731 configurations: 0 disagreements`; `asc nodes by dimension = sweep's phase-1 nodes + sweep's legal spans: 0 disagreements over 3731`. Given the core files (not here), also `114702 distinct cores (asc) ...; asc and desc sets differ on 0 []; asc file count != R-line count on 0` (one truncated record, from a killed run, skipped) | seconds | `logs/rerun/seven/s7-tally.log` |
| `python src/seven_check_numpy.py search data/configs7.txt I` for 18 configurations: the first of each type (0, 359, 691, 996, 2371, 3648), 3, 9, 41, 47 and 49 (the last four with cores of dimension 5 or 7) and seven drawn at random (1326, 617, 1617, 2666, 197, 296, 3363); then `python src/seven_check_tally.py numpy logs` | per configuration `config I: nodes [...] p1 [...] triesD [...] cores [...] fp F full 0 allcovered-below-8 0`; then `numpy: 18 configurations compared, 0 differ, 0 waiting for asc; tries in D over those compared 2,037,841,792` | 4.9 h of running time, six at a time | `logs/20261031-063551-check7-numpy.log`, `logs/check7-numpy-output.txt`; the tally `logs/rerun/seven/s7-tally.log` |
| `python src/seven_check_numpy.py control apn_8bit.txt 1 20 7`, `... control new_apns.txt 1 10 7` (planted: a known function's rank-4 radicals plus random solids, no classification) | each control `T re-found: YES`; `controls: 20 run, 20 re-found, 0 skipped ...`, `controls: 10 run, 10 re-found, 0 skipped ...` | 5.8 h for the two, one after the other | `logs/20261031-065338-check7-npctl.log`, `logs/check7-numpy-ctl-bl22-output.txt`, `logs/check7-numpy-ctl-new-output.txt` |
| `seven_check search8 data/configs8.txt s8-bB.prog asc B B+23` for B = 4728, 7104, 9240, 14808, 31824, 35928, 38808, 52656, 57288, 63984 (`seven_check_slices.py crosssample 10 7`: 10 random blocks of 24 consecutive inputs of the 9,372 × 8 deletions, 30 of eight_check's configurations with each of their 8 members removed, taken as they are, no canonical form) | `solutions: 0` in each; over the blocks (`python src/seven_check_tally.py prog - 240 logs/check7-s8-R.txt`, whose `missing` count is against indices 0 to 239, not the blocks) nodes by dim `240 6720 188160 3401106 20184575 12679118 125400 229 0`, cores `97` / `5529` / `190` at dims 5 / 6 / 7 (5,816), 2,274,531,950 tries at points of D, bent extension 0, full spans 0 | 8 min | `logs/20261031-171708-check7-s8sample.log`, `logs/check7-s8-output.txt`, `logs/check7-s8-R.txt` |
| `seven_check ctlbatch data/configs7.txt apn_8bit.txt 2 1000000 3 ctl-ind-bl22.txt nond` (from the second build: Build, above) | the D phase exercised from a hyperplane of the function: `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `started from S0 of dimension: 7:1440`; `nonext: 1397 controls; re-found 1397 ...`, `7:1397`; 43 `k = 4` functions skipped on the nonext path (`no member 4-tuple of the class`) | under 27 min (the job was cancelled after this batch, to change the others' design) | `logs/20261031-062245-check7-controls.log`, `logs/check7-ctl-ind-bl22-output.txt`, `logs/check7-ctl-ind-bl22.txt` |
| `seven_check ctlbatch data/configs7.txt LUT 2 MAXN CH ctl-NAME.txt OPT` for (NAME, LUT, CH, OPT) = (`ind2-bl22`, `apn_8bit.txt`, 3, `nond2+guided`), (`bl22`, `apn_8bit.txt`, 3, `guided`), (`blp22`, `blp22_8bit.txt`, 6, `guided`), (`yuperrin`, `yuperrin_qam.txt`, 6, `guided`), (`yuwangli`, `yuwangli_qam.txt`, 6, `guided`), (`new`, `new_apns.txt`, 1, `guided`), (`ind2-new`, `new_apns.txt`, 1, `nond2+guided`), (`swind-bl22`, `apn_8bit.txt`, 1, `nond+dsweep+guided`), (`swind2-bl22`, `apn_8bit.txt`, 1, `nond2+dsweep+guided`), MAXN 1000000; and (`full-bl22`, `apn_8bit.txt`, 1, none) with MAXN 20 | `ind2-bl22` and `bl22`: `1440` of 1440 and `1397` of 1397 (the same 43 skipped), `reps hit: 305` / `1389`, two-step starts `6:1440` / `6:1397`; `blp22` `222` and `222`; `yuperrin` `24` and `24`; `yuwangli` `72` of 78 on each path with 6 setup failures (the list's line 3, APN but not quadratic; the job exits 3 because of these); `new` and `ind2-new` `7199` and `7199`, `reps hit: 341` / `2974`; `swind-bl22` `479` and `468` (2 skipped: `S0 already the whole space`, 11 `k = 4` without a member of the class); `swind2-bl22` `471` and `459` (19 and 11 skipped likewise); `full-bl22` `20` and `20`; in every batch `NOT found 0`, `other spaces found 0` | 8.5 h of running time for the job, which then ran the two steps below | `logs/20261031-065328-check7-controls2.log`, `logs/check7-ctl-*-output.txt`, `logs/check7-ctl-*.txt` |
| `python src/seven_check_tally.py reps logs data/configs7.txt` (the job ran it on `runs\check7`) | `configurations reached by the control batches: 3553 of 3731; unreached 178 {'spr': 16, 'ext13': 25, 'ext12': 15, 'ext11': 61, 'ext10': 59, 'ext9': 2}`; here also `the same list as check7-reps-unreached.txt: True` | seconds | `logs/check7-reps-unreached.txt`; the re-run `logs/rerun/seven/s7-tally.log` |
| `seven_check ctlrot data/configs7.txt apn_8bit.txt 1 40 REPS guided` on those 178 (one control per configuration: the first 40 functions of the list with k ≥ 3 in rotation) | `rep N (ext)` or `(nonext): line L, re-found, 0 other spaces` for each of the 178 | (same job) | `logs/check7-ctlrot-output.txt` |
| `seven_check mapa data/configs7.txt data/p7.txt` | `mapa: 3731 of A's configurations, 0 on none of ours, 3731 of our 3731 configurations hit, 0 hit more than once`, with `canon aut` equal to search A's tag on every line (search A's list mapped onto search B's, after both verdicts) | 20 s | `logs/20261101-064335-mapa7.log`, `logs/check7-mapa-output.txt` |

The near misses, compared by `python src/seven_check_vs_a.py` from the logs above (`seven-check.md` §7;
`logs/check7-vs-a-output.txt`): both searches log one at the same moment here, when a span first holds all seven
rank-4 forms with holes still uncovered (search B calls these spans cores, and both searches cover every point off the
holes before any hole), so no containment step is needed as at eight and nine. Through `mapa`, the classifications are
a bijection with the same `|Aut|` and with w = 2, 3, 4, 5, 6, 10 matching ext9, ext10, ext11, ext12, ext13, spr one to
one; search A's 114,702 `NEAR` lines (no configuration and basis repeated) and search B's cores agree in number on every
(configuration, dimension) pair, in `asc` and in `desc`: `(configuration, dim) pairs 3063, unequal 0`, the same 156
configurations with dimension-5 near misses, 2,785 with dimension 6, 121 with dimension 7, and the dimension-4 one on
search A's `p7-2222-a168-w4`, search B's 1527.

The hole sets against Langevin and Hou's count of partial spread functions in eight variables: download `fix-7.txt`
from https://langevin.univ-tln.fr/project/spread/psp.html; both `python src/seven_check_mass.py fix-7.txt` (search B's
hole sets, from `data/configs7.txt`) and `python src/holes_tally.py logs/holes7.txt data/p7.txt fix-7.txt 3726 0`
(search A's, from its `holes` run) give 3,726 classes with the same number for each of Langevin and Hou's 27 stabiliser
orders, and `17833337132662061531136` point sets, their total (`logs/rerun/seven/s7-tally.log`).

### Six (`i = 11`)

Kurz's table gives 341 hole sets for six solids, none maximal (every 6-partial spread extends), so both searches classify
the 6-partial spreads again. Times are running time between the runner's pauses, on a box often shared with other work;
a program's own clock counts the pauses too. The two searches finished in the other order from the cases above: search
B's verdict came first (2026-11-04 09:10) and search A's final run after it, since at six the hole stage of `cover`
became the whole cost.

Search A, `pspread3.cpp` (the same program as for eight and seven) then `cover.cpp`. A first run with the `cover4` build
(split at dimension 5 from the start) finished two of its eight chunks in 27 hours and was cancelled 49 hours into the
third; Rowen then rewrote the hole stage (Build, above) and ran all eight chunks again. That run, on `cover5t`, is the
result, and the `cover4` run stays as a check on the rewrite.

**How independent the rewrite is.** Six is the first case where the two searches aren't independent in their main idea, and the result should be weighed
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

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread3 enum 6 p6.txt --ckpt p6.ckpt --wcheck 1000` | `w=4: reduced 993, weighted 806400, orbits 1, orbit sum 806400 MATCH`; `w=5: reduced 142071, weighted 155232000, orbits 68, ... MATCH`; `w=6: reduced 97245, weighted 104630400, orbits 48, ... MATCH`; `w=7: reduced 185718, weighted 188596800, orbits 92, ... MATCH`; `w=11: reduced 162593, weighted 156917824, orbits 132, ... MATCH` (no w below 4: every 6-partial spread extends to at least ten solids); `orbits 341`; `w spot checks: 551, mismatches 0`; writes `p6.txt` (= `data/p6.txt`, tags `p6-<n>-a<\|Aut P\|>-w<w>`) | 20 s | `logs/20261101-065050-p3-k6.log`, `logs/pspread6-output.txt` |
| `pspread3 selftest 6 data/p6.txt` | `selftest: 341 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 24 s | `logs/rerun/six/p3-selftest6.log` (run for this folder) |
| `pspread3 holes 6 data/p6.txt holes6.txt`, then `python src/holes_tally.py logs/holes6.txt data/p6.txt fix-6.txt 341 0` (Langevin and Hou's `fix-6.txt`, downloaded from the page below) | `341 configurations, keys distinct`; `holes: 341 done this run, 0 with a cover off the list`; then `configurations 341, hole-set classes 341 (Kurz 341)`, `maximal (hole set contains no solid, w = 0): 0`, `configurations per hole set: {1: 341}`, `\|Aut H\| distribution vs Langevin-Hou: identical (17 orders)`, `point sets: 1339989812392369324032 vs Langevin-Hou 1339989812392369324032: MATCH`, `problems: 0` | under a second; the tally a second | `logs/20261101-065145-holes6.log`, `logs/holes6-output.txt`, `logs/holes6.txt`; the tally `logs/rerun/six/s6-tally.log` |
| `cover_chunks_split_exe.cmd DIR 5 cover4` (the first run, superseded), with `DIR\c0.txt` … `c7.txt` the lines of `data/p6.txt` in order, 43 to a chunk (the last 40), each run as `cover4 DIR\cK.txt -1 5 > DIR\cK.out` (the job used `cover_chunks_split.cmd`, the same script with `cover4` fixed) | `c0` `done: covers 43, tasks 1204, nodes 18427765, solutions 0, ...` (nodes by dimension and `NEAR` lines as `cover5t`'s `c0`, below), `c1` `nodes 17691942` likewise; `c2` cancelled at frontier item 17,689 of 72,221 after its solid stage had finished, so with all its `NEAR` lines and no `done:` line | c0 15.7 h and c1 11.5 h by their clocks; c2 49 h on the clock when cancelled, shared with the other jobs on the box | `logs/20261101-065248-cover6.log`, `logs/six-cover4-output-c0.txt`, `-c1.txt`, `-c2-cancelled.txt`, the spans in `logs/six-near.txt` (label `A-cover4`) |
| `cover_chunks_split_exe.cmd DIR 5 cover5t` over the same chunks: **the result** | each chunk ends `done: covers 43, tasks 1204, ...` (`c7`: `covers 40, tasks 1120`), `solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`, a `nodes by span dimension` line, a `near misses` line and `hole-stage candidates ..., passing the filter (sent to push) ...` (here the candidates are the tree's nodes: over the chunks 27,883,874,524 of them, and 29,898,477 forms sent to push) | 3.2 h on a mostly open box, 11,531 s by the chunks' clocks (c0 421 s, c1 325 s, c2 2,056 s, c3 401 s, c4 406 s, c5 492 s, c6 1,200 s, c7 6,230 s) | `logs/20261104-120910-cover6t.log`, `logs/six-output-c0.txt` … `-c7.txt` (without their `NEAR` lines: Left out, below) |
| `python src/cover_tally.py runs/i11/cover5/c?.out` (the eight outputs whole, `NEAR` lines included: 790 MB, not in this folder; on `logs/six-output-c*.txt` it gives the same but for the near misses) | `8/8 outputs done, 11531 s`; `covers: 341`, `tasks: 9,548`, `nodes: 169,008,542`, `solutions: 0`, `full (dim < 8): 0`, `dim 8 non-APN: 0`, `bad seeds: 0`; `nodes by span dimension 1..8: 9,548 / 267,344 / 5,307,371 / 45,933,458 / 80,431,617 / 36,891,380 / 167,824 / 0`; `near misses 7941871 by dim {4: 7, 5: 113279, 6: 7828585}` on all 341 configurations, none reaching dim 7 or 8; `lowest near-miss dim 4` on `p6-101`, `p6-114`, `p6-304` (2) and `p6-310` (3); the most, 115,660, on `p6-340-a2160-w11` | seconds | `logs/six-tally-output.txt` |

The rewritten hole stage was checked against the old one three ways (`logs/rerun/six/compare.log` repeats all three from
this folder's source, and adds positive controls through it: Re-run, below): whole runs on eleven's and ten's configurations (`cover5t` on `data/p11.txt` and `data/p10.txt`,
nodes at every dimension and `NEAR` lines equal to `logs/eleven-output.txt` and `logs/ten-output.txt`:
`logs/20261104-114824-cover5t-reg67.log`, `logs/cover5t-regress-p11-output.txt`, `-p10-output.txt`, the filter sending
1 and 173 candidates to `push`); forty frontier items of `p6-101-a24-w11` under the three modes of the hole stage (9,347
nodes, by dimension `0 0 0 0 40 9277 30 0` in all three, and 44,545 `NEAR` lines with one digest; candidates 4,128,555,200
with `--old`, all sent to `push`; 19,602,070,397 Gray steps and 6,487,289 tree nodes, each sending exactly 9,601 to
`push`; times 3,334 s, 92 s and 8.3 s by the programs' clocks, `--old` on a held box, so about 1,400 s of running:
`logs/cover5-prof-old40-output.txt`, `-gray40-`, `-tree40-`); and six's own chunks, where `c0` and `c1` of the `cover4`
run and of the `cover5t` run are identical in every count with the same `NEAR` lines as sets, and the cancelled `c2`'s
1,164,968 `NEAR` lines (its solid stage had finished) are `cover5t`'s (the digests in `logs/six-near.txt`, per
configuration).

On the 255 configurations outside `c0` and `c1`, a sample run under both hole stages: `src/cover_sample.cmd`, on `cover5s`
(`src/cover.cpp` with `src/cover5s.patch`, which adds `--fsample N SEED`: N frontier items per cover by a fixed seed, and an
`ITEM` line for each with its nodes by dimension, solutions, full spans and non-APN spaces). Eight dimension-6 items on
every configuration, one dimension-5 item on each of the 210 that have any. `python src/cover_sample_cmp.py
logs/cover-sample s6 s5` gives `s6: tree 2040 items, old 2040 items, on 255 configurations; lines only in tree 0, only in
old 0`, `s5: tree 210 items, old 210 items, on 210 configurations; lines only in tree 0, only in old 0`, `identical`
(`logs/cover-sample/`, `logs/20261104-175003-cover-sample.log`). `src/cover_sample_low.cmd` tried the same on the seven
dimension-4 items, all of them, which those splits can't draw. Its tree half finished (`logs/cover-sample/tree-s4.out`:
14,691,503 nodes, 11,866,974,280 hole-stage candidates, 26,597 s, 0 solutions). Its filter-free half was cancelled after
2.5 h, with nothing to compare (`logs/20261104-220921-cover-sample-low.log`): by the ratios of `s5` and `s6` it had about
7 x 10^12 candidates to run at 5 to 6 µs a thread, some two months on the job's seven threads (one per item), or two
weeks at the 1.4 µs of the `d5` controls under `--old`. The
verdict on those seven doesn't wait on it: they sit on `p6-101`, `p6-114`, `p6-304` and `p6-310`, search B's configurations
31, 24, 103 and 130 (`logs/check6-mapa-output.txt`), and the point branching searched those whole, with no filter.

Known-function controls through the tree, separate from the re-run's: `src/cover_ctl.py gen` cuts a known function down to
a D-dimensional subspace holding all its rank-4 forms, and `cover_ctl.py check` looks for the whole space among the
`SOLUTION` lines (`logs/cover-ctl/`). Each function has four rank-4 forms and runs on its own four solids (195 holes), so these
aren't six's configurations; they enter the hole stage by the real path, `seed` then `dfs` into `holes_tree`, nothing relaxed.
`d6-tree.out`: 20 of 20 and no others. `d5-tree.out` and `d5-old.out` (`--old`, job cover5t-ctl-d5old): 8 of 8 and 2 other
solutions, both APN, on `ctl5-104` and `ctl5-177`, with nodes `0 0 0 0 8 7480 1450 10` and the same `SOLUTION` lines under both
stages (candidates 4,109,334,208 under `--old`, all sent to `push`; 5,716,174 under the tree, 9,458 sent), in 742 s and 3 s. From
4 dimensions (job cover5t-ctl-d4, split at 5, replacing an unsplit job that was cancelled at the runner's restart): `d4-split.out`, 8 of 8 and 217 other
solutions, all APN and with the same four-form rank profile, 225 distinct spaces in all (48 each through `ctl4-73`, `ctl4-104` and `ctl4-169`, 20 each through
`ctl4-53`, `ctl4-65`, `ctl4-177` and `ctl4-204`, the planted one alone through `ctl4-6`), with nodes `0 0 0 8 571421 303339066 1328301 225`, 223,274,124,767
candidates and 310,407,448 sent to `push`, in 56,961 s by the program's clock, which counts paused time too. There is no filter-free run from 4 dimensions.

Search B, `six_check` (`six-check.md`), its jobs run by `six_check_slices.py` (the job logs give its commands), one
process per configuration in the searches, so the outputs here keep each configuration's line. Two searches run over the
same tree: the point branching, which is seven's search with one member fewer, ran its off-D phase in full (both orders)
and its whole search on a sample, since in full it would have taken more than ten days; and `lin`, written for six with
its own primitives, ran in full in both orders:

| Command | Expected | Time | Log |
|---|---|---|---|
| `six_check count3` | `count3 (a): k-cliques of V5824 by plain enumeration: 1:5824 2:4583936 3:606183424 4:14122766336`: every clique once, no symmetry used; `count3 (b): triangles 606183424, 4-cliques 14122766336 (by edges)`; `count3: 4-cliques = seven-check's 14122766336 OK; (a) and (b) agree OK` | 30 s | `logs/20261103-224209-check6-count3.log`, `logs/check6-count3-output.txt` |
| `six_check parta C parta-cC-s0.txt maxc-cC-s0.txt 0 1` for class C = 0 … 4 | per class `...: R reduced 3-cliques, 0 maximal`, R = 316,416, 161,408, 80,276, 28,600, 1,920 (588,620 in all), every `maxc` file empty, each output ending in `END`; writes `parta-cC-s0.txt` (= `data/parta6-cC-s0.txt`) | under a minute | `logs/20261103-224209-check6-class.log`, `logs/check6-parta-c*-s0-output.txt` |
| `six_check merge orbits.txt data/configs7.txt data/parta6-c0-s0.txt ... data/parta6-c4-s0.txt` | `seven-check's configurations: 3731; sum 210*20160/\|Aut\| = 14122766336 (its count4: 14122766336)`; `reduced 3-cliques 588620, of them maximal (Q empty) 0; weighted 3-cliques 606183424`; `sum over 3-cliques of \|Q\| = 56491065344; 4 * (weighted 4-cliques, seven-check) = 56491065344  OK`; `orbits 341; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 4:1 5:68 6:48 7:92 11:132`; `orbits by type: spr 132 ext13 92 ext12 48 ext11 68 ext10 1 ...`; `... C(14,3) * 1128864 = 410906496  OK`; `per 7-orbit O: ... fails for 0 of 3731`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits6.txt`) | seconds | same job log, `logs/check6-merge-output.txt` |
| `six_check configs orbits.txt data/configs7.txt data/configs8.txt data/spreads.txt configs6.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; 341 hole sets, none without a solid (Kurz's table); `hole sets with a cover mismatch: 0; covers in all 344`; `orbits reached by no cover: 0`; outside sources `seven-check's 3731 x 7: 26117 configurations; on no listed orbit 0; ...` and `eight-check's 9372 x 28: 262416 ...; on no listed orbit 0; ...`, each reaching all 341 orbits, and `fourteen-check's spread reps x 3003: 219219 ...; on no listed orbit 0; ... spr 132/132`; `orbits reached by no deletion from seven-check's 3731: 0`; writes `configs6.txt` (= `data/configs6.txt`; spr 0..131, ext13 132..223, ext12 224..271, ext11 272..339, ext10 340) | seconds | same job log, `logs/check6-configs-output.txt` |
| `python src/six_check_mass.py fix-6.txt` | `(a) hole classes 341, containing no solid 0   [Kurz r = 6: 341, 0]  OK`; `Langevin-Hou: 341 classes over 17 fixator orders; ours: 341 over 17`; `(b) classes per \|Aut H\| against Langevin-Hou: 0 orders differ []; their qte = classes * \|GL\| / order: OK`; `sum_H \|GL\|/\|Aut H\| = 1339989812392369324032; Langevin-Hou Total = 1339989812392369324032  OK`; `(c) labelled 6-partial spreads: from covers 1340074242473170305024; from configurations 1340074242473170305024; from Part A's W3 (606183424) 1340074242473170305024; count3's plain 3-clique count 606183424 = W3  OK` | seconds | `logs/rerun/six/s6-tally.log` |
| `six_check selftest data/configs6.txt` | `selftest: ND4 28 (28), ND6 13888 (13888), GL4 20160 (20160), V5824 5824`; `lines in D over the configurations: 2860:341`; `brute force: ...: 0 disagreements`; `selftest: 341 configurations, 119216328 generated forms checked, 0 failures` | 3 min | `logs/20261103-224438-check6-selftest.log`, `logs/check6-selftest-output.txt` |
| `six_check linselftest data/configs6.txt 16` | `K_a: 0 failures`; a line per configuration sampled, `200 cores sampled, N accepted in full scans, failures so far 0` (N a running total, 3,597 from the eighth configuration on: past it, the 200 cores sampled on each had no legal extension for the full scan to accept, so the scan's acceptances were checked on the first eight; the comparisons with `L_legal_add` run on all 16); `linselftest: 209920397 checks (3641 legal), 0 failures` | 36 min | `logs/20261103-230201-check6-linselftest.log`, `logs/check6-linselftest-output.txt` |
| `six_check search data/configs6.txt p1asc.prog asc I I phase1 cores p1asc.cores` for I = 0 … 340, and the same with `desc`: the part off the holes alone, every configuration's cores | `python src/six_check_tally.py prog data/configs6.txt logs/check6-p1asc-R.txt`: `341 of 341`, nodes by dim `341 9548 267344 5305886 45791519 80157984 7828585 0 0`, cores by dim `0 0 0 0 7 113279 7828585 0 0` (7,941,871 in all, every configuration has some); `desc`: nodes `341 9548 267344 5292965 45576666 79748390 7828585 0 0`, the same cores; `six_check_tally.py cores` on the two core files: `7941871 cores on 341 configurations, 0 repeated`, `0 configurations with cores missing, 0 extra cores with a seventh solid, 0 other extra cores` | 2 h each, side by side; the tally 7 min | `logs/20261103-224914-check6-p1asc.log`, `logs/20261103-225721-check6-p1desc.log`, `logs/check6-p1asc-output.txt`, `logs/check6-p1desc-output.txt`, `logs/check6-p1asc-R.txt`, `logs/check6-p1desc-R.txt`, `logs/20261104-005608-check6-tally-p1.log`, `logs/check6-tally-p1-output.txt`, `logs/check6-tally-p1desc-prog-output.txt` |
| `six_check lin data/configs6.txt lin.prog I I cores lin.cores` for I = 0 … 340: **`lin`, asc** | `L: 341 of 341 configurations`; nodes by dim `341 9548 267344 5308580 45963730 81635576 40014483 340912 0`, of them with a seventh solid (x4) `0 0 0 2694 172211 1318748 3134339 173214 0`; candidates tried `9548 267344 7485632 148640240 1599459715 1176846915427 18588809554755 425655971658 0` (20.19 × 10¹²); cores `8` / `125504` / `8689010` at dims 4 / 5 / 6, x4 `1` / `12225` / `860425`; `solutions 0`; no node of dimension 8. `six_check_tally.py cores` against `p1asc.cores`: `8814522 cores on 341 configurations, 0 repeated`, `0 configurations with cores missing, 872651 extra cores with a seventh solid, 0 other extra cores` | 4.1 h; the tally 44 min | `logs/20261103-225938-check6-lin.log`, `logs/check6-lin-output.txt`, `logs/check6-lin-prog.txt`, `logs/20261104-030844-check6-tally-lin.log`, `logs/check6-tally-lin-output.txt` |
| `six_check lin data/configs6.txt lindesc.prog I I desc cores lindesc.cores` for I = 0 … 340: **`lin`, desc** | `L: 341 of 341 configurations`; nodes by dim `341 9548 267344 5296100 45777356 81273352 40362043 352939 0`, x4 `0 0 0 3135 200690 1361065 3193017 179843 0`; 20.19 × 10¹² candidates; cores, x4 cores and `lin`'s fingerprint equal to `asc`'s on all 341; `solutions 0`; no node of dimension 8 | 3.2 h | `logs/20261104-030955-check6-lindesc.log`, `logs/check6-lindesc-output.txt`, `logs/check6-lindesc-prog.txt` |
| `six_check search data/configs6.txt asc.prog asc I I` for the 16 of `logs/check6-asc-sample.txt`, and on configuration 340 as a timing run (`logs/check6-asc340-*`): the point branching's whole search; `python src/six_check_tally.py nodes logs/check6-asc-R.txt logs/check6-lin-prog.txt` | `nodes by dimension, point branching (R) vs lin less x4 (L): 16 configurations in both, 0 differ`; `prog`: `16 of 341 configurations`, nodes by dim `16 448 12544 249619 2185056 4174853 18023357 108524 0`, cores `7` / `23020` / `517801` at dims 4 / 5 / 6, 636,498,734,410 tries at points of D, `bentext 0 full 0`. Configuration 340: nodes `1 28 784 15568 131388 230894 75824 23 0`, equal to `lin`'s less x4 (`1 configurations in both, 0 differ`), `triesD 10052894181`. | 30.6 h (829 to 33,470 s a configuration, by the program's clock; 340: 2,586 s, against 112 s for `lin`) | `logs/20261103-234500-check6-asc.log`, `logs/check6-asc-output.txt`, `logs/check6-asc-R.txt`, `logs/20261103-225232-check6-asc340.log`, `logs/check6-asc340-output.txt`, `logs/check6-asc340-R.txt` |
| `six_check plantfind LUT KMIN 100000000 CODIM plants.txt` over the known lists (Data, below), then `six_check ctlplant data/configs6.txt LUT logs/check6-plants-NAME.txt ctlplant-RUN.txt I 1 [lin]` per planting: the main positive controls, a known function's rank-4 radicals plus solids that the off-D phase can't complete, and the configuration's whole search unguided | codimension 2: 5 plantings, all from `apn_8bit.txt` (`d0 = 6`), `FOUND ...: 1 spaces (1 APN)` for each by `lin` and by the point branching (`c2bl-lin`, `c2bl-pb`; 33.5 × 10⁹ tries of a line form in D in the first); codimension 1: 40 of `apn_8bit.txt`'s 464 and all 44 of `new_apns.txt`'s (from its 7,199 with k ≥ 2), by `lin` (`c1bl-lin`, `c1new-lin`); `blp22_8bit.txt`, `yuperrin_qam.txt` and `yuwangli_qam.txt` have no planting; 94 controls, 94 re-found, 0 other spaces | 2.2 h for the two `apn_8bit.txt` `lin` batches (with `lfull`'s first try, below), 29 min for `new_apns.txt`, 6.8 h by their clocks for the five by point branching; `plantfind` over all lists 80 min | `logs/20261103-225721-check6-plantfind.log`, `logs/20261103-233727-check6-plantfind2.log`, `logs/20261104-000612-check6-plantfind3.log`, `logs/check6-plantfind-*-output.txt`, `logs/check6-plants-*.txt`; `logs/20261103-234511-check6-ctlplant-lin.log`, `logs/20261104-012507-check6-ctlplant-new.log`, `logs/20261104-000959-check6-ctl-pb.log`; `logs/check6-ctlplant-*-output.txt`, `logs/check6-ctlplant-*.txt` |
| `six_check ctlbatch data/configs6.txt apn_8bit.txt 2 MAXN CH ctl-NAME.txt OPT` for (NAME, CH, OPT) = (`lfull-bl22`, 1, `lin`) with MAXN 20, the unguided full controls; (`lind-bl22`, 3, `nond+lin`), (`lind2-bl22`, 3, `nond2+lin`), (`lind3-bl22`, 3, `nond3+lin`), (`ind2-bl22`, 3, `nond2`), (`ind3-bl22`, 3, `nond3`) with MAXN 1000000, the D phase from a subspace of the function of codimension 1, 2 or 3 | `lfull`: `ext controls 20 refound 20`, `nonext controls 20 refound 20`, 2.16 × 10¹² candidates over every other core; the in-D batches each `ext controls 1440 refound 1440`, `nonext controls 1383 refound 1383` (`reps 129` and `209`; 57 `k = 4` functions skipped on the nonext path, `no member 4-tuple of the class`); `notfound 0 fail 0 extra 0 nonapn 0` in every line. `ctl-lfull-bl22-nofix` is the first `lfull` batch, 0 of 40, which caught a bug in `lin` (a core that is already a whole space wasn't recorded, which arises only in controls: `six-check.md` §6); the job holding it failed with exit 3 for that reason | 23 min for the five in-D batches; `lfull` 46 min | `logs/20261103-233738-check6-ctl-ind.log`, `logs/20261104-020445-check6-lfull.log`, `logs/check6-ctl-*-output.txt`, `logs/check6-ctl-*.txt` |
| `six_check ctlrot data/configs6.txt apn_8bit.txt 1 40 REPS lin` in chunks of 31 over every configuration (one unguided control per configuration, the first 40 functions of the list with k ≥ 3 in rotation) | `rep N (ext)` or `(nonext): line L, re-found, 0 other spaces; ...` for each of the 341 | 11 h on the shared E-cores | `logs/20261104-035403-check6-ctlrot-lin.log`, `logs/check6-ctlrot-lin-output.txt` |
| `six_check ctlbatch data/configs6.txt apn_8bit.txt 2 5 1 ctl-full-bl22.txt` (the point branching's unguided full controls) | `ext controls 5 refound 5 notfound 0 fail 0 extra 0 nonapn 0 reps 4`, `nonext controls 5 refound 5 notfound 0 fail 0 extra 0 nonapn 0 reps 5`, each over the configuration's whole search, tries at points of D `1678434240` / `42200996663` / `9323128` and `2320948672` / `44106529070` / `7405291` at dims 5 / 6 / 7 | 10 h of the 16.7 h job (which ran the five planted ones first), mostly on the shared E-cores | `logs/20261104-000959-check6-ctl-pb.log`, `logs/check6-ctl-full-bl22-output.txt`, `logs/check6-ctl-full-bl22.txt` |
| `six_check mapa data/configs6.txt data/p6.txt` | `mapa: 341 of A's configurations, 0 on none of ours, 341 of our 341 configurations hit, 0 hit more than once`, with `canon aut` equal to search A's tag on every line (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/20261104-152337-check6-mapa.log`, `logs/check6-mapa-output.txt` |

The two searches compared, by `python src/six_check_vs_a.py` from the logs above (`six-check.md` §7;
`logs/check6-vs-a-output.txt` is the original's output, which read the spans themselves): through `mapa`, the
classifications are a bijection with the same `|Aut|` and with w = 4, 5, 6, 7, 11 matching ext10, ext11, ext12, ext13,
spr one to one. Search A logs a `NEAR` line when a span first holds all six rank-4 forms with holes uncovered, the event
search B's cores record, so its 7,941,871 near misses and search B's cores agree in number on every (configuration,
dimension) pair, in `p1asc` and in `p1desc`: `(configuration, dim) pairs 627, unequal 0`, the same 4 configurations with
dimension-4 near misses, 282 with dimension 5, all 341 with dimension 6. Where a configuration's representatives share a
frame, 225 of the 341, search A's spans and search B's are the same sets (5,520,490 spans); on the other 116 they're
disjoint, being in different frames. `logs/six-near.txt` has these per configuration, for each run a count by dimension
and a digest of the set of spans (`src/six_near.py`), so that the set comparison can be repeated here without the spans.
Their node totals agree through dimension 2 and then differ. Per configuration (`python src/cover_percfg_cmp.py logs`, its output
`logs/cover-percfg-cmp-output.txt`, from `logs/cover-percfg-output.txt`: `src/cover_percfg.py`, the published `cover5t` run one
configuration at a time, job cover5t-percfg, its totals the published ones): on the 225 shared frames search A's counts equal
`lin`'s less x4 at dimensions 1-3 and 6-8 on every configuration, and differ at 4 and 5 on 175 and 91, which is the first phase
and not the holes, since the hole stage has only the near misses there and, on the four configurations with 4-dimensional near misses, their
first children. On the other 116 they differ at dimensions 3 to 7. Against the point branching, on the 13 of its 17 with shared
frames, they are equal at 6-8, and on those four from 5 up.

The hole sets against Langevin and Hou's count of partial spread functions in eight variables: download `fix-6.txt`
from https://langevin.univ-tln.fr/project/spread/psp.html; both `python src/six_check_mass.py fix-6.txt` (search B's
hole sets, from `data/configs6.txt`) and `python src/holes_tally.py logs/holes6.txt data/p6.txt fix-6.txt 341 0`
(search A's, from its `holes` run) give 341 classes with the same number for each of Langevin and Hou's 17 stabiliser
orders, and `1339989812392369324032` point sets, their total (`logs/rerun/six/s6-tally.log`).

## Re-run for this folder

On 2026-10-25, from a fresh MSVC build of `src/` with the Build line above, on an idle machine (output in
`logs/rerun/`):

- `sixteen_check selftest`: every line `ok`.
- `fourteen_check selftest data/spreads.txt`, `search data/spreads.txt asc14.prog asc`, `holes 11`, `spreads`: the
  logged results, the search's node and try counts identical at every dimension (5 s); the regenerated `spreads.txt`,
  `spreads.txt.all` and `holes11.txt` are identical to the originals (line endings aside).
- `python src/fourteen_check_levels.py data/spreads.txt 0 4`: `[1, 28, 784, 10234, 15917]`.
- `cover data/i3-reps.txt`: `nodes 2284642, solutions 0` (3 s); `src/spreads_to_cover.py` reproduces `data/i3-reps.txt`.
- `thirteen_check selftest data/configs.txt` and `search data/configs.txt asc13.prog asc`: the logged results,
  identical counts and the same 4 near-miss nodes (29 s).
- `thirteen_check holes`: the logged counts. This build also splits the non-extendable hole invariants by plane counts
  (three lines, 25,081,056 + 33,153,120 + 4,324,320 = 62,558,496), which is why its `holes10.txt` is larger.
- `thirteen_check configs`: log identical apart from the output path, `configs.txt` identical to `data/configs.txt`
  (2 min). `thirteen_check covers`: `covers.txt` identical to search B's.
- `python src/thirteen_check_levels.py data/configs.txt 0 3` / `170 3`: `[1, 28, 784, 10573]` / `[1, 28, 784, 10288]`.
- `python src/i4_holes.py data/kurz13.txt`: as in the table; `python src/cover_prep.py max`: `max.txt` identical to the
  one cover.cpp's `i4-max` run read.
- `python src/cover_prep.py ctl ... apn_8bit.txt 12631 1 KEEP` for KEEP 4, 3, 2: identical to `data/ctl-12631-k*.txt`;
  `cover` on each: the known space re-found as the only solution (4, 5 and 13 nodes).

Not re-run: the hour-long searches, the control batches, `cover_prep.py ext` (22 min) and the 15-hour `sixteen i3`.

For twelve, later the same day:

- `python src/holes_iso12.py data/p12.txt`: the 5 hole classes with `|Aut H|` and cover sums as in the table.
- `python src/twelve_check_levels.py data/configs12.txt 0 4` / `280 4` / `364 4`: the logged counts.
- On 2026-10-26, from a fresh build of `src/pspread.cpp`, `cover.cpp` and `twelve_check.cpp` (logs in `logs/rerun/`,
  times in `times-twelve.txt`): `cover data/p12.txt` and `cover data/i3-reps.txt` (`cover-i3-reps-2.log`), the three
  `pspread selftest`s, `pspread enum 17`, and `twelve_check` `configs`, `selftest`, `search` asc / desc / rand 7,
  `search13`, `searchsp` (summary lines only), `maprowen`, `ctlreps` on representatives 56 and 279, and the `nond`
  trial `2 10 1`: every result as in the tables, node counts identical at every dimension; the regenerated
  `configs12.txt` identical to `data/configs12.txt`, the `nond` trial's output identical apart from the list's path.

Not re-run: the three `pspread enum` runs for K = 12, 13, 14 (an hour or more each), `twelve_check parta` (26 min),
`searchmax` and the control batches.

For eleven, on 2026-10-27, from a fresh build of `src/cover.cpp`, `pspread.cpp` and `eleven_check.cpp` with the Build
line (the first two identical to the binaries that ran, by disassembly), as one job beside a long CPU job of
Rowen's (logs in `logs/rerun/`, times in `times-eleven.txt`): `cover data/p11.txt` (`cover-p11.log`,
100 s), `pspread selftest 11 data/p11.txt`, `eleven_check parta 4` (2.3 min), `merge`, `configs`, `selftest`, `search`
asc / desc / rand 7 (2 min each), the four `eleven_check_levels.py` recounts, `mapa`, `searchmax data/maxc11-o798.txt`,
and `ctlreps` on representatives 93 and 452: every result as in the tables. Node, try and near-miss counts are
identical at every dimension and per configuration; the `NEAR` lines of `cover` and the `P1` spans of the three
searches are the same sets as the originals'; the regenerated `parta-c4.txt` (apart from its timing), `orbits.txt` and
`configs11.txt` are identical to `data/parta11-c4.txt`, `data/orbits11.txt` and `data/configs11.txt`, and the other
logs to the originals apart from timings and output paths.

Not re-run: `pspread enum 11` (about 2 h), `eleven_check parta` for classes 0 to 3 (up to 2 h), the cross-runs
`search12`, `search13`, `searchsp` and `searchmax` on the sample (50 min to 4.5 h), and the control batches.

For ten, on 2026-10-28, from a fresh build of `src/cover.cpp`, `pspread3.cpp` and `ten_check.cpp` with the Build
line, as two jobs on the E-cores beside other work (logs in `logs/rerun/`, times in `times-ten.txt`): `cover
data/p10.txt` (`cover-p10.log`, 5.5 min), `pspread3 selftest 10 data/p10.txt`, `ten_check parta` for class 4 and for
class 3 slice 1 (`parta 3 ... 1 2`), `merge`, `configs`, the five `ten_check_levels.py` recounts, `mapa`, the `nond`
control batch (`ctl-ind-bl22`, 41 s), `ctlreps` on representatives 131 and 611, the first 500 configurations of
`search11` and `searchmax data/maxc10-s400.txt` and the first 800 of `searchsp` (configurations 0 to 499 or 799, given
as the last two arguments), `search` asc and desc (7 and 10 min) and `selftest` (13 min); and `ten_check_tally.py
reps` and `maxall` and `near_misses10.py` on the files here (`t10-tally.log`): every result as in the tables. Node,
try and near-miss counts are identical at every dimension and per configuration (in the slices, configuration by
configuration); the `NEAR` lines of `cover` and the `P1` spans of the two searches are the same sets as the
originals'; the regenerated `parta-c4-s0.txt`, `parta-c3-s1.txt`, `orbits.txt` and `configs10.txt` are identical to
`data/parta10-c4-s0.txt`, `data/parta10-c3-s1.txt` (apart from the timing on their `CLASS` lines),
`data/orbits10.txt` and `data/configs10.txt`, and the other logs to the originals apart from timings and
output paths. The two regenerated lists of maximal cliques hold all 57 lines of `data/maxc10-s400.txt` that fall in
their stretch of the concatenation (lines 1, 401, 801, … counted over the 17 files in class and slice order).

Not re-run: `pspread3 enum 10` (about 5 h), `ten_check parta` for the other slices (44 min on an open box for all 17),
the rest of `search11`, `search13`, `searchsp` and `searchmax` (47 min to 6.5 h), `search` rand 7, and the control
batches other than `nond`.

For nine, on 2026-10-29, from a fresh build of `src/cover.cpp`, `pspread3.cpp` and `nine_check.cpp` with the Build
line, as two jobs, mostly on the E-cores beside other work (logs in `logs/rerun/`, times in `times-nine.txt`): `cover
data/p9.txt` (`cover-p9.log`, 13.5 min), `pspread3 selftest 9 data/p9.txt`, `nine_check parta` for class 4 (both
slices) and class 3 slice 3, `merge` (on the 34 `data/parta9-*.txt`), `configs`, the
seven `nine_check_levels.py` recounts, `mapa`, the `nond` control batch (`ctl-ind-bl22`), `ctlreps` on
representatives 123 and 5655 (the first and last of the 851), the first 500 configurations of `search10`, `search13`
and `searchmax data/maxc9-s1000.txt` and the first 800 of `searchsp`, the first 1,000 of `search` desc and rand 7,
`search` asc in full (32 min) and `selftest` (36 min); `orbitcov` on `data/maxc9-s1000.txt` (new: `16493 cliques, 0 on no
listed orbit, 0 on a non-max orbit; 1335 orbits hit of 1348 max orbits`) and on the 6 cliques of
`logs/check9-maxs1000-near.txt` (all on orbit 4776); and the tallies (`nine_check_tally.py prog` on the three main
runs here and on the cross-runs' `.prog` files, `ten_check_tally.py reps` and `maxall`, `nine_check_compare.py`;
`t9-tally.log`): every result as in the tables. Node, try and near-miss counts are identical at every dimension and
per configuration (the searches' `R` lines, configuration by configuration, against the originals': the main runs'
in `logs/`, the cross-runs' in their `.prog` files, which aren't in this folder); the `NEAR` lines
of `cover` and the `P1` spans of `asc` are the same sets as the originals'; the regenerated `parta-c4-s0.txt`,
`parta-c4-s1.txt` and `parta-c3-s3.txt` are identical to `data/parta9-c4-s0.txt`, `data/parta9-c4-s1.txt` and
`data/parta9-c3-s3.txt` (apart from the timing on their `CLASS` lines), `orbits.txt` and `configs9.txt` to
`data/orbits9.txt` and `data/configs9.txt`, and the other logs to the originals apart from timings and output paths.

Not re-run: `pspread3 enum 9` (about 3.8 h of running time), `nine_check parta` for the other 31 slices (2.2 h on the
E-cores for all 34), the rest of `search10`, `search13`, `searchsp` and `searchmax` (50 min to about 17 h), `orbitcov`
over all of Part A's maximal cliques, and the control batches other than `nond`.

For eight, on 2026-10-30, from a fresh build of `src/cover.cpp`, `pspread3.cpp` and `eight_check.cpp` with the Build
line, as one job on the E-cores beside other work (logs in `logs/rerun/`): `pspread3 selftest 8 data/p8.txt` (as in
the table), `pspread3 holes 8 data/p8.txt` (its `holes8.txt` the same set of lines as `logs/holes8.txt`; the program
writes them in the order its threads finish), `pspread3 enum 17` (a check that the `holes` edit left the enumeration
alone: `p17.txt` identical to `logs/rerun/p17.txt`, and the log the same as that of the binary that ran eight's
enumeration, `pspread3-p17-runbinary.log`, apart from timings), `cover` on the first chunk of `data/p8.txt` (lines
1 to 1,172: `done: covers 1172, tasks 32816, nodes 97450561, solutions 0`, node counts identical at every dimension, and
the same 305 `NEAR` lines as `logs/eight-output-c0.txt`, as a set); and the tallies (`nine_check_tally.py prog` on the
three main runs, `ten_check_tally.py reps` and `maxall`, `eight_check_compare.py orders` and `numpy`,
`eight_check_vs_a.py`, `eight_check_mass.py`, `holes_tally.py` for eight and nine, `cover_tally.py`;
`e8-tally.log`) and `eight_check_percfg.py` (`e8-percfg.log`): every result as in the tables, the numpy comparison and
the masses identical to the originals' output. The same job's remaining steps (logs, outputs, the job script `rerun.py` and the comparison `compare.py` with its output `compare.log` in `logs/rerun/eight/`; compared by Rowen): `cover` on both the first and the last chunk (lines 8,205 to 9,372: node counts identical at every dimension, the same 4,213 `NEAR` lines as a set); `eight_check` `parta` for class 4 (both slices) and class 3 slice 3, each the same lines as `data/` apart from the timing field, and the same maximal cliques; `merge` from the 34 files in `data/`, its `orbits.txt` identical to `data/orbits8.txt` and its output line for line; `configs`, its list identical to `data/configs8.txt` and its output line for line apart from the path it names; `mapa` line for line; `orbitcov` on the sample; the `nond` batch (1,440 and 1,421 re-found, the summaries as in the original apart from a path's case); `ctlrot` on 91, 92, 9364, 9365 (the same four lines); `search` asc on 8,204 to 9,371, desc and rand 7 on 2,344 to 2,843, `searchmax` on 0 to 499 and `search9` on 0 to 99, every `R` line identical to the original's (1,168, 500, 500, 500 and 100); numpy on configuration 9191 and two planted controls, identical apart from timings.

Not re-run: `pspread3 enum 8` (87 min), `pspread3 holes 9` (25 min), `eight_check count5` (12 min of running time),
`selftest` (56 min), the other Part A slices, the rest of the three main searches and of the cross-runs, and the control
batches other than `nond`.

For seven, on 2026-11-01, from a fresh build of `src/cover.cpp`, `pspread3.cpp` and `seven_check.cpp` with the Build
line (each giving the run binaries' object files, Build above), as three jobs on the E-cores beside other work, by
`logs/rerun/seven/rerun.py` from this folder's files (logs, `times.txt`, and the comparison `compare.py` with its
output `compare.log` in `logs/rerun/seven/`; the `cover` logs there without their `NEAR` lines, and the `configs` and
`mapa` logs left out, since `compare.log` finds them equal to the originals in `logs/`), and the tallies on the files
here (`s7-tally.log`):

- `pspread3 enum 7`: its `p7.txt` byte for byte `data/p7.txt`, and every class total, `w=` line and the orbit count as
  in `logs/pspread7-output.txt`; the w spot checks, every 1,000th clique of each thread and so dependent on the thread
  count, were 15,238 here against 15,245, 0 mismatches in both (52 min). `pspread3 selftest 7` as in the table;
  `pspread3 holes 7`: its `holes7.txt` the same set of lines as `logs/holes7.txt`, its output the original's.
- `cover` (the source as it then was, `cover4`'s) on the third chunk, lines 935 to 1,401 of `data/p7.txt`: unsplit, and split
  at dimension 5, where its 56 dimension-5 near misses, in 6 configurations, become 56 frontier nodes searched in the
  second phase (`phase 1 done 836 s: frontier 56 items at dim 5`). Both give `cover3`'s original `done:` line,
  `70571705` nodes, identical at every dimension, and the same 3,841 `NEAR` lines as a set (21 and 17 min).
- `seven_check` `count4` line for line; `parta` for all 37 slices, each the original file apart from the timing field,
  every `maxc` file empty (34 min); `merge` from them, its `orbits.txt` identical to `data/orbits7.txt` and its output
  line for line; `configs`, its list identical to `data/configs7.txt` and its output line for line apart from the
  file name it writes; `mapa` byte for byte.
- `search` `asc` on configurations 0 to 10 (run by the first build) and 1,000 to 1,099 (the second), both with
  `SEVEN_NOCACHE` set, and on 3,631 to 3,730 (the seventh, with the caches); `desc` on 1,000 to 1,099 (`SEVEN_NOCACHE`);
  the sweep on 0 to 8 (the third build) and 2,000 to 2,099 (the seventh); `phase1` on 1,000 to 1,099; and the
  `search8` block at 4,728 (24 inputs): every `R` line identical to the original's (211, 100, 109, 100 and 24), except
  the sweep's count of element checks on configuration 1, whose two dimension-5 cores the sweep's cache now serves
  (920,973,542 checks then, 242,798,460 now; its spans, cores and every other field the same); and the cores of `asc`
  and `desc`, written again, the same sets as the originals' on every configuration (2,339 and 904 cores).
- The `ind-bl22` control batch (`nond`, run by the second build): its `ctl-ind-bl22.txt` identical, its output line
  for line apart from one line the later source prints first (`ctlbatch options: start from a hyperplane of T (in-D
  test)`); `ctlrot ... guided` on 99, 102, 3707 and 3729: the same four lines.
- `seven_check_numpy.py` on configuration 3 and the first planted control of `control apn_8bit.txt 1 20 7`: the same
  lines apart from the times (7.5 min each).
- `selftest data/configs7.txt` line for line (`896507276 generated forms checked, 0 failures`; 27 min).
- The tallies and comparisons on the files here (`s7-tally.log`): `holes_tally.py` and `seven_check_mass.py` against
  `fix-7.txt`, `cover_tally.py` over the nine outputs, `seven_check_tally.py` `prog` on the four runs and the cross-run,
  `orders`, `numpy` and `reps`, and `seven_check_vs_a.py` (its output the original's line for line): every result as
  in the tables.

So the current sources reproduce the records of every build that ran seven, the three overwritten ones included,
apart from that one count of the sweep's.

Not re-run: the rest of the searches (each 9 to 11 h of running time, and 10 h for search A's cover chunks), the
control batches other than `ind-bl22` (8.5 h), the numpy sample beyond configuration 3 (4.9 h) and the planted
controls beyond the first.

For six, on 2026-11-04, from a fresh build of `src/cover.cpp`, `pspread3.cpp` and `six_check.cpp` with the Build line
(each giving the object file of the run binary built from that source, Build above; the earlier builds of `six_check` are
the same source with less in it), as three jobs on an open box beside search B's last jobs, by `logs/rerun/six/rerun.py`
from this folder's files (logs, `times.txt`, and the comparison `compare.py` with its output `compare.log` in
`logs/rerun/six/`; the `cover` logs there without their `NEAR` lines, and the `configs` and `mapa` logs left out, since
`compare.log` finds them equal to the originals), and the tallies on the files here (`s6-tally.log`):

- `pspread3 enum 6`: its `p6.txt` byte for byte `data/p6.txt`, and every class total, `w=` line and the orbit count as in
  `logs/pspread6-output.txt`; the w spot checks, every 1,000th clique of each thread and so dependent on the thread count,
  were 544 here against 551, 0 mismatches in both. `pspread3 selftest 6` as in the table; `pspread3 holes 6`: its
  `holes6.txt` the same set of lines as `logs/holes6.txt`, its output the original's.
- `cover` (the current source, `cover5t`'s) on the first three chunks of `data/p6.txt`, split at dimension 5 (12, 16 and
  55 min): each gives `cover5t`'s `done:` line, nodes identical at every dimension, the same filter counts and the same
  `NEAR` lines as a set (940,270, 906,078 and 1,164,968). `c2` holds two of the configurations with a dimension-4 near
  miss. Against the old hole stage: `c0` and `c1` equal `cover4`'s in the same ways, and `c2`'s `NEAR` lines are the
  cancelled run's.
- `cover` on `data/p11.txt` and `data/p10.txt` (2 and 5 min): nodes at every dimension, `done:` lines and `NEAR` lines as a
  set equal to `logs/eleven-output.txt` and `logs/ten-output.txt`, which `cover3` wrote, and to `cover5t`'s regressions.
- `cover` on the forty frontier items of `p6-101-a24-w11` (`--only p6-101-a24-w11 --frange 0 40`), by default and with
  `--gray` (3 and 38 s; `--old`, nearly an hour by its clock, wasn't re-run): every count and the `NEAR` lines as
  set equal to `cover5t`'s run and to `cover5`'s (the Gray build), filter counts included, so `--gray` is that build's
  filter; and nodes and `NEAR` lines equal to the `--old` run's.
- The positive controls through the rewritten hole stage, new here: the `SOLUTION` and `done:` lines of `cover` on
  `data/ctl-12631-k4.txt`, `-k3` and `-k2` are the old hole stage's (they leave it little to do), and `cover` on
  `python src/cover_prep.py ctl ... apn_8bit.txt 12631 1 1` (one completing form kept: the case cancelled as too slow on
  the old hole stage, Left out below) leaves it most of the search: `done: covers 1, tasks 1, nodes 499550, solutions 1, ...`,
  `nodes by span dimension 1..8: 0 1 15 150 773 495431 3179 1`, 773 near misses, 389,475,073 tree nodes, 506,819 forms sent
  to `push`, and the one `SOLUTION` is the known function (41 min, `cover-ctl-k1.summary.log`). With none kept the search
  is too large for this folder and was cancelled.
- `six_check` `count3`, `merge` (its `orbits.txt` identical to `data/orbits6.txt`), `configs` (its list identical to
  `data/configs6.txt`), `selftest` (4 min) and `linselftest` (10 min) line for line, apart from timings and the file name
  `configs` writes; `parta` for the 5 classes, each the original file apart from the timing field, every `maxc` file empty,
  and its output the original's apart from timings and the running clique counts in its progress lines, which depend on
  the threads; `mapa` byte for byte.
- `plantfind` on `apn_8bit.txt` (seconds): codimension 1
  writes `logs/check6-plants-bl22-c1.txt` line for line, its output the original's. Codimension 2 with k ≥ 4, the command
  in the job log, gives that log's 3 plantings; the list search B's controls read, `logs/check6-plants-bl22-c2.txt`, has 5,
  with k ≥ 3, from a later run that left no output (its file is newer than the job's end). k ≥ 3 here writes those 5
  line for line (`172 functions with 3 <= k <= 6, 5 planted`), which as far as can be told is that run.
- `search phase1 cores` on configurations 0 to 40 and 300 to 340 in `asc` and 100 to 140 in `desc` (about 90 s each), one
  process per range where the originals ran one per configuration; `lin` on 0, 24, 57, 190, 248, 340 (1 to 3 min each), and
  `lin desc` on 0 and 340; `search asc` on 248, the point branching's whole search (15 min): every `R` and `L` line
  identical to the original's (82 + 41, 6 + 2, 1), and the core files the same sets as the originals' on every
  configuration (3,206,897 off-D cores; 228,535 of `lin`'s).
- The controls: `lind-bl22` (44 s), its `ctl-lind-bl22.txt` identical apart from a path and its output line for line;
  the planted `c2bl-lin` plantings 0 and 3, `FOUND` lines identical apart from the time; `ctlrot ... lin` on 0 and 340,
  the same two lines.
- The tallies and comparisons on the files here (`s6-tally.log`): `holes_tally.py` and `six_check_mass.py` against
  `fix-6.txt`, `cover_tally.py` over the eight outputs (the near misses aside, which it reads from the `NEAR` lines:
  `logs/six-tally-output.txt` has its run on the whole outputs), `six_check_tally.py` `prog` on the five runs and `nodes`,
  `six_near.py compare` (`A-cover4` equal to `A` on its 129 configurations, `p1desc` to `p1asc` on 341, the three
  profile runs to each other) and `six_check_vs_a.py` (the original's output line for line, apart from the two lines
  that say how many files it read and how it compared the sets): every result as in the tables.

So the current sources reproduce the records of the builds that ran six wherever they were re-run; `lfull`'s build alone
had no run here.

Not re-run: search A's other five chunks (about 2.4 h on an open box, 1.7 h of it `c7`), `--old` on the forty items, the
rest of the off-D phases (2 h each), of `lin` (3 to 4 h each) and of the point branching's sample (up to 4.5 h a
configuration), `six_check_tally.py cores` over the whole core files (7 and 44 min; they aren't in this folder), and the
control batches other than `lind-bl22` and those few of the planted and rotated ones (from 3 min to 10 h a batch).

## Data

`data/`:
- `kurz13.txt`: Kurz's three 8 × 60 generator matrices (arXiv 2112.11763, §9, after Table 9.3), cut from the paper's
  text, not retyped. Search B's `kurz.txt` was the same file.
- `spreads.txt`: search B's 73 `i = 3` configurations (output of `fourteen_check spreads`).
- `configs.txt`: search B's 190 `i = 4` configurations (output of `thirteen_check configs`).
- `i3-reps.txt`: `spreads.txt` as cover.cpp input (`python src/spreads_to_cover.py data/spreads.txt`).
- `ctl-12631-k4.txt`, `-k3`, `-k2`: cover.cpp's positive controls (`python src/cover_prep.py ctl ... 12631 1 KEEP`).
- `p12.txt`: search A's 376 `i = 5` configurations (output of `pspread enum 12`; cover.cpp input, tags
  `p12-<n>-a<|Aut P|>-w<w>`). `p13.txt`, `p14.txt`: `pspread enum 13` / `14` (190 and 73 configurations, the checks
  against the cases above; inputs to `pspread selftest`).
- `orbits12.txt`: search B's 376 orbits with their invariants (output of `twelve_check parta`; input to `configs`).
- `configs12.txt`: search B's 376 `i = 5` configurations (output of `twelve_check configs`).
- `p11.txt`: search A's 817 `i = 6` configurations (output of `pspread enum 11`; cover.cpp input, tags
  `p11-<n>-a<|Aut P|>-w<w>`).
- `parta11-c0.txt` … `parta11-c4.txt`: search B's Part A for `i = 6`, one file per minimal class (outputs of
  `eleven_check parta`; inputs to `merge`).
- `orbits11.txt`: search B's 817 orbits with their invariants (output of `eleven_check merge`; input to `configs`).
- `configs11.txt`: search B's 817 `i = 6` configurations (output of `eleven_check configs`).
- `maxc11-o798.txt`: the 85 maximal 8-cliques on orbit 798 (output of `eleven_check orbitcov`; input to the last
  `searchmax`).
- `p10.txt`: search A's 1,777 `i = 7` configurations (output of `pspread3 enum 10`; cover.cpp input, tags
  `p10-<n>-a<|Aut P|>-w<w>`).
- `parta10-c0-s0.txt` … `parta10-c4-s0.txt`: search B's Part A for `i = 7`, one file per slice of a minimal class (17;
  outputs of `ten_check parta`; inputs to `merge`).
- `orbits10.txt`: search B's 1,777 orbits with their invariants (output of `ten_check merge`; input to `configs`).
- `configs10.txt`: search B's 1,777 `i = 7` configurations (output of `ten_check configs`).
- `maxc10-s400.txt`: every 400th of Part A's 6,042,938 maximal 7-cliques (input to `searchmax` and `orbitcov`).
  Included because it can't be cut again exactly: `parta` writes its maximal cliques in the order its threads finish.
- `p9.txt`: search A's 5,663 `i = 8` configurations (output of `pspread3 enum 9`; cover.cpp input, tags
  `p9-<n>-a<|Aut P|>-w<w>`).
- `parta9-c0-s0.txt` … `parta9-c4-s1.txt`: search B's Part A for `i = 8`, one file per slice of a minimal class (34;
  outputs of `nine_check parta`; inputs to `merge`).
- `orbits9.txt`: search B's 5,663 orbits with their invariants (output of `nine_check merge`; input to `configs`).
- `configs9.txt`: search B's 5,663 `i = 8` configurations (output of `nine_check configs`).
- `maxc9-s1000.txt`: every 1000th of Part A's 16,492,152 maximal 6-cliques (input to `searchmax`), included for the
  same reason as `maxc10-s400.txt`.
- `maxc9-missed.txt`: every raw maximal 6-clique of Part A on the 13 max orbits the sample misses (14,646, each
  labelled with its orbit; input to `searchmax`).
- `p8.txt`: search A's 9,372 `i = 9` configurations (output of `pspread3 enum 8`; cover.cpp input, tags
  `p8-<n>-a<|Aut P|>-w<w>`; input to `pspread3 holes`).
- `parta8-c0-s0.txt` … `parta8-c4-s1.txt`: search B's Part A for `i = 9`, one file per slice of a minimal class (34;
  outputs of `eight_check parta`; inputs to `merge`).
- `orbits8.txt`: search B's 9,372 orbits with their invariants (output of `eight_check merge`; input to `configs`).
- `configs8.txt`: search B's 9,372 `i = 9` configurations (output of `eight_check configs`).
- `maxc8-s100.txt`: every 100th of Part A's 1,410,915 maximal 5-cliques (input to `searchmax`), included for the
  same reason as `maxc10-s400.txt`.
- `p7.txt`: search A's 3,731 `i = 10` configurations (output of `pspread3 enum 7`; cover.cpp input, tags
  `p7-<n>-a<|Aut P|>-w<w>`; input to `pspread3 holes` and `seven_check mapa`).
- `orbits7.txt`: search B's 3,731 orbits with their invariants (output of `seven_check merge`; input to `configs`).
- `configs7.txt`: search B's 3,731 `i = 10` configurations (output of `seven_check configs`).
- `p6.txt`: search A's 341 `i = 11` configurations (output of `pspread3 enum 6`; cover.cpp input, tags
  `p6-<n>-a<|Aut P|>-w<w>`; input to `pspread3 holes` and `six_check mapa`).
- `parta6-c0-s0.txt` … `parta6-c4-s0.txt`: search B's Part A for `i = 11`, one file per minimal class (outputs of
  `six_check parta`; inputs to `merge`).
- `orbits6.txt`: search B's 341 orbits with their invariants (output of `six_check merge`; input to `configs`).
- `configs6.txt`: search B's 341 `i = 11` configurations (output of `six_check configs`).

Not included, all regenerated by the commands above: `spreads.txt.all` (280 KB, 40 s), `holes11.txt`, `holes10.txt`
(16 to 22 MB depending on the build, 8 s), search B's `covers.txt` (1.9 MB, seconds), cover.cpp's `max.txt` (1.9 MB)
and `ext.txt` (4.5 MB, 22 min), `twelve_check parta`'s `maxc.txt` (5.6 MB, 26 min), `eleven_check parta`'s
`maxc-c0.txt` … `maxc-c3.txt` (36 MB together, `maxc-c4.txt` is empty; 3.9 h of running time on the E-cores for the
five classes, 2 h of it class 0) and the 1-in-50 sample `maxc-s50.txt` taken from them (710 KB; how in the table),
`ten_check parta`'s `maxc-c*-s*.txt` (about 290 MB together; 44 min on an open box for all 17 slices; deleted after
the `orbitcov` tally, so `maxc10-s400.txt` is in `data/`), `nine_check parta`'s `maxc-c*-s*.txt` (about 645 MB; 2.2 h
on the E-cores for all 34 slices, one after another, so `maxc9-s1000.txt` is in `data/`), `eight_check parta`'s
`maxc-c*-s*.txt` (about 51 MB; 70 min of running time for all 34 slices, one after another, so `maxc8-s100.txt` is in
`data/`), and the `.prog` resume files (search B's near-miss spans, its `P1` lines, are in `logs/check11-*-p1.txt` and
`logs/check10-*-p1.txt`; for nine the three main runs' `.prog` files are in `logs/` whole, as
`logs/check9-{asc,desc,rand7}-prog.txt`, since their `R` lines are the per-configuration record those runs left, and
for eight likewise as `logs/check8-{asc,desc,rand7}-prog.txt`, each order's 8 ranges concatenated; for seven their
`R` lines are in `logs/check7-{asc,desc,sweep,phase1,s8}-R.txt`). Search A's cover chunks for eight (`c0.txt` …
`c7.txt`) are `data/p8.txt` cut into consecutive runs of 1,172 lines, and for seven `data/p7.txt` cut into runs of 467
(`c4x.txt` is `c4.txt` without the line `p7-2222-a168-w4`, and `p2222.txt` that line alone). For seven, `seven_check
parta`'s 37 slices `parta-c*-s*.txt` (7.6 MB) aren't included either: they take 9 minutes on an open box, one after
another, and the ones made for this folder are identical to the originals apart from the timing field (Re-run, above);
its `maxc-c*-s*.txt` are all empty (no 7-partial spread is maximal). For six, search A's chunks are `data/p6.txt` cut into
runs of 43 and `six_check parta`'s `maxc-c*-s0.txt` are empty too; its core files (`p1asc.cores`, `p1desc.cores`,
`lin.cores`, `lindesc.cores`: 457 MB to 973 MB each) and `.prog` files' `P1` lines aren't here either (Left out, below).

The known-function lists used by the controls aren't included either: `apn_8bit.txt` (Beierle and Leander, Zenodo
4738942), `blp22_8bit.txt` (Beierle, Leander and Perrin, Zenodo 5821018), `new_apns.txt` (Beierle, Langevin, Leander,
Polujan and Rasoolzadeh, Zenodo 16752428, 3.4 GB), and `yuwangli_qam.txt` / `yuperrin_qam.txt` (the QAM lists of Yu,
Wang and Li with Edel and Pott's, and of Yu and Perrin, from sboxU's `apnDB`, one look-up table per line).

## Left out of `logs/`

- The full log of the `searchmem` run over the 11,004 Kurz covers (1.5 MB, one line per cover): the copy here keeps
  its header and summary lines.
- Runs superseded or repeated by the ones above: `fourteen_check search ... 0 0` on spread 0 alone (the same numbers
  as spread 0 in the full run); a first `yuwangli_qam.txt` control batch that crashed on the list's line 3 (an older
  build, before an APN filter for list lines was added; the batch above is its re-run); a cancelled cover.cpp control
  with only 1 completing form (KEEP 1); and `cover` on a 50-cover sample of `max.txt`.
- For twelve: `cover` on all 376 and on a 9-configuration sample with the build before the per-dimension counts were
  added (the same 12,303,628 nodes, 0 solutions; the run above is its repeat); two first submissions of the
  `pspread enum` runs that exited at once (a malformed command line); search B's per-configuration lines of the
  `searchsp` run (961 KB; the copy here keeps its header and summary); and the same for the `searchmax` run (100,783 lines,
  ~14 MB).
- For eleven: search B's per-configuration lines of the `search12`, `search13`, `searchsp` and `searchmax` (sample)
  runs (4,512, 2,418, 26,572 and 14,208 lines, 7 MB together; the `*.summary.txt` copies keep the headers and summary
  lines); a test build's trial run before the real ones: Part A cut to the first few branches of each class (445
orbits), with `configs` and a control batch on that partial list (whose `FAIL` lines are functions landing on
configurations the partial list doesn't hold).
- For ten: search B's per-configuration lines of the `search11`, `search13`, `searchsp` and `searchmax` (sample) runs
  (8,987, 8,866, 73,073 and 15,108 lines, 15 MB together; the `*.summary.txt` copies keep the headers and summary
  lines). `logs/pspread10-output.txt` has both of the K = 10 run's starts (the program appends to it), and its job log
  both commands.
- For nine: search B's per-configuration lines of the `search10`, `search13`, `searchsp` (8 ranges) and `searchmax`
  (sample) runs (17,770, 22,165, 146,146 and 16,493 lines, 30 MB together) and those runs' `.prog` files (44 MB; their
  totals, retallied from them, are in `logs/rerun/t9-tally.log`), and the per-configuration lines of the three main
  runs, which their `-prog.txt` files hold as `R` lines; the `*.summary.txt` copies keep the headers and summary lines.
- For eight: the per-configuration lines of the three main runs (their `-prog.txt` files hold them as `R` lines) and of
  the `searchmax` sample run (14,110 lines; its `.prog` files too, 1.8 MB), each range's header and summary kept in the
  `*.summary.txt` copies; likewise the `search9` cross-run's per-configuration lines (50,967, 7.8 MB) and `.prog` files
  (15 MB), each range's header and summary kept in `logs/check8-s9-output.summary.txt`; a `search13` cross-run over thirteen_check's 31 maximal 13-configurations × 1,287
  (39,897), cancelled after 2,705 configurations with 0 solutions to free the box (every configuration it reaches is
  in the `search9` family's; `eight-check.md`, State); a trial of the numpy `search9` mode on two configurations
  (`check8-pytest`); and a third `ctlreps` job, cancelled before it ran a control (`check8-ctlreps1`).
- For seven: the per-configuration lines of the searches' `.out` files (their `R` lines hold the same record; each
  range's header and summary are kept in the `*.summary.txt` copies); the `.prog` files' `P1` lines (the distinct spans
  with every off-D point covered, up to 200,000 per range: 1.1 million lines, 70 MB); search B's core files
  `asc-r*.cores` and `desc-r*.cores` (114,702 cores in each order, 13 MB; their counts and fingerprints, per
  configuration, are in the `R` lines, and `orders` compares the files as sets when given them); search A's first run
  of chunk `c4`, cancelled after about 17 hours with its last task unfinished and so without a `done:` line or a usable
  count (2.8 MB, 27,693 `NEAR` lines; redone as `c4x` and `p2222` above); the first `selftest` job, cancelled after 10
  minutes with the first `asc` job and run again; a cost probe of configuration 1527 (`SEVEN_DMAX=6`, cancelled after
  3 minutes); and first submissions of the split regression and of `mapa` that ran nothing (their command lines were
  mangled before they reached `cmd`).
- For six: search A's `NEAR` lines, 7,941,871 in the final run (790 MB, with its progress lines; the copies
  `logs/six-output-c*.txt` keep everything else, as do `logs/six-cover4-output-c*.txt` for the first run's 3,011,316, and
  `logs/cover5-prof-*-output.txt` for the 44,545 of each profile run), and search B's `P1` lines and core files
  (`p1asc.prog` and `p1desc.prog` 507 MB each, `asc.prog` over 78 MB; the `.cores` files 457 MB to 973 MB), together 3.7 GB.
  `logs/six-near.txt` keeps each run's count by dimension and a digest of its set of spans, per configuration; the `R`
  lines (`logs/check6-{p1asc,p1desc,asc,asc340}-R.txt`) carry the cores' counts and fingerprints, and `lin` and `lindesc`'s
  `.prog` files are here whole (`logs/check6-{lin,lindesc}-prog.txt`), `L` and `T` lines. From the searches' `.out`
  files (one process per configuration, each repeating a header and totals) `logs/check6-{p1asc,p1desc,asc,lin,lindesc}-output.txt`
  keep the header once, every `config` line, and, for `lin`, the last process's running totals. Also left out: two first
  submissions that ran nothing (`check6-selftest` and `cover5-reg6`, their command lines mangled before they reached
  `cmd`); `six_check`'s `tmp/` trial runs (other than the runs on configuration 340 here); the empty lists `plants-blp-*`,
  `plants-new-c2`, `plants-yp-c2`, `plants-ywl-c2` (no planting: `check6-plantfind-*-output.txt` say so); and, of `plantfind`'s
  runs, the output of the one that wrote `logs/check6-plants-bl22-c2.txt` (Re-run, above).
- Local paths in the logs were shortened to paths relative to the working tree.

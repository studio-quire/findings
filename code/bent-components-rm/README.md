# Checks for "The non-bent components of a quadratic function form a Reed–Muller codeword"

`checks.py` (Python 3 with numpy, a few minutes) runs eight checks; `expected-output.txt` is its output.

| Check | What it tests | Expected |
|---|---|---|
| 1 | Lemma 1: the bent indicator of a quadratic function has degree at most n/2 (Gold functions in dimensions 6, 8, 10; random quadratic (n, m)-functions) | degrees n/2 for the Gold functions, never above n/2 |
| 2 | Kasami–Tokura's weights in [d, 2d) against all codewords of RM(2,5) and RM(2,6) | {8, 12} and {16, 24, 28}, as predicted |
| 3 | n = 4, every subspace of Alt(F₂⁴): smallest non-MNBC number of non-bent components for m = 4..7 | 5, 11, 23, 47 = Theorem 2's bound |
| 4 | a quadratic (6,5)-function with 5 non-bent components, and the same with one or two zero coordinates appended | 5, 11, 23 = Theorem 2's bound |
| 5 | Theorem 4 on the APN functions of check 1 | holds |
| 6 | the (8,6)-function (xy, Tr(αx³), Tr(α²y³)) on F₁₆ × F₁₆ and its paddings to (8,7), (8,8); which components are not bent; every (xy, q₁(x), q₂(y)) with q₁, q₂ quadratic forms on F₁₆ | 5, 11, 23 = Theorem 2's bound; Tr(αx³), Tr(α²y³) and the three with a³ = α³; 50 choices, all Tr(λx³), Tr(μy³) with λ, μ not cubes and λμ a cube |
| 7 | Theorem 5: the inequality of its steps 1 and 6 for n = 8, 12, …, 48; for n = 8, step 4's condition by brute force over every set T of size ≥ 12 in F₂⁴ ∖ 0 and every word of RM(2,4) | contradiction for every n; at most 4 points of any T meet the condition, where all of T would have to |
| 8 | Theorem 6: Kasami, Tokura and Azumi's weights of RM(r, 2r) in [2d, 2.5d) that are 2 mod 4, for r = 5..12; their Table I entries of degree r evaluated directly for r = 5, 6, 7; the affine stabilizer of x₁x₂x₃x₄ + x₅x₆(x₇x₈ + x₃x₄) in AGL(8,2) | 2.5d − 6 and 2.5d − 2 only; every entry as the table states, the 2 mod 4 ones (24) and (26) (and (18) for r = 5); order 663,552, orbit = their count of weight-34 words of RM(4,8) |

`python checks.py 1 --luts FILE` also runs check 1 on every 8-bit LUT in FILE (one LUT per line, any separators), for
instance Beierle and Leander's `apn_8bit.txt` (Zenodo 4738942).

`python checks.py 7 --luts FILE` also checks, on every 8-bit APN function in FILE, the radical partition and, on every pair of rank-4 components, the two facts of Theorem 5's steps 3 and 4 (their sum is bent; the Pfaffian of each form on the other's radical is 1). On 521 functions from Beierle and Leander's and Beierle, Leander and Perrin's lists, 517 of them with two to four rank-4 components: 1,878 pairs, all as stated.

# Checks for "The non-bent components of a quadratic function form a Reed–Muller codeword"

`checks.py` (Python 3 with numpy, a few minutes) runs six checks; `expected-output.txt` is its output.

| Check | What it tests | Expected |
|---|---|---|
| 1 | Lemma 1: the bent indicator of a quadratic function has degree at most n/2 (Gold functions in dimensions 6, 8, 10; random quadratic (n, m)-functions) | degrees n/2 for the Gold functions, never above n/2 |
| 2 | Kasami–Tokura's weights in [d, 2d) against all codewords of RM(2,5) and RM(2,6) | {8, 12} and {16, 24, 28}, as predicted |
| 3 | n = 4, every subspace of Alt(F₂⁴): smallest non-MNBC number of non-bent components for m = 4..7 | 5, 11, 23, 47 = Theorem 2's bound |
| 4 | a quadratic (6,5)-function with 5 non-bent components, and the same with one or two zero coordinates appended | 5, 11, 23 = Theorem 2's bound |
| 5 | Theorem 4 on the APN functions of check 1 | holds |
| 6 | the (8,6)-function (xy, Tr(αx³), Tr(α²y³)) on F₁₆ × F₁₆ and its paddings to (8,7), (8,8); which components are not bent; every (xy, q₁(x), q₂(y)) with q₁, q₂ quadratic forms on F₁₆ | 5, 11, 23 = Theorem 2's bound; Tr(αx³), Tr(α²y³) and the three with a³ = α³; 50 choices, all Tr(λx³), Tr(μy³) with λ, μ not cubes and λμ a cube |

`python checks.py 1 --luts FILE` also runs check 1 on every 8-bit LUT in FILE (one LUT per line, any separators), for
instance Beierle and Leander's `apn_8bit.txt` (Zenodo 4738942).

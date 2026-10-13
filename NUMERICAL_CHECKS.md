# Numerical checks

Checked on 12 October 2026 with Wolfram 14.3.0 for Windows (64-bit), dated
8 July 2025. The official Wolfram/AgentTools MCP package 2.2.0 was used during
development and to run the test suite. MCP is not a release dependency.

## Data preservation and main results

All 13 retained array archives, containing 269 arrays, were read back
through Mathematica and verified exactly. `figure_data/array_conversion.json`
records the archive hashes, converted paths and shapes.

The 16 tests in `Tests/Core.wlt` passed through the Wolfram MCP TestReport tool.
They include fresh pair evolution, principal reconstruction, response-kernel
regeneration, canonical energies, six Kerr scaling exponents and GR residue
checks. The evolution, reconstruction and regeneration tests verify that fresh
runs of the released code agree with the retained data. The numerical checks below are distinct from the physical
approximation errors quoted in the paper.

| Quantity | Mathematica result |
| --- | --- |
| Infinity mismatch to the midpoint, 0-400M | 9.89378159961143e-7 |
| Initial canonical energies | 0.3284468105600703; 0.3515542044861051 |
| Peak cubic differences / common peak | 8.10947781837%; 11.1231719292% |
| Native complex L2 reconstruction errors | 1.0210317770e-5; 7.3019181685e-6 |
| Maximum isolated-pole factorization discrepancy | 2.5840794423e-10 |
| Kerr pole factor at chi = 0.999999 | 1.5072533194642 |
| Smooth-packet response-row vs direct evolution | Relative difference below 4.7e-13 |
| Causal full-fit coefficient matrix recovery | Maximum relative error 1.17e-14 |

## Notes on the native runs

**Two-pole integration.** The native 40-digit complex-ray integration gives
bounded closure **0.259555041620**. The native radius sequence is in
`results/two_pole.json`; the paper's matrices remain in
`source/output/schwarzschild_direct_exchange_rank_audit/`. The causal workflow
uses the retained paper matrix by default. To use the newly integrated matrix,
pass `QCSTwoPole[]["matrices"]` as the second argument of
`QCSCausal[100000, matrixAssociation]`.

**RWZ benchmark.** At N = 80 with a 12-term tail, native direct evaluation of
the retained Wolfram source gives radial ratio
`0.573212911473 + 0.037401087153 i` and normalized TT coefficient
`0.918741150276 - 0.062665238942 i`. The source, N/tail scan and reference
values are retained; no higher-order convergence claim is added.

**Pole scattering.** The twelve native refinements give maximum relative error
**1.1915457365e-7**. This rounds to the paper's **1.2e-7** and exceeds the
historical **1e-7** uniform target. The JSON result reports
`original_uniform_target_passed: false`. The workflow returns its diagnostics
normally; this is not a claim that the target passed.

**Monte Carlo forecasts.** Native forecasts use 100,000 complex-Gaussian draws
with Wolfram's MersenneTwister and seed 20260814; the contact ensemble uses
seed 20260816. The reported binomial standard error conditions on the
estimated threshold; independent estimation of the 99th-percentile null
threshold adds uncertainty.

**Continuum scattering.** All eight native window/cutoff/refinement cases ran.
For the 409.6M record, cutoff 16, step 0.01 case, the largest of the five
spatial point errors is **1.1969308354e-6** over 0.5-70M. This is the maximum
point error, not an aggregate norm over all points. The full results are in
`results/scattering.json`.

## What was run

Fresh native runs covered `check`, `pair`, `pair-from-operators`, `gr`,
`flux`, `wavepacket`, `two-pole`, `causal`, `kerr`, `rwz`, `scattering`,
`pole-scattering` and the three numerical figures. Response kernels at n=32,
48 and 64 were regenerated; the n=32 regeneration is also exercised in the
test suite. The five figure PDFs are retained under `figures/`; the
Mathematica plots reproduce the same data.

This release is a Mathematica implementation with preserved reference data.

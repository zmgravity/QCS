# Upstream numerical resources

## Kerr frequencies and excitation factors

Hayato Motohashi, *Kerr quasinormal mode frequencies and excitation factors*,
version 0.2.0, 19 February 2025, [Zenodo](https://zenodo.org/records/14380191),
version DOI [10.5281/zenodo.14380191](https://doi.org/10.5281/zenodo.14380191).
The upstream archive is `KerrQNMEFs-2.zip`; its published MD5 is
`4f14846fa2205f3409bf6f1e1ec3120c`. This snapshot retains 26 extracted mode
files, identified individually by SHA-256 in MANIFEST.sha256.

The nine table columns are a/M, Re(2M omega), Im(2M omega), Re(A), Im(A),
Re(nu), Im(nu), Re((2M)^(-2s) B), Im((2M)^(-2s) B). The calculations divide
frequencies by 2 and gravitational excitation factors by 16 to use M=1.

The dataset requests citation of both its record and the associated article:
H. Motohashi, *Resonant Excitation of Quasinormal Modes of Black Holes*,
Phys. Rev. Lett. 134, 141401 (2025),
[10.1103/PhysRevLett.134.141401](https://doi.org/10.1103/PhysRevLett.134.141401).

The upstream API metadata, checked on 2 October 2026, specifies
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
The retained tables are unmodified; mode selection and unit conversion are
documented above. The metadata response is preserved in
`source/external_inputs/KerrQNMEFs-2/zenodo_metadata.json`. Citation, attribution, and a link to the
license accompany these data.

## Quadratic GR sources and response tables

The [QuadraticQNM repository](https://github.com/akuntz00/QuadraticQNM) accompanies
the work of Bruno Bucciotti, Leonardo Juliano, Adrien Kuntz, and Enrico Trincherini
on amplitudes and polarizations of quadratic Schwarzschild QNMs. The associated
QCS paper, [arXiv:2610.04595v1](https://arxiv.org/abs/2610.04595v1), supplies the
research-paper citations; its bibliography is available with the paper.

The GR package records commit `374d354ce276f49bcadec8c07dd657e5e6c54d54`.
That provenance applies to its vendored inputs under
`source/gr_exchange/investigation/public_inputs/`. The broader local copy under
`source/external_inputs/QuadraticQNM/` is identified by this archive's file hashes;
no separate commit was recorded for that copy. The retained inputs are
identified by this release's MANIFEST.sha256.

The original MIT license (copyright 2024 Adrien Kuntz) is included in both input
directories. The specialized expressions in source/analysis and the exported
GR operators derive from these upstream resources; their generating programs
and retained samples are included.

## Original materials

Historical implementation attribution for the original calculation programs,
preparation coefficients, numerical results, and packaging documentation remains
Mian Zhu. The associated arXiv paper lists Mian Zhu, Jun Tao and Peng Wang.
Original code is distributed under MIT and original data, documentation and
figure artwork under CC BY 4.0, as specified in LICENSE.md. Upstream material
and derived expressions retain their existing notices and terms. Wolfram Mathematica is installed separately and is governed by its own license.
This release was checked with version 14.3.

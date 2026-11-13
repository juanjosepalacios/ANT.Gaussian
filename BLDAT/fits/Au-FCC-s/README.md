How the `FCC-s` Au Bethe-lattice set in `BLDAT/BL079.dat` was made: a CRYSTAL23 calculation of bulk
fcc Au with one s orbital and one electron per atom, and a fit of a single-orbital Bethe lattice to
its DOS. Use the set with `BLPAR1/2 = FCC-s` for electrode atoms that carry only the s basis below,
as the outer atoms of `examples/Au-nanocontacts/Au4-Au-Au--Au4.gjf`.

## Basis and ECP

The s-only Au atoms of the examples: Gaussian's `SHF` ECP, which leaves one valence electron (6s),
and one s shell of three primitives (exponents 0.4409, 0.2626, 0.0617; coefficients -0.7478, 0.8176,
0.8113). `crystal/INPUT` gives both in CRYSTAL's format (the ECP as `INPUT`, Z = 1, with its local
term and its s, p and d semilocal terms).

## CRYSTAL23 (`crystal/`)

fcc Au, a = 4.078 A (experiment), SVWN, not spin-polarized, `SHRINK 24 24`, `SMEAR 0.01`,
`FMIXING 80`, `TOLDEE 8`. It converges in 6 cycles to E = -0.0835723 Ha with E_F = 0.26327 Ha.

`d3w` writes the total DOS (the only orbital) from E_F - 0.6 to E_F + 2.0 Ha (`DOSS.DAT`, states per
Hartree per cell, both spins). The band runs from -14.73 to +16.79 eV around E_F and holds 1.995
states; 0.998 are below E_F, so the band is half filled. Per spin, the DOS is nearly flat at about
0.035 states/eV from the bottom to E_F (0.0359 at E_F), peaks at 0.067 at +6.7 eV and has a tail up
to +16.8 eV: a nearly-free-electron band, not symmetric about E_F.

## Fit (`fit/`)

A single-orbital, orthogonal Bethe lattice with 12 directions (fcc, <110>) has a symmetric DOS of
width 4 sqrt(11) |sss| centred on es, so at half filling E_F = es, and only |sss| sets its shape.
`bethe1.py` gives its DOS per spin (broadening ETA = 0.1 eV, as in the other fits); `fits.py` compares
three choices of |sss| against the CRYSTAL23 DOS per spin (relative L2 misfit from -15 to +4 eV):

| Choice | \|sss\| (eV) | Band (eV) | DOS(E_F) per spin (states/eV) | Misfit |
|---|---|---|---|---|
| A: DOS at E_F | 2.4412 | +-16.19 | 0.0359 | 0.141 |
| B: occupied band width (bottom at -14.73 eV) | 2.2206 | +-14.73 | 0.0394 | 0.166 |
| C: least squares | 2.3974 | +-15.90 | 0.0365 | 0.140 |

`FCC-s` is A: it reproduces the DOS at E_F, which is what transport sees, with a misfit as small as
the least-squares one. sss is taken negative (an s bond; the sign does not change the DOS of a Bethe
lattice) and es = 0 (ANT shifts all on-site energies to put E_F at zero). The earlier 1-orbital set
`lsda-6s` (removed from BL079.dat after this fit) had sss = -0.0814 Ha = -2.21 eV, close to choice
B. `dos_fit.dat` has the CRYSTAL23 and `FCC-s` DOS per spin on the CRYSTAL23 grid.

The Bethe lattice cannot reproduce the asymmetry of the band: it puts the same weight above and below
E_F, so the CRYSTAL23 peak at +6.7 eV is missing and the band bottom is 1.5 eV lower. An overlap
`S_sss` would skew the band, but non-orthogonal Bethe-lattice electrodes have not worked well in ANT
(see the `Papacon-nonortho` notes in `examples/Al-nanocontacts`), so it is not fitted.

## Check in ANT

ANT 3.0.0 with `FCC-s` on both electrodes of Au4-Au-Au--Au4 (`LEADDOS`, 13 November 2026) reads the
set, finds 12 lattice directions and a Bethe-lattice Fermi energy of exactly 0 eV. Its `Lead1DOS.dat`
(both spins) divided by 2 agrees with `bethe1.py` at ETA = 0.001 eV to 1e-6 states/eV over -20 to
+20 eV. The contact converges in 14 cycles to E = -272.526621 a.u., device E_F = -2.5120 eV and
T(E_F) = 1.578 e^2/h (1.553 with the earlier `lsda-6s`).

## Reproducing

- Run CRYSTAL23 with `crystal/INPUT`, then `properties` with `crystal/d3w` (writes `DOSS.DAT`);
  `crystal/out` is the CRYSTAL23 output.
- `cd fit; python3 fits.py` prints the table above and writes `dos_fit.dat` (needs numpy and scipy).

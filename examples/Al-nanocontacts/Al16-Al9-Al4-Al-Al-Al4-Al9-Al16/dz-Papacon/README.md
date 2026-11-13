Case dz-Papacon of the [Al16-Al9-Al4-Al-Al-Al4-Al9-Al16 example](../README.md): LANL2DZ basis and ECP
on the inner atoms 17-44, minimal sp basis and ECP10CE on the outer 16-atom layers, orthogonal
Bethe-lattice set `Papacon`.

Run `Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.gjf` with `Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.ant`, with
`Scf=(conver=6,dsymm,nodamp)`. The `.gjf` reads the basis and ECP from
`Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.gbsdz` (352 basis functions). `.ant` settings: `PACC = 1.0d-6`,
`CHARGEACC` and `FERMIACC` 1.0d-4, `NPULAY = 5`, `SL = 0.001`, LDOS on the apex atom 30, energy
window `EW1 = -5` to `EW2 = 5` eV.

## Outputs

Results and the timing study are in the [structure README](../README.md#results).

- `Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.log`: the 8-core run
- `T.Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.dat`: transmission from E - E_F = -5 to +5 eV
- `DOS.Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.dat`: energy, total DOS, LDOS of atom 30, then its
  orbital-resolved LDOS
- `Al16-Al9-Al4-Al-Al-Al4-Al9-Al16.xyz`: the geometry written by ANT; the 82 `Ar` atoms mark the
  Bethe-lattice attachment sites

## Notes

- The `.gjf` uses `%NProcShared=8` (135 s); change it to the number of cores available. It
  substitutes links from `../../../../bin`.

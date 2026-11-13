Case dzp of the [Al9-Al4-Al-Al-Al4-Al9 examples](../README.md): LANL2DZ plus one d shell on the
central atoms 10-19, minimal sp basis and ECP10CE on the outer layers, orthogonal Bethe-lattice set
`Papacon`.

Run `Al9-Al4-Al-Al-Al4-Al9.gjf` with `Al9-Al4-Al-Al-Al4-Al9.ant` (`CHARGEACC` and `FERMIACC` 1.0d-5,
`FULLACC`, `DEVBRENT`). The `.gjf` reads the basis and ECP from `Al9-Al4-Al-Al-Al4-Al9.gbsdzp` (202
basis functions).

`P.Al9-Al4-Al-Al-Al4-Al9.dat` is the converged density matrix of this run; the bias case starts from
it. ANT reads it at start-up, so delete it to run this case from scratch.

## Outputs

Results are in the [structure README](../README.md#results).

- `Al9-Al4-Al-Al-Al4-Al9.log`
- `T.Al9-Al4-Al-Al-Al4-Al9.dat`: transmission from E - E_F = -5 to +5 eV
- `DOS.Al9-Al4-Al-Al-Al4-Al9.dat`: energy, total DOS, LDOS of atom 14, then its orbital-resolved LDOS
- `Al9-Al4-Al-Al-Al4-Al9.xyz`: the geometry written by ANT; the `Ar` atoms mark the Bethe-lattice
  attachment sites

## Notes

- The `.gjf` uses `%NProcShared=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.

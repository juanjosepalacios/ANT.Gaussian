Case dzp-bias of the [Al9-Al4-Al-Al-Al4-Al9 examples](../README.md): dzp-Papacon at a bias voltage of
1 V (`BIASVOLT = 1.0`), restarted from its zero-bias density matrix.

Copy `P.Al9-Al4-Al-Al-Al4-Al9.dat` from [dzp-Papacon](../dzp-Papacon/README.md) into this directory,
then run `Al9-Al4-Al-Al-Al4-Al9.gjf` with `Al9-Al4-Al-Al-Al4-Al9.ant`. ANT reads the density matrix at
start-up and connects the leads from the first cycle. The `.gjf` has `nosymm` in place of `dsymm`
in `Scf=(...)`, as required for finite bias, and reads the basis from `Al9-Al4-Al-Al-Al4-Al9.gbsdzp`
(202 basis functions).

## Outputs

Results are in the [structure README](../README.md#results).

- `Al9-Al4-Al-Al-Al4-Al9.log`
- `T.Al9-Al4-Al-Al-Al4-Al9.dat`: transmission from E - E_F = -5 to +5 eV
- `DOS.Al9-Al4-Al-Al-Al4-Al9.dat`: energy, total DOS, LDOS of atom 14, then its orbital-resolved LDOS
- `Al9-Al4-Al-Al-Al4-Al9.xyz`: the geometry written by ANT; the `Ar` atoms mark the Bethe-lattice
  attachment sites

## Notes

- `nosymm` inside `Scf=(...)` does not stop Gaussian from reorienting the contact; it uses the same
  standard frame as dzp-Papacon, so the restart density matrix matches.
- The `.gjf` uses `%NProcShared=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.

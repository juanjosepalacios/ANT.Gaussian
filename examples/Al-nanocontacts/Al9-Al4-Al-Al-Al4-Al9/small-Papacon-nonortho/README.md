Case small-nonortho of the [Al9-Al4-Al-Al-Al4-Al9 examples](../README.md): small-Papacon with the
non-orthogonal Bethe-lattice set `Papacon-nonortho` at `OVERLAP = 0.5`.

Run `Al9-Al4-Al-Al-Al4-Al9.gjf` with `Al9-Al4-Al-Al-Al4-Al9.ant` (`BLPAR1 = BLPAR2 = Papacon-nonortho`,
`OVERLAP = 0.5`) and `Scf=(conver=6,...)`. The `.gjf` reads the basis and ECP from
`Al9-Al4-Al-Al-Al4-Al9.gbssmall` (112 basis functions).

## Outputs

Results are in the [structure README](../README.md#results).

- `Al9-Al4-Al-Al-Al4-Al9.log`
- `T.Al9-Al4-Al-Al-Al4-Al9.dat`: transmission from E - E_F = -5 to +5 eV
- `DOS.Al9-Al4-Al-Al-Al4-Al9.dat`: energy, total DOS, LDOS of atom 14, then its orbital-resolved LDOS
- `Al9-Al4-Al-Al-Al4-Al9.xyz`: the geometry written by ANT; the `Ar` atoms mark the Bethe-lattice
  attachment sites

## Notes

- This SCF does not reach a physical solution: the device Fermi energy ends at +13.5 eV, 20 eV above
  small-Papacon, and the energy is 2 a.u. higher. Gaussian stops on `Density matrix is not changing`.
  The case documents the current behaviour of non-orthogonal Bethe-lattice electrodes; do not use
  its transmission.
- `OVERLAP = 0.5` halves the book overlaps. With `OVERLAP = 1.0` (the book values) ANT spent more
  than 7 minutes in the Bethe-lattice energy-window search and was stopped.
- The `.gjf` uses `%NProcShared=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.

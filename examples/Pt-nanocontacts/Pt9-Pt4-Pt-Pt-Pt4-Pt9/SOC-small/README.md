Case 1 of the [Pt9-Pt4-Pt-Pt-Pt4-Pt9 examples](../README.md): SOC with the small basis set.

Run `Pt9-Pt4-Pt-Pt-Pt4-Pt9.gjf` with `Pt9-Pt4-Pt-Pt-Pt4-Pt9.ant` (`SOC`, `LEADDOS`). The `.gjf`
reads the basis from `Pt9-Pt4-Pt-Pt-Pt4-Pt9.gbssmall` (252 basis functions).

`P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat` is the density matrix of this run. ANT reads it at start-up as the
initial density matrix, so delete it to run from scratch.

## Outputs

These come from an earlier version of ANT (10 April 2024).

- `Pt9-Pt4-Pt-Pt-Pt4-Pt9.log`
- `T.Pt9-Pt4-Pt-Pt-Pt4-Pt9.SOC.dat`: total transmission in units of e^2/h
- `Lead1DOS.dat`, `Lead2DOS.dat`: electrode DOS written with `LEADDOS`
- `P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat`: density matrix
- `Pt9-Pt4-Pt-Pt-Pt4-Pt9.xyz`: the geometry written by ANT; the 56 `Ar` atoms mark the Bethe-lattice attachment sites
- `name`: the job name file written by ANT

## Notes

- The `.gjf` uses `%nproc=8`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.

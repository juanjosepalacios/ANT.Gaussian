Case 2 of the [Pt9-Pt4-Pt-Pt-Pt4-Pt9 examples](../README.md): SOC with the large basis set and
modified values of the SOC integrals.

Run `Pt9-Pt4-Pt-Pt-Pt4-Pt9.gjf` with `Pt9-Pt4-Pt-Pt-Pt4-Pt9.ant` (`SOC`, `LEADDOS`). The `.gjf` is
the one of [case 1](../SOC-small/README.md) with the basis read from
`Pt9-Pt4-Pt-Pt-Pt4-Pt9.gbslarge` (442 basis functions).

## Outputs

These come from the rerun of 11 November 2026 with the corrected Pt d on-site energies; results are
in the [structure README](../README.md#soc-results-cases-1-and-2).

- `Pt9-Pt4-Pt-Pt-Pt4-Pt9.log`
- `T.Pt9-Pt4-Pt-Pt-Pt4-Pt9.SOC.dat`: total transmission in units of e^2/h

## Notes

- The `.gjf` uses `%nproc=8`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.

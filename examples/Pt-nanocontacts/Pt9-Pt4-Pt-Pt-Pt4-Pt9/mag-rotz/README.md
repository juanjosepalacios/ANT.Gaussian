Case 4 of the [Pt9-Pt4-Pt-Pt-Pt4-Pt9 examples](../README.md): case 3 restarted with every atom
listed in `SPINROTATOM` at theta = 0, phi = 0, which rotates nothing but goes through the same
spinor code as case 5.

Copy `../mag/P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat` here, then run `Pt9-Pt4-Pt-Pt-Pt4-Pt9.gjf` with
`Pt9-Pt4-Pt-Pt-Pt4-Pt9.ant`. The `.gjf` and basis are those of [case 3](../mag/README.md).

## Outputs

- `Pt9-Pt4-Pt-Pt-Pt4-Pt9.log`
- `T.Pt9-Pt4-Pt-Pt-Pt4-Pt9.SOC.dat`: total transmission in units of e^2/h

Results and further notes are in the [structure README](../README.md#spin-rotation-test-cases-3-5).

## Notes

- The `.gjf` uses `%nproc=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.

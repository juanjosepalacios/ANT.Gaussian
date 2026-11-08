Case 3 of the [Pt nanocontact examples](../README.md): collinear spin-polarized calculation of a Pt
contact forced to be magnetic, the reference of the spin rotation test.

Run `Pt9-Pt4-Pt-Pt-Pt4-Pt9.gjf` with `Pt9-Pt4-Pt-Pt-Pt4-Pt9.ant`: charge 0, multiplicity 15,
`nosymm`, `NSPINLOCK = 1000`, no SOC, small basis (`Pt9-Pt4-Pt-Pt-Pt4-Pt9.gbssmall`).

`P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat` is the density matrix of this run, which cases 4 and 5 restart from.
ANT reads it at start-up as the initial density matrix, so delete it before repeating this case from
scratch.

## Outputs

- `Pt9-Pt4-Pt-Pt-Pt4-Pt9.log`
- `T.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat`: transmission in units of e^2/h, spin up and then spin down, separated by a blank line
- `P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat`: density matrix

Results and notes are in the [parent README](../README.md#spin-rotation-test-cases-3-5).

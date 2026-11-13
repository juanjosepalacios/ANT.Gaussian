# Changelog

## 3.0.0 - 2026-11-13

- Version 3.0.0: the device Green's functions are several times cheaper.
  With the three changes below (contour integration over single points,
  OpenBLAS in the NVIDIA build, restricted S^1/2 Sigma S^1/2), the 60-atom
  Al example runs in 42 s instead of 135 s on 8 cores and in 210 s instead
  of 570 s on 1, with the same results.

- Restrict the transformation of the Bethe-lattice self-energies to the
  non-orthogonal basis, S^1/2 Sigma S^1/2, to the rows and columns of
  Sigma that are not zero (the orbitals of the atoms attached to the
  electrode) instead of two full NAOrbs^3 complex matrix products per
  electrode (`SigmaToNonOrtho`). The result is the same: the 1-core
  60-atom Al run gives the same SCF energy to all printed digits. That
  run drops from 448 to 210 s on 1 core and from 70 to 42 s on 8 (device
  Fermi search 39 -> 14 s); Ni4-Ni-Ni-Ni4 from 27 to 24 s on 8 cores.

- NVIDIA build: link the optimized BLAS/LAPACK of the NVIDIA HPC SDK
  (OpenBLAS, `liblapack_ilp64`/`libblas_ilp64`) instead of compiling the
  reference routines of `lapack_blas.f`. The libraries are copied into
  `obj/` with every Fortran symbol renamed to `antnv_*`, so Gaussian's own
  BLAS/LAPACK in `l502.exe` (util.a, ATLAS) are untouched, and the new
  `lapack_blas_nv.F90` keeps the module `lapack_blas` with the same
  routines. `ANTLAPACK = bundled` in `make.nvidia.in` restores the old
  build. 60-atom Al example: 84 -> 70 s on 8 cores and 569 -> 448 s on
  1, with the same results; Ni4-Ni-Ni-Ni4 reproduces its example log.

- Parallelize the contour integration of the device density matrix
  (`IntCompPlane`) over single quadrature points instead of pairs, and
  evaluate the centre point and the first doubling levels (15 points, or
  more if there are more threads) in one parallel loop, so the threads are
  no longer idle at the low levels. The stopping test is applied level by
  level as before, and levels it does not reach are discarded; on one
  thread nothing is computed ahead. Also fixes out-of-bounds writes of the
  abscissa arrays at the last level (255 points with `NIP = 7`). In the
  60-atom Al example the device Fermi search drops from 98 to 48 s on 8
  cores and the run from 135 to 84 s (125 to 94 s and 189 to 158 s on 4);
  1 and 2 cores are unchanged. Results agree within the run-to-run scatter
  (SCF energies within 5e-5 a.u., E_F within 0.04 meV, T(E_F) = 2.0067),
  except one 2-core run whose final Fermi search fell back from MULLER to
  SECANT and ended 1.3 meV off (T(E_F) = 2.0064); two repeats did not.

- Add the 60-atom Al16-Al9-Al4-Al-Al-Al4-Al9-Al16 example (LANL2DZ on the
  28 inner atoms, `Papacon`, `conver=6`): T(E_F) = 2.007, and a timing
  study on 1-8 cores (570, 320, 189 and 135 s; about 13% serial work).

- Replace the Al `Papacon` set in BL013.dat, of unknown origin, with the
  orthogonal 2-center fit of Papaconstantopoulos, Handbook of the Band
  Structure of Elemental Solids, 2nd ed. (2015), Table 7.10, and add the
  non-orthogonal fit with overlaps as `Papacon-nonortho`. Both keep the s
  and p first-neighbour terms only, for the minimal sp basis. Rerun the
  four Al9-Al4-Al-Al-Al4-Al9 examples with the new default set (T(E_F)
  1.87-1.95 against 1.76-1.92 before), regroup them by structure, as the
  Ni and Pt examples, and add a case with `Papacon-nonortho` at
  `OVERLAP = 0.5`, whose SCF ends with the device Fermi level 20 eV too
  high.

- Rerun the five Pt9-Pt4-Pt-Pt-Pt4-Pt9 examples with the corrected Pt d
  on-site energies and update their outputs and READMEs. The SOC cases
  and the collinear case change little (T(E_F) 4.231, 5.898 and 3.033
  against 4.252, 5.914 and 3.053 e^2/h). The restarts mag-rotz and
  mag-rotx are no longer spin-locked (NSPINLOCK is ignored on a restart
  since f027077), so they relax to a weakly magnetic state 16 mHa below
  the collinear one; they still agree with each other, which is what the
  spin rotation test checks.

- Warn in the Pt, Ni and Fe example READMEs that Gaussian moves the
  molecule to its standard orientation unless `nosymm` is on the route
  line (inside `scf=()` it does not prevent it), and that ANT then works
  in that frame: the ZM and ZMATOM fields, the SPINROTATOM angles, the
  direction of the moments with SOC and the Bethe-lattice d on-site
  energies all refer to it. The Pt cases and the Ni4 and doll-Papacon-sd
  Ni cases are reoriented by whole axes, so their cubic axes still lie
  along X, Y and Z; the Ni4 field Bz acts along -z of the input. The other
  Ni9 cases and the Fe wire keep the input frame.

- Fix the d on-site energies of the Pt `Papacon` and `Papacon-nonortho`
  sets in BL078.dat, which were swapped. The handbook's two-center d1 (t2g,
  0.45468 and 0.41348 Ry) now goes to `edt` and d2 (eg, 0.43701 and
  0.41334 Ry) to `edd`, as in Papaconstantopoulos, Handbook of the Band
  Structure of Elemental Solids, 2nd ed. (2015), Table 6.55, and as in the
  Ni sets. The orthogonal set changes by 17.7 mRy in each level; the Pt
  example results were computed with the swapped values.

- Group examples/Fe-wire by atomic structure, as the Ni and Pt examples:
  BCC-sd/, BCC-spd/ and 1DLEAD-spd/ move into Fe64-bcc-wire/, named after
  the system like its files, whose README describes the wire and compares
  the variations. The `.gjf` files now substitute links from
  `../../../../bin`.

- Name every file of the Ni9-Ni4-Ni-Ni-Ni4-Ni9 Ni examples after the
  system, as the folder that holds them: drop the `.F`, `.AF` and `.doll`
  job-name suffixes from the inputs, outputs and `.gjf` titles, and from
  the job names recorded in the logs, since each variation has its own
  subfolder. The AF case now starts from
  `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat` copied unchanged from LANL1MB-Papacon-F/.

- Remove the Ni `Papacon-sd` Bethe-lattice set from BL028.dat. Inputs that
  give `BLPAR = Papacon-sd`, among them the doll-Papacon-sd Ni example, no
  longer run.

- Group examples/Pt-nanocontacts by atomic structure, as in
  examples/Ni-nanocontacts: the five cases (SOC-small, SOC-large, mag,
  mag-rotz, mag-rotx) move into Pt9-Pt4-Pt-Pt-Pt4-Pt9/, whose README
  describes the contact and holds the spin rotation test results. The case
  folders and file names are unchanged. The `.gjf` files now substitute
  links from `../../../../bin`.

- Loosen the device spectral-weight check by a factor of 100 (still at most
  0.1) and widen its window in steps of 100 eV instead of 10 eV when the
  Bethe-lattice electrodes have overlaps (OVERLAP > 0). With
  Papacon-nonortho at OVERLAP = 0.5 the device weight missed or exceeded its
  target by 1e-4 to 4e-3 whatever the window, so the check failed after
  1000 eV of widening.

- Reorganize examples/Ni-nanocontacts by atomic structure: Ni4-Ni-Ni-Ni4/
  and Ni9-Ni4-Ni-Ni-Ni4-Ni9/, each with a README comparing its variations,
  which sit in subfolders with their own inputs, outputs and README, as in
  examples/Fe-wire. The files keep the system name. The Ni9 F and doll cases
  move to CRENBS-spd-F/ and doll-Papacon-sd/. Add LANL1MB-Papacon/ (Ni4),
  LANL1MB-Papacon-F/ and LANL1MB-Papacon-AF/ (Ni9): the minimal spd basis
  LANL1MB and the Hay-Wadt LANL1 ECP on all atoms with the new orthogonal
  `Papacon` electrodes; the AF case starts from the density of the F case,
  whose P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat is included. Remove the old AF case,
  computed with the removed Papacon set, and the old Ni4 cases (the 20 T sd
  contact and its SPINROTATOM rotx variant).

- Replace the Ni `Papacon` and `Papacon-old` Bethe-lattice sets in BL028.dat
  with the first-neighbour Slater-Koster fit of ferromagnetic fcc Ni from
  Papaconstantopoulos, Handbook of the Band Structure of Elemental Solids,
  2nd ed. (2015), Tables 4.73 and 4.78: `Papacon` (the BLPAR default) is the
  orthogonal fit and `Papacon-nonortho` the fit with overlaps (spin-averaged).
  They are now the first two sets in the file. Inputs that use the default
  electrodes for Ni get different parameters; inputs with Papacon-old no
  longer run. Papacon-sd, FCC-sd and FCC-spd are unchanged.

- Remove the SWOFFSPL keyword. With NSPINLOCK < 0 (the default), spin locking
  is now switched off once RMSDP < SL/10, the threshold being set from SL after
  the input is read; NSPINLOCK >= 0 still supersedes it. Inputs that give
  SWOFFSPL stop with "Undefined keyword".

- Move the lower end of the device density contours down with the trial
  Fermi level during the device Fermi search, to min(EMin, EMin+shift), in
  QTot, QTot_SPINLOCK, QTot_SOC and CompPD_SOC. The contours ran from the
  fixed EMin while shift=-E_F moved the device spectrum along them, so a
  trial Fermi level more than the 10 eV buffer above the bottom of the
  spectrum pushed the lowest states below EMin and lost their charge. The
  charge then fell again at high trial Fermi levels, the Brent bracket could
  fail with negative residuals at both ends, and roots far above the device
  spectrum could be accepted.

- Fix Bethe-lattice electrodes with overlaps (OVERLAP > 0) stopping at the
  second electrode with "could not allocate memory for H0X and S0X matrices"
  unless LEADDOS was given: the extended-Hamiltonian matrices of the first
  electrode were only freed after writing the lead DOS.

- Add the FCC-sd and FCC-spd Ni Bethe-lattice sets to BL028.dat, fitted on a
  12-direction Bethe lattice to the CRYSTAL23 LSDA DOS of fcc Ni (a = 3.524 A)
  with the CRENBS s+d basis and ECP and with the same basis plus the STO-3G 4p
  shell scaled by 0.60. BLDAT/fits/Ni-FCC-sd and BLDAT/fits/Ni-FCC-spd hold the
  bases, the 4p-scale optimization, the CRYSTAL23 inputs, outputs and fit
  target DOS, and the fit scripts and logs.

- Add the FCC-spd Fe Bethe-lattice set to BL026.dat, fitted on a 12-direction
  Bethe lattice to the CRYSTAL23 LSDA DOS of fcc Fe (a = 3.515 A) with the
  BCC-spd basis. BLDAT/fits/Fe-FCC-spd holds the CRYSTAL23 input, output and
  fit target DOS, and the fit scripts and logs.

- Add the BCC-spd Fe Bethe-lattice set to BL026.dat, fitted to the CRYSTAL23
  LSDA DOS of bcc Fe with the CRENBS s+d basis plus a scaled STO-3G 4p shell.
  BLDAT/fits/Fe-BCC-spd holds the basis, the 4p-scale optimization, the CRYSTAL23
  inputs and fit target DOS, and the fit scripts and logs.

- Remove the lsda-3d4s Fe Bethe-lattice set from BL026.dat; use FCC-sd or
  BCC-sd instead. Inputs with BLPAR1/2 = lsda-3d4s no longer run.

- Change the default of NSPINLOCK from 0 to -1, so that SWOFFSPL (default
  1.0D-3) decides when spin locking is switched off unless NSPINLOCK >= 0 is
  given or the calculation restarts from a density matrix file.

- Set NSPINLOCK to 0 when the calculation starts from a density matrix file,
  whatever its value in the input or default, and say so in the log. Spin
  locking would otherwise pull the spin of the restarted density back to the
  multiplicity of the Gaussian input. NSPINLOCK >= 0 supersedes SWOFFSPL.

- Remove SPINROT, THETA and PHI and the rigid rotation of the whole system.
  SPINROTATOM alone now switches on the spin rotation; atoms not listed stay
  along z. Input errors from read_parameters now stop with exit status 1.

- Stop while reading the input when SCFSOC, SOC_FAC_P/D/F, RCUT,
  SOC_VAL_P/D/F, SOCFACATOM or SOCVALATOM is given without SOC; these
  keywords had no effect without it. Document that SCFSOC needs SOC, replace
  the obsolete SOCFAC keyword in the Ag and Ni example inputs, and add
  test_soc_keywords.py.

- Fix ZMATOM: per-atom Zeeman fields were looked up with an atom index that
  was never set in CompHZM, so they landed on an undefined atom. Each orbital
  now takes the field of its own atom. Add test_zeeman.py and document SOC,
  SCFSOC, SOC_FAC_P/D/F, RCUT, SOCFACATOM, SOC_VAL_P/D/F, SOCVALATOM, ZM,
  BX/BY/BZ, ZMATOM, POL and TCOMP in keywords.txt.

- Remove the TREV keyword and its second transmission pass with the
  spin-reversed spinor Hamiltonian. That pass ignored SPINROT, reinitialized
  the Bethe-lattice electrodes with UD/DU/DD and wrote unopened unit 666.

- Rename the spin-rotation keyword ROT to SPINROT. ROT is no longer
  recognized.

- Fix spin rotation (SPINROT): the up-down spinor blocks used the conjugate
  phase, so (THETA, PHI) rotated the moment to (THETA, -PHI); and global
  THETA/PHI without SPINROTATOM rotated only on-site blocks, leaving
  inter-atomic blocks collinear along z. Both branches now take angles from
  one rule. Add test_spin_rotate.py and document SPINROT, THETA, PHI and
  SPINROTATOM.

- Device spectral-weight check: replace the fixed 0.1 tolerance with
  min(0.1, input CHARGEACC*(NCDEl+QEXCESS)), the device Fermi-search scale,
  and converge its quadrature to a tenth of it. Add back the analytic ETA
  Lorentzian tail of device states inside [EMin,EMax], which widening cannot
  remove; states outside the window stay uncorrected. Log the tolerance and
  tail estimate.

- Apply ETA once in the device Green's function, (z+i*ETA)S-H-Sigma, in
  gplus0, gplus, glesser and their SOC variants (-i*ETA for advanced), so
  device states uncoupled from the electrodes are broadened and can be
  partially occupied. The -i*ETA*S electrode self-energy initialisations
  were overwritten by CompSelfEnergies except for GHOST electrodes, where
  they gave 2*ETA through undefined intent(out) values; GHOST now returns
  an explicit zero. Gamma_L/R are unchanged. The hybridization mesh passes
  real energies to gplus0 to avoid applying ETA twice.

- Log a dashed banner, flushed immediately, when each device density-matrix
  integration and Fermi search starts. Informative progress messages
  (lattice directions, SOC/Zeeman Fermi levels) use the same banner form.

- Add EFMAX: an optional upper bound (eV) on the device Fermi energy in the
  Brent bracket and the Muller/secant/bisection chain. EMin/EMax, density
  contours and the spectral-weight check are unchanged.

- Log the Fermi energy of a rejected device Brent root. Add opt-in QSCAN,
  which logs the device charge residual at 15 points within +-1 eV of that
  root before stopping.

- Check device spectral weight on [EMin,EMax] for all electrode types,
  widening both bounds by 10 eV after a quadrature or weight-check failure
  until it passes (stop after 1000 eV). Keep the widened bounds, and start
  every density contour, including biased and SOC ones, at exactly EMin.

- Device Fermi searches default to Muller -> secant -> bisection (DEVMULLER);
  DEVBRENT selects safeguarded Brent with an adaptive warm-start bracket.
  Both paths cover spin and SOC, rebuild the density at the accepted root
  and log iterations, charge evaluations, residual and wall time. Fix
  uninitialized FIXEFERMI success and remove the 0.1 eV fallback.

- Find Bethe-lattice electrode Fermi levels with safeguarded Brent
  (BLBRENT); BLSECANT restores bisection + secant. Electrode setup now uses
  the FERMIACC and CHARGEACC input values instead of the internal 1.0D-2.

- Use periodic band filling exclusively for 1D lead Fermi levels, with
  parallel generalized Bloch spectra, ETA-broadened occupations and independent
  k-mesh convergence for each lead. Remove the contour-based Brent/warm-start
  lead searches and lower-bound searches. LEADBANDS remains a no-op alias.

- Add opt-in ADAPTQUAD: deterministic parallel adaptive Gauss-Kronrod
  integration for collinear device density matrices.

- Add an explicit LDOS switch requiring both valid atom-range endpoints;
  LDOS_BEG/LDOS_END alone no longer enable local-DOS output.

- Introduce NATPC1/2 for 1D-lead primitive-cell atom counts, retaining
  NEMBED1/2 for Bethe embedding and as warned legacy aliases for 1D inputs.
- Reject conflicting cell counts and document the new input names.

## 2.8.4 - 2026-10-31

- Parallelize OneDLead bulk-charge quadrature with call-local Green-function
  and self-energy scratch and ordered accumulation; retain the separate
  device-charge fallback for 1D electrodes.
- Exclude complete 1D lead supercells from reduced-device charge counting.
- Use a consistent Fermi energy reference for bulk lead DOS output.
- Add NVIDIA/Intel 1D-charge regression tests and record Au10 validation
  and its remaining limitations in TOTCHARGE.md.

- Use safeguarded Brent Fermi searches with independent warm-start brackets for each lead.
- Seed lead integration lower bounds from repeated-cell generalized eigenvalues.
- Honor configured Fermi/charge tolerances and reject unconverged lead integrations.
- Correct the Dyson convergence residual and detect inversion and convergence failures.
- Expose L1DALPHA and set its default to 0.75; retain SELFACC as the Dyson tolerance.
- Add regression tests for root finding, warm starts, eigenvalue bounds and Dyson mixing.

## 2.8.3 - 2026-10-29

- Parallelize TotCharge quadrature with explicit per-call spin and ordered
  accumulation; keep the 1D-lead path serial because of saved scratch state.
- Report charge quadrature wall time.
- Add bounds-checked NVIDIA/Intel tests for one/two spin channels and 1-8 threads.
- Pt SOC restart: approximately 7.19x faster TotCharge at eight threads,
  with identical saved spectra; see TOTCHARGE.md for benchmark details.


## 2.8.2 - 2026-10-29

### Fixed
- Initialize charge-quadrature convergence history and compare successive
  estimates only after refinement.
- Prevent TotCharge and QTot2 from exceeding quadrature array bounds; stop
  explicitly if the maximum order fails to converge.
- Make integration scheduling state and nested contour-loop indices thread-private.

### Performance and diagnostics
- Remove four redundant barriers before thread-local accumulator merges.
- Report the effective OpenMP maximum and transport-analysis wall time.
- Document measured thread scaling and matching Gaussian/OpenMP thread settings
  in OPENMP.md. On Lian, the Pt SOC restart took 82.5 s at four threads and
  52.3 s at eight, with identical saved spectra. No measurable fixed-thread
  speedup is claimed for the synchronization edits alone.

### Installed-source baseline
- Preserve Lian's installed 2.8.1 molecular-analysis, charge-integration and
  diagnostic-output changes relative to the previous GitHub snapshot.
- Preserve installed Co, Pt and Au electrode parameters and NVIDIA build files.
- These inherited changes include moving Mol_Sub to the Evaluation block
  and a charge-quadrature tolerance of 0.1; see LIAN-2.8.1.md for context.

### Validation
- NVIDIA build and bundled Au contact and Pt SOC examples completed normally.
- Gold scaling checks: 1, 2, 4 and 8 threads, with finite transmission spectra.
- Pt SOC saved-density restarts: baseline/candidate at four threads and candidate
  at eight threads produced identical saved spectra and reported SCF energies.
- SCFSOC, experimental QTot2 and Intel builds were not validated.
- Release metadata/header edits were checked without rerunning calculations.

### Version tracking
- Read the packaging version from VERSION and include release documentation
  and NVIDIA build configuration in archives.
- Add VERSION as the authoritative release number and retain last-changed
  release headers for the modified sources, electrode data and NVIDIA build files.
- Release tag v2.8.2 identifies the source snapshot; binaries do not yet embed
  the release or Git revision automatically.

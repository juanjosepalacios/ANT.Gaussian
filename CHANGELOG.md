# Changelog

## Unreleased

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

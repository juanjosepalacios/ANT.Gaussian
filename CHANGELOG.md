# Changelog

## Unreleased

- Parallelize OneDLead bulk-charge quadrature with call-local Green-function
  and self-energy scratch and ordered accumulation; retain the separate
  device-charge fallback for 1D electrodes.
- Exclude complete 1D lead supercells from reduced-device charge counting.
- Use a consistent Fermi energy reference for bulk lead DOS output.
- Add NVIDIA/Intel 1D-charge regression tests and record Au10 validation
  and its remaining limitations in TOTCHARGE.md.

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

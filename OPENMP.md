# OpenMP tuning on Lian

Branch: perf/openmp-integration, based on f8ebfc7. Tested NVIDIA HPC SDK 23.11.

## Runtime configuration

Set Gaussian %NProcShared and OMP_NUM_THREADS to the same requested thread count.
In these runs, Gaussian's %NProcShared controlled the effective OpenMP maximum.
Changing OMP_NUM_THREADS alone while %NProcShared=1 did not increase the effective team size.
Keep MKL_NUM_THREADS=1 and OMP_DYNAMIC=FALSE to avoid unintended nested BLAS threading.
Use OMP_STACKSIZE=1G and ulimit -s unlimited for these test cases.
CPU affinity was not tuned; this hybrid CPU warrants a separate placement benchmark.

## Source changes

- Private chunk in four integration parallel regions; eliminates concurrent shared writes.
- Private nested k loop index in the two complex-contour integration regions.
- END DO NOWAIT before the four critical accumulator merges. Each thread merges its own
  completed accumulator; END PARALLEL still synchronizes before the result is consumed.
- Log effective OpenMP maximum and elapsed transport-analysis wall time.

## Gold contact results

Same 601-point Au4-Au-Au-Au4 example; Gaussian memory 1 GB.
One timed run per table entry: small differences between builds are timing noise,
not evidence of a source-level speedup. Transport timing includes final charge work.

| Build | Threads | Full run (s) | Transport stage (s) | Max absolute delta T vs serial baseline |
|---|---:|---:|---:|---:|
| baseline | 1 | 12.886 | not instrumented | 0 |
| baseline | 2 | 9.147 | not instrumented | 0.0001 |
| baseline | 4 | 7.145 | not instrumented | 0.0001 |
| baseline | 8 | 6.456 | not instrumented | 0.0002 |
| candidate | 1 | 13.035 | 3.0456 | 0 |
| candidate | 2 | 9.165 | 1.6555 | 0.0001 |
| candidate | 4 | 7.259 | 0.9120 | 0.0001 |
| candidate | 8 | 6.463 | 0.5535 | 0.0001 |

All runs terminated normally with 601 finite samples. Differences are measured
at the precision of the saved transmission file, not full internal precision.
At the same thread count, candidate and baseline gold spectra agree within 0.0001.

## Further optimization opportunities

The transmission loops hold a large CRITICAL region around DOS, eigenchannel analysis,
polarization and output. Move independent work into thread-local buffers and write
energy-ordered results after the loop in a subsequent change. Channel tracking
has cross-energy state and needs special handling before parallelizing that work.
Do not simply enable the charge-quadrature OpenMP directives for PGI: DDOS writes
the module-level ispin and requires a thread-safety refactor first.

Raw benchmark records and logs: /home/palacios/ANT-runs/openmp-integration/.
Existing installations and the prior fixed build are preserved. Nothing is pushed.

## Platinum SOC saved-density restart

28 atoms, 252 basis functions, 601 energies, 4 GB Gaussian memory.
All runs begin with the same saved density, converge in two SCF cycles,
and produce SCF energy -770.753121135 hartree.

| Build | Threads | Full run (s) | Max absolute delta T |
|---|---:|---:|---:|
| baseline | 4 | 82.561 | 0 |
| candidate | 4 | 82.496 | 0 |
| candidate | 8 | 52.308 | 0 |

All saved SOC transmission samples agree exactly at output precision.
Candidate transport-stage timings: 53.3128 s at four threads, 31.1716 s at eight.
Observed full-run reduction from four to eight threads: approximately 37%.
These are individual timings, not statistical confidence intervals.
The source edits alone have not shown measurable speedup at fixed thread count.

Recommended starting point for these examples on Lian: %NProcShared=8,
OMP_NUM_THREADS=8, MKL_NUM_THREADS=1, OMP_DYNAMIC=FALSE.
Do not extrapolate this to 16 or 32 threads without benchmarking.

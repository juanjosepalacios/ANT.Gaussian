Nanographene flake (170 C atoms) with a carbon Bethe-lattice electrode around its edge, spin-polarized
LSDA, and the local DOS of an atom at the centre of the flake.

The structure and inputs are those of `examples/graphene/nanographene.gjf` and `.ant`, with
`SL = 0.1` instead of `1.0d-2`. The first electrode is a Bethe lattice with the default carbon set
`Papacon` of `BLDAT/BL006.dat` (sp, in eV: es = -8, ep = 0, sss = -6, sps = 5, pps = 5,
ppp = -3); the second is a ghost (`TYPE2 = GHOST`), so there is no transmission. ANT attaches the
Bethe lattice at 38 sites all around the perimeter (the `Ar` atoms of `nanographene.xyz`).

## Running the example

Two steps in the same directory:

1. `nanographene.gjf` with `nanographene.ant`: the SCF from scratch (`guess=mix`, charge 0,
   multiplicity 3). It writes the density matrix `P.nanographene.dat` (25 MB, not stored here).
2. Copy `P.nanographene.dat` to `P.nanographene.ldos.dat`, then run `nanographene.ldos.gjf` with
   `nanographene.ldos.ant`. ANT starts from that density matrix; the route adds `nodamp` to
   `Scf=(...)`, and the `.ant` adds `LDOS` with `LDOS_BEG = LDOS_END = 81` and the default energy
   window (`EW1 = -3`, `EW2 = 3`, `ESTEP = 0.01` eV around E_F).

`.ant` settings of both steps: `ALPHA = 0.02`, `NPULAY = 4`, `CHARGEACC` and `FERMIACC` 1.0d-4,
`PACC = 1.0d-6`, `SELFACC = 1.0d-5`, `SL = 0.1`, `MAX = 10`, `NSPINLOCK = 0`, `OVERLAP = 0`,
`LEADDOS`. The `.gjf` files use `%NProcShared=4`; change it to the number of cores available. They
substitute links from `../../../../bin`.

## Outputs

- `nanographene.log.gz`, `nanographene.ldos.log.gz`: the logs of steps 1 and 2, compressed with
  gzip (10.7 MB each uncompressed, mostly the full Mulliken matrix printed by `pop=reg`). Step 2 was
  run under the job name `nanographene`, which appears in its log.
- `DOS.nanographene.ldos.dat`: total DOS and the DOS of atom 81 for E - E_F = -3 to +3 eV in 0.01 eV
  steps, in states/eV. Columns: E - E_F, total DOS, atom 81, and the s, p1, p2 and p3 orbitals of
  atom 81. Two blocks separated by a blank line, spin up and then spin down (spin down negative).
- `Lead1DOS.dat`: bulk DOS of the carbon Bethe lattice (`LEADDOS`), the same in both steps.
- `nanographene.xyz`: the geometry written by ANT; the `Ar` atoms mark the Bethe-lattice attachment
  sites.

## Results

ANT.Gaussian 3.0.0, Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026. Step 1 used
the `l502.exe` of commit 7709e49, step 2 that of commit a0e09f7; they differ only in how
Bethe-lattice parameter files are read, which does not affect this set.

| | Step 1 (from scratch) | Step 2 (restart, LDOS) |
|---|---|---|
| Basis functions | 680 | 680 |
| SCF energy (a.u.) | -950.603267 | -950.587449 |
| SCF cycles (device Fermi searches) | 23 (23) | 5 (6) |
| Bethe-lattice Fermi energy (eV) | 0.0000 | 0.0000 |
| Device Fermi energy (eV) | -7.1107 | -7.0650 |
| Largest atomic spin (ANT Mulliken) | 0.0047 | 0.0009 |
| Total spin of the flake | 0.021 | 0.0028 |
| Unsuccessful density integrations | 5 | 0 |
| Wall time | 889 s | 302 s |

LDOS of atom 81 (step 2), the same for both spins to within 1e-4 states/eV:

| | States/eV per spin |
|---|---|
| Total DOS at E_F | 1.850 |
| Atom 81 at E_F | 0.0051 |
| Minimum, at E - E_F = -0.29 eV | 0.0046 |
| Largest peak, at -2.80 eV | 0.111 |
| Peak above E_F, at +2.30 eV | 0.079 |

Both runs terminate normally. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log).

## Notes

- The flake is essentially non-magnetic: every atomic spin is below 0.005 in step 1 and 0.001 in
  step 2, largest on the zigzag edges and alternating between the two sublattices. The Bethe lattice
  attached all around the perimeter broadens the edge states that would carry edge magnetism in an
  isolated flake. The `<S**2>` of about 17 printed by Gaussian refers to its Kohn-Sham orbitals, not
  to ANT's density matrix; Gaussian's own Mulliken spins are also below 0.004 per atom in step 2.
- Atom 81 is one of the six atoms around the hexagon at the centre of the flake. Its DOS near E_F is
  all in the pi orbital, which appears as p1 (px) because the route has no `nosymm`: Gaussian moves
  the flake to its standard orientation, perpendicular to X, and ANT works in that frame.
- The restart does not end where step 1 did: its energy is 15.8 mHa higher, E_F 46 meV higher and the
  largest atomic spin smaller (0.0009 against 0.0047). Both states are nearly non-magnetic.

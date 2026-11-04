Five cases for Ni nanocontacts with Bethe-lattice electrodes and spin-polarized LSDA:

1) A small contact (Ni4-Ni-Ni-Ni4, 10 atoms) in a Zeeman field of 20 T along z (`Ni4-Ni-Ni-Ni4`)
2) A larger contact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) with parallel (ferromagnetic) electrodes (`.F`)
3) The same contact with antiparallel electrodes, a domain wall in the constriction (`.AF`),
   started from the density matrix of case 2
4) The same contact with an all-electron basis on the central atoms 10-19 (`.doll`)
5) Case 1 with the spins of the two central atoms rotated to the x direction (`.rotx`)

Each case has its own `.gjf`, `.ant` and `.xyz`, so its outputs carry the job name and don't overwrite
each other.

## Running the cases

Run cases 1, 2, 4 and 5 in any order. Case 3 must run after case 2:

1. Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.gjf`. It writes `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat`.
2. Copy it to `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.dat`, then run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.gjf`.

ANT reads that file at start-up as the initial density matrix. `MRSTART = 15` reverses the spin of
atoms 15-28 in it, and `UD` makes the right Bethe-lattice electrode spin-down. Without the copy, case 3
does not start from the converged parallel state. `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat` from the run
reported below is included, so case 3 can be run on its own. Because of that file, running case 2 in
this folder restarts from the converged density; delete or move it to repeat case 2 from scratch.

| Case | Basis | Electrodes (`BLPAR`) | `.ant` settings |
|---|---|---|---|
| Ni4 | Minimal sd + ECP on all atoms | `Papacon-sd` | `ZM`, `Bz = 20.0`; `PACC` 1.0d-8, `CHARGEACC`/`FERMIACC` 1.0d-5 |
| F | Minimal spd + ECP on all atoms | `Papacon` (default) | `MAX = 10`; `PACC` 1.0d-7, `CHARGEACC`/`FERMIACC` 1.0d-5 |
| AF | As F | `Papacon`, right electrode reversed (`UD`) | `MRSTART = 15`, `LEADDOS`; `PACC` 1.0d-8 |
| doll | sd + ECP on atoms 1-9 and 20-28; all-electron on 10-19 | `Papacon-sd` | `NSPINLOCK = 5`, `NIP = 7`; `PACC` 1.0d-6, `CHARGEACC`/`FERMIACC` 1.0d-4 |
| rotx | As Ni4 | `Papacon-sd` | Ni4 settings plus `SPINROTATOM` (atoms 5 and 6 at theta = 90, phi = 0) |

The Zeeman field is given in tesla (`BX`, `BY`, `BZ`), with E = muB * B * sigma and
muB = 5.788e-5 eV/T, so 20 T splits the spin levels by 2.3 meV. The commented `ZMATOM`, `SPINROTATOM` and
`SPINEDIT` blocks in the Ni4 `.ant` show per-atom fields, spin rotations and spin edits. Spin
rotations use `SPINROTATOM` (see `doc/keywords.txt`), as in case 5:

```
SPINROTATOM
2
5, 90.0, 0.0
6, 90.0, 0.0
```

Each line gives an atom, theta and phi in degrees; the majority spin of that atom points along
(sin theta cos phi, sin theta sin phi, cos theta), here +x. Atoms not listed stay along z. The
rotation is applied once, after self-consistency: ANT converges the collinear solution, rebuilds the
Hamiltonian in the spinor basis with the rotated blocks, finds the device Fermi energy again and then
computes the populations and transmission.

**Warning:** `SPINROTATOM` rotates the spins but is not self-consistent. The rotated moments don't relax
after the rotation, so case 5 is a non-self-consistent estimate of a non-collinear spin texture. A
true spin flip is better done with `SPINEDIT`, which reverses the spin of selected atoms in the
starting density matrix and lets the SCF converge the new texture, as `MRSTART` does in case 3.

The `doll` case was set up with spin-orbit coupling (`SOC`, `SOCFACATOM`). With SOC and the default
`NIP = 6`, its density integration failed 249 times at the 127-point limit in the first 18 ANT
cycles. The run below therefore uses `NIP = 7` (up to 255 points) without SOC. The SOC lines remain
in the `.ant`, commented out.

Outputs of each case: `<job>.log`, `T.<job>.dat` (`T.<job>.SOC.dat` for case 5) and, for case 3, `Lead1DOS.dat` and `Lead2DOS.dat`
(electrode DOS from `LEADDOS`). The transmission covers E - E_F = -3 to +3 eV in 0.01 eV steps. For
cases 2-4 the file has two blocks separated by a blank line, spin up and then spin down, each in units
of e^2/h. Cases 1 and 5 are computed in the spinor basis because of the Zeeman term and the
rotation, so their files have one block with the total transmission in units of e^2/h. In all cases G = T_total * e^2/h = (T_total / 2) G0.

## Results

ANT.Gaussian 2.8.4, Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads. Cases 1-4: commit 1e4a005
sources, 3 November 2026. Case 5: commit 1442c1d plus the removal of `SPINROT`, `THETA` and `PHI`, so that
`SPINROTATOM` alone switches on the rotation, 4 November 2026.

| | Ni4 (20 T) | F | AF | doll | rotx |
|---|---|---|---|---|---|
| Atoms / basis functions | 10 / 60 | 28 / 252 | 28 / 252 | 28 / 418 | 10 / 60 |
| SCF energy (a.u.) | -402.179706 | -1126.758239 | -1126.768283 | -15786.572053 | -402.179887 |
| SCF cycles (device Fermi searches) | 37 (26) | 111 (99) | 19 (21) | 41 (35) | 37 (26) |
| Device Fermi energy (eV) | -2.6383 | -4.3635 | -4.2798 | -6.4278 | -2.6671 |
| T(E_F) up / down | spinor | 0.837 / 2.864 | 1.202 / 1.227 | 0.871 / 2.548 | spinor |
| T(E_F) total (e^2/h) | 2.101 | 3.700 | 2.429 | 3.418 | 2.706 |
| Conductance G (G0) | 1.05 | 1.85 | 1.21 | 1.71 | 1.35 |
| Spin polarization of T(E_F), (up - down)/(up + down) | | -0.55 | -0.01 | -0.49 | |
| T total at E_F - 1 eV / E_F + 1 eV | 2.512 / 0.389 | 2.419 / 3.880 | 0.900 / 2.938 | 4.365 / 3.478 | 2.282 / 0.380 |
| Spin of the two central atoms | 1.15 / 1.15 | 0.87 / 0.85 | 0.87 / -0.86 | 0.37 / 0.37 | <Sx> 1.02, <Sz> 0.21 each |
| Electrons in left / right half | 50.00 / 50.00 | 139.73 / 140.27 | 139.76 / 140.24 | 230.00 / 230.00 | 50.00 / 50.00 |
| Spin of left / right half | 5.77 / 5.77 | 12.96 / 12.50 | 13.16 / -12.78 | 11.47 / 11.47 | 5.71 / 5.71 |
| Unsuccessful density integrations (point limit) | 45 (127) | 0 | 88 (127) | 10 (255) | 45 (127) |
| Wall time | 15 s | 1935 s | 376 s | 4375 s | 14 s |

All runs terminate normally, pass the spectral-weight check without widening the energy window, and
report no Fermi-search errors.

Charges and spins come from ANT's Mulliken population analysis (`--- Mulliken population analysis ---`
in the log), not from Gaussian's population analysis at the end of the log. ANT splits the cluster
into two halves, one per electrode: atoms 1-5 / 6-10 for Ni4 and 1-14 / 15-28 for the larger contact.
The central atoms are 5 and 6 for Ni4 and 14 and 15 for the larger contact. Every case is neutral:
100, 280, 280, 460 and 100 electrons in total. Cases 1 and 5 are computed in the spinor basis; the
spins of case 1 are <Sz>, and case 5 gives the components <Sx>, <Sy>, <Sz> of each atom.

## Notes

- The antiparallel configuration is lower in energy than the parallel one by 0.0100 a.u. (0.27 eV)
  with these electrode parameters. Its moments are nearly antisymmetric across the contact (+0.87
  and -0.86 on atoms 14 and 15), so the domain wall sits between the two central atoms.
- Parallel electrodes give 1.85 G0, carried mostly by the spin-down channel, which has partly filled
  d states at E_F (T down / T up = 3.4). In the antiparallel configuration both spin channels carry
  the same transmission (1.20 and 1.23) and the conductance drops to 1.21 G0. The magnetoresistance
  (G_F - G_AF) / G_AF is 52%.
- Case 2 needs many more cycles (111) than case 3 (19), which starts from the converged parallel
  density with the spins of one half reversed.
- The all-electron basis on the central atoms (case 4) changes the conductance from 1.85 to 1.71 G0
  and keeps the spin-down dominance (polarization -0.49 against -0.55). The device Fermi energy and
  SCF energy are not comparable with cases 2-3, since the basis, ECPs and Bethe-lattice parameters
  differ.
- With the all-electron basis the moments of the central atoms are smaller: 0.58 on atoms 10-13 and
  16-19 and 0.37 on atoms 14 and 15, against 0.72-0.87 in case 2. Most outer atoms keep about
  1.0 each, so each half carries a spin of 11.47 against 12.5-13.0 in case 2.
- Cases 1 and 3 report unsuccessful density integrations at the 127-point limit (`NIP = 6`). Both
  converge normally and pass the spectral-weight check; `NIP = 7` would allow 255 points at a higher
  cost, as in case 4.
- At 20 T the Zeeman splitting (2.3 meV) is small compared with the exchange splitting of Ni, so case
  1 mainly fixes the magnetization along z. The original input used `Bz = 20000.0`, a field of 20 kT.
- Case 5 shows that the rotation works. Atoms 5 and 6 end up with <Sx> = 1.02 and <Sz> = 0.21, with
  no y component, while the other atoms stay along z (<Sz> = 1.17, <Sx> = 0.05). The moments of the
  rotated atoms keep a small z component and shrink from 1.15 to 1.04, and they take up 0.10 e,
  because the rotation is not followed by self-consistency. Restarting from one converged density and
  changing only theta of atoms 5 and 6 gives:

  | theta | <Sx> / <Sz> of atoms 5, 6 | Electrons on atoms 5, 6 | Device E_F (eV) | T(E_F) (e^2/h) |
  |---|---|---|---|---|
  | 0 | 0.00 / 1.15 | 10.10 | -2.638 | 2.101 |
  | 90 | 1.02 / 0.21 | 10.19 | -2.666 | 2.708 |
  | 180 | 0.00 / -0.83 | 10.34 | -2.707 | 2.972 |

  theta = 0 reproduces case 1 exactly, and theta = 180 reverses the two moments. The rise in
  conductance with the rotation angle comes from a non-self-consistent spin texture and should not
  be read as the conductance of a domain wall.
- Repeated runs of the collinear Ni4 SCF, which cases 1 and 5 share, ended between -402.17969 and
  -402.18001 a.u. in five runs, since the 45 unsuccessful integrations leave the convergence loose.
  Across four runs of case 5, T(E_F) agrees to 0.001 and the device Fermi energy to 3 meV.
- The `.gjf` files substitute links from `../../bin`. `Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.gjf` asks for
  `%nproc=16`; lower it to the number of cores available.

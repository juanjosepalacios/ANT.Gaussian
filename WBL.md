# WBL: transmission in the wide-band limit of the molecule's embedding

## What it is for

The Seebeck coefficient follows from the slope of T(E) at E_F. In a metal /
molecule / metal junction that slope has two sources. One is the molecular
levels. The other is the electrodes themselves: both the Bethe-lattice
self-energies and the metal atoms of the cluster have an energy-dependent DOS
near E_F. WBL separates the two. In the same run, from the same converged H
and S, it writes an extra transmission in which everything outside the
molecule is frozen at E_F. The full matrix structure of the coupling is kept,
so the approximation is a wide-band limit in energy only, and no Mott-type
additivity is assumed. The difference between S from T and S from T_WBL is
the electrode contribution.

## Usage

Add the keyword to the `.ant` file:

    WBL

ANT then writes `T.<job>.WBL.dat` next to the normal `T.<job>.dat`, which
does not change. The WBL file has the format of the T file:

- E - E_F (eV);
- T, summed over spin when NSpin = 1, or one block per spin separated by a
  blank line when NSpin = 2;
- the NCHANNELS largest eigenvalues of the transmission matrix, sorted at each
  energy and not followed across energies.

Compute S from it exactly as from T.

The molecule is ANT's own partition of the cluster, the one printed in the log
under "Atoms in left electrode" and "Atoms in right electrode": the atoms
between the two electrodes. Set NATOMEL1 and NATOMEL2 to move the boundary,
for example to count anchor atoms as part of the electrode. If either
electrode is empty, or there is no molecule, WBL is skipped with a warning.

Not implemented for SOC, SCFSOC, SPINROTATOM or ZM. With any of those, the log
gets a warning and no WBL file is written.

## Formulation

ANT's device Green function at z = E - E_F is

    G(z) = [ A(z) - Sigma_L(z) - Sigma_R(z) ]^-1,
    A(z) = (z - shift + i eta) S - H,

with the Bethe-lattice self-energies (times GLUE, in the non-orthogonal
basis). The transmission is T = Tr[Gamma_L G^a Gamma_R G^r], or its HTRANSM
form.

**The basis.** In the non-orthogonal AO basis the energy enters through z S,
including the molecule-electrode block z S_ME, so a "molecule block" does not
separate the molecule from the electrodes. WBL is therefore defined after
Loewdin orthogonalization, X' = S^-1/2 X S^-1/2. There

    A'(z) = (z - shift + i eta) I - H',

the energy enters as z times the identity, and the molecule-electrode
coupling H'_ME does not depend on energy at all. The molecule M is the set of
Loewdin orbitals of the molecule atoms.

**The approximation.** With GInv0' = A'(0) - Sigma'_L(0) - Sigma'_R(0), the
inverse at E_F, WBL keeps the energy only on the molecule orbitals:

    G'_WBL(z) = [ GInv0' + z P_M ]^-1,
    T_WBL(z)  = Tr[ Gamma'_L(0) G'_WBL^a Gamma'_R(0) G'_WBL^r ].

By the Schur complement, the molecule block is

    G'_WBL,MM(z) = [ z - shift + i eta - H'_MM - Sigma_emb(0) ]^-1,
    Sigma_emb(0) = Sigma'_MM(0) + V_ME g'_E(0) V_EM,
    V_ME = H'_ME + Sigma'_ME(0),
    g'_E(0) = [ A'_EE(0) - Sigma'_EE(0) ]^-1,

with Sigma' = Sigma'_L + Sigma'_R. Sigma' acts on the electrode orbitals; its
small parts on M come from the Loewdin tails. The molecule keeps its own
Hamiltonian, and its embedding self-energy (the cluster electrode atoms and
the Bethe lattice together) is frozen at E_F. That is the wide-band limit, and
the only approximation.

**Checks and properties.**

- T does not depend on the basis, so at E_F T_WBL equals the normal T exactly.
  The log checks this to rounding (`WBL check` lines).
- T_WBL depends only on the orthogonal model, not on how much overlap the AO
  basis carries.
- At every energy T_WBL is the transmission of a physical system (a Hermitian
  Hamiltonian with energy-independent self-energies), so it never exceeds the
  number of open channels.

**Why not the AO basis.** A first version used [GInv0 + z S_MM]^-1 in the AO
basis. That freezes the overlap part z S_ME of the coupling as well. With
diffuse metal basis functions and short contacts (LANL2DZ, C-Au near 2.3 A),
S_ME is large. Compressed helicene junctions then showed series of
unphysically sharp peaks above E_F that did not move with UPLUS, while the
normal T and the DOS were smooth. The Loewdin form removes that freezing.

**Why the full matrices and not a molecule-only formula.** One could write
T_WBL on the molecule alone, as Tr[Gamma_eff,L G_MM^a Gamma_eff,R G_MM^r]
with Gamma_eff,X = V_ME g'_E Gamma'_X g'_E^+ V_ME^+. That form is exact only
when the left and right electrodes do not couple to each other directly. With
diffuse metal basis functions across a molecule a few Angstrom long, they do.
In one tested Au junction the molecule-only form overestimated T(E_F) by more
than an order of magnitude. Keeping the full device matrices includes the
direct electrode-electrode path exactly, also frozen at E_F. To show how large
that path is, the log prints three numbers:

- the molecule-only value;
- the direct transmission with the molecule removed,
  Tr[Gamma'_L g'_E^a Gamma'_R g'_E^r];
- the largest |H| and |S| between the two electrodes, in the AO basis.

## Range of validity

The approximation is the wide-band limit itself: the embedding self-energy is
taken at E_F. Use T_WBL for the slope at E_F (the Seebeck coefficient), not
for the shape of T(E) over eV, where the real embedding changes.

## Cost

S^-1/2 is computed once, and GInv0, Gamma_L and Gamma_R are transformed once
per spin. Each energy point then costs one NAOrbs inversion, three or four
ZGEMMs and, with NCHANNELS, one diagonalization. That is close to a normal
point without the Bethe-lattice solve. WBLTransmission runs after
transmission, with its own OpenMP energy loop and no critical section. The
SCF, the density and every other output are untouched.

## Limits

- Collinear runs only.
- NSpin = 2 is implemented per spin but not yet checked on a production run.
- Loewdin orbitals on the molecule atoms have small tails on the electrode
  atoms; that is the price of an orthogonal partition.

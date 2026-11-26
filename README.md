# Code and data for Paper 1: NS-MFP for a sharp-interface solver

This repository contains the code and the data used in the paper *A Jacobian-free global stability method for sharp-interface multiphase solvers: extending NS-MFP to capillary-driven flows* (A. Okojunu, M. Sussman and U. S. Nair). Every table and every data figure of the paper can be recomputed from it, and every simulation can be repeated with it.

The repository depends on one other code, the two-phase solver `amrex_implicit_interfaces` of M. Sussman, at a fixed commit. The analysis code is held in the folder `LSA_Rayleigh_Fortran/`, which is placed inside the solver tree to run the simulations. Levels 1 and 2 below need neither the solver nor a cluster.

## Contents

| Path | Contents |
|---|---|
| `LSA_Rayleigh_Fortran/src`, `tests`, `Makefile` | the analysis code (Fortran 2008, LAPACK): the twin differences, the interface-amplitude fit, the dynamic mode decomposition (DMD), the dispersion relations and the input decks of the solver, with its test suite |
| `LSA_Rayleigh_Fortran/scripts` | the scripts that run the simulations under Slurm and analyse them |
| `LSA_Rayleigh_Fortran/plotting` | the Python scripts that print the tables and draw the figures; `plotting/paper1/` holds the superposition test and the two schematic figures |
| `LSA_Rayleigh_Fortran/results` | the stored output of every run used in the paper: amplitude records, decomposition outputs, studies and run records |
| `LSA_Rayleigh_Fortran/*.txt` | the case files, one line per simulation |
| `LSA_Rayleigh_Fortran/reproduce_paper1.sh` | the commands of the three levels below |
| `solver_patches/` | the change to the solver needed by the superposition test (Level 3 only) |

## Level 1. Check the code (about one minute)

On Linux, or WSL on Windows, with gfortran, make and LAPACK:

```bash
git clone https://github.com/abelokoj/lsa-nsmfp-sharp-interface.git
cd lsa-nsmfp-sharp-interface/LSA_Rayleigh_Fortran
./reproduce_paper1.sh 1
```

The test suite must end with `ALL TEST PROGRAMS PASSED`: 367 checks in eight programs (Section 5.1 of the paper). The script then prints the growth rates of Eq. (17) and of the confined two-fluid relation with the wall at 4 R0, at the wavenumbers of Table 3; at k R0 = 0.7 they are 0.329367 and 0.328617.

## Level 2. Recompute the tables and figures from the stored data (about one minute)

```bash
python3 -m pip install -r plotting/requirements.txt
./reproduce_paper1.sh 2
```

The tables are written to `figures/paper1_tables.txt` and the figures to `figures/`. The output is labeled with the table numbers of the dissertation on which the paper is based. The correspondence is as follows.

| In the paper | Printed or drawn by Level 2 as |
|---|---|
| Table 2 and Section 5.4, superposition | "Superposition of two perturbations" |
| Table 3, growth rates of the 16-cell sweep | Table 5.3, and Table 6.1 for the column of the confined relation |
| Table 4, the three grids | Table 5.4 |
| Table 5, sampling stride | Table 5.6 |
| Section 5.7, the stable wavenumber | Table 5.5 |
| Section 5.2, the base state | "drift of runs/sw070_16.base" |
| Section 5.3, proportionality on 8, 16 and 32 cells | `linearity_08`, `_16`, `_32` and the field comparison |
| Section 5.10, the mode shapes | "rms with wall" for each wavenumber |
| Fig. 1 and Fig. 3 | `figures/nsmfp_variants`, `figures/ch3_column_schematic` |
| Fig. 4 | `figures/ch5_proportionality` |
| Fig. 5 | `figures/campaign_dispersion` |
| Figs. 6 and 7 | `figures/ch5_dmd_evidence`, `figures/ch5_mode_shapes` |

Fig. 2 is a drawing of Algorithm 1 and is not computed from data. Level 2 also draws some figures of the dissertation that are not used in the paper.

Level 2 refits every interface growth rate from its amplitude record and recomputes the tables. It does not repeat the decompositions: these are read from their stored outputs, because their inputs, the twin-difference fields, are too large for this repository (about 1 GB). They are repeated in Level 3.

For the superposition test, Level 2 refits the four interface growth rates of Section 5.4 from their records. The departures of Table 2 and the decomposition of each harmonic need the six twin-difference fields of these runs (215 MB), which are archived on Zenodo [TODO: add the DOI of the archive]. Without them, Level 2 prints their stored output. With them, it recomputes both:

```bash
tar xf p1_superposition_exports.tar -C results/superposition_2026-11-25
./reproduce_paper1.sh 2
```

## Level 3. Run the simulations again

The runs were made on the FSU Research Computing Center cluster with `gnu/13` and `openmpi/4.1.6`, under Slurm. They used commit `dddb827c` of a fork of the solver, whose solver source is identical to commit `0b8df42` of the original repository; the run records in `results/*/run_records/` give the commit and the checksum of the executable for every run. The superposition runs used the same source with the patch in `solver_patches/`.

```bash
export AII=$HOME/amrex_implicit_interfaces
git clone https://github.com/msussman42/amrex_implicit_interfaces.git $AII
cd $AII && git checkout 0b8df42
git apply /path/to/lsa-nsmfp-sharp-interface/solver_patches/*.patch
cp -r /path/to/lsa-nsmfp-sharp-interface/LSA_Rayleigh_Fortran $AII/NS_fluids_lib/
module purge; module load gnu/13; module load openmpi
cd $AII/NS_fluids_lib && make -j16              # writes amr2d.gnu.MPI.ex
cd LSA_Rayleigh_Fortran
ACCOUNT=<your Slurm account> ./reproduce_paper1.sh 3
```

The patch changes only the initial interface in `NS_fluids_lib/MITSUHIRO_PIPE.F90`, so that a sum of axial harmonics can be imposed, and it applies to commit `0b8df42` without conflict. It is used only for the superposition runs; with `ns.radblob3 = 0`, the default of every other deck, the patched solver computes what the unpatched one does. Use commit `0b8df42` and not the current head: at some later commits a calculation with fluids only stops at start-up with `expecting num_FSI_outer_sweeps>=2` (`LSA_Rayleigh_Fortran/docs/solver_notes.md`).

Level 3 submits six case files as Slurm job arrays. For the first five, the number of ranks and the wall-clock limit of each task are those of the original runs. The number of ranks of the superposition runs is not recorded in this repository, and they are run with the default of `scripts/submit_campaign.sh`:

| Case file | Runs | Ranks | Limit per task | In the paper |
|---|---|---|---|---|
| `campaign_16_wide.txt` | sweep at k R0 = 0.2, 0.3, 0.4 | 32 | 13 h | Tables 3 and 5 |
| `campaign_16.txt` | sweep at k R0 = 0.5 to 0.9; ε = 5 × 10⁻³ and 10⁻²; time step halved; wall at 6 R0; gas density and viscosity divided by ten; k R0 = 1.1 | 16 | 6 h | Tables 3 and 5, Sections 5.3 and 5.5 to 5.7 |
| `rerun_meas.txt` | k R0 = 0.7 on 8, 16 and 32 cells | 16 | 30 h | Table 4, Figs. 6 and 7 |
| `campaign_amp08.txt`, `campaign_amp32.txt` | three amplitudes on 8 and on 32 cells | 8, 16 | 4 h, 36 h | Section 5.3 |
| `campaign_p1_superposition.txt` | the three superposition pairs | 16 | 12 h | Table 2 |

Each perturbed run has an unperturbed twin, and the two must run on the same number of ranks; the case files and the scripts ensure this. A task that reaches its limit is not lost: submitting the same command again skips the finished cases and restarts the others from their last checkpoint. When every job has ended, `./reproduce_paper1.sh 3-analyse` analyses the new runs and compares them with the stored ones. Bitwise reproducibility has not been tested. Another compiler, MPI library or rank count changes the round-off, so the interface rates are expected to agree closely with the stored ones, and the decomposition rates, which come from the difference of two runs, less well.

## Origin

The analysis code is that of the dissertation of A. Okojunu (Florida State University, 2026), published at https://github.com/abelokoj/LSA_Rayleigh_Fortran, with the superposition test added for this paper. The dissertation results that the paper does not use (the Chapter 6 maps, the restarting variant and the September 2026 runs) are not included in this repository.

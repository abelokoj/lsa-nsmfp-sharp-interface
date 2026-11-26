#!/usr/bin/env bash
# Reproduces the results of Paper 1 (the NS-MFP method for a sharp-interface solver) at three levels.
#   ./reproduce_paper1.sh 1          build the analysis code, run its tests and print the theoretical growth rates
#   ./reproduce_paper1.sh 2          recompute every table and figure of the paper from the stored data (no solver)
#   ./reproduce_paper1.sh 3          submit the simulations (this folder must be inside amrex_implicit_interfaces/NS_fluids_lib)
#   ./reproduce_paper1.sh 3-analyse  analyse the new runs once every job has ended, and compare them with the stored data
# Run it from this folder. Level 2 writes its figures to figures/ and its tables to figures/paper1_tables.txt.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
SNAP=results/followup_2026-10-22
SP=results/superposition_2026-11-25

level1 () {
  make -j8
  make test | tail -12
  # Eq. (17) (Chandrasekhar) and the confined two-fluid relation with the wall at 4 R0, at the wavenumbers of Table 3
  ./bin/run_lsa theory --kr0 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 --r-wall 4
}

level2 () {
  [[ -x bin/run_lsa ]] || level1
  "$PY" -c "import numpy, scipy, matplotlib" 2>/dev/null || { echo "install the Python packages first: $PY -m pip install -r plotting/requirements.txt"; exit 1; }
  mkdir -p figures
  {
    echo "== Interface growth rates refitted from the amplitude records (Section 5.5)"
    scripts/refit_records.sh | tail -3
    echo
    echo "== Tables 3, 4 and 5 and the stable mode (dissertation Tables 5.3 to 5.6 and 6.1; see README.md for the mapping)"
    "$PY" plotting/thesis/print_thesis_tables.py
    echo
    echo "== The base state (Section 5.2): drift of the unperturbed twin at k R0 = 0.7"
    cat $SNAP/analysis/studies/drift_sw070_16.txt
    echo
    echo "== Proportionality on 8, 16 and 32 cells per radius (Section 5.3)"
    for g in 08 16 32; do echo "-- $g cells"; cat $SNAP/analysis/linearity_$g.txt; done
    echo "-- the scaled fields compared in direction:"; sed -n '1,60p' results/linearity_fields_2026-10-23/README.md
    echo
    echo "== Mode shapes against the confined and unbounded eigenfunctions (Section 5.10)"
    "$PY" plotting/check_mode_shapes_confined.py $SNAP/analysis
    echo
    echo "== Superposition of two perturbations (Section 5.4, Table 2)"
    superposition
  } 2>&1 | tee figures/paper1_tables.txt
  "$PY" plotting/plot_campaign_results.py results/campaign_2026-10-09 --outdir figures                    # Fig. 5
  "$PY" plotting/thesis/make_extra_figures.py results/campaign_2026-10-09 \
        results/estimators_2026-10-22/weighting_test_output_2026-10-22.txt plotting figures              # Fig. 4
  "$PY" plotting/thesis/make_dmd_evidence_figures.py $SNAP/analysis plotting figures                      # Figs. 6 and 7
  "$PY" plotting/paper1/plot_schematics.py                                                                # Figs. 1 and 3
  echo "Tables: figures/paper1_tables.txt. Figures: figures/ (see README.md for the file of each figure)."
}

superposition () {
  local a=$SP/analysis h=$SP/records/spsum_16_harmonics tmax
  echo "-- interface rates of A (k R0 = 0.35) and B (k R0 = 0.7) imposed alone, refitted from their records"
  scripts/refit_records.sh $SP | grep -E "^sp0|^case"
  # the window of C closes when the sum of the two harmonic amplitudes first exceeds 0.1 R0 (as in scripts/analyse_broadband.sh)
  tmax=$(paste $h/h01.dat $h/h02.dat | awk '!/^#/ {s=0; for (i=2; i<=NF; i+=6) if ($i>0) s+=$i; if (s>0.1) {print t; exit} t=$1}')
  echo "-- interface rates of the two harmonics of C, up to t = $tmax"
  ./bin/run_lsa record --record $h/h01.dat --kr0 0.35 --t-max "$tmax" | grep -E "window|growth rate"
  ./bin/run_lsa record --record $h/h02.dat --kr0 0.7 --t-max "$tmax" | grep -E "window|growth rate"
  if [[ -f $a/spsum_16/export/perturbation.bin ]]; then
    echo "-- departure from superposition of the fields (Table 2)"
    "$PY" plotting/paper1/superposition_check.py $a/sp035_16 $a/sp070_16 $a/spsum_16 figures/superposition.csv
    echo "-- the decomposition of each harmonic of C"
    mkdir -p figures/superposition
    # the level-set export of C in spsum_16/broadband/ was identical to spsum_16/export_ls/, so only the latter is archived
    "$PY" plotting/broadband_spectrum.py --archive $a/spsum_16/broadband/dmd_ls --records $a/spsum_16/broadband/records \
        --export $a/spsum_16/export_ls --t-max "$tmax" --zlen 17.951958 --nharm 2 --outdir figures/superposition
    cat figures/superposition/spectrum.csv
  else
    echo "-- Table 2 and the decomposition of each harmonic need the twin-difference fields (README.md, Level 2); their stored output is:"
    cat $SP/superposition_summary.txt
    cat $a/spsum_16/broadband/spectrum.csv
  fi
}

level3 () {
  [[ -x ../amr2d.gnu.MPI.ex ]] || { echo "build the patched solver first: ../amr2d.gnu.MPI.ex not found (see README.md, Level 3)"; exit 1; }
  make -j8 && make test | tail -1
  ACCOUNT="${ACCOUNT:?set ACCOUNT to your Slurm account}" scripts/submit_campaign.sh \
      campaign_16_wide.txt campaign_16.txt rerun_meas.txt campaign_amp08.txt campaign_amp32.txt campaign_p1_superposition.txt
  echo "When every job has ended (squeue --me), run: ./reproduce_paper1.sh 3-analyse"
}

level3_analyse () {
  local acc="${ACCOUNT:?set ACCOUNT to your Slurm account}"
  sbatch -A "$acc" -n 1 -c 4 --mem=16G -t 04:00:00 -o runs/array/studies-%j.out scripts/run_studies.sh all
  sbatch -A "$acc" -n 1 -c 2 --mem=8G -t 02:00:00 -o runs/array/linearity-%j.out scripts/linearity_08_32.sh 08 32
  sbatch -A "$acc" -n 1 -c 4 --mem=32G -t 04:00:00 -o runs/array/linfields-%j.out scripts/linearity_fields.sh
  "$PY" plotting/paper1/superposition_check.py runs/analysis/sp035_16 runs/analysis/sp070_16 runs/analysis/spsum_16 runs/analysis/superposition.csv
  "$PY" plotting/thesis/compare_summaries.py runs/analysis/summary.csv
  scripts/refit_records.sh runs
  echo "When the three batch jobs have ended: $PY plotting/thesis/print_thesis_tables.py --snapshot runs"
}

case "${1:-}" in
  1) level1 ;;
  2) level2 ;;
  3) level3 ;;
  3-analyse) level3_analyse ;;
  *) sed -n '2,8p' "$0"; exit 1 ;;
esac

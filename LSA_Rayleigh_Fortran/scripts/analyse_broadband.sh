#!/bin/bash
# analyse_broadband.sh -- several growth rates from one broadband pair of runs (campaign_broadband.txt).
#
# For a run whose deck has ns.radblob3 = M > 0 (M axial harmonics of the box length L = ns.yblob), and its twin <label>.base:
#   1. interface records: run_all.sh has written runs/records/<label>_harmonics/hNN.dat, the amplitude of harmonic n; each unstable harmonic (k_n R0 = 2 pi n R0 / L < 1) is fitted with run_lsa record at k r0 = 2 pi n / L, up to the time at which the total amplitude sum_n a_n first exceeds 0.1 R0 (beyond it the harmonics interact);
#   2. DMD: the level-set twin difference is exported (run_lsa export --fields L0101) and decomposed in the band |r - r0| < 2 dr over the same time window, keeping all modes in the archive;
#   3. plotting/broadband_spectrum.py identifies the axial harmonic of each DMD mode and tabulates DMD rate, interface rate and theory (Chandrasekhar, two-fluid viscous) per harmonic, with a figure.
#
# Usage (from $LSA): scripts/analyse_broadband.sh bb72_08 [more labels]
#   sbatch -A genacc_q -n 1 -c 4 --mem=32G -t 04:00:00 -o runs/array/broadband-%j.out scripts/analyse_broadband.sh bb72_08
# Outputs: runs/analysis/<label>/broadband/ (records/, export_ls/, dmd_ls/, spectrum.csv, figures). FORCE=1 redoes every step.

set -u
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RUNROOT="${RUNROOT:-runs}"
OUT="${RUNROOT}/analysis"
FORCE="${FORCE:-0}"
need () { [[ "$FORCE" == 1 || ! -e "$1" ]]; }
[[ -x bin/run_lsa ]] || { echo "bin/run_lsa missing; run make first" >&2; exit 1; }
[[ $# -ge 1 ]] || { echo "usage: $0 LABEL [LABEL...]" >&2; exit 1; }

deckval () { awk -v k="$1" '$1==k {print $3; exit}' "$2"; }
props_of () {
  awk '/^ns\.denconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],d," "); rg=d[2]}
       /^ns\.viscconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],v," "); ml=v[1]; mg=v[2]}
       END{if (ml!="") printf "--mu-l %s ", ml; if (rg!="") printf "--rho-g %s ", rg; if (mg!="") printf "--mu-g %s", mg}' "$1"
}

for label in "$@"; do
  deck="$RUNROOT/$label/inputs"
  [[ -f "$deck" ]] || { echo "$label: no deck $deck" >&2; continue; }
  M=$(deckval ns.radblob3 "$deck"); M=${M%%.*}
  L=$(deckval ns.yblob "$deck")
  [[ -n "$M" && "$M" -ge 1 ]] || { echo "$label: not a broadband run (ns.radblob3 = ${M:-unset})" >&2; continue; }
  hdir="$RUNROOT/records/${label}_harmonics"
  dest="$OUT/$label/broadband"
  mkdir -p "$dest/records"
  # shellcheck disable=SC2046
  props=$(props_of "$deck")
  echo "== $label: $M harmonics of L = $L"

  # --- 1. interface records per harmonic, up to the end of the linear range of the sum ---
  tmax=$(paste $(for n in $(seq -f %02g 1 "$M"); do echo "$hdir/h$n.dat"; done) 2>/dev/null |
    awk '!/^#/ {s=0; for (i=2; i<=NF; i+=6) if ($i>0) s+=$i; if (s>0.1) {print t; exit} t=$1} END{if (t!="") print t}' | head -1)
  echo "  linear range of the sum ends at t = ${tmax:-end of run}"
  for n in $(seq 1 "$M"); do
    nn=$(printf %02d "$n")
    k=$(awk -v n="$n" -v L="$L" 'BEGIN{printf "%.10f", 2*3.141592653589793*n/L}')
    unstable=$(awk -v k="$k" 'BEGIN{print (k<1)?1:0}')
    [[ "$unstable" == 1 && -f "$hdir/h$nn.dat" ]] || continue
    if need "$dest/records/h$nn.csv"; then
      ./bin/run_lsa record --record "$hdir/h$nn.dat" --kr0 "$k" $props ${tmax:+--t-max $tmax} \
        --csv "$dest/records/h$nn.csv" > "$dest/records/h$nn.txt" 2>&1 || true
    fi
  done

  # --- 2. level-set DMD of the twin difference over the same window ---
  if [[ -d "$RUNROOT/$label.base" ]]; then
    if need "$dest/export_ls/perturbation.bin"; then
      ./bin/run_lsa export --perturbed "$RUNROOT/$label" --base "$RUNROOT/$label.base" --fields L0101 \
        --export "$dest/export_ls" > "$dest/export_ls.log" 2>&1 || true
    fi
    if [[ -f "$dest/export_ls/perturbation.bin" ]] && need "$dest/dmd_ls/eigenvalues.txt"; then
      nmax=$(awk -v t="${tmax:-1e30}" '$1+0 <= t+1e-9 {n++} END{print n+0}' "$dest/export_ls/times.txt")
      skip=$(( nmax / 10 ))          # drop the first tenth of the window: the start-up transient of every harmonic
      ./bin/run_lsa analyse --from-export "$dest/export_ls" --kr0 0.7 $props --weight volume --band 2 \
        --skip "$skip" --max-snapshots "$(( nmax - skip ))" --n-keep 64 --archive "$dest/dmd_ls" \
        > "$dest/dmd_ls.log" 2>&1 || true
    fi
  else
    echo "  no twin $RUNROOT/$label.base: DMD skipped"
  fi

  # --- 3. spectrum table and figure ---
  if command -v python3 >/dev/null 2>&1; then
    python3 plotting/broadband_spectrum.py --archive "$dest/dmd_ls" --records "$dest/records" \
      --export "$dest/export_ls" ${tmax:+--t-max $tmax} --zlen "$L" --nharm "$M" $props --outdir "$dest" || true
  else
    echo "  python3 not found: run plotting/broadband_spectrum.py on the workstation"
  fi
done

#!/usr/bin/env bash

set -u
shopt -s nullglob

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <baseline_raw_file>"
    exit 1
fi

baseline="$1"

base="$PWD"
orig_dir="${base}/orig"
tsss_dir="${base}/tsss_tans_mvcomp"

mkdir -p "$orig_dir" "$tsss_dir"

# Move raw files into orig/ if they are still in the current folder
rawfiles=( *raw*fif* )
if (( ${#rawfiles[@]} > 0 )); then
    mv "${rawfiles[@]}" "$orig_dir"/
fi

cd "$orig_dir" || {
    echo "ERROR: Cannot access $orig_dir"
    exit 1
}

baseline_file="${orig_dir}/${baseline}"

if [[ ! -e "$baseline_file" ]]; then
    echo "ERROR: baseline file not found: $baseline_file"
    exit 1
fi

echo "Using baseline for -trans: $baseline_file"

runs=( *raw*fif* )

if (( ${#runs[@]} == 0 )); then
    echo "No raw FIF files found in $orig_dir"
    exit 1
fi

for run in "${runs[@]}"; do
    echo "Running: $run"
    stem="${run%.fif}"
    hpfile="${tsss_dir}/${stem}_headpos.txt"

 #   /neuro/bin/util/maxfilter \
 #       -f "$run" \
 #       -o "${tsss_dir}/${run}" \
 #       -ctc /neuro/databases/ctc/ct_sparse.fif \
 #       -cal /neuro/databases/sss/sss_cal.dat \
 #       -autobad off \
 #       -st 10 \
 #       -corr .9 \
 #       -trans "$baseline_file" \
 #       -movecomp \
 #       -hpisubt amp \
 #       -hp "$hpfile" \
 #       -force


 #    /neuro/bin/util/maxfilter \
 #       -f "$run" \
 #       -o "${tsss_dir}/${run}" \
 #       -ctc /neuro/databases/ctc/ct_sparse.fif \
 #       -cal /neuro/databases/sss/sss_cal.dat \
 #       -autobad off \
 #       -st 10 \
 #       -corr .9 \
 #       -movecomp \
 #       -hpisubt amp \
 #       -hp "$hpfile" \
 #       -force

#/neuro/bin/util/maxfilter \
#    -f "$run" \
#    -o "${tsss_dir}/${run}" \
#    -ctc /neuro/databases/ctc/ct_sparse.fif \
#    -cal /neuro/databases/sss/sss_cal.dat \
#    -autobad off \
#    -st 10 \
#    -corr .9 \
#    -trans "$baseline_file" \
#    -movecomp inter \
#    -hpistep 200 \
#    -hpicons \
#    -hpisubt amp \
#    -hp "$hpfile" \
#    -force

  
/neuro/bin/util/maxfilter \
    -f "$run" \
    -o "${tsss_dir}/${run}" \
    -ctc /neuro/databases/ctc/ct_sparse.fif \
    -cal /neuro/databases/sss/sss_cal.dat \
    -autobad off \
    -st 10 \
    -corr .9 \
    -force
    
done

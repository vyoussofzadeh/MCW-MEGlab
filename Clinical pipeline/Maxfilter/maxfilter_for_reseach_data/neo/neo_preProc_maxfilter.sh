#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

base="$PWD"
orig_dir="${base}/orig"
tsss_dir="${base}/tsss"

cd -- "$orig_dir" || {
    printf 'ERROR: Cannot access %s\n' "$orig_dir" >&2
    exit 1
}

runs=()
for file in *.fif; do
    if [[ -f "$file" && "$file" =~ _raw(_part[0-9]+)?\.fif$ ]]; then
        runs+=("$file")
    fi
done

if (( ${#runs[@]} == 0 )); then
    printf 'ERROR: No raw FIF files found in %s\n' "$orig_dir" >&2
    exit 1
fi

mkdir -p -- "$tsss_dir"

printf 'Input folder:  %s\n' "$orig_dir"
printf 'Output folder: %s\n' "$tsss_dir"

for run in "${runs[@]}"; do
    output="${tsss_dir}/${run}"

    if [[ -f "$output" ]]; then
        printf 'Skipping existing output: %s\n' "$run"
        continue
    fi

    printf 'Running tSSS: %s\n' "$run"

    /neuro/bin/util/maxfilter \
        -f "$run" \
        -o "$output" \
        -ctc /neuro/databases/ctc/ct_sparse.fif \
        -cal /neuro/databases/sss/sss_cal.dat \
        -autobad off \
        -st 10 \
        -corr .9 \
        -force
done

printf 'tSSS processing completed.\n'
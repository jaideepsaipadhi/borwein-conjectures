#!/bin/bash
# The two jobs the in-session validation (22 Sept 2026) could not finish in its 10-minute windows.
# Run from the repository root:  bash validation/run_remaining.sh
# Needs ~8 GB RAM (the exact n up to 820 segment) and about 1-1.5 h total on 2 cores.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/validation/third" "$ROOT/validation/cubic"

echo "== Third Borwein: exact recomputation n = 580..820 (every coefficient, all classes)"
( cd "$ROOT/third/certificates" && python3 centre3_exact.py 580 820 > "$ROOT/validation/third/rerun_exact_580_820.log" 2>&1 )
tail -1 "$ROOT/validation/third/rerun_exact_580_820.log"      # must print: DONE n=580..820 total violations 0

echo "== Cubic: expansion corroboration check (not part of the proof; ~13 min)"
( cd "$ROOT/cubic/certificates" && python3 expansion_check.py 36000 > "$ROOT/validation/cubic/rerun_expansion_check.log" 2>&1 )
tail -1 "$ROOT/validation/cubic/rerun_expansion_check.log"    # shipped log: max |rest|/bound = 0.8660

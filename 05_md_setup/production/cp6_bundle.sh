#!/usr/bin/env bash
# Build the GPU transfer bundle for all four GO19 v2 systems.
#
# The rented GPU runs the ENTIRE chain from the built system:
#   complex.gro -> min1 -> min2 -> nvt -> npt1..npt4 -> production
# so the required start point is complex.gro, not min2.gro.  Yoltla's
# min2.gro and nvt.gro are carried as *_yoltla.gro reference files, for
# cross-checking the GPU minimisation against the CPU one (same GROMACS
# 2024.3, different hardware -- they should agree).
#
# Runs on the Yoltla login node, seconds.  No `set -e`: an unguarded
# failure must not skip the gates.
set -uo pipefail

IAN="${IAN:-/LUSTRE/home/lancad/2026/l.100073/alumnos/ian}"
OUT="$IAN/gpu_bundle"
SYSTEMS="go19_v2 go19_v2_R53A go19_v2_R272L go19_v2_R53A_R272L"
FAIL=0

echo "=============================================================="
echo " CP6 GPU bundle   (full chain from complex.gro)"
echo " source: $IAN"
echo "=============================================================="

rm -rf "$OUT"; mkdir -p "$OUT" || { echo "cannot create $OUT"; exit 1; }

for D in $SYSTEMS; do
    echo ""
    echo "--- $D"
    mkdir -p "$OUT/$D"

    # ---- required: everything the GPU needs to run the chain from scratch
    for f in gromacs/complex.gro gromacs/complex.top gromacs/index.ndx \
             gromacs/posre_system1.itp gromacs/posre_FAD.itp gromacs/posre_LIG.itp; do
        if [[ -f "$IAN/$D/$f" ]]; then
            cp "$IAN/$D/$f" "$OUT/$D/$(basename "$f")"
        else
            echo "  MISSING (required): $D/$f"; FAIL=1
        fi
    done

    # ---- reference only: Yoltla results, for comparison not for input
    [[ -f "$IAN/$D/min/min2.gro"  ]] && cp "$IAN/$D/min/min2.gro"  "$OUT/$D/min2_yoltla.gro"
    [[ -f "$IAN/$D/equil/nvt.gro" ]] && cp "$IAN/$D/equil/nvt.gro" "$OUT/$D/nvt_yoltla.gro"

    T="$OUT/$D/complex.top"
    G="$OUT/$D/complex.gro"
    n="-"; other="-"

    # ---- gate 1: posres include blocks present (fix_posres.py applied)
    if [[ -f "$T" ]]; then
        n=$(grep -c '^#ifdef POSRES' "$T")
        [[ "$n" == "3" ]] || { echo "  GATE FAIL: $n '#ifdef POSRES' blocks, expect 3"; FAIL=1; }

        # ---- gate 2: self-contained (no includes outside the bundle)
        other=$(grep '^#include' "$T" | grep -vc 'posre_')
        if [[ "$other" != "0" ]]; then
            echo "  GATE FAIL: $other #include beyond posre_*.itp:"
            grep '^#include' "$T" | grep -v 'posre_' | sed 's/^/      /'
            FAIL=1
        fi
    fi

    # ---- gate 3: the corrected ligand, in the file the GPU will actually read
    if [[ -f "$G" ]]; then
        tot=$(sed -n '2p' "$G" | tr -d ' ')
        lig=$(awk  'NR>2 && substr($0,6,5) ~ /LIG/ {c++} END{print c+0}' "$G")
        ligH=$(awk 'NR>2 && substr($0,6,5) ~ /LIG/ && substr($0,11,5) ~ /H/ {c++} END{print c+0}' "$G")
        fad=$(awk  'NR>2 && substr($0,6,5) ~ /FAD/ {c++} END{print c+0}' "$G")
        [[ "$lig"  == "16" ]] || { echo "  GATE FAIL: $lig LIG atoms, expect 16"; FAIL=1; }
        [[ "$ligH" == "6"  ]] || { echo "  GATE FAIL: LIG has $ligH H, expect 6 (corrected -2 ligand)"; FAIL=1; }
        [[ "$fad"  == "84" ]] || { echo "  GATE FAIL: $fad FAD atoms, expect 84"; FAIL=1; }
        echo "  POSRES=$n  foreign-includes=$other  atoms=$tot  LIG=$lig  LIG-H=$ligH  FAD=$fad"
    fi

    # ---- gate 4: index.ndx carries the Solute group production needs
    if [[ -f "$OUT/$D/index.ndx" ]]; then
        if grep -qE '^\[ *Solute *\]' "$OUT/$D/index.ndx"; then
            echo "  index.ndx: Solute group present"
        else
            echo "  GATE FAIL: index.ndx has no [ Solute ] group (compressed-x-grps needs it)"
            grep -oE '^\[ *[^]]+ *\]' "$OUT/$D/index.ndx" | tr -d '[] ' | tr '\n' ' ' | sed 's/^/      groups: /'; echo
            FAIL=1
        fi
    fi

    ls "$OUT/$D" | tr '\n' ' ' | sed 's/^/  files: /'; echo
    du -sh "$OUT/$D" | awk '{print "  size: "$1}'
done

echo ""
echo "=============================================================="
if [[ "$FAIL" != "0" ]]; then
    echo " GATES FAILED — bundle NOT archived. Fix the above first."
    exit 1
fi

cd "$IAN" || exit 1
tar czf gpu_bundle.tar.gz gpu_bundle || { echo "tar failed"; exit 1; }
echo " OK — all gates passed"
ls -lh "$IAN/gpu_bundle.tar.gz" | awk '{print "  "$5"  "$9}'
echo ""
echo " md5 (compare after upload to the GPU):"
md5sum "$IAN/gpu_bundle.tar.gz" | sed 's/^/   /'
echo "=============================================================="

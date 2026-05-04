#!/bin/bash
set -euo pipefail

make

echo "| Lattice | Mode | Connective Constant |"
echo "| :--- | :--- | :--- |"

declare -A expected=(
    ["2:standard"]="2.286331"
    ["2:hamiltonian"]="1.553873"
    ["2:2sap"]="1.529273"
    ["2:2sap_ham"]="1.329970"
    ["3:standard"]="2.572051"
    ["3:hamiltonian"]="1.629231"
    ["3:2sap"]="1.907603"
    ["3:2sap_ham"]="1.490649"
)

tolerance="0.000010"

for L in 2 3; do
    declare -A modes=( ["standard"]=0 ["hamiltonian"]=1 ["2sap"]=2 ["2sap_ham"]=3 )
    
    for mode_name in "standard" "hamiltonian" "2sap" "2sap_ham"; do
        mode_id=${modes[$mode_name]}
        output=$(./tm_master -L "$L" -M 1 -m "$mode_id")
        result=$(awk '/Connective Constant/ {print $4}' <<< "$output")

        if [[ -z "$result" ]]; then
            echo "Failed to parse connective constant for ${L}x1 ${mode_name}" >&2
            exit 1
        fi

        target=${expected["$L:$mode_name"]}
        awk -v got="$result" -v want="$target" -v tol="$tolerance" '
            BEGIN {
                diff = got - want
                if (diff < 0) diff = -diff
                if (diff > tol) {
                    printf("Regression failure: got %.12f, expected %.12f, tolerance %.12f\n", got, want, tol) > "/dev/stderr"
                    exit 1
                }
            }
        '
        
        echo "| ${L}x1 | ${mode_name} | ${result} |"
    done
done

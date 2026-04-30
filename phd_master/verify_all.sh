#!/bin/bash

# Ensure engine is compiled
make > /dev/null 2>&1

echo "| Lattice | Mode | Connective Constant |"
echo "| :--- | :--- | :--- |"

for L in 2 3; do
    # Map mode names to integers
    declare -A modes=( ["standard"]=0 ["hamiltonian"]=1 ["2sap"]=2 ["2sap_ham"]=3 )
    
    for mode_name in "standard" "hamiltonian" "2sap" "2sap_ham"; do
        mode_id=${modes[$mode_name]}
        
        # Run simulation and capture result
        RESULT=$(./tm_master -L $L -M 1 -m $mode_id | grep "Connective Constant" | awk '{print $4}')
        
        # Format for table
        echo "| ${L}x1 | ${mode_name} | ${RESULT} |"
    done
done

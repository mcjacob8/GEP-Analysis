#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=3000
#SBATCH --output=/volatile/halla/sbs/mcjacob/Fit_optics_%j.out
#SBATCH --error=/volatile/halla/sbs/mcjacob/Fit_optics_%j.err

if [ $# -ne 3 ]; then
    echo "Usage: sbatch $0 <inputlist> <configfile> <skiptree>" >&2
    exit 1
fi

inputlist=$1
configfile=$2
skiptree=$3

export G4SBS=/work/halla/sbs/mcjacob/g4sbs/install
source $G4SBS/bin/g4sbs.sh
export ROOT_INCLUDE_PATH=$G4SBS/include${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}
export LD_LIBRARY_PATH=$G4SBS/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

g4sbsdir=/work/halla/sbs/mcjacob/g4sbs/g4sbs/root_macros

echo "=== Step 1: MakeGEP_optics_tree ==="
if [ "$skiptree" == "1" ]; then
    echo "Skipping Step 1: MakeGEP_optics_tree"
else
    analyzer -b -q "$g4sbsdir/MakeGEP_optics_tree.C+(\"$inputlist\")" || { echo "Step 1 failed" >&2; exit 1; }
fi

echo "=== Step 2: Fit_SBS_optics_and_spin ==="
analyzer -b -q "$g4sbsdir/Fit_SBS_optics_and_spin.C+(\"$configfile\")" || { echo "Step 2 failed" >&2; exit 1; }
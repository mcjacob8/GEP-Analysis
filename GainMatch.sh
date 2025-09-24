#!/bin/sh

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=3000
#SBATCH --output=/volatile/halla/sbs/mcjacob/GainMatch_%j.out
#SBATCH --error=/volatile/halla/sbs/mcjacob/GainMatch_%j.err

configfile=$1

echo "Input file = $configfile"
analyzer -b -q 'GEM_GainMatch_gep.C+("'$configfile'")'

#!/bin/bash

## ------------------------------------------------------------------------ ##
## This shell script is to streamline the physics analysis process. It is   ##
## designed to run all the individual analysis scripts in successive order  ##
## and ultimately produce the final form factor ratio and errors.           ##
##                                                                          ##
## This script should be continously modified overtime as the analysis.     ##
## procedures become more mature.                                           ##
##                                                                          ##
## ---------                                                                ##
##  Jacob McMurtry, rby2vw@virginia.edu CREATED 11-19-2025                  ##
## ---------                                                                ##
## ** Do not tamper with this sticker! Log any updates to the script above. ##
## ------------------------------------------------------------------------ ##

# List of arguments
configfile=$1
run_on_ifarm=$2

#Description of how to run if the user puts in a wrong input
if [ "$#" -ne 2 ]; then
    echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo -e "This script expects the format:\n"
    echo -e "./CalculateRatio.sh <configfile> <run_on_ifarm>"
    echo -e " "
    exit;
fi

if [ "$run_on_ifarm" -eq 1 ]; then
    echo "Submitting jobs to batch farm..."
    # sbatch commands
else
    echo "Running locally..."
    analyzer -b -q 'GetAsymmetry.C+("'$configfile'")'
    analyzer -b -q 'RatioPlot.C+()'
fi
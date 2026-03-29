#!/bin/bash

args=$#
if (( args <= 3 )); then
    failwith "Usage: ./shrinkdcd.sh (psf) (dcd) (#FRAMES) (#FRAMES TO LEAVE)"
fi

PSF=$1
DCD=$2
FRAMES=$3
FRAMES_TO_LEAVE=$4


(( FRAMESTRIDE = FRAMES / FRAMES_TO_LEAVE ))

DCD_WRITE=${DCD%%.dcd}_sampled.dcd

# make smaller trajectory by removing 1 in $FRAMESTRIDE snapshots
if [[ X$FRAMELEAVE != X ]]; then
    vmd -e <(cat - <<"EOF"
mol load psf $PSF
animate read dcd $DCD skip $FRAMESTRIDE
animate write dcd $DCD_WRITE
quit
EOF
    ) -dispdev text
fi

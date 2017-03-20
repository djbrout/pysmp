#!/usr/bin/env bash
awk 'FNR > 50 { nextfile }; /SN-X3/ { print FILENAME }' /pnfs/des/scratch/pysmp/imglistfake/all/des_fake_*.dat | uniq > data/x3lightcurves.txt
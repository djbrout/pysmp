#!/usr/bin/env bash
awk 'FNR > 50 { nextfile }; /SN-S1/ { print FILENAME }' /pnfs/des/scratch/pysmp/imglistfake/des_fake_*.dat | uniq > data/s1lightcurves.txt
#!/usr/bin/env bash
awk 'FNR > 50 { nextfile }; /SN-S2/ { print FILENAME }' /pnfs/des/scratch/pysmp/imglistfake/des_fake_*.dat | uniq > data/s2lightcurves.txt
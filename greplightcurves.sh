#!/usr/bin/env bash
LC_ALL=C fgrep -m 1 -l "SN-S1" /pnfs/des/scratch/pysmp/imglistfake/*.dat > data/s1lightcurves.txt
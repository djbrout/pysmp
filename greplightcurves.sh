#!/usr/bin/env bash
awk 'FNR > 50 { nextfile }; /SN-C2/ { print FILENAME }' /project/projectdirs/dessn/dbrout/imgList/all/des_fake_*.dat | uniq > data/c2lightcurves.txt
#!/usr/bin/env bash
awk 'FNR > 50 { nextfile }; /SN-X2/ { print FILENAME }' /project/projectdirs/dessn/dbrout/imgList/all/des_fake_*.dat | uniq > data/x2lightcurves.txt
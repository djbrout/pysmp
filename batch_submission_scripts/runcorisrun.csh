#!/bin/csh

cd /global/u1/d/dbrout/projects/pysmp
 
#diffimg
source /global/u1/d/dbrout/projects/.diffimg_setup_sextractor
source ~/.cshrc.ext

python smp.py --index=$1 -f $2
###python makeplotsmy.py --index=$1 -f $2



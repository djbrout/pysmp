import os
from subprocess import *
import numpy as np
import time

allindexes = range(0,40)
filts = ['g','r','i','z']
walltime= '10:00:00'

doskipping = True
snfilelist = 'badinputs.txt'
outdir = '/project/projectdirs/dessn/dbrout/simdummytest3'
snfiles = open(snfilelist).read()
snfiles = snfiles.split('.smp')

text = '#!/bin/bash -l\n#SBATCH --partition=shared\n' \
       '#SBATCH -N 1\n#SBATCH -C knl,quad,flat\n' \
       '#SBATCH -A m2875\n' \
       '#SBATCH --time=' + walltime + '\n' + \
       '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/knldummy.log\n' + \
       '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/knldummy.log\n' + \
       '#SBATCH --job-name=preps1\n' + \
       '#SBATCH --mail-type=NONE\n' + \
       '#SBATCH --mail-user=bdrizzle@yahoo.com\n' + \
       '#SBATCH --gres=craynetwork:1\n' + \
       '\n' + \
       'cd /project/projectdirs/des/djbrout/pysmp/\n' + \
       'source setup_scripts/setupcori2.sh\n' + \
       'echo "RUNNING NOW"\n' + \
       'export WALLTIME=' + walltime.split(':')[0] + '\n'

for i in allindexes:
    if snfiles[i] == '\n': continue
    filt = snfiles[i].split('_')[-1]
    snfiles[i] = snfiles[i][:-2]+'.dat'
    text += 'python smptest.py -f ' + filt +\
            ' -o '+outdir+' -s /project/projectdirs/des/djbrout/pysmp/imglist/all/'+snfiles[i]+' '\
            '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ & \n'\
            '\n'\

text += 'wait \n \n \n'

script = '/global/cscratch1/sd/dbrout/logs/knlprep.sh'
f = open(script, 'w')
f.write(text)
f.close()
output = Popen(["sbatch", script], stdout=PIPE).communicate()
print output[0]

import os
from subprocess import *
import numpy as np
import time

allindexes = range(75,90)
filts = ['g','r','i','z']
walltime= '00:30:00'

doskipping = True
snfilelist = 'data/s1lightcurves.txt'
outdir = '/project/projectdirs/dessn/dbrout/simdummytest3'
snfiles = open(snfilelist).readlines()
#snfiles = snfiles.split('.smp')

text = '#!/bin/bash -l\n#SBATCH --partition=debug\n' \
       '#SBATCH -N 1\n#SBATCH -C haswell\n' \
       '#SBATCH -A des\n' \
       '#SBATCH --time=' + walltime + '\n' + \
        '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/debugdummy4.log\n' + \
       '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/debugdummy4.log\n' + \
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
    #if snfiles[i] == '\n': continue
    for filt in filts:
    #filt = snfiles[i].split('_')[-1]
    #snfiles[i] = snfiles[i][:-2]+'.dat'
        snfile = snfiles[i].split('/')[-1].strip()
        text += 'python smptest.py -f ' + filt +\
            ' -o '+outdir+' -s /project/projectdirs/des/djbrout/pysmp/imglist/all/'+snfile+' '\
            '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ & \n'\
            '\n'\

text += 'wait \n \n \n'

script = '/global/cscratch1/sd/dbrout/logs/knlprep.sh'
f = open(script, 'w')
f.write(text)
f.close()

output = Popen(["sbatch", script], stdout=PIPE).communicate()
print output[0]
#print text
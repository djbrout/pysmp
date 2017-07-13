import os
from subprocess import *
import numpy as np
import time

allindexes = range(0,1800)
#allindexes = [100,107,113,120,13,178,214,269,278,40,60,80,92]
filts = ['g','r','i','z']
cntr = 0
#filts = ['r']
walltime= '48:00:00'
#np.random.shuffle(allindexes)

doskipping = True
snfilelist = 'data/x3lightcurves.txt'
outdir = '/project/projectdirs/des/djbrout/simtestx3'
snfiles = open(snfilelist).readlines()

script = '/global/cscratch1/sd/dbrout/logs/scavenger.sh'
f = open(script, 'w')
strg = '#!/bin/bash -l\n' +\
            '#SBATCH --partition=scavenger\n' +\
            '#SBATCH -N 230\n' +\
            '#SBATCH -C haswell\n'+\
            '#SBATCH --time='+walltime+'\n' +\
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/scavenger.log\n' +\
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/scavenger.log\n' +\
            '#SBATCH --job-name=scav \n' +\
            '#SBATCH --mail-type=NONE\n' +\
            '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +\
            '#SBATCH --gres=craynetwork:1\n' +\
            '\n' +\
            'cd /project/projectdirs/des/djbrout/pysmp/\n' +\
            'source setup_scripts/setupcori2.sh\n'+\
            'echo "RUNNING NOW"\n'+\
            'export WALLTIME='+walltime.split(':')[0]+'\n'



for i in allindexes:
    for filt in filts:
        if doskipping:
            sn = snfiles[i].split('/')[-1].split('.')[0]
            if os.path.exists(outdir+'/lightcurves/'+sn+'_'+filt+'.smp'):
                print 'skipping ',outdir+'/lightcurves/'+sn+'_'+filt+'.smp  because already exists a good fit...'
                continue
            # else:
            #     print 'nope',outdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #     continue


            strg += 'srun -n 1 python smpall.py --index=' + str(i) + ' -f ' + filt +\
            ' -o '+outdir+' --snfilelist='+snfilelist+' --usefake --skipdone '+\
            '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ & \n'


            #'python smpshift.py --index=' + str(i) + ' -f ' + filt + ' --nozpt \n'
            # 'python smpshift.py --index=' + str(i) + ' -f '+filt+' --nozpt --snfilelist=data/x3lightcurves.txt '
            #                                                       '-o /project/projectdirs/des/djbrout/116simdeep '
            #                                                      '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ \n' +
            '\n'
f.write(strg)
f.close()
output = Popen(["sbatch", script], stdout=PIPE).communicate()
print output[0]
        #time.sleep(1)
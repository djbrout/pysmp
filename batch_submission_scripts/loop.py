import os
from subprocess import *
import numpy as np
import time

allindexes = range(0,1800)
filts = ['g','r','i','z']
#filts = ['g']
#np.random.shuffle(allindexes)

for i in allindexes:
    for filt in filts:
        print i
        script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=shared\n' +
            '#SBATCH -n 1\n' +
            '#SBATCH -c 1\n'+
            '#SBATCH -C haswell\n'+
            '#SBATCH -A dessn\n' +
            '#SBATCH --time=42:00:00\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_'+filt+'spc.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_'+filt+'spc.log\n' +
            '#SBATCH --job-name=sm_'+filt+'_' + str(i) + '\n' +
            '#SBATCH --mail-type=NONE\n' +
            #'#SBATCH --qos=premium\n'+
            '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
            '#SBATCH --gres=craynetwork:1\n' +
            '\n' +
            'cd /project/projectdirs/des/djbrout/pysmp/\n' +
            'source setup_scripts/setupcori2.sh\n'+
            #'source /scratch3/scratchdirs/masao/setup_DiffImg.sh\n'
            'echo "RUNNING NOW"\n'+
            #'python test.py\n'
            #'cd /global/u1/d/dbrout/SEaR/\n' +
            #'echo "--start='+str(i*nproc)+' --stop='+str((i+1)*nproc)+'" \n'+
            #'python mpp.py --start='+str(i*nproc)+' --stop='+str((i+1)*nproc)+' \n'
            #'python mpp.py --start=' + str(i * nproc) + ' --stop=' + str((i + 1) * nproc) + ' \n'

            'python smpall.py --index=' + str(i) + ' -f ' + filt +
            ' -o /project/projectdirs/des/djbrout/simv2 --snfilelist=data/s2lightcurves.txt --usefake '+
            '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/'


            #'python smpshift.py --index=' + str(i) + ' -f ' + filt + ' --nozpt \n'
            # 'python smpshift.py --index=' + str(i) + ' -f '+filt+' --nozpt --snfilelist=data/x3lightcurves.txt '
            #                                                       '-o /project/projectdirs/des/djbrout/116simdeep '
            #                                                      '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ \n' +
            '\n'
        )
        f.close()
        output = Popen(["sbatch", script], stdout=PIPE).communicate()
        print output[0]
        #time.sleep(1)
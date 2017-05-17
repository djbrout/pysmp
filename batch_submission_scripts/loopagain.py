import os
from subprocess import *
import numpy as np
import time
outdir = '/project/projectdirs/des/djbrout/spec_shallow_v2/'
redofiles = open(outdir+'/missing.txt','r')

#allindexes = range(250,500)
filts = ['g','r','i','z']
#filts = ['g']
#np.random.shuffle(allindexes)

for i in redofiles:

    if True:
    #for filt in filts:
        print i
        fl = i.split()[0]
        filt = i.split()[1]
        script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=shared\n' +
            '#SBATCH -n 1\n' +
            '#SBATCH -c 1\n'+
            '#SBATCH -C haswell\n'+
            '#SBATCH -A dessn\n' +
            '#SBATCH --time=35:00:00\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(fl) + '_'+filt+'.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(fl) + '_'+filt+'.log\n' +
            '#SBATCH --job-name=spec'+filt+'_' + str(fl) + '\n' +
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

            #'python smpshift.py --index=' + str(i) + ' -f ' + filt + ' --nozpt \n'
            'python smpshift.py  -f '+filt+' -o '+outdir+
                                                ' -s /project/projectdirs/des/djbrout/pysmp/imglist/spec/'+fl+'.dat'+
                                                ' \n' +
            '\n'
        )
        f.close()
        output = Popen(["sbatch", script], stdout=PIPE).communicate()
        print output[0]
        #time.sleep(1)
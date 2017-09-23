import os
from subprocess import *
import numpy as np
import time

#allindexes = range(250,500)
#filts = ['g','r','i','z']
#filts = ['g']
#np.random.shuffle(allindexes)

#fields = ['S1','S2','X1','X2','X3','C1','C2','C3']#,'E1','E2'
ccdnums = range(0,63)
#for i in ccdnums:
if True:
    for field in fields:
        #print i
        script = '/global/cscratch1/sd/dbrout/logs/sm_' +field+ '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=xfer\n' +
            '#SBATCH -M escori\n'+
            '#SBATCH -A dessn\n' +
            '#SBATCH --time=00:30:00\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/'+field+'getim.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' +field+'getim.log\n' +
            '#SBATCH --job-name='+field+' \n' +
            '#SBATCH --mail-type=NONE\n' +
            #'#SBATCH --qos=premium\n'+
            '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
            '\n' +
            'cd /project/projectdirs/des/djbrout/pysmp/\n' +
            'source  ' + field +
            '\n'
        )
        f.close()
        output = os.system("sbatch --array=1-63%1 "+ script)
        print 'Field Submitted',field
        #print output[0]
        #time.sleep(1)
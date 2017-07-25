import os
from subprocess import *
import numpy as np
import time

allindexes = range(0,500)
#allindexes = [100,107,113,120,13,178,214,269,278,40,60,80,92]
filts = ['g','r','i','z']
#filts = ['r']
walltime= '02:00:00'
#np.random.shuffle(allindexes)

doskipping = False
#snfilelist = 'badinputs.txt'
snfilelist = 'data/s1lightcurves.txt'
outdir = '/project/projectdirs/des/djbrout/simdummytest3'
snfiles = open(snfilelist).read()
#snfiles = snfiles.split('.smp')

for i in allindexes:
    for filt in filts:

        if doskipping:
            print snfiles[i]
            sn = snfiles[i].split('/')[-1].split('.')[0]
            if os.path.exists(outdir+'/lightcurves/'+sn+'_'+filt+'.smp'):
                print 'skipping ',outdir+'/lightcurves/'+sn+'_'+filt+'.smp  because already exists a good fit...'
                continue
            if os.path.exists('/global/cscratch1/sd/dbrout/specnpzfiles/'+sn+'_'+filt+'.mcmcinput.npz'):
                print 'skipping ', outdir + '/lightcurves/' + sn + '_' + filt + '.smp  because already exists a good fit...'
                continue
            # else:
            #     print 'nope',outdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #     continue
        print i
        script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=shared\n' +
            '#SBATCH -n 1\n' +
            '#SBATCH -c 1\n'+
            '#SBATCH -C haswell\n'+
            '#SBATCH -A des\n' +
            '#SBATCH --time='+walltime+'\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_'+filt+'reallysim2.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_'+filt+'reallysim2.log\n' +
            '#SBATCH --job-name=preps_'+filt+'' + str(i) + '\n' +
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
            'export WALLTIME='+walltime.split(':')[0]+'\n'+
            'python smptest.py --index=' + str(i) + ' -f ' + filt +
            ' -o '+outdir+' --snfilelist='+snfilelist+' --usefake '+
            '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ \n'


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
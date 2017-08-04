import os
from subprocess import *
import numpy as np
import time

allindexes = range(420,421)
#allindexes = [100,107,113,120,13,178,214,269,278,40,60,80,92]
#filts = ['g','r','i','z']
#filts = ['r']
walltime= '15:00:00'
#np.random.shuffle(allindexes)

doskipping = False
#snfilelist = 'badinputs.txt'
#snfilelist = 'data/s2lightcurves.txt'
outdir = '/project/projectdirs/dessn/dbrout/simtestdummy/lightcurves/'
npzdir = '/global/cscratch1/sd/dbrout/simnpzfiles/'
#snfiles = open(snfilelist).readlines()
#snfiles = snfiles.split('.smp')

for i in allindexes:
    #for filt in filts:
    if True:
        if doskipping:
            print snfiles[i]
            sn = snfiles[i].split('/')[-1].split('.')[0]
            #if os.path.exists(outdir+'/lightcurves/'+sn+'_'+filt+'.smp'):
            #    print 'skipping ',outdir+'/lightcurves/'+sn+'_'+filt+'.smp  because already exists a good fit...'
            #    continue
            if os.path.exists(npzdir+'/'+sn+'_'+filt+'.mcmcinput.npz'):
                print 'skipping ', outdir + '/lightcurves/' + sn + '_' + filt + '.smp  because already exists a good fit...'
                continue
            # else:
            #     print 'nope',outdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #     continue
        #print i,'submitted'
        #continue
        script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=shared\n' +
            '#SBATCH -n 1\n' +
            '#SBATCH -c 1\n'+
            #'#SBATCH -C haswell\n'+
            '#SBATCH -A des\n' +
            '#SBATCH --time='+walltime+'\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_mcmcsim.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_mcmcsim.log\n' +
            '#SBATCH --job-name=sim_' + str(i) + '\n' +
            '#SBATCH --mail-type=NONE\n' +
            #'#SBATCH --qos=premium\n'+
            '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
            '#SBATCH --gres=craynetwork:0\n' +
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
            'python mcmc_manager.py --index=' + str(i) + ' --outpath='+outdir+' --npzfolder='+npzdir+' '+
            ' \n'


            #'python smpshift.py --index=' + str(i) + ' -f ' + filt + ' --nozpt \n'
            # 'python smpshift.py --index=' + str(i) + ' -f '+filt+' --nozpt --snfilelist=data/x3lightcurves.txt '
            #                                                       '-o /project/projectdirs/des/djbrout/116simdeep '
            #                                                      '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ \n' +
            '\n'
        )
        f.close()
        output = Popen(["sbatch", script], stdout=PIPE).communicate()
        print output[0]
        #print open(script).read()
        #time.sleep(1)
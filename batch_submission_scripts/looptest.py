import os
from subprocess import *
import numpy as np
import time

#allindexes = np.arange(0,2500)
#18820
allindexes = np.arange(7000,8000)
#edison needs to run 3002 - 4000
#now up to 7000

#allindexes = np.arange(362,4089,1)
#allindexes = [100,107,113,120,13,178,214,269,278,40,60,80,92]
filts = ['g','r','i','z']
#filts = ['z']
#'g','r',
walltime= '05:00:00'
#np.random.shuffle(allindexes)

doskipping = True
#snfilelist = 'badinputs.txt'

snfilelist = 'data/alllightcurves.txt'
#snfilelist = 'data/speclist.txt'


outdir = '/project/projectdirs/dessn/dbrout/simv2.0/'
npzdir = '/global/cscratch1/sd/dbrout/simnpzfilesv2.1/'

outdir = '/project/projectdirs/dessn/dbrout/simv3.0/'
npzdir = '/global/cscratch1/sd/dbrout/simnpzfilesv3.0/'


#outdir = '/project/projectdirs/dessn/dbrout/specv2.0/'
#npzdir = '/global/cscratch1/sd/dbrout/specnpzfilesv2.0/'

snfiles = open(snfilelist).readlines()
#snfiles = snfiles.split('.smp')
count = 0
tot = 0
for ii in allindexes:
    i = int(round(ii))

    for filt in filts:
        tot += 1
        #if count > 10: continue
        if doskipping:
            #print snfiles[i]
            sn = snfiles[int(i)].split('/')[-1].split('.')[0]
            #if os.path.exists(outdir+'/lightcurves/'+sn+'_'+filt+'.smp'):
            #    print 'skipping ',outdir+'/lightcurves/'+sn+'_'+filt+'.smp  because already exists a good fit...'
            #    continue
            print npzdir+'/'+sn+'_'+filt+'.mcmcinput.npz'
            if os.path.exists(npzdir+'/'+sn+'_'+filt+'.mcmcinput.npz'):
                print 'skipping ', outdir + '/lightcurves/' + sn + '_' + filt + '.mcmcinput.npz  because already exists a good fit...'
                print count, tot

                continue
            # else:
            #     print 'nope',outdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #     continue
        if not os.path.exists(snfiles[int(i)].strip):
            print  snfiles[int(i)].strip(),'PASSED'
            continue
        print snfiles[int(i)],'submitted'
        count += 1
        #if count < 78: continue
        print count, tot
        #continue
        script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=shared\n' +
            '#SBATCH -n 1\n' +
            '#SBATCH -c 1\n'+
            '#SBATCH -C haswell\n'+
            '#SBATCH -A dessn\n' +
            '#SBATCH --time='+walltime+'\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_'+filt+'simv3.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_'+filt+'simv3.log\n' +
            '#SBATCH --job-name=sim-'+filt+'' + str(i) + '\n' +
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

            'python smptest.py --usefake --index=' + str(int(i)) + '  -f  ' + filt +
            ' -o '+outdir+' --snfilelist='+snfilelist+' --savenpzfilesdir='+npzdir+' '+
            ' --snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ \n' +
            #' --snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/v4/real/ \n' +





            #
            # 'python smptest.py --usefake --index=' + str(int(i*2+1)) + '  -f  ' + filt +
            # ' -o '+outdir+' --snfilelist='+snfilelist+' --savenpzfilesdir='+npzdir+' '+
            # ' --snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ & \n \n' +
            #
            # 'wait\n'+
            
            #' --snfilepath=/project/projectdirs/dessn/dbrout/imgList/all/ \n'


            #'python smpshift.py --index=' + str(i) + ' -f ' + filt + ' --nozpt \n'
            # 'python smpshift.py --index=' + str(i) + ' -f '+filt+' --nozpt --snfilelist=data/x3lightcurves.txt '
            #                                                       '-o /project/projectdirs/des/djbrout/116simdeep '
            #                                                      '--snfilepath=/project/projectdirs/des/djbrout/pysmp/imglist/all/ \n' +
            '\n'
        )
        f.close()
        #if count >= 269: continue
        output = Popen(["sbatch", script], stdout=PIPE).communicate()
        print output[0]
        #print script
        #raw_input()

        #raw_input('stopppp')
        # print open(script).read()
        # raw_input('stopppp')
        #time.sleep(1)


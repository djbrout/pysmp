import os
from subprocess import *
import numpy as np
import time

allsn = ["des_real_01248677_r.mcmcinput.npz",
"des_real_01248677_i.mcmcinput.npz",
"des_real_01248677_z.mcmcinput.npz",
"des_real_01253039_g.mcmcinput.npz",
"des_real_01253039_r.mcmcinput.npz",
"des_real_01253920_g.mcmcinput.npz",
"des_real_01253920_r.mcmcinput.npz",
"des_real_01250017_g.mcmcinput.npz",
"des_real_01259412_g.mcmcinput.npz",
"des_real_01259412_i.mcmcinput.npz",
"des_real_01259412_z.mcmcinput.npz",
"des_real_01257366_g.mcmcinput.npz",
"des_real_01257366_r.mcmcinput.npz",
"des_real_01258940_g.mcmcinput.npz",
"des_real_01261086_g.mcmcinput.npz",
"des_real_01261086_r.mcmcinput.npz",
"des_real_01261086_i.mcmcinput.npz",
"des_real_01261086_z.mcmcinput.npz",
"des_real_01263715_g.mcmcinput.npz",
"des_real_01256823_z.mcmcinput.npz",
"des_real_01275946_g.mcmcinput.npz",
"des_real_01277168_g.mcmcinput.npz",
"des_real_01266657_g.mcmcinput.npz",
"des_real_01277971_g.mcmcinput.npz",
"des_real_01277971_r.mcmcinput.npz",
"des_real_01277971_i.mcmcinput.npz",
"des_real_01277971_z.mcmcinput.npz",
"des_real_01280201_r.mcmcinput.npz",
"des_real_01280201_i.mcmcinput.npz",
"des_real_01280201_z.mcmcinput.npz",
"des_real_01257695_g.mcmcinput.npz",
"des_real_01262214_g.mcmcinput.npz",
"des_real_01262214_r.mcmcinput.npz",
"des_real_01262214_i.mcmcinput.npz",
"des_real_01262214_z.mcmcinput.npz",
"des_real_01300516_g.mcmcinput.npz",
"des_real_01289664_g.mcmcinput.npz",
"des_real_01285160_g.mcmcinput.npz",
"des_real_01285160_r.mcmcinput.npz",
"des_real_01285160_i.mcmcinput.npz",
"des_real_01285160_z.mcmcinput.npz",
"des_real_01297026_r.mcmcinput.npz",
"des_real_01297026_i.mcmcinput.npz",
"des_real_01297026_z.mcmcinput.npz",
"des_real_01290127_i.mcmcinput.npz",
"des_real_01291090_g.mcmcinput.npz",
"des_real_01294014_g.mcmcinput.npz",
"des_real_01293758_g.mcmcinput.npz",
"des_real_01293758_r.mcmcinput.npz",
"des_real_01299775_g.mcmcinput.npz",
"des_real_01292336_g.mcmcinput.npz",
"des_real_01281668_g.mcmcinput.npz",
"des_real_01283373_g.mcmcinput.npz",
"des_real_01281886_g.mcmcinput.npz",
"des_real_01302117_g.mcmcinput.npz",
"des_real_01292332_g.mcmcinput.npz",
"des_real_01295256_g.mcmcinput.npz",
"des_real_01292145_g.mcmcinput.npz",
"des_real_01291639_g.mcmcinput.npz",
"des_real_01291639_r.mcmcinput.npz",
"des_real_01291639_i.mcmcinput.npz",
"des_real_01291639_z.mcmcinput.npz",
"des_real_01294743_i.mcmcinput.npz",
"des_real_01294743_z.mcmcinput.npz",
"des_real_01286398_r.mcmcinput.npz",
"des_real_01286398_i.mcmcinput.npz",
"des_real_01286398_z.mcmcinput.npz",
"des_real_01283936_g.mcmcinput.npz",
"des_real_01296321_g.mcmcinput.npz",
"des_real_01283923_g.mcmcinput.npz",
"des_real_01283923_z.mcmcinput.npz",
"des_real_01289656_g.mcmcinput.npz",
"des_real_01289656_r.mcmcinput.npz",
"des_real_01291957_g.mcmcinput.npz",
"des_real_01302187_g.mcmcinput.npz",
"des_real_01302187_r.mcmcinput.npz",
"des_real_01302187_i.mcmcinput.npz",
"des_real_01302187_z.mcmcinput.npz",
"des_real_01289288_g.mcmcinput.npz",
"des_real_01295027_g.mcmcinput.npz",
"des_real_01285317_g.mcmcinput.npz",
"des_real_01279500_r.mcmcinput.npz",
"des_real_01279500_i.mcmcinput.npz",
"des_real_01279500_z.mcmcinput.npz",
"des_real_01287626_g.mcmcinput.npz",
"des_real_01283878_g.mcmcinput.npz",
"des_real_01283878_r.mcmcinput.npz",
"des_real_01283878_i.mcmcinput.npz",
"des_real_01283878_z.mcmcinput.npz",
"des_real_01313594_g.mcmcinput.npz",
"des_real_01308884_g.mcmcinput.npz",
"des_real_01308884_r.mcmcinput.npz",
"des_real_01308884_i.mcmcinput.npz",
"des_real_01308884_z.mcmcinput.npz",
"des_real_01306626_g.mcmcinput.npz",
"des_real_01308751_g.mcmcinput.npz",
"des_real_01306299_g.mcmcinput.npz",
"des_real_01308582_g.mcmcinput.npz",
"des_real_01308582_r.mcmcinput.npz",
"des_real_01308582_i.mcmcinput.npz",
"des_real_01308582_z.mcmcinput.npz",
"des_real_01306360_g.mcmcinput.npz",
"des_real_01302648_g.mcmcinput.npz",
"des_real_01314897_i.mcmcinput.npz",
"des_real_01314897_z.mcmcinput.npz",
"des_real_01310338_g.mcmcinput.npz",
"des_real_01316431_g.mcmcinput.npz",
"des_real_01262715_g.mcmcinput.npz",
"des_real_01262715_r.mcmcinput.npz",
"des_real_01262715_i.mcmcinput.npz",
"des_real_01262715_z.mcmcinput.npz",
"des_real_01306073_z.mcmcinput.npz",
"des_real_01303496_g.mcmcinput.npz",
"des_real_01305504_r.mcmcinput.npz",
"des_real_01305504_z.mcmcinput.npz",
"des_real_01307830_g.mcmcinput.npz",
"des_real_01307830_r.mcmcinput.npz",
"des_real_01307830_i.mcmcinput.npz",
"des_real_01307830_z.mcmcinput.npz",
"des_real_01303952_g.mcmcinput.npz",
"des_real_01258906_g.mcmcinput.npz",
"des_real_01305626_g.mcmcinput.npz",
"des_real_01305626_r.mcmcinput.npz",
"des_real_01305626_i.mcmcinput.npz",
"des_real_01305626_z.mcmcinput.npz",
"des_real_01301933_g.mcmcinput.npz",
"des_real_01301933_i.mcmcinput.npz",
"des_real_01301933_z.mcmcinput.npz",
"des_real_01315259_g.mcmcinput.npz",
"des_real_01306537_g.mcmcinput.npz",
"des_real_01303004_g.mcmcinput.npz",
"des_real_01316385_g.mcmcinput.npz",
"des_real_01316385_r.mcmcinput.npz",
"des_real_01316385_i.mcmcinput.npz",
"des_real_01316385_z.mcmcinput.npz",
"des_real_01317164_g.mcmcinput.npz",
"des_real_01316465_g.mcmcinput.npz",
"des_real_01316465_r.mcmcinput.npz",
"des_real_01316465_i.mcmcinput.npz",
"des_real_01316465_z.mcmcinput.npz",
"des_real_01329196_g.mcmcinput.npz",
"des_real_01329196_r.mcmcinput.npz",
"des_real_01329196_i.mcmcinput.npz",
"des_real_01329196_z.mcmcinput.npz",
"des_real_01327978_g.mcmcinput.npz",
"des_real_01318142_g.mcmcinput.npz",
"des_real_01318142_r.mcmcinput.npz",
"des_real_01318142_i.mcmcinput.npz",
"des_real_01318142_z.mcmcinput.npz",
"des_real_01319366_g.mcmcinput.npz",
"des_real_01319821_g.mcmcinput.npz",
"des_real_01329312_g.mcmcinput.npz",
"des_real_01329312_r.mcmcinput.npz",
"des_real_01329312_i.mcmcinput.npz",
"des_real_01329312_z.mcmcinput.npz",
"des_real_01322979_g.mcmcinput.npz",
"des_real_01322979_r.mcmcinput.npz",
"des_real_01322979_i.mcmcinput.npz",
"des_real_01322979_z.mcmcinput.npz",
"des_real_01325358_r.mcmcinput.npz",
"des_real_01325358_i.mcmcinput.npz",
"des_real_01325358_z.mcmcinput.npz",
"des_real_01330426_g.mcmcinput.npz",
"des_real_01329166_g.mcmcinput.npz",
"des_real_01329166_r.mcmcinput.npz",
"des_real_01329166_i.mcmcinput.npz",
"des_real_01329166_z.mcmcinput.npz",
"des_real_01320166_g.mcmcinput.npz",
"des_real_01297465_g.mcmcinput.npz",
"des_real_01297465_r.mcmcinput.npz",
"des_real_01297465_i.mcmcinput.npz",
"des_real_01297465_z.mcmcinput.npz",
"des_real_01291080_g.mcmcinput.npz",
"des_real_01291080_i.mcmcinput.npz",
"des_real_01291080_z.mcmcinput.npz",
"des_real_01330044_g.mcmcinput.npz",
"des_real_01330031_g.mcmcinput.npz",
"des_real_01295305_i.mcmcinput.npz",
"des_real_01295305_z.mcmcinput.npz",
"des_real_01317666_r.mcmcinput.npz",
"des_real_01317666_i.mcmcinput.npz",
"des_real_01317666_z.mcmcinput.npz",
"des_real_01339450_g.mcmcinput.npz",
"des_real_01334858_g.mcmcinput.npz",
"des_real_01334858_r.mcmcinput.npz",
"des_real_01334858_i.mcmcinput.npz",
"des_real_01334858_z.mcmcinput.npz",
"des_real_01337703_g.mcmcinput.npz",
"des_real_01337703_r.mcmcinput.npz",
"des_real_01336453_g.mcmcinput.npz",
"des_real_01334879_r.mcmcinput.npz",
"des_real_01334879_i.mcmcinput.npz",
"des_real_01334879_z.mcmcinput.npz",
"des_real_01331993_i.mcmcinput.npz",
"des_real_01331993_z.mcmcinput.npz",
"des_real_01334597_g.mcmcinput.npz",
"des_real_01334597_r.mcmcinput.npz",
"des_real_01334597_i.mcmcinput.npz",
"des_real_01334597_z.mcmcinput.npz",
"des_real_01335472_g.mcmcinput.npz",
"des_real_01334087_g.mcmcinput.npz",
"des_real_01334423_g.mcmcinput.npz",
"des_real_01334707_r.mcmcinput.npz",
"des_real_01334707_i.mcmcinput.npz",
"des_real_01334707_z.mcmcinput.npz",
"des_real_01333438_g.mcmcinput.npz",
"des_real_01333438_r.mcmcinput.npz",
"des_real_01333438_i.mcmcinput.npz",
"des_real_01333438_z.mcmcinput.npz",
"des_real_01338649_g.mcmcinput.npz",
"des_real_01336975_g.mcmcinput.npz",
"des_real_01338675_g.mcmcinput.npz",
"des_real_01338675_r.mcmcinput.npz",
"des_real_01338675_i.mcmcinput.npz",
"des_real_01338675_z.mcmcinput.npz",
"des_real_01341370_g.mcmcinput.npz",
"des_real_01338843_g.mcmcinput.npz",
"des_real_01338843_r.mcmcinput.npz",
"des_real_01338843_i.mcmcinput.npz",
"des_real_01299643_g.mcmcinput.npz",
"des_real_01300912_g.mcmcinput.npz",
"des_real_01300912_r.mcmcinput.npz",
"des_real_01300912_i.mcmcinput.npz",
"des_real_01300912_z.mcmcinput.npz",
"des_real_01324542_g.mcmcinput.npz",
"des_real_01334302_g.mcmcinput.npz",
"des_real_01335694_g.mcmcinput.npz",
"des_real_01333483_g.mcmcinput.npz",
"des_real_01333483_r.mcmcinput.npz",
"des_real_01333483_i.mcmcinput.npz",
"des_real_01333483_z.mcmcinput.npz",
"des_real_01336009_g.mcmcinput.npz",
"des_real_01335717_g.mcmcinput.npz",
"des_real_01335717_r.mcmcinput.npz",
"des_real_01335717_i.mcmcinput.npz",
"des_real_01335717_z.mcmcinput.npz",
"des_real_01335868_g.mcmcinput.npz",
"des_real_01336480_g.mcmcinput.npz",
"des_real_01336008_g.mcmcinput.npz",
"des_real_01337117_g.mcmcinput.npz",
"des_real_01330642_g.mcmcinput.npz",
"des_real_01330642_i.mcmcinput.npz",
"des_real_01330642_z.mcmcinput.npz",
"des_real_01337687_g.mcmcinput.npz",
"des_real_01338278_g.mcmcinput.npz",
"des_real_01332059_g.mcmcinput.npz",
"des_real_01339392_g.mcmcinput.npz",
"des_real_01335718_r.mcmcinput.npz",
"des_real_01335718_i.mcmcinput.npz",
"des_real_01335718_z.mcmcinput.npz",
"des_real_01336662_g.mcmcinput.npz",
"des_real_01336662_r.mcmcinput.npz",
"des_real_01336662_i.mcmcinput.npz",
"des_real_01336662_z.mcmcinput.npz",
"des_real_01337649_g.mcmcinput.npz",
"des_real_01337228_g.mcmcinput.npz",
"des_real_01336002_g.mcmcinput.npz",
"des_real_01343208_r.mcmcinput.npz",
"des_real_01343208_i.mcmcinput.npz",
"des_real_01343208_z.mcmcinput.npz",
"des_real_01337655_g.mcmcinput.npz",
"des_real_01337655_i.mcmcinput.npz",
"des_real_01337655_z.mcmcinput.npz",
"des_real_01346387_r.mcmcinput.npz",
"des_real_01346387_i.mcmcinput.npz",
"des_real_01346387_z.mcmcinput.npz",
"des_real_01344274_g.mcmcinput.npz",
"des_real_01346956_g.mcmcinput.npz",
"des_real_01346956_r.mcmcinput.npz",
"des_real_01343871_g.mcmcinput.npz",
"des_real_01343759_i.mcmcinput.npz",
"des_real_01343759_z.mcmcinput.npz",
"des_real_01304442_g.mcmcinput.npz",
"des_real_01304442_r.mcmcinput.npz",
"des_real_01304442_z.mcmcinput.npz",
"des_real_01248907_g.mcmcinput.npz",
"des_real_01248907_r.mcmcinput.npz",
"des_real_01345899_g.mcmcinput.npz",
"des_real_01345899_r.mcmcinput.npz",
"des_real_01345899_i.mcmcinput.npz",
"des_real_01345899_z.mcmcinput.npz",
"des_real_01297501_g.mcmcinput.npz",
"des_real_01297501_r.mcmcinput.npz",
"des_real_01297501_i.mcmcinput.npz",
"des_real_01297501_z.mcmcinput.npz",
"des_real_01292195_g.mcmcinput.npz",
"des_real_01292195_r.mcmcinput.npz",
"des_real_01292195_i.mcmcinput.npz",
"des_real_01292195_z.mcmcinput.npz",
"des_real_01284587_g.mcmcinput.npz",
"des_real_01289600_g.mcmcinput.npz",
"des_real_01289600_r.mcmcinput.npz",
"des_real_01289600_i.mcmcinput.npz",
"des_real_01289600_z.mcmcinput.npz",
"des_real_01290568_g.mcmcinput.npz",
"des_real_01290568_r.mcmcinput.npz",
"des_real_01290568_i.mcmcinput.npz",
"des_real_01290568_z.mcmcinput.npz",
"des_real_01290816_g.mcmcinput.npz",
"des_real_01290816_r.mcmcinput.npz",
"des_real_01302058_g.mcmcinput.npz",
"des_real_01303279_g.mcmcinput.npz",
"des_real_01303279_r.mcmcinput.npz",
"des_real_01304678_g.mcmcinput.npz",
"des_real_01304678_r.mcmcinput.npz",
"des_real_01308314_g.mcmcinput.npz",
"des_real_01309288_g.mcmcinput.npz",
"des_real_01309288_r.mcmcinput.npz",
"des_real_01309492_g.mcmcinput.npz",
"des_real_01309492_r.mcmcinput.npz",
"des_real_01309492_i.mcmcinput.npz",
"des_real_01309492_z.mcmcinput.npz",
"des_real_01309749_g.mcmcinput.npz",
"des_real_01309749_r.mcmcinput.npz",
"des_real_01309749_i.mcmcinput.npz",
"des_real_01309749_z.mcmcinput.npz",
"des_real_01315192_g.mcmcinput.npz",
"des_real_01322229_g.mcmcinput.npz",
"des_real_01322229_r.mcmcinput.npz",
"des_real_01322229_i.mcmcinput.npz",
"des_real_01322229_z.mcmcinput.npz",
"des_real_01340454_g.mcmcinput.npz",
"des_real_01261579_g.mcmcinput.npz",
"des_real_01261579_r.mcmcinput.npz",
"des_real_01261579_i.mcmcinput.npz",
"des_real_01261579_z.mcmcinput.npz",
"des_real_01347120_g.mcmcinput.npz",
"des_real_01347120_r.mcmcinput.npz",
"des_real_01347120_i.mcmcinput.npz",
"des_real_01347120_z.mcmcinput.npz",
"des_real_01343337_g.mcmcinput.npz",
"des_real_01291794_g.mcmcinput.npz",
"des_real_01345798_g.mcmcinput.npz",
"des_real_01345798_r.mcmcinput.npz",
"des_real_01345798_i.mcmcinput.npz",
"des_real_01345798_z.mcmcinput.npz",
"des_real_01302523_r.mcmcinput.npz",
"des_real_01302523_i.mcmcinput.npz",
"des_real_01302523_z.mcmcinput.npz",
"des_real_01346137_g.mcmcinput.npz",
"des_real_01253101_g.mcmcinput.npz",
"des_real_01292315_g.mcmcinput.npz",
"des_real_01292315_r.mcmcinput.npz",
"des_real_01292315_i.mcmcinput.npz",
"des_real_01292315_z.mcmcinput.npz",
"des_real_01249851_g.mcmcinput.npz",
"des_real_01249851_r.mcmcinput.npz",
"des_real_01336687_g.mcmcinput.npz",
"des_real_01303883_g.mcmcinput.npz"]
#filts = ['g','r','i','z']
#filts = ['r']
walltime= '00:30:00'
#np.random.shuffle(allindexes)

doskipping = False
#snfilelist = 'badinputs.txt'
#snfilelist = 'data/s2lightcurves.txt'
outdir = '/project/projectdirs/dessn/dbrout/specv1_1/lightcurvestest/'
npzdir = '/global/cscratch1/sd/dbrout/specnpzfiles/'
#npzdir = '/project/projectdirs/dessn/dbrout/specv1_1/npzfiles/'
#snfiles = open(snfilelist).readlines()
#snfiles = snfiles.split('.smp')

#for i in allsn[::-1]:
if True:
    #for filt in filts:
    if True:
        if doskipping:
            pass
            # print snfiles[i]
            # sn = snfiles[i].split('/')[-1].split('.')[0]
            # if os.path.exists(outdir+'/'+i.split()+'.smp'):
            #     print 'skipping ',outdir+'/lightcurves/'+sn+'_'+filt+'.smp  because already exists a good fit...'
            #     continue
            # if os.path.exists(npzdir+'/'+sn+'_'+filt+'.mcmcinput.npz'):
            #     print 'skipping ', outdir + '/lightcurves/' + sn + '_' + filt + '.smp  because already exists a good fit...'
            #     continue
            # else:
            #     print 'nope',outdir+'/lightcurves/'+sn+'_'+filt+'.smp'
            #     continue
        # print i,'submitted'
        # continue
        script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
        f = open(script, 'w')
        f.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --partition=shared\n' +
            '#SBATCH -n 1\n' +
            '#SBATCH -c 1\n'+
            #'#SBATCH -C haswell\n'+
            '#SBATCH --array=1-50\n'
            '#SBATCH -A des\n' +
            '#SBATCH --time='+walltime+'\n' +
            #'#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_mcmcspec.log\n' +
            #'#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_mcmcspec.log\n' +
            '#SBATCH --output=/project/projectdirs/des/djbrout/pysmp/logs/${SLURM_ARRAY_TASK_ID}_mcmcspec.log\n'+
            '#SBATCH --error=/project/projectdirs/des/djbrout/pysmp/logs/${SLURM_ARRAY_TASK_ID}_mcmcspec.log\n'+
            '#SBATCH --job-name=' + str(i)[8:] + '\n' +
            '#SBATCH --mail-type=NONE\n' +
            #'#SBATCH --qos=premium\n'+
            #'#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
            '#SBATCH -L SCRATCH,project,cscratch1'+
            #'#SBATCH --gres=craynetwork:0\n' +
            '\n' +
            'cd /project/projectdirs/des/djbrout/pysmp/\n' +
            'source setup_scripts/setupcori2.sh\n'+
            #'source /scratch3/scratchdirs/masao/setup_DiffImg.sh\n'
            'echo "RUNNING NOW"\n'+
            'echo $HOSTNAME\n '
            #'python test.py\n'
            #'cd /global/u1/d/dbrout/SEaR/\n' +
            #'echo "--start='+str(i*nproc)+' --stop='+str((i+1)*nproc)+'" \n'+
            #'python mpp.py --start='+str(i*nproc)+' --stop='+str((i+1)*nproc)+' \n'
            #'python mpp.py --start=' + str(i * nproc) + ' --stop=' + str((i + 1) * nproc) + ' \n'
            'export WALLTIME='+walltime.split(':')[0]+'\n'+
            'python mcmc_manager.py --index=0 --sn=${SLURM_ARRAY_TASK_ID} --outpath='+outdir+' --npzfolder='+npzdir+' '+
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
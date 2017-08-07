import os
from subprocess import *
import numpy as np
import time

allindexes = range(1350,1351)
allindexes = [27
,29
,147
,121
,241
,242
,258
,260
,173
,172
,174
,336
,545
,508
,346
,348
,320
,620
,619
,621
,250
,252
,316
,771
,772
,774
,763
,765
,764
,766
,808
,810
,179
,180
,182
,679
,680
,745
,747
,746
,748
,683
,682
,684
,599
,832
,831
,833
,842
,844
,843
,845
,934
,936
,944
,903
,905
,904
,906
,915
,914
,916
,930
,932
,550
,551
,553
,524
,876
,875
,877
,1058
,1060
,1062
,1064
,974
,976
,1036
,1054
,1056
,999
,998
,1000
,1219
,591
,593
,1001
,1003
,1002
,1004
,1086
,1088
,964
,1090
,1092
,1130
,1129
,1131
,1270
,1167
,1330
,1332
,1288
,664
,1322
,1324
,556
,467
,466
,468
,371
,370
,372
,396
,783
,787
,789
,788
,790
,899
,901
,900
,902
,166
,1342
,1344
,1320
,625]
#filts = ['g','r','i','z']
#filts = ['r']
walltime= '48:00:00'
#np.random.shuffle(allindexes)

doskipping = False
#snfilelist = 'badinputs.txt'
#snfilelist = 'data/s2lightcurves.txt'
outdir = '/project/projectdirs/dessn/dbrout/specv1_1/lightcurves/'
npzdir = '/global/cscratch1/sd/dbrout/specnpzfiles/'
#npzdir = '/project/projectdirs/dessn/dbrout/specv1_1/npzfiles/'
#snfiles = open(snfilelist).readlines()
#snfiles = snfiles.split('.smp')

for i in allindexes:
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
            '#SBATCH -C haswell\n'+
            '#SBATCH -A des\n' +
            '#SBATCH --time='+walltime+'\n' +
            '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_mcmcspec.log\n' +
            '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_mcmcspec.log\n' +
            '#SBATCH --job-name=spec_' + str(i) + '\n' +
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
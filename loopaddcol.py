from subprocess import *

allindexes = range(0,1)

for i in allindexes:

    print i
    script = '/global/cscratch1/sd/dbrout/logs/sm_' + str(i) + '.sh'
    f = open(script, 'w')
    f.write(
        '#!/bin/bash -l\n' +
        '#SBATCH --partition=shared\n' +
        '#SBATCH -n 1\n' +
        '#SBATCH -c 1\n'+
        '#SBATCH --array=1-255 \n'+
        '#SBATCH -C haswell\n'+
        '#SBATCH -A dessn\n' +
        '#SBATCH --time=00:19:00\n' +
        '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_sim.log\n' +
        '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_sim.log\n' +
        '#SBATCH --job-name=ss_' + str(i) + '\n' +
        '#SBATCH --mail-type=NONE\n' +
        #'#SBATCH --qos=premium\n'+
        '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
        '#SBATCH --gres=craynetwork:1\n' +
        '\n' +
        'cd /project/projectdirs/des/djbrout/pysmp/\n' +
        'source setup_scripts/setupcori2.sh\n'+

        'python  addcoltoDESlightcurve.py --index=$SLURM_JOBID \n'
        #SIM
        #'python addcoltoDESlightcurve.py --index=$SLURM_JOB_ID --savelcdir=/project/projectdirs/des/djbrout/specv1_3/SMP_SPEC_v1_3/ '
        #                                                     '--resultsdir=/project/projectdirs/des/djbrout/specv1_3//  \n' +
        # 'python addcoltoDESlightcurve.py --index=' + str(
        #     i) + ' --savelcdir=/project/projectdirs/des/djbrout/spec_v7/SMP_RAW_specv1 '
        #          '--resultsdir=/project/projectdirs/des/djbrout/spec_v7/ \n' +


        '\n'
    )
    f.close()
    output = Popen(["sbatch", script], stdout=PIPE).communicate()
    print output[0]
    #print script

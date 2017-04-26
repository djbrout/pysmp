from subprocess import *

allindexes = range(0,1000)

for i in allindexes:

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
        '#SBATCH --time=0:10:00\n' +
        '#SBATCH --output=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_dsim.log\n' +
        '#SBATCH --error=/global/cscratch1/sd/dbrout/logs/' + str(i) + '_dsim.log\n' +
        '#SBATCH --job-name=addcol_' + str(i) + '\n' +
        '#SBATCH --mail-type=NONE\n' +
        #'#SBATCH --qos=premium\n'+
        '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
        '#SBATCH --gres=craynetwork:1\n' +
        '\n' +
        'cd /project/projectdirs/des/djbrout/pysmp/\n' +
        'source setup_scripts/setupcori2.sh\n'+
        'python addcoltoDESlightcurve.py --index=' + str(i) + ' --dofakes --savelcdir=/project/projectdirs/des/djbrout/114sim/SMP_RAW_SIM_v1_13 \n' +
        '\n'
    )
    f.close()
    output = Popen(["sbatch", script], stdout=PIPE).communicate()
    print output[0]
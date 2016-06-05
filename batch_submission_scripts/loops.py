





























import os
from subprocess import *	
import numpy as np
for i in np.arange(0,500):
	print i
	script = '/global/u1/d/dbrout/PySMP/submission_scripts/smp_'+str(i)+'.sh'
	f = open(script,'w')
	f.write(
		'#!/bin/csh\n' + 
		'#SBATCH --partition=shared\n' +
		'#SBATCH --nodes=1\n' +
		#'#SBATCH -c 8\n' + 
		'#SBATCH -A des\n' +
		'#SBATCH --time=20:00:00\n' +
		'#SBATCH --output=/global/cscratch1/sd/dbrout/v3/batchout/smp_'+str(i)+'_v3_64newsky_floatpos_galsim_galsimzpt_skyerr_2cores_r.log\n'+
		'#SBATCH --error=/global/cscratch1/sd/dbrout/v3/batchout/smp_'+str(i)+'_v3_64newsky_floatpos_galsim_galsimzpt_skyerr_2cores_r.err\n'+
		'#SBATCH --job-name=r_'+str(i)+'\n' +
		'#SBATCH --mail-type=All\n' +
		'#SBATCH --mail-user=djbrout@gmail.com\n' +
		'#SBATCH --gres=craynetwork:1\n' +
		'\n'+
		'cd  /global/u1/d/dbrout/PySMP\n'+
		'srun -n 1 /global/u1/d/dbrout/projects/pysmp/runcorisrun.csh '+str(i)+' r \n' +
		'\n'
		)
	f.close()
	#batcmd='sbatch '+script
	output = Popen(["sbatch", script], stdout=PIPE).communicate()
	jobid = output[0].strip().split(' ')[3]
	#os.system('cp '+script+' /global/u1/d/dbrout/PySMP/submission_scripts/smp_jobid_'+str(jobid)+'.sh')
	print output[0]


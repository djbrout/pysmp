#!/bin/bash -l                                                                                                                                                                 
#SBATCH --array=52-53
#SBATCH -M escori
#SBATCH --qos=xfer                                                             
#SBATCH -t 12:00:00                                                      
#SBATCH -J hstar_job                                                               
#SBATCH -o hstar_job.o%j                                                         
#SBATCH --licenses=SCRATCH      

version=10
cd /global/cscratch1/sd/dbrout/images_v"$version"/
aaa=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-S1_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-S2_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar


#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-C1_CCD"$aaa"_v"$version".tar
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-C2_CCD"$aaa"_v"$version".tar
htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-C3_CCD"$aaa"_v"$version".tar


#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-X1_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-X2_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-X3_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-E1_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar
#htar -xf /home/projects/dessn/diffim/FinalPhoto/v"$version"/SN-E2_CCD"$SLURM_ARRAY_TASK_ID"_v"$version".tar

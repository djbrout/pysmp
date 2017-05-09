resultsdir = '/project/projectdirs/des/djbrout/114sim/'
savelcdir = resultsdir+'/SMP_RAW_SIMshallow_v3'


import os


readmetext = '## Example output lightcurve file\n\
Contains the original forced photometry lightcurve file structure with the addition of several columns\n\
```\n\
SMP_FLUX - fit flux (at zpt of 31)\n\
SMP_FLUXERR - fit flux error (at zpt of 31)\n\
SMP_FLUX_ZPT - zeropoint of smp flux (always 31)\n\
SMP_FIT_ZPT - smp fit zeropoint which is used to scale all fluxes to 31\n\
SMP_FIT_ZPT_STD - uncertainty in smp fit zeropoint which is used to scale all fluxes to 31\n\
SMP_CHISQ - (sim_stamp - image_stamp)^2/(err)^2\n\
SMP_SKY - fit value for the sky\n\
SMP_SKYERR - uncertainty in the sky\n\
SMP_FIX - 1 means SN flux was fixed to zero in smp fit\n\
SMP_FLAG - 1 means this epoch was flagged for some reason inside smp pipeline\n\
```\n\
(just to be clear: ALL FLUXES, SKYS, SKYERRS, ETC... ARE REPORTED AT A ZEROPOINT OF 31)\n\
-999 means missing data\n\
\n\
## Flag Bit Definitions\n\
```\n\
CHISQ_FLAG = 2\n\
PIPELINE_FLAG = 1\n\
BADSKY_FLAG = 4\n\
BADSKYERR_FLAG = 8\n\
BADZPT_FLAG = 16\n\
BADZPTERR_FLAG = 32\n\
DONTFIT_FLAG = 65536\n\
```\n'


readme = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.README', 'w')
readme.write(readmetext)
readme.close()

os.popen('cd '+savelcdir+'\n ls *.dat > '+savelcdir.split('/')[-1]+'.LIST')

os.popen('cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/')
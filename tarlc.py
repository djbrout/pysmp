# resultsdir = '/project/projectdirs/des/djbrout/simv1'
# savelcdir = '/project/projectdirs/des/djbrout/SMP_SIM_v1'
resultsdir = '/project/projectdirs/des/djbrout/simv1'
savelcdir = '/project/projectdirs/des/djbrout/SMP_SIM_v1'


import os


readmetext = '## ' \
             'This is a spectroscopically confirmed DES-SN sample\n\
   for seasons Y1-Y3. Based on fakes, we expect .034% 8 sigma lightcurve fit outliers.\n ' \
             'We find very little hostgal sb fudge needed for shallow fields based on fakes. Deep fields we need to run more bright hostgalaxy sb fakes in order to better determine fudge\n\
And also include the TYPE definitions from the DIFFIMG README.\n\
 Spectroscopic Tags and Statistics (from SNCAND table)\n\
         SNTYPE                Ncand\n\
         ------------------------------\n\
           0 (NO_SPEC     )       0\n\
           1 (SNIa        )       239\n\
         101 (SNIa?       )       0\n\
           3 (SNIax       )       0\n\
           4 (SNIa-pec    )       0\n\
          20 (SNIIP       )       0\n\
          21 (SNIIL       )       0\n\
          22 (SNIIn       )       0\n\
          23 (SNIIb       )       0\n\
          29 (SNII        )       0\n\
         122 (SNIIn?      )       0\n\
         129 (SNII?       )       0\n\
          32 (SNIb        )       0\n\
          33 (SNIc        )       0\n\
          39 (SNIbc       )       0\n\
         139 (SNIbc?      )       0\n\
          41 (SLSN-I      )       0\n\
          42 (SLSN-II     )       0\n\
          43 (SLSN-R      )       0\n\
         141 (SLSN-I?     )       0\n\n\n\
             Example output lightcurve file\n\
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
Column Definition\n\
PHOTPROB\n\
This is the reduced chi squared of the fit. Should be close to 1. We need to understand where we should be making the cut on photprob before performing light curve fits\n\
Flag Bit Definitions\n\
CHISQ_FLAG = 2 #none have been set\n\
BADSKY_FLAG = 4 # FITSKY > 999.\n\
BADSKYERR_FLAG = 8 # FIT_SKYERR > .5\n\
BADZPT_FLAG = 16 # 27 > FITZPT or FITZPT > 35\n\
BADSCATTER_FLAG = 32# the scatter about the zeropoint greater that .2 mags\n\
SMP_PIPELINE_FLAG = 32768 # A flag was set inside the SMP pipeline (ex. bad PSF, too much masking, etc...)\n\
DIDNT_ENTER_SMP_FLAG = 65536 #These epochs never even made it into the scene modeling pipeline (for example the image never existed on NERSC)\n\
Issues\n\
We\'re seeing relatively lower fitprob values compared to DIFFIMG\n\
Bright galaxy fitprob issue persists\n\
To Do List\n\
Re-Fit with Y4 images included for better galaxy models (and to resolve 1 case SNID 01315841 band r\n\
Provide x and y chi squared stamp slopes for further fit quality assessment\n\
```\n'


readme = open(savelcdir + '/' + savelcdir.split('/')[-1] + '.README', 'w')
readme.write(readmetext)
readme.close()

os.popen('cd '+savelcdir+'\n ls *.dat > '+savelcdir.split('/')[-1]+'.LIST \n cd .. \n rm '+savelcdir.split('/')[-1]+'.tar.gz \n tar -zcvf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+' \n')

#os.popen('cd '+savelcdir+'\n ls *.dat > '+savelcdir.split('/')[-1]+'.LIST \n cd .. \n rm '+savelcdir.split('/')[-1]+'.tar.gz \n tar -xf '+savelcdir.split('/')[-1]+'.tar.gz ./'+savelcdir.split('/')[-1]+'\n')
#os.popen('rm '+savelcdir.split('/')[-1]+'.tar.gz ')
#os.popen('cd '+resultsdir+'\n tar -zcf '+savelcdir.split('/')[-1]+'.tar.gz '+savelcdir.split('/')[-1]+'/')

# print ('cd '+resultsdir+'/lightcurves/ \n ls *.pdf > pdflist \n tar -zcf '+savelcdir.split('/')[-1]+'_stamps.tar.gz '
#                                                     '-I pdflist \n mv'
#                                                     ' '+savelcdir.split('/')[-1]+'_stamps.tar.gz'
#                                                 ' --exclude=\'*.npz\' --exclude=\'*.smp\' --exclude=\'*.png\' --exclude=\'*.gz\' *' \
#                                                 ''+savelcdir+'/'+savelcdir.split('/')[-1]+'_stamps.tar.gz \n')
#os.popen('rm '+savelcdir.split('/')[-1]+'_stamps.tar.gz ')
raw_input('done with first tar')
os.popen('cd '+resultsdir+'/lightcurves/ \n ls *.pdf > pdflist \n tar -zcvf '
                                                    ' '+savelcdir.split('/')[-1]+'_stamps.tar.gz '
                                                    '  --exclude=\'*.npz\' --exclude=\'*.smp\' --exclude=\'*.png\' --exclude=\'*.gz\' * \n')
                                                    # '\n mv'
                                                    # ' '+savelcdir.split('/')[-1]+'_stamps.tar.gz'
                                                    # ' '+savelcdir+'/'+savelcdir.split('/')[-1]+'_stamps.tar.gz \n')

os.popen('cd '+resultsdir+'/lightcurves/ \n ls *.smp > pdflist \n tar -zcvf '
                                                    ' '+savelcdir.split('/')[-1]+'_smpfiles.tar.gz '
                                                    '  --exclude=\'*.npz\' --exclude=\'*.pdf\' --exclude=\'*.png\' --exclude=\'*.gz\' * \n')
                                                    # '\n mv'
                                                    # ' '+savelcdir.split('/')[-1]+'_smpfiles.tar.gz'
                                                    # ' '+savelcdir+'/'+savelcdir.split('/')[-1]+'_smpfiles.tar.gz \n')
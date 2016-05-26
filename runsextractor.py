'''

Dillon Brout
dbrout@physics.upenn.edu


Python function to grab
sextractor sky and skyerr
values from any image


'''


#import os
import subprocess

def getsky_and_skyerr(imagefilename):

    proc = subprocess.Popen("sex "+imagefilename+" -c defaults/default.sex", stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    #print "program output:", out
    #proc = subprocess.Popen(['python', 'fake_utility.py'], stdout=subprocess.PIPE)
    while True:
        line = proc.stdout.readline()
        if line != '':
            # the real code does filtering here
            print "test:", line.rstrip()
        else:
            break

            #os.system('sex '+imagefilename+' -c defaults/default.sex')


im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
getsky_and_skyerr(im)
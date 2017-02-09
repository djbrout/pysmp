import os
import numpy as np
import time

maxlightcurves = 1000
hostname = os.popen('hostname').read()
if 'dsp053' in hostname:
    isdonedir = '/home/dscolnic/testdir/isdone'
    logdir = '/home/dscolnic/testdir/smplogs'

if 'dsp057' in hostname:
    isdonedir = '/home/dscolnic/testdir/isdone'
    logdir = '/home/dscolnic/testdir/smplogs'

lcf = open('/home/dscolnic/testdir/runfile.txt','r')
smpdir = '/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST20'


runninglist = np.chararray(50,itemsize=300)

lightcurves = lcf.readlines()
lcf.close()

corelist = np.arange(24)

cntr = -1
i = -1
offset = 150
while i < 23:
    cntr += 1
    if lightcurves[cntr+offset][:11]+lightcurves[cntr+offset][9:].split('.')[0]+'.smp' in os.listdir(os.path.join(smpdir,'lightcurves')):
        print lightcurves[cntr+offset].split('.')[0]+'.smp', 'already exists'
        cntr += 1
    else:
        i += 1
        print ''
        #print 'taskset -c '+str(int(i))+' python smp.py --nozpt --dontglobalstar --index='+str(int(i))+' > ' \
        #    ''+ os.path.join(logdir,lightcurves[i].split('.')[0]+'.log')+' &'
        os.popen('taskset -c '+str(int(i))+' python smpnsc.py --nozpt --index='+str(int(cntr+offset))+' 1>& '+
                 os.path.join(logdir,lightcurves[cntr+offset].split('.')[0]+'.log')+' &')
        print lightcurves[cntr+offset].strip(),'Submitted to SMP. Core #'+str(int(i))
        print 'See log file here',os.path.join(logdir,lightcurves[cntr+offset].split('.')[0]+'.log')
        print ''
        runninglist[i] = lightcurves[cntr+offset]


j=cntr+1
while j <= maxlightcurves:
    donefiles = os.listdir(isdonedir)
    for core in corelist:
        if runninglist[core].split('.')[0]+'.done' in donefiles:
            kg = True
            while kg:
                if lightcurves[j+offset][:11]+lightcurves[j+offset][9:].split('.')[0]+'.smp' in os.listdir(os.path.join(smpdir, 'lightcurves')):
                    print lightcurves[j+offset].split('.')[0] + '.smp','already exists'
                    j += 1
                else:
                    kg = False
            print ''
            print 'Running SN '+str(int(j+offset))+'/'+str(int(maxlightcurves))
            print runninglist[core],'Has finished photometry on core',int(core)
            #print 'taskset -c ' + str(int(core)) + ' python smp.py --nozpt --dontglobalstar --index=' + str(int(j)) + ' &'
            os.popen('taskset -c ' + str(int(core)) + ' python smpnsc.py --nozpt --index=' + str(int(j+offset)) +' 1>& '+
                os.path.join(logdir,lightcurves[j+offset].split('.')[0]+'.log')+' &')
            print lightcurves[j+offset].strip(),'Submitted to SMP. Core #'+str(int(core))
            print 'See log file here', os.path.join(logdir, lightcurves[j+offset].split('.')[0] + '.log')

            runninglist[core] = lightcurves[j+offset]
            j += 1
            print ''

    time.sleep(10)
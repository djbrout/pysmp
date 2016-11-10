import os
import numpy as np
import time

maxlightcurves = 419
isdonedir = '/home/dscolnic/testdir/isdone'
lcf = open('/home/dscolnic/testdir/runfile.txt','r')
logdir = '/home/dscolnic/testdir/smplogs'
smpdir = '/export/scratch0/ps1sn1/data/v10.0/GPC1v3/eventsv1/smpworkspace/PS_TEST1'


runninglist = np.chararray(24,itemsize=300)

lightcurves = lcf.readlines()
lcf.close()

corelist = np.arange(24)

cntr = -1
i = -1
while i < 23:
    cntr += 1
    #print os.listdir(os.path.join(smpdir,'lightcurves'))
    print lightcurves[cntr][:11]+lightcurves[cntr][9:].split('.')[0]+'.smp'
    raw_input()
    if lightcurves[cntr][:9]+lightcurves[cntr][11:].split('.')[0]+'.smp' in os.listdir(os.path.join(smpdir,'lightcurves')):
        print lightcurves[cntr].split('.')[0]+'.smp', 'already exists'
        cntr += 1
    else:
        i += 1
        print ''
        #print 'taskset -c '+str(int(i))+' python smp.py --nozpt --dontglobalstar --index='+str(int(i))+' > ' \
        #    ''+ os.path.join(logdir,lightcurves[i].split('.')[0]+'.log')+' &'
        os.popen('taskset -c '+str(int(i))+' python smpnsc.py --index='+str(int(cntr))+' 1>& '+
                 os.path.join(logdir,lightcurves[cntr].split('.')[0]+'.log')+' &')
        print lightcurves[cntr].strip(),'Submitted to SMP. Core #'+str(int(i))
        print 'See log file here',os.path.join(logdir,lightcurves[cntr].split('.')[0]+'.log')
        print ''
        runninglist[i] = lightcurves[cntr]


j=cntr+1
while j <= maxlightcurves:
    donefiles = os.listdir(isdonedir)
    for core in corelist:
        if runninglist[core].split('.')[0]+'.done' in donefiles:
            kg = True
            while kg:
                if lightcurves[j].split('.')[0] + '.smp' in os.listdir(os.path.join(smpdir, 'lightcurves')):
                    print lightcurves[j].split('.')[0] + '.smp','already exists'
                    j += 1
                else:
                    kg = False
            print ''
            print 'Running SN '+str(int(j))+'/'+str(int(maxlightcurves))
            print runninglist[core],'Has finished photometry on core',int(core)
            #print 'taskset -c ' + str(int(core)) + ' python smp.py --nozpt --dontglobalstar --index=' + str(int(j)) + ' &'
            os.popen('taskset -c ' + str(int(core)) + ' python smpnsc.py --index=' + str(int(j)) +' 1>& '+
                os.path.join(logdir,lightcurves[j].split('.')[0]+'.log')+' &')
            print lightcurves[j].strip(),'Submitted to SMP. Core #'+str(int(core))
            print 'See log file here', os.path.join(logdir, lightcurves[j].split('.')[0] + '.log')

            runninglist[core] = lightcurves[j]
            j += 1
            print ''

    time.sleep(10)
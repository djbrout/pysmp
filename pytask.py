import os
import numpy as np
import time

maxlightcurves = 419
isdonedir = '/home/dscolnic/testdir/isdone'
lcf = open('/home/dscolnic/testdir/runfile.txt','r')
logdir = '/home/dscolnic/testdir/smplogs'



runninglist = np.chararray(24,itemsize=300)

lightcurves = lcf.readlines()
lcf.close()

corelist = np.arange(24)

for i in range(24):
    print ''
    #print 'taskset -c '+str(int(i))+' python smp.py --nozpt --dontglobalstar --index='+str(int(i))+' > ' \
    #    ''+ os.path.join(logdir,lightcurves[i].split('.')[0]+'.log')+' &'
    os.popen('taskset -c '+str(int(i))+' python smp.py --index='+str(int(i))+' 1>& '+
             os.path.join(logdir,lightcurves[i].split('.')[0]+'.log')+' &')
    print lightcurves[i].strip(),'Submitted to SMP. Core #'+str(int(i))
    print 'See log file here',os.path.join(logdir,lightcurves[i].split('.')[0]+'.log')
    print ''
    runninglist[i] = lightcurves[i]


j=24
while j <= maxlightcurves:
    donefiles = os.listdir(isdonedir)
    for core in corelist:
        if runninglist[core].split('.')[0]+'.done' in donefiles:
            print ''
            print 'Running SN '+str(int(j))+'/'+str(int(maxlightcurves))
            print runninglist[core],'Has finished photometry on core',int(core)
            #print 'taskset -c ' + str(int(core)) + ' python smp.py --nozpt --dontglobalstar --index=' + str(int(j)) + ' &'
            os.popen('taskset -c ' + str(int(core)) + ' python smp.py --index=' + str(int(j)) +' 1>& '+
                os.path.join(logdir,lightcurves[j].split('.')[0]+'.log')+' &')
            print lightcurves[j].strip(),'Submitted to SMP. Core #'+str(int(core))
            print 'See log file here', os.path.join(logdir, lightcurves[j].split('.')[0] + '.log')

            runninglist[core] = lightcurves[j]
            j += 1
            print ''

    time.sleep(10)
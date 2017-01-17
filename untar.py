import os
import sys


f = open('alreadyuntarredv2.txt','r')
afiles = f.read().split()
f.close()
tarfiles = os.listdir('/pnfs/des/persistent/smp/v2')

for tfile in tarfiles[::-1]:
    if not tfile.split('.')[-1] == 'tar':
        continue
    if not 'S1' in tfile:
        if not 'S2' in tfile:
            continue
    if tfile.split('/')[-1] in afiles:
        print 'already untarred',tfile
    else:
        print 'untarring ',tfile
        f = open('alreadyuntarredv2.txt', 'a')
        f.write(tfile.split('/')[-1] + ' ')
        f.close()
        out = os.popen('tar -xvf /pnfs/des/persistent/smp/v2/'+tfile.split('/')[-1]+' -C /pnfs/des/persistent/smp/v2/').read()
        #print out
        #if not 'Exiting' in out:
        #    f = open('alreadyuntarredv6.txt','a')
        #    f.write(tfile.split('/')[-1]+' ')
        #    f.close()

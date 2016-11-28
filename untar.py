import os
import sys


f = open('alreadyuntarredv6.txt','r')
afiles = f.read().split()
f.close()
tarfiles = os.listdir('/pnfs/des/persistent/smp/v6')

for tfile in tarfiles:
    if not tfile.split('.')[-1] == 'tar':
        continue
    if tfile.split('/')[-1] in afiles:
        print 'already untarred',tfile
    else:
        print 'untarring ',tfile
        out = os.popen('tar -xvf /pnfs/des/persistent/smp/v6/'+tfile.split('/')[-1]+' -C /pnfs/des/persistent/smp/v6/').read()
        print out
        if not 'Exiting' in out:
            f = open('alreadyuntarredv6.txt','a')
            f.write(tfile.split('/')[-1]+' ')
            f.close()

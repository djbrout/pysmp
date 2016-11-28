import os
import sys


f = open('alreadyuntarred.txt','r')
afiles = f.read().split()
f.close()
tarfiles = os.listdir('/pnfs/des/persistent/smp/v2')

for tfile in tarfiles[::-1]:
    if not tfile.split('.')[-1] == 'tar':
        continue
    if tfile.split('/')[-1] in afiles:
        print 'already untarred',tfile
    else:
        print 'untarring ',tfile
        out = os.popen('tar -xvf /pnfs/des/persistent/smp/v2/'+tfile.split('/')[-1]+' -C /pnfs/des/persistent/smp/v2/')
        print out
        f = open('alreadyuntarred.txt','a')
        f.write(tfile.split('/')[-1]+' ')
        f.close()

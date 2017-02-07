import numpy as np
import os

v6dir = '/pnfs/des/scratch/pysmp/smp_v622/np_data/r/'
v4dir = '/pnfs/des/scratch/pysmp/smp_v42/np_data/r/'

f6 = os.listdir(v6dir)
ff6 = []
for f in f6:
    if 'withSn' in f:
        ff6.append(f)

f4 = os.listdir(v4dir)
ff4 = []
for f in f4:
    if 'withSn' in f:
        ff4.append(f)


commonfiles = []
for f in ff6:
    if f in ff4:
        commonfiles.append(f)

bigv6stamps = []
bigv4stamps = []
bigfakemags = []
bigfakefluxs = []
for f in commonfiles:
    v6dat = np.load(v6dir+f)
    v4dat = np.load(v4dir + f)
    print v6dat.keys()
    raw_input()
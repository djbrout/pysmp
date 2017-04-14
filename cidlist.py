import numpy as np
import os


outfile = open('CID_FAKES.LIST','w')
outfile.write('# Variables extracted from\n\
#FILE:  OUT_SNFIT_SMP.HBOOK\n\
#TABLE: 7788\n\
#\n\
NVAR: 2\n\
VARNAMES:  CID zHD\n')

for i in os.listdir():
    outfile.write('SN:  '+str(int(cid))+' '+str(round(float(z),6))+'\n')

outfile.close()
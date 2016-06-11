import os
import numpy as np

class tmpwriter():
    # tempdir = location to write files
    # tmp_index = index for parallel computation to avoid over-writing files
    def __init__(self, tempdir,tmp_index=0):
        self.tmpdir = tempdir
        self.tmp_index = str(round(tmp_index))
    def writefile(self,text,filename):
        tempfile = os.path.join(self.tmpir, 'tmp_' + self.tmp_index + '.txt')
        a = open(tempfile,'w')
        a.write(text)
        a.close()
        os.system('cp ' + tempfile + ' ' + filename)

    def appendfile(self,text,filename):
        tempfile  = os.path.join(self.tmpir, 'tmp_' + self.tmp_index + '.txt')
        os.system('cp ' + filename + ' ' + tempfile)
        a = open(tempfile,'a')
        a.write(text)
        a.close()
        os.system('cp ' + tempfile + ' ' + filename)

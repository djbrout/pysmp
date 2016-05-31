import os
import matplotlib as m

m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(threshold=np.nan)
import sys
import dilltools as dt

default_checkstar_file = '/Volumes/ExtraSpace/pysmp_downloads/des_fake_00229567_r_standardstarfits.txt'


def checkstars(checkstarfile=default_checkstar_file):
    cols = dt.readcol(checkstarfile,delim='\t')

    print cols.keys()
    return cols


if __name__ == '__main__':
    a = checkstars()

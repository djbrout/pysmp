import numpy as np
import pyfits as pf
import os
import scipy.signal
import scipy.ndimage as nd

#hello from fermilab2
# Returns xvals, medians, mads
def bindata(x, y, bins, returnn=False):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)

    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs) & (x < bf)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere) & ~np.isnan(yhere)]
        ss = [abs(yhere) < 3. * np.std(yhere)]
        try:
            nums[i] = len(yhere[ss])
            medians[i] = np.median(yhere[ss])
            mads[i] = 1.48 * np.median(abs(yhere[ss] - medians[i])) * 1 / np.sqrt(len(yhere[ss]))
        except IndexError:
            print 'excepted'
            nums[i] = 0.
            medians[i] = np.nan
            mads[i] = np.nan
    xvals = (bins[1:] + bins[:-1]) / 2.
    if returnn:
        return xvals, medians, mads, nums
    return xvals, medians, mads

# Takes in Filename, reads file columnwise, and returns dictionary such that:
# import rdcol
# a = rdcol.read('filename',headline,datastartline)
# a["Column_Name"] -> returns the list for that column
#
# headline and datastartline are always > 0
#
# By Dillon Brout
# dbrout@physics.upenn.edu
def read(filename, headline, startline, delim=' '):
    linenum = 0
    go = 0
    column_list = []
    return_cols = {}
    inf = open(filename)
    for line in inf:
        line = line.replace('#', '')
        line = line.strip()
        cols = line.split(delim)
        cols[:] = (value for value in cols if value != '')
        if linenum == headline - 1:
            for col in cols:
                return_cols[col.strip()] = []
                column_list.append(col.strip())
                go += 1
        if linenum >= startline - 1:
            index = 0
            for col in cols:
                try:
                    return_cols[column_list[index]].append(float(col.strip()))
                except:
                    return_cols[column_list[index]].append(col.strip())
                index += 1
        linenum += 1
    inf.close()
    return return_cols


def save_fits_image(image,filename):
    hdu = pf.PrimaryHDU(image)
    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)

    return


def psfphotometry(im, psf, sky, weight, gal, guess_scale):
    chisqvec = []
    fluxvec = []

    galconv = scipy.signal.fftconvolve(gal, psf, mode='same')

    radius = 12
    substamp = galconv.shape[0]
    # Make a mask with radius
    fitrad = np.zeros([substamp, substamp])
    for x in np.arange(substamp):
        for y in np.arange(substamp):
            if np.sqrt((substamp / 2. - x) ** 2 + (substamp / 2. - y) ** 2) < radius:
                fitrad[int(x), int(y)] = 1.

    if guess_scale is None:
        for i in np.arange(-10000, 200000, 5):
            sim = galconv + sky + i * psf
            chisqvec.append(np.sum((im - sim) ** 2 * weight * fitrad))
            fluxvec.append(i)
    else:
        for i in np.arange(guess_scale - 2000, guess_scale + 2000, 1):
            sim = galconv + sky + i * psf
            chisqvec.append(np.sum((im - sim) ** 2 * weight * fitrad))
            fluxvec.append(i)

    ii = fitrad.ravel()
    i = ii[ii != 0]

    ndof = len(i) + 1

    fluxvec = np.array(fluxvec)
    chisqvec = np.array(chisqvec)
    hh = chisqvec * 0 + min(chisqvec)
    mchisq = min(chisqvec)
    idx = np.isclose(chisqvec, hh, atol=1.)

    sim = galconv + sky + fluxvec[chisqvec == min(chisqvec)] * psf
    sum_data_minus_sim = np.sum(im - sim)
    return fluxvec[chisqvec == min(chisqvec)], fluxvec[chisqvec == min(chisqvec)] - fluxvec[idx][
        0], mchisq / ndof, sum_data_minus_sim


# Takes in Filename, reads file columnwise, and returns dictionary such that:
# import rdcol
# a = rdcol.read('filename',headline,datastartline)
# a["Column_Name"] -> returns the list for that column
#
# headline and datastartline are always > 0
#
# By Dillon Brout
# dbrout@physics.upenn.edu

def readcol(filename,headline=1,startline=2,delim=' '):
    linenum = 0
    go = 0
    column_list = []
    return_cols = {}
    inf = open(filename)
    for line in inf:
        line = line.replace('#', '')
        line = line.strip()
        cols = line.split(delim)
        cols[:] = (value for value in cols if value != '')
        if linenum == headline - 1:
            for col in cols:
                return_cols[col.strip()] = []
                column_list.append(col.strip())
                go += 1
        if linenum >= startline - 1:
            index = 0
            for col in cols:
                try:
                    return_cols[column_list[index]].append(float(col.strip()))
                except:
                    return_cols[column_list[index]].append(col.strip())
                index += 1
        linenum += 1
    inf.close()
    for k in return_cols.keys():
        return_cols[k] = np.array(return_cols[k])
    return return_cols


def pixelate(matrix, pixelation_factor):
    zmatrix = nd.interpolation.zoom(matrix, 1. / float(pixelation_factor))
    return zmatrix

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
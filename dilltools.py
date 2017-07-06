import numpy as np
import pyfits as pf
import os
import scipy.signal
import scipy.ndimage as nd
import scipy.stats
import sys

#hello from fermilab2
# Returns xvals, medians, mads
def bindata(x, y, bins, returnn=False, window=0.,dontrootn=False):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)

    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs-window) & (x < bf+window)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere) & ~np.isnan(yhere)]
        ss = [abs(yhere) < 3. * np.std(yhere)]
        try:
            nums[i] = len(yhere[ss])
            medians[i] = np.median(yhere[ss])
            if dontrootn:
                mads[i] = 1.48 * np.median(abs(yhere[ss] - medians[i]))
            else:
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




def fitprobfromchisq(x,df):
    return scipy.stats.chisquare(x,ddof=df)

def binrms(x, y, bins,rad):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)
    rms = np.zeros(len(bins) - 1)
    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs - rad) & (x < bf + rad)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere) & ~np.isnan(yhere)]
        ss = [abs(yhere) < 3. * np.std(yhere)]
        try:
            nums[i] = len(yhere[ss])
            d = yhere[ss]
            dc = d[abs(d) < 3]

            rms[i] = np.sqrt(np.nanmean(np.square(d[d < 3.])))
        except IndexError:
            print 'excepted'
            nums[i] = 0.
            rms[i] = np.nan
    xvals = (bins[1:] + bins[:-1]) / 2.
    return xvals, rms


def binstd(x, y, bins,rad):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)
    std = np.zeros(len(bins) - 1)
    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs - rad) & (x < bf + rad)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere) & ~np.isnan(yhere)]
        ss = [abs(yhere) < 3. * np.std(yhere)]
        try:
            nums[i] = len(yhere[ss])
            d = yhere[ss]
            dc = d[abs(d) < 3]
            std[i] = np.std(dc)
        except IndexError:
            print 'excepted'
            nums[i] = 0.
            std[i] = np.nan
    xvals = (bins[1:] + bins[:-1]) / 2.
    return xvals, std

def iterstat(d, startMedian=False, sigmaclip=3.0,iter=6):
    """Get the sigma-clipped mean of
    a distribution, d.
    Usage: mean,stdev = iterstat.iterstat
    Input:
    d:           the data
    Optional Inputs:
    sigmaclip:   number of standard deviations to clip
    startMedian: if True, begin with the median of the distribution
    iter:        number of iterations
    """


    clip = sigmaclip
    img = d.astype('float64')
    if startMedian:
        md = np.median(img)
    else:
        md = np.mean(img)
    n = float(len(img))
    std = np.sqrt(np.sum((img - md) ** 2.) / (n - 1))

    for ii in range(iter):
        gd = np.where((img < md + clip * std) &
                      (img > md - clip * std))

        md = np.mean(img[gd])
        num = len(img[gd])
        n = float(len(gd[0]))
        std = np.sqrt(np.sum((img[gd] - md) ** 2.) / (n - 1.))

    return (md, std, num)


def savefits(data, filename, fermigrid=True):
    if os.path.isfile(filename):
        os.remove(filename)
    if fermigrid:
        tempfile = 'tmp.fits'
        # print 'saving to temporary file',tempfile
        #return
        if os.path.isfile(tempfile):
            try:
                os.remove(tempfile)
            except:
                pass
        #if os.path.isfile(filename):
        #    os.remove(filename)
        save_fits_image(data, tempfile,go=True)
        # print 'ifdh cp to ',filename

        os.system('mv ' + tempfile + ' ' + filename)
    else:
        save_fits_image(data, filename, go=True)
    print 'saved', filename

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


def save_fits_image(image,filename,go=False):
    try:
        os.remove(filename)
    except:
        pass
    if go:
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
    numcols = None
    inf = open(filename)
    for line in inf:
        #if linenum == 0:
        #    line = line[:191] + ' '+ line[191:]
        #    #print line
        #    #sys.exit()
        #print line
        line = line.replace('#', '')
        line = line.strip()
        cols = line.split(delim)
        #print len(cols)
        #sys.exit()
        cols[:] = (value for value in cols if value != '')
        if linenum == headline - 1:
            numcols = len(cols)
            #print numcols
            for col in cols:
                return_cols[col.strip()] = []
                column_list.append(col.strip())
                go += 1
        if linenum >= startline - 1:
            index = 0
            if len(cols) != numcols:
                #print 'passinge',len(cols),numcols
                pass
                 #print 'WARNING: Could not read line ' + str(linenum + 1) + ' of ' + filename
            else:
                for col in cols:
                    try:
                        return_cols[column_list[index]].append(float(col.strip()))
                    except:
                        try:
                            return_cols[column_list[index]].append(col.strip())
                        except:
                            pass
                            #print 'WARNING2: Could not read line '+str(linenum+1)+' of '+filename
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
    def __init__(self, tempdir='./tmp/',tmp_subscript=0,usedccp=False,useifdh=False):
        self.tmpdir = tempdir
        self.tmp_index = str(tmp_subscript)
        self.usedccp = usedccp
        self.useifdh = useifdh
    def writefile(self,text,filename):
        tempfile = os.path.join(self.tmpdir, 'tmp_' + self.tmp_index + '.txt')
        if os.path.isfile(tempfile):
            os.remove(tempfile)
        try:
            if os.path.isfile(filename):
                os.remove(filename)
        except:
            print 'could not remove file'
        a = open(tempfile,'w')
        a.write(text)
        a.close()
        if self.usedccp:
            os.system('dccp ' + tempfile + ' ' + filename)
        elif self.useifdh:
            os.system('ifdh cp ' + tempfile + ' ' + filename)
        else:
            os.system('mv ' + tempfile + ' ' + filename)

        print 'saved', filename

    def appendfile(self,texts,filename):
        tempfile  = os.path.join(self.tmpdir, 'tmp_' + self.tmp_index + '.txt')
        if os.path.isfile(tempfile):
            os.remove(tempfile)

        if self.usedccp:
            os.system('dccp ' + filename + ' ' + tempfile)
        elif self.useifdh:
            #print 'ifdh cp ' + filename + ' ' + tempfile
            os.system('ifdh cp ' + filename + ' ' + tempfile)
        else:
            os.system('mv ' + filename + ' ' + tempfile)

        a = open(tempfile, 'a')
        for text in texts:
            a.write(text)
        a.close()

        if self.usedccp:
            os.system('dccp ' + tempfile + ' ' + filename)
        elif self.useifdh:
            os.system('ifdh rm ' + filename)
            os.system('ifdh cp ' + tempfile + ' ' + filename)
        else:
            if os.path.isfile(filename):
                os.remove(filename)
            os.system('mv ' + tempfile + ' ' + filename)


    def savez(self,filename,**kwargs):
        tempfile  = os.path.join(self.tmpdir, 'tmp_' + self.tmp_index + '.npz')
        if os.path.isfile(tempfile):
            os.remove(tempfile)
        if os.path.isfile(filename):
            os.remove(filename)
        np.savez(tempfile,**kwargs)
        if self.usedccp:
            os.system('dccp ' + tempfile + ' ' + filename)
        elif self.useifdh:
            print 'ifdh cp ' + tempfile + ' ' + filename
            os.system('ifdh rm '+filename)
            os.system('ifdh cp ' + tempfile + ' ' + filename)
            os.popen('rm '+tempfile)
        else:
            os.system('mv ' + tempfile + ' ' + filename)
        print 'saved',filename


    def savefits(self,data,filename):

        tempfile = 'tmp.fits'
        #print 'saving to temporary file',tempfile
        return
        if os.path.isfile(tempfile):
            os.remove(tempfile)
        if os.path.isfile(filename):
            os.remove(filename)
        save_fits_image(data,tempfile)
        #print 'ifdh cp to ',filename
        if self.usedccp:
            os.system('dccp ' + tempfile + ' ' + filename)
        elif self.useifdh:
            os.system('ifdh cp ' + tempfile + ' ' + filename)
        else:
            os.system('mv ' + tempfile + ' ' + filename)
        print 'saved', filename

    def cp(self,src,dst):
        if os.path.isfile(dst):
            os.remove(dst)
        # print 'ifdh cp to ',filename
        if self.usedccp:
            os.system('dccp ' + src + ' ' + dst)
        elif self.useifdh:
            ls = os.popen('ifdh ls '+dst).read()
            if len(ls) > 0:
                print os.popen('ifdh rm ' + dst).read()
            print os.popen('ifdh cp ' + src + ' ' + dst).read()
        else:
            os.system('mv ' + src + ' ' + dst)
        print 'saved', dst

    #def doesfileexist(self):

import numpy as np
import pyfits as pf


def bindata(x, y, bins, returnn=False):
    medians = np.zeros(len(bins) - 1)
    mads = np.zeros(len(bins) - 1)
    nums = np.zeros(len(bins) - 1)

    for i in np.arange(len(bins) - 1):
        bs = bins[i]
        bf = bins[i + 1]
        ww = [(x > bs) & (x < bf)]
        yhere = y[ww]
        yhere = yhere[np.isfinite(yhere)]
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
import numpy as np

def bindata(x,y,bins,returnn=False):
        medians = np.zeros(len(bins)-1)
        mads = np.zeros(len(bins)-1)
        nums = np.zeros(len(bins)-1)

        for i in np.arange(len(bins)-1):
                bs = bins[i]
                bf = bins[i+1]
                ww = [(x>bs)&(x<bf)]
                yhere = y[ww]
                yhere = yhere[np.isfinite(yhere)]
                ss = [abs(yhere) < 3.*np.std(yhere)]
                try:
                        nums[i] = len(yhere[ss])
                        medians[i] = np.median(yhere[ss])
                        mads[i] = 1.48*np.median(abs(yhere[ss]-medians[i]))*1/np.sqrt(len(yhere[ss]))
                except IndexError:
                        print 'excepted'
                        nums[i] = 0.
                        medians[i] = np.nan
                        mads[i] = np.nan
        xvals = (bins[1:] + bins[:-1])/2.
        if returnn:
                return xvals,medians,mads,nums
        return xvals,medians,mads

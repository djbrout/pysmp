# pysmp v1.0

## Knobs

### Data Ingestion Knobs
#### Inside default.config
* -r /path/to/topdir/of/images
* --filter=g: Each filter is fit individually so you must specify your filter of interest
* -p /path/to/default.param: discussed below
* -o /path/to/outdir/for/smp
* --snfilepath=/path/to/forced/photometry/smp/input/lightcurves: SMP requires forced photometry lightcurves as an input.
* --nozpt: will always calculate the zeropoint of every image containin the SN ra and dec even if they already been calculated on a past run.
* --usefake: if fakes are present in the image given by Kessler pipeline then this will handle them accordingling and propagate data for later analysis.
* --bigstarcatalog=/path/to/big/star.cat
* --snfilelist=/path/to/file/containing/all/forced/photometry/input/files.txt
#### Inside default.param
* mjdplus: (mjd-peakmjd) larger than this value are fixed to zero SN flux (necessary for SMP)
* mjdminus: (mjd-peakmjd) smaller than this value are fixed to zero SN flux (necessary for SMP)
* HDR_AIRMASS_NAME: AIRMASS
* HDR_PLATESCALE_NAME: PIXSCAL1
* HDR_PSF_FWHM: PSF_FWHM
* FORCERADEC: False
* FRA: 40.514263
* FDEC: -1.630960
### MCMC Knobs
#### Inside default.param
* SN_PLUS_GALMODEL_STEPS: 500000
* SN_SHIFT_STD: 0.0001
* GALMODEL_PROPOSAL_STD: .6


## Doing a Fit

```
python smpshift.py --index=$1 --nozpt
```

### Submitting multiple fits

Submit to your favorite batch grid system with each job containing a different --index=$1 .


Afterwards you will run the following commands to summarize your jobs and regenerate lightcurves with SMP fluxes.

```
python addcoltoDESlightcurve.py --lcdir=/path/to/new/lc/files --resultsdir=/path/to/smp/results
```

## Example output lightcurve file

Will contain the original lightcurve file structure with the addition of several columns

* asdf
* asdf

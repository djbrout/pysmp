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
* SN_FLUX_PROPOSAL_MULT: 2. # Proposal STD is sqrt(diffim_flux)*SN_FLUX_PROPOSAL_DIV
* GALMODEL_PROPOSAL_STD: .6


## Doing a Fit

```
python smpshift.py --index=$1 --nozpt
```

Submit to your favorite batch grid system with each job containing a different --index=$1 .


Afterwards you will run the following command to regenerate lightcurves with SMP fluxes in a new directory! (dont overwrite old lightcurves!)

```
python addcoltoDESlightcurve.py --lcdir=/path/to/new/lc/files --resultsdir=/path/to/smp/results
```

## Example output lightcurve file

Contains the original forced photometry lightcurve file structure with the addition of several columns
```
SMP_FLUX - fit flux (at zpt of 31)
SMP_FLUXERR - fit flux error (at zpt of 31)
SMP_FLUX_ZPT - zeropoint of smp flux (always 31)
SMP_FIT_ZPT - smp fit zeropoint which is used to scale all fluxes to 31
SMP_FIT_ZPT_STD - uncertainty in smp fit zeropoint which is used to scale all fluxes to 31
SMP_CHISQ - (sim_stamp - image_stamp)^2/(err)^2
SMP_SKY - fit value for the sky
SMP_SKYERR - uncertainty in the sky
SMP_FIX - 1 means SN flux was fixed to zero in smp fit
SMP_FLAG - 1 means this epoch was flagged for some reason inside smp pipeline
```
(just to be clear: ALL FLUXES, SKYS, SKYERRS, ETC... ARE REPORTED AT A ZEROPOINT OF 31)

-999 means missing data

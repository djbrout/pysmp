# pysmp v1.0

## Knobs

### Data Ingestion Knobs
#### Inside default.config
* asdf
* asdf
* asdf
#### Inside default.param
* asdf
* asdf
* asdf

### MCMC Knobs
#### Inside default.config
* asdf
* asdf
* asdf
#### Inside default.param
* asdf
* asdf
* asdf


### Doing a Fit

```
python smpshift.py --index=$1 --nozpt
```

### Submitting multiple fits

Submit to your favorite batch grid system with each job containing a different --index=$1
Afterwards you will run the following commands to summarize your jobs and regenerate lightcurves

```
python addcoltoDESlightcurve.py --lcdir=/path/to/new/lc/files --resultsdir=/path/to/smp/results
```


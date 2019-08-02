# Chasing Sea Ice
#### Alek Petty, Linette Boisvert, Jeremy Harbeck

*Original code written by Jeremy Harbeck and Linette Boisvert, adapted to Python and modularized by Alek Petty*

Python code to correct IceBridge flight paths to account for ice drift when underflying ICESat-2
 
## Introduction

### Conda installation

The code was written in Python 3 (3.6) but appears to also work in Python 2.7.

If you're having problems try using the included conda environment file - chasingseaice.yml - to ensure consistency in the Python environment

```
conda env create -f chasingseaice.yml
```

Alternatively you can try generating your own conda environment using the following packages

```
conda create -n chasingseaice python=3.6 numpy scipy matplotlib h5py basemap

```
The conda Python environment can be activated with 

```
source activate chasingseaice
```
### Running the drift correction code

Create a new folder for a given flight date and copy in the sequence file if you wish to produce a drift corrected version.

The code can be run with the following command

```
python driftcorrect.py
```
All the input data are hard-coded into the python file itself, including the measurement time, the IS-2 time, the plane position and some extra options.

Currently it's just setup to do a single drift correction based on the wind measurement time, the plane position and the IS-2 cross-over time. The script should return single values for the drift correction (distance and the direction the ice is coming from) as well as an updated sequence file (ending in .driftcorrected) if the OUT_SEQUENCE flag is True.







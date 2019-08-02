# Chasing Sea Ice
### Python code to correct IceBridge flight paths to account for ice drift when underflying ICESat-2
 

### Alek Petty, Linette Boisvert, Jeremy Harbeck

Original code written by Jeremy Harbeck and Linette Boisvert, adapted to Python and modularized by Alek Petty

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

Create a folder and place in your sequence file 
The code can be run with the following command

```
python driftcorrect.py
```
All the input data are hard-coded into the python file itself.

Currently it's just setup to do a single drift correction based on the wind measurement time, the plane position and the IS-2 cross-over time.







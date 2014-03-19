#!/bin/bash

##NB: Custom pythonpath is due to local build on a cluster,
##    this would be unnecessary with a normal installation

cd `pwd`

## location of GMPY2
export PYTHONPATH=/home/jmblab/statonse/apps/gmpy2/src/gmpy2-2.0.0b3/build/lib.linux-x86_64-2.7
## location of Ecolopy
export PYTHONPATH=$PYTHONPATH:/home/jmblab/statonse/apps/ecolopy/build/lib

python2.7 ecolopy_noopt_harg.py

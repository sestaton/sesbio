#!/bin/bash

cd apps
mkdir gmpy2
cd gmpy2

wget https://ftp.gnu.org/gnu/gmp/gmp-5.1.3.tar.bz2 && tar xjf gmp-5.1.3.tar.bz2
cd gmp-5.1.3
./configure --prefix=/home/jmblab/statonse/apps/gmpy2
make
make check
make install
cd ..

## Download mpir, mpfr, mpc, gmpy2 full source from https://code.google.com/p/gmpy/downloads/list
unzip full-src-mpir-mpfr-mpc-gmpy2-2.0.2.zip
cd src/mpfr-3.1.1
chmod +x configure
./configure --prefix=/home/jmblab/statons/apps/gmpy2 --with-gmp=/home/jmblab/statonse/apps/gmpy2
make
make check
make install

cd ../mpc
chmod +x configure
./configure --prefix=/home/jmblab/statons/apps/gmpy2 --with-gmp=/home/jmblab/statons/apps/gmpy2 --with-mpfr=/home/jmblab/statons/apps/gmpy2
make
make check
make install

cd ../gmpy2-2.0.0b3
/usr/local/python/2.7.2/bin/python setup.py build_ext -Ddir=/home/jmblab/statonse/apps/gmpy2

export PYTHONPATH=/home/jmblab/statonse/apps/gmpy2/build/lib.linux-x86_64-2.7

cd ~/apps
git clone https://github.com/fransua/ecolopy.git
cd ecolopy
/usr/local/python/2.7.2/bin/python setup.py build
export PYTHONPATH=$PYTHONPATH:/home/jmblab/statonse/apps/ecolopy/build/lib

## To test, $DISPLAY needs to be set, so login with `ssh xxxx@zcluster.rcc.uga.egu -X`
/usr/local/python/2.7.2/bin/python examples/little_tour_with_bci_dataset.py
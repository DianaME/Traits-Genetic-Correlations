#!/bin/bash -l
#  FILENAME: gibbs_run.sh
cd /home/descamil/Multivariate/BLUPF90/opt1/

echo -e 'renf90.par\n 1200000\n 100000\n 20'|./gibbs3f90

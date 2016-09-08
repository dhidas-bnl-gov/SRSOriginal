#!/bin/bash

PREFIX=Test_qsub
LOG=$PREFIX.log

qsub -m abe -M dhidas@bnl.gov -v OUTPREFIX=$PREFIX -N $LOG -o $PWD -e $PWD -t 1-40 GridUtilities/run_grid_qsub.sh


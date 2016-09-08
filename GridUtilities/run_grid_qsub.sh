#!/bin/bash
SGE_ROOT=/nsls2/apps/sge62
export SGE_ROOT
. $SGE_ROOT/default/common/settings.sh


export PATH=/nsls2/projects/bldev/SoftDev/dhidas/Python-2.7.10/bin:$PATH
export PYTHONDIR=/nsls2/projects/bldev/SoftDev/dhidas/Python-2.7.10
export LD_LIBRARY_PATH=$PYTHONDIR/lib:$LD_LIBRARY_PATH

SECTION=$SGE_TASK_ID

cd $DIR


Example_Grid_MultiElectron_Flux.py $SECTION $OUTPREFIX

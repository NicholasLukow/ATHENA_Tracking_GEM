#!/bin/env bash

export EIC_LEVEL=dev
source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_bash.sh


#EDIT the paths below to point to your own installs
singularity shell -B /cvmfs:/cvmfs /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.simg
source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/eic_setup.sh -n
source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/setup_local.sh /eic/u/lukown/myinstall
export LD_LIBRARY_PATH=/eic/u/lukown/EicToyModel/build/lib:${OPT_SPHENIX}/vgm/lib64:${LD_LIBRARY_PATH}
source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/setup_local.sh /eic/u/lukown/EicToyModel/fun4all_with_eicroot
export ROOT_INCLUDE_PATH=/eic/u/lukown/EicToyModel/build/include/etm:${ROOT_INCLUDE_PATH}
source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/setup_local.sh /eic/u/lukown/myinstall


echo START
date

EVENTS=$1
PMIN=$2
PMAX=$3
ETAMIN=$4
ETAMAX=$5
GEN=$6
BFIELD=$7
NAME=$8
number=$9
BASEDIR=${10}
SOURCE_DIR=${11}

WRKDIR=${BASEDIR}

# change to Working Directory
cd $WRKDIR

# Run Simulation
root -b -l -q ${BASEDIR}/Fun4All_G4_HybridGEM.C\(${EVENTS},${PMIN},${PMAX},${ETAMIN},${ETAMAX},${GEN},${BFIELD},\"${NAME}${number}\"\)

date
echo DONE

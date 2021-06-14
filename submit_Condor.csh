#! /bin/env csh

# Output Directory
#currently have this set as the working directory since I dont think I can write to /gpfs02/eic/ from within the container
#set OUTPUT_DIR = /gpfs02/eic/${USER}/scratch/
set OUTPUT_DIR = /eic/u/lukown/Simulations/WorkingDirectory

if ($#argv > 0) then 
    set OUTPUT_DIR = $1
endif

set SOURCE_DIR = `pwd`
if ($#argv > 1) then 
    set SOURCE_DIR = $2
endif
set EXECUTABLE   = "SimuCondor.sh"

####### Initialize condor file
echo ""  > CondorFile
echo "Universe     = vanilla" >> CondorFile
echo "Executable   = ${EXECUTABLE}" >> CondorFile
echo "getenv = true" >> CondorFile

# Make sure structure exists
mkdir -p ${OUTPUT_DIR}/condor || exit

#JOBS
set NUMBER = 1
set LIMIT = 100

#Simulation Parameters
set EVENTS = 5000


set PMIN = 1
set PMAX = 30
set ETAMIN = -3.5
set ETAMAX = 3.5
set BFIELD = 5
set NAME = "4BARR"

while ( "$NUMBER" <= "$LIMIT" )

    set LogFile = ${OUTPUT_DIR}/condor/${NAME}.P=${PMIN}-${PMAX}.5kevents.BFIELD=${BFIELD}.ETA=${ETAMIN}-${ETAMAX}_${NUMBER}.out
    set ErrFile = ${OUTPUT_DIR}/condor/${NAME}.P=${PMIN}-${PMAX}.5kevents.BFIELD=${BFIELD}.ETA=${ETAMIN}-${ETAMAX}_${NUMBER}.err

   
    set Args = ( $EVENTS $PMIN $PMAX $ETAMIN $ETAMAX $BFIELD $NAME $NUMBER $OUTPUT_DIR $SOURCE_DIR)
    echo "" >> CondorFile
    echo "Output       = ${LogFile}" >> CondorFile
    echo "Error        = ${ErrFile}" >> CondorFile
    echo "Arguments    = ${Args}" >> CondorFile
    echo "Queue" >> CondorFile   

    echo Submitting:
    echo $EXECUTABLE $Args
    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile
    echo

    @ NUMBER++

end

condor_submit CondorFile



This contains the code for a hybrid tracking detector with GEMs. 
The current model uses the 2m4o BMT from Saclay, GEM Disks from the EIC Toy Model, and modified simple Si tracker (vertexing, barrel, and disks) based on Rey's code using parameters from Table 11.12 in the Yellow Report.
There is also a simple mid-rapidity Beryllium beam pipe as well as optional forward and backward beampipe extensions.

In order to make use of the new field map (option 5 in the .C script), you must retrieve it from here: https://github.com/eic/BeastMagneticField

In order to make use of the GEMs, you must first install the EicToyModel. This can be found here: https://github.com/eic/EicToyModel
Follow the install instructions there.

First open a shell in the singularity container using:
>singularity shell -B /cvmfs:/cvmfs -B /direct/eic+u/USERNAME/scratch:/scratch/ /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.simg


You must also install the MicroMegas simulation: https://github.com/hqh0127/EIC_MMStripCZ
Install this by going into the source directory, making a build folder, running the autogen.sh script and installing
> cd source

> mkdir build

> cd build

> ../autogen.sh --prefix=$HOME/myinstall

> make install

Once you have installed the necessary items, each time you open a new singularity shell you will need to execute a series of source statements (There may be redundant steps, but this sequence gets things working properly for me):
> source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/eic_setup.sh -n

> source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/setup_local.sh $HOME/myinstall

> export LD_LIBRARY_PATH=$HOME/EicToyModel/build/lib:${OPT_SPHENIX}/vgm/lib64:${LD_LIBRARY_PATH}

> source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/setup_local.sh $HOME/EicToyModel/fun4all_with_eicroot

> export ROOT_INCLUDE_PATH=$HOME/EicToyModel/build/include/etm:${ROOT_INCLUDE_PATH}

> source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/setup_local.sh $HOME/myinstall


A couple paths in the Fun4All_G4_HyrbidGem.C macro will need to be updated:
line ~176 where "EicGeoParData::ImportMediaFile" is found. Edit the path so it points to the EicToyModel you installed
line ~199 where the "ATHENA" magnetic field map is set. Edit the path so it points to where you have downloaded the the file from https://github.com/eic/BeastMagneticField

Then to test things are working, from the directory with the Fun4All_G4_HybridGem.C script:
> root

> root[]  .x Fun4All_G4_HybridGem.C(-1)

> root[] .L DisplayOn.C

> root[] PHG4Reco \*g4 = QTGui();

And you should see the visualization of the detector setup with a cutaway.

To add or remove parts of the detector, open the detector_setup.h file and comment out the corresponding definition.

There are two scripts to allow for job submission:

SimuCondor.sh:

	-requires a one-time change: update the "source" lines to point to your installation

submit_condor.csh:

	-change the OUTPUT_DIR to your working directory where you have the Fun4All_G4_HybridGEM.C script and the other .C scripts

	-LIMIT and EVENTS are used to set the number of jobs and number of events per job

	-NAME is the name you wish to give the output .root files (the job number will be appended to this)

	-the rest are self explanatory inputs to the simulation

If you wish to submit jobs, you do not need to open a singularity shell and perform the source commands as written (this is taken care of inside the job).
Just set the parameters, edit the Fun4All_G4_HybridGEM.C and detector_setup.h scripts for the detector configuration you want, and the submit the jobs via:

> csh submit_condor.csh

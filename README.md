In order to make use of the GEMs, you must first install the EicToyModel. This can be found here: https://github.com/eic/EicToyModel
Follow the install instructions there.


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

Then to test things are working, from the directory with the Fun4All_G4_HybridGem.C script:
> root
[] .x Fun4All_G4_HybridGem.C(-1)
[] .L DisplayOn.C
[] PHG4Reco *g4 = QTGui();

And you should see the visualization of the detector setup with a cutaway.

To add or remove parts of the detector, open the detector_setup.h file and comment out the corresponding definition.

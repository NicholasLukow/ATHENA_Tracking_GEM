#ifndef MACRO_G4_DUMMYMRICH_C
#define MACRO_G4_DUMMYMRICH_C

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4main/PHG4Reco.h>
#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

// ======================================================================================================================
void add_dummymRICH(PHG4Reco *g4Reco, double R, double zpos){

	double thinness = 15./100.*8.897;  // black hole thickness, needs to be taken into account for the z positions
	double hole = 0;	// If beampipe needs to be fed through a hole
	PHG4CylinderSubsystem * cyl_mRICH;
	cyl_mRICH = new PHG4CylinderSubsystem("DUMMYMRICH", 1);
	cyl_mRICH->set_double_param("length"   , thinness );
	cyl_mRICH->set_double_param("radius"   , hole          );
	cyl_mRICH->set_double_param("thickness", R - hole   );
	cyl_mRICH->set_string_param("material" , "G4_Al"   );
	cyl_mRICH->set_double_param("place_z"  , zpos + thinness/2);
	cyl_mRICH->SuperDetector("DUMMYMRICH");
	cyl_mRICH->set_color(0.9,0.5,1,0.5);
	//cyl_mRICH->BlackHole();
	g4Reco->registerSubsystem(cyl_mRICH);

}
#endif

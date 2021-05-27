#ifndef MACRO_G4GEM_C
#define MACRO_G4GEM_C

#include "GlobalVariables.C"
#include <Geant4/G4Material.hh>
#include <Geant4/G4Element.hh>
#include <g4detectors/PHG4SectorSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

	int make_GEM_station(string name, PHG4Reco *g4Reco, double zpos, double etamin, double etamax, const int N_Sector = 8);
	void AddLayers_MiniTPCDrift(PHG4SectorSubsystem *gem);
	void AddLayers_Tracker(PHG4SectorSubsystem *gem);

	namespace Enable
{
	bool EGEM = true;
	bool FGEM = true;
}  // namespace Enable
// ======================================================================================================================
void EGEM_Init(){
	BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 80.);
	// extends only to -z
	BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -160.);
}
// ======================================================================================================================
void FGEM_Init(){
	BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 150.);
	BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 282.);
}
// ======================================================================================================================
void EGEMSetup(PHG4Reco *g4Reco){
//make_GEM_station("", , z, eta_min, eta_max)
	/* Careful with dimensions! If GEM station volumes overlap, e.g. with TPC volume, they will be drawn in event display but will NOT register any hits.
	 * Geometric constraints:
	 * TPC length = 211 cm --> from z = -105.5 to z = +105.5
	 */
	float thickness = 3.;	
	make_GEM_station("EGEM1", g4Reco, -180.0 , -1.15 , -3.9,12);
	//make_GEM_station("EGEM2", g4Reco, -180.0 , -1.15 , -3.9,12);
	//make_GEM_station("EGEM3", g4Reco, -180.0 , -1.15 , -3.9,12);
}
// ======================================================================================================================
void FGEMSetup(PHG4Reco *g4Reco){
	float thickness = 3.;
	make_GEM_station("FGEM1", g4Reco, 180.0 , 1.15, 3.9,12);
	make_GEM_station("FGEM2", g4Reco, 187.0 , 1.15, 3.9,12);
	make_GEM_station("FGEM3", g4Reco, 194.0 , 1.15, 3.9,12);
	
	make_GEM_station("FFGEM1", g4Reco, 307.0 , 1.15, 3.9,12);
	make_GEM_station("FFGEM2", g4Reco, 314.0 , 1.15, 3.9,12);
	make_GEM_station("FFGEM3", g4Reco, 321.0 , 1.15, 3.9,12);
}
// ======================================================================================================================
int make_GEM_station(string name, PHG4Reco *g4Reco, double zpos, double etamin,
		double etamax, const int N_Sector = 8)
{

	double polar_angle = 0;

	if (zpos < 0){
		zpos = -zpos;
		polar_angle = M_PI;
	}
	if (etamax < etamin){
		double t = etamax;
		etamax = etamin;
		etamin = t;
	}

	PHG4SectorSubsystem *gem;
	gem = new PHG4SectorSubsystem(name);

	gem->SuperDetector(name);

	gem->get_geometry().set_normal_polar_angle(polar_angle);
	gem->get_geometry().set_normal_start(zpos * PHG4Sector::Sector_Geometry::Unit_cm());
	gem->get_geometry().set_min_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamax));
	gem->get_geometry().set_max_polar_angle(PHG4Sector::Sector_Geometry::eta_to_polar_angle(etamin));
	gem->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	gem->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	gem->get_geometry().set_N_Sector(N_Sector);
	gem->get_geometry().set_material("G4_AIR");
	gem->OverlapCheck(Enable::OVERLAPCHECK);

	AddLayers_Tracker(gem);
	g4Reco->registerSubsystem(gem);

	return 0;
}
// ======================================================================================================================
// Component Material X0 (cm) Thickness (cm) Area (%) Rad. Length (%)
void AddLayers_Tracker(PHG4SectorSubsystem *gem)
{
	assert(gem);

	const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
	const double mm = 0.1 * cm;
	const double um = 1e-3 * mm;

        const int n_GEM_layers = 3;
        const string gem_gas = "G4_Ar";
        const double frame_fill_area = 6.5; //assume frame covers 6.5% of GEM area. 
        const double gem_fill_area = 80.0; //assume 20% of GEM material area is missing (holes). 
      
        //Front
	gem->get_geometry().AddLayer("Frame", "G10", 1.5 * mm, false, frame_fill_area);
	gem->get_geometry().AddLayer("EntranceWindow", "G4_MYLAR", 25 * um, false, 100);
	gem->get_geometry().AddLayer("Frame", "G10", 1.5 * mm, false, frame_fill_area);
	gem->get_geometry().AddLayer("WindowGas", gem_gas, 1.5 * mm, false, 100);
	gem->get_geometry().AddLayer("Cathode", "G4_KAPTON", 50 * um, false, 100);
	gem->get_geometry().AddLayer("Frame", "G10", 3 * mm, false, frame_fill_area);
	gem->get_geometry().AddLayer("DriftGas", gem_gas, 3 * mm, false, 100); //drift gap

        for(int igem = 1; igem <= n_GEM_layers; igem++)
          {
	    gem->get_geometry().AddLayer(Form("GEMFrontCu%d",igem), "G4_Cu",5 * um, false, gem_fill_area);
	    gem->get_geometry().AddLayer(Form("GEMKapton%d",igem), "G4_KAPTON",50 * um, false, gem_fill_area);
	    gem->get_geometry().AddLayer(Form("GEMBackCu%d",igem), "G4_Cu",5 * um, false, gem_fill_area);
	    gem->get_geometry().AddLayer(Form("Frame%d",igem), "G10", 2 * mm, false, frame_fill_area);
	    gem->get_geometry().AddLayer(Form("GEMGas%d",igem), gem_gas, 2 * mm, false, 100); //1st, 2nd, and induction gaps
          }

        //readout
        //

	gem->get_geometry().AddLayer("ReadoutKapton", "G4_KAPTON", 50 * um, false, 100);
	gem->get_geometry().AddLayer("ReadoutCu", "G4_Cu", 5 * um, false, 100);
	gem->get_geometry().AddLayer("ExitWindowGas", gem_gas, 1.5 * mm, false, 100);
	gem->get_geometry().AddLayer("Frame", "G10", 1.5 * mm, false, frame_fill_area);
	gem->get_geometry().AddLayer("ExitWindow", "G4_MYLAR", 25 * um, false, 100);
	gem->get_geometry().AddLayer("Frame", "G10", 1.5 * mm, false, frame_fill_area);
}
#endif

//For tracking (end cap resolutions)
/*
 kalman->add_phg4hits(
 "G4HIT_EGEM_1",                //      const std::string& phg4hitsNames,
  PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
  400e-4,                      //       radial-resolution [cm]   // 400 um based on eRD6 EIC prototype, may get down to 250 um
  100e-4,                       //        azimuthal-resolution [cm] //100 um based on eRD6 EIC prototype, may get down to 50 um
  999.0,                           //      z-resolution [cm] //not used in vert. plane
  1,                           //      efficiency,
  0                            //      noise hits
 );
*/

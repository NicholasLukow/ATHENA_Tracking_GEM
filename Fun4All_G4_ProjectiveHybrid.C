/*
================================================================================================================

================================================================================================================
*/
#pragma once
#include "/eic/u/lukown/Simulations/WorkingDirectory/detector_setup.h"
#include <phgenfit/Track.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4histos/G4HitNtuple.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/HepMCNodeReader.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>
#include <g4lblvtx/PHG4ParticleGenerator_flat_pT.h>
//#include <g4lblvtx/AllSi_Al_support_Subsystem.h>
//#include <g4lblvtx/EicFRichSubsystem.h>
#include <g4_modifieddrich/dRICH_Subsystem.h>
#include "G4_BlackHole.C"


#include <g4alsupp/Al_support_Subsystem.h>
#include <gem_g4alsupp/GEM_support_Subsystem.h>



#include <g4mrich/PHG4mRICHSubsystem.h>


#include "G4_Pipe_EIC.C"

#include "G4_DIRC_SMALL.C"

#include "G4_GridPix.C"

#ifdef _EICTOYVST_
#include <EicRootVstSubsystem.h>
#include <EtmOrphans.h>
#endif


//gdml
#include <gdmlimporter/GdmlImportDetectorSubsystem.h>
#include <gdmlimporter/SimpleNtuple.h>
#include <gdmlimporter/TrackFastSimEval.h>
#include <g4detectors/PHG4GDMLSubsystem.h>

#include "G4_GEM_EIC_v1.C"
//to include gems                                                                                                                                     
#include <EicToyModelSubsystem.h>
#include <EicRootGemSubsystem.h>
#include <EtmOrphans.h>

//For MPGD from QH                                                                                                                                    
#include <g4exampledetector/PHG4CylinderStripSubsystem.h>
#include <g4exampledetector/CreateCZHitContainer.h>

#include <TrackFastSimEval.h>

R__LOAD_LIBRARY(libeicdetectors.so)

R__LOAD_LIBRARY(libgdmlimportdetector.so)

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)


R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libg4example01detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libPHPythia6.so)


//Function to make GEM disk
void MakeGEM(array<double,6> Params, EicRootGemSubsystem *&fgt)
{
	    //auto sbs = ModifiedGEM();// creates GEM Module with modified material
	    auto sbs = new GemModule();// creates GEM Module with modified material
	    sbs->SetDoubleVariable("mDriftFoilCopperThickness", 5 * etm::um);
	    sbs->SetDoubleVariable("mGemFoilCopperThickness", 5 * etm::um);
	    sbs->SetDoubleVariable("mGemFoilKaptonThickness", 50 * etm::um);
	    sbs->SetDoubleVariable("mReadoutSupportThickness", 0 * etm::um);
	    sbs->SetDoubleVariable("mReadoutKaptonThickness", 50 * etm::um);
	    sbs->SetDoubleVariable("mFrameThickness", 17 * etm::mm);
	    sbs->SetDoubleVariable("mFrameBottomEdgeWidth", 10 * etm::mm); //30
	    sbs->SetDoubleVariable("mFrameTopEdgeWidth", 50 * etm::mm);
	    sbs->SetDoubleVariable("mFrameSideEdgeWidth", 10 * etm::mm); //15 
            sbs->SetDoubleVariable("mEntranceWindowThickness", 50 * etm::um); //25 um for both the entrance and exit
	    sbs->SetDoubleVariable("mActiveWindowBottomWidth", Params[3] * etm::mm);
	    sbs->SetDoubleVariable("mActiveWindowTopWidth", Params[2] * etm::mm);
	    sbs->SetDoubleVariable("mActiveWindowHeight", Params[0] * etm::mm);
	    fgt->AddWheel(sbs, Params[5], Params[1] * etm::mm, Params[4] * etm::mm, 0);
}

//Function to calculate the parameters for GEM disk geometry given the Z position, minimum eta covered, inner radius clearance, and the number of modules
array<double,6> FullGEMParameters(double Z, double EtaMin, double InnerRadius, double NModules)
{
	    double Height = -50 + TMath::Abs(Z)*TMath::Tan(2*TMath::ATan(TMath::Exp(-1*EtaMin))); 
	    double ActiveHeight = Height - InnerRadius; //NB: this may not account for the frame material thickness at the outer edge of the module, so it may not precisely get the full eta coverage
	    double CenterRadius = 0.5*ActiveHeight + InnerRadius;
	    double TopWidth = 2*Height*TMath::Tan(TMath::Pi()/NModules);
	    double BottomWidth = (TopWidth/Height)*InnerRadius;

	    array<double,6 > Params = {ActiveHeight, CenterRadius, TopWidth, BottomWidth, Z, NModules};
	    return Params;
}

double RadiusFromZEta(double Z, double Eta)
{
	double Radius = TMath::Abs(Z)*TMath::Tan(2*TMath::ATan(TMath::Exp(-1*Eta))); 
	return Radius;
}

void Fun4All_G4_ProjectiveHybrid(
			int nEvents = -1,			// number of events
			double pmin = 1., 			// GeV/c
			double pmax = 30., 			// GeV/c
			double etamin = -3.7,
			double etamax = 3.7,
			int generatorVersion = 1, 		// Generator setting
			int magnetic_field = 6, 		// Magnetic field setting
			TString out_name = "TEST")	// output filename
{	
	// ======================================================================================================
  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors");
  gSystem->Load("libg4testbench");
  gSystem->Load("libg4histos");
  gSystem->Load("libg4example01detector.so");
  gSystem->Load("libg4trackfastsim.so");
	// Input from the user
	const int particle_gen = generatorVersion;     // 1 = particle generator, 2 = particle gun, 3 = simple event generator, 4 = pythia8 e+p collision, 5 = particle generator flat in pT
	double pix_size_vtx = 10.; // um - size of pixels in vertexing layers
	double pix_size_bar = 10.; // um - size of pixels in barrel layers
	double pix_size_dis = 10.; // um - size of pixels in disk layers
	const int nDisks_per_side = 5;
	const int do_projections =  1;
	// Parameters for projections
	string projname1   = "DIRC";            // Cylindrical surface object name
	double projradius1 = 95;// 100;//112;// 80.;               // [cm] 
	//NOTE: these surfaces are black holes. Care must be taken in the choice of dimensions as to not absorb particles within the acceptance of other detectors
	double length1     = 310; //200.;              // [cm]
	// ---
	double thinness    = 0.01;               // black hole thickness, needs to be taken into account for the z positions
	// ---
	string projname2   = "FOR";             // Forward plane object name
	double projzpos2   = 154+thinness/2;//315+thinness/2.;   // [cm]
	double projradius2 = 110;//210.;               // [cm]
	// ---
	string projname3   = "BACK";            // Backward plane object name
	double projzpos3   = -(136+thinness/2.);// [cm]
	double projradius3 = 95.;               // [cm]
	// ---
	string projname4   = "FOREXIT";            // Backward plane object name
	double projzpos4   = 342;// [cm]
	double projradius4 = 200.;               // [cm]
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// If you want to fix the random seed for reproducibility
	 //recoConsts *rc = recoConsts::instance();
	 //rc->set_IntFlag("RANDOMSEED", 12345);
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string("pi-"));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);			// Vertex generation range
	gen->set_mom_range(pmin,pmax);		// Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(etamin,etamax);//4.0
	gen->set_phi_range(0,2.*TMath::Pi());
	// --------------------------------------------------------------------------------------
	// Particle generator flat in pT
	PHG4ParticleGenerator_flat_pT *gen_pT = new PHG4ParticleGenerator_flat_pT();
	gen_pT->set_name(std::string("pi-"));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen_pT->set_vtx(0,0,0);                    // Vertex generation range
	gen_pT->set_pT_range(pmin,pmax);         // Momentum generation range in GeV/c
	gen_pT->set_z_range(0.,0.);
	gen_pT->set_eta_range(etamin, etamax);               // Detector coverage corresponds to |η|< 4
	gen_pT->set_phi_range(0.,2.*TMath::Pi());
	// ======================================================================================================
	if     (particle_gen==1){se->registerSubsystem(  gen); cout << "Using particle generator"     << endl;}
	else if(particle_gen==5){se->registerSubsystem(gen_pT); cout << "Using particle generator flat in pT"  << endl;}
	else{ cout << "Particle generator option requested has not been implemented. Bailing out!" << endl; exit(0); }
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	//g4Reco->SetWorldMaterial("G4_Galactic");	
        EicGeoParData::ImportMediaFile("/eic/u/lukown/EicToyModel/examples/eicroot/media.geo");
	// ======================================================================================================
	// Magnetic field setting
	TString B_label;
	if(magnetic_field==1){          // uniform 1.5T
		B_label = "_B_1.5T";
		g4Reco->set_field(1.5);
	}
	else if(magnetic_field==2){     // uniform 3.0T
		B_label = "_B_3.0T";
		g4Reco->set_field(3.0);
	}
	else if(magnetic_field==3){     // sPHENIX 1.4T map
		B_label = "_sPHENIX";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
		g4Reco->set_field_rescale(-1.4/1.5);
	}
	else if(magnetic_field==4){     // Beast 3.0T map
		B_label = "_Beast";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);
	}
	else if(magnetic_field==5){     // ATHENA 3.0T map (updated BeAST Field map with dual coil configuration)
		B_label = "_ATHENA0528";
		g4Reco->set_field_map( string("../BeastMagneticField/data/EIC_v.0.1.0_Magnetic_Field_Map_2021_05_28_radial_coords_[cm]_[T].401301.line.Bmap"), PHFieldConfig::kFieldBeast);
	}
	else if(magnetic_field==6){     // updated BeAST 3.0T map
		B_label = "_ATHENA0507";
		g4Reco->set_field_map( string("../BeastMagneticField/data/EIC_Magnetic_Field_Map_2021_05_07_radial_coords_[cm]_[T].120000.lines.Bmap"), PHFieldConfig::kFieldBeast);
	}
	else{                           // The user did not provide a valid B field setting
		cout << "User did not provide a valid magnetic field setting. Set 'magnetic_field'. Bailing out!" << endl;
	}	








	// ======================================================================================================
	// Detector setup
	PHG4CylinderSubsystem *cyl;
	
	if(do_projections){
		PHG4CylinderSubsystem *cyl;
		cyl = new PHG4CylinderSubsystem(projname1);
		cyl->set_double_param("length", length1);
		cyl->set_double_param("radius", projradius1); // dirc radius
		cyl->set_double_param("thickness", 0.01); // needs some thickness
		cyl->set_string_param("material", "G4_AIR");
		//cyl->set_string_param("material", "G4_Galactic");
		cyl->SetActive();
		cyl->SaveAllHits();
		cyl->SuperDetector(projname1);
		cyl->set_color(1,0,0,0.7); //reddish
		g4Reco->registerSubsystem(cyl);

		cyl = new PHG4CylinderSubsystem(projname2);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 10); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius2-10); // 
		cyl->set_string_param("material", "G4_AIR");
		//cyl->set_string_param("material", "G4_Galactic");
		cyl->set_double_param("place_z", projzpos2);
		cyl->SetActive();
		cyl->SaveAllHits();
		cyl->SuperDetector(projname2);
		cyl->set_color(0,1,1,0.3); 
		g4Reco->registerSubsystem(cyl);

		cyl = new PHG4CylinderSubsystem(projname3);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 10); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius3-10); // 
		cyl->set_string_param("material", "G4_AIR");
		//cyl->set_string_param("material", "G4_Galactic");
		cyl->set_double_param("place_z", projzpos3);
		cyl->SetActive();
		cyl->SaveAllHits();
		cyl->SuperDetector(projname3);
		cyl->set_color(0,1,1,0.3);
		g4Reco->registerSubsystem(cyl);
		
		cyl = new PHG4CylinderSubsystem(projname4);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 15); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius4-15); // 
		cyl->set_string_param("material", "G4_AIR");
		//cyl->set_string_param("material", "G4_Galactic");
		cyl->set_double_param("place_z", projzpos4);
		cyl->SetActive();
		cyl->SaveAllHits();
		cyl->SuperDetector(projname4);
		cyl->set_color(0,1,1,0.3); 
		g4Reco->registerSubsystem(cyl);

	}

	#ifdef _DIRC_
	
	DIRCSetup(g4Reco);

	#endif

	#ifdef _TPC_
	GridPixSetup(g4Reco);
	#endif


	#ifdef _SIVTX_
	//---------------------------
	// Vertexing
	double si_vtx_r_pos[] = {3.64,4.45,5.26};
	const int nVtxLayers = sizeof(si_vtx_r_pos)/sizeof(*si_vtx_r_pos);
	double si_z_vtxlength[] = {42.0, 42.0, 42.0};
	//double si_thick_vtx = vtx_matBud/100.*9.37;
	double si_thick_vtx = 0.05/100.*9.37;

	for (int ilayer = 0; ilayer < nVtxLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
		cyl->set_string_param("material" , "G4_Si"               );
		cyl->set_double_param("radius"   , si_vtx_r_pos[ilayer]  );
		cyl->set_double_param("thickness", si_thick_vtx          );
		cyl->set_double_param("place_z"  , 0                     );
		cyl->set_double_param("length"   , si_z_vtxlength[ilayer]);
		cyl->SetActive();
		cyl->SuperDetector("SVTX");
		cyl->set_color(0,0.8,0.1);
		g4Reco->registerSubsystem(cyl);
	}
	#endif

	#ifdef _SIBARR_
	//---------------------------
	// Barrel
	double si_r_pos[] = {13.38, 18.0};
	const int nTrckLayers = sizeof(si_r_pos)/sizeof(*si_r_pos);
        
	//projective Si Barrel
	double barr_prapidity = 1.1;
	double si_barr_length_1 = (1-exp(-2*barr_prapidity))/exp(-barr_prapidity)*si_r_pos[0];
	double si_barr_length_2 = (1-exp(-2*barr_prapidity))/exp(-barr_prapidity)*si_r_pos[1];
	double si_z_length[] = {si_barr_length_1, si_barr_length_2};

	cout << "LENGTHS OF SI BARR: " << si_barr_length_1 << " , " << si_barr_length_2 << endl;
	
	//double si_z_length[] = {84.0, 84.0};
	//double si_thick_bar = barr_matBud/100.*9.37;
	double si_thick_bar = 0.55/100.*9.37;

	for (int ilayer = 0; ilayer < nTrckLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("BARR", ilayer);
		cyl->set_string_param("material" , "G4_Si"            );
		cyl->set_double_param("radius"   , si_r_pos[ilayer]   );
		cyl->set_double_param("thickness", si_thick_bar       );
		cyl->set_double_param("place_z"  , 0                  );
		cyl->set_double_param("length"   , si_z_length[ilayer] - 4);
		cyl->SetActive();
		cyl->SuperDetector("BARR");
		cyl->set_color(0,0.5,1);
		g4Reco->registerSubsystem(cyl);	
	}
	#endif

	#ifdef _SIDISKS_
	//---------------------------
	// Disks
/*	
	double si_z_pos[] = {-121.,-105.4,-89.8,-74.2,-58.6,-43.0,-22.0,22.0,43.0,58.6,74.2,89.8,105.4, 121.};
 	const int nDisks = sizeof(si_z_pos)/sizeof(*si_z_pos);
 	double si_r_max[] = {19.0, 19.0, 19.0, 19.0, 19.0, 13.94, 7.13, 7.13, 13.94, 19.0, 19.0, 19.0, 19.0, 19.0};
 	double si_r_min[] = {9.93, 8.35, 6.67, 4.99, 3.64, 3.64, 3.64, 3.64, 3.64, 3.64, 4.99, 6.67, 8.35, 9.93};
 */
  	//double si_z_pos[] = {-121., -96.25, -71.5, -46.75, -22.0, 22.0, 46.75, 71.5, 96.25, 121.};
  	
	//expanded spacing 
	//double si_z_pos[] = {-143., -115.5, -88., -60.5, -33, 33, 60.5, 88, 115.5, 143 };
	
	//Asymmetric expanded spacing 
	double si_z_pos[] = {-135., -109.5, -84., -58.5, -33, 33, 63., 93., 123., 153 };
	
	//original projective
	//double si_z_pos[] = {-121., -99., -77., -55., -33.0, 33.0, 55., 77., 99., 121.};
  	const int nDisks = sizeof(si_z_pos)/sizeof(*si_z_pos);
  	//double si_r_max[] = {19.0, 19.0, 19.0, 19.0, 7.13, 7.13, 19.0, 19.0, 19.0, 19.0};
  	double si_r_max[nDisks];
	for (int iDisk = 0 ; iDisk < nDisks; iDisk++)
	{
		si_r_max[iDisk] = RadiusFromZEta(TMath::Abs(si_z_pos[iDisk]), 2.5);
		cout << "r max for disk " << iDisk + 1 << " : " << si_r_max[iDisk] << endl;
		//if (iDisk == 4 || iDisk == 5) si_r_max[iDisk] = 7.13;
	}

	//projective Si Disk (only first disk changes)
	//double si_r_max[] = {19.0, 19.0, 19.0, 19.0, 16.47, 16.47, 19.0, 19.0, 19.0, 19.0};
  	
	double si_r_min[] = {9.93, 7.25, 4.65, 3.64, 3.64, 3.64, 3.64, 4.65, 7.25, 9.93};
  	double si_thick_disk = 0.24/100.*9.37;

 	for (int ilayer = 0; ilayer < nDisks ; ilayer++){
 		cyl = new PHG4CylinderSubsystem("FBVS", ilayer);
 		cyl->set_string_param("material" , "G4_Si"         );
 		cyl->set_double_param("radius"   , si_r_min[ilayer]);
 		cyl->set_double_param("thickness", si_r_max[ilayer]-si_r_min[ilayer]);
 		cyl->set_double_param("place_z"  , si_z_pos[ilayer]);
 		cyl->set_double_param("length"   , si_thick_disk   );
 		cyl->SetActive();
		cyl->SuperDetector("FBST");
		cyl->set_color(1,0,0);
		g4Reco->registerSubsystem(cyl);
	
	}
	#endif
		
	#ifdef _ALSUPP_
	// Al Support Structure
	Al_support_Subsystem *Al_supp = new Al_support_Subsystem("Al_supp");
	g4Reco->registerSubsystem(Al_supp);	
	//Experimental GEM Support
	GEM_support_Subsystem *GEM_supp = new GEM_support_Subsystem("GEM_supp");
	g4Reco->registerSubsystem(GEM_supp);	
	#endif
	// ------------	

	#ifdef _SIMPLEVST_
	//---------------------------
	// Barrel

	double si_r_pos[] = { 7.0, 14.0,  21, 28.0}; //Simplified Toy Model
	const int nTrckLayers = sizeof(si_r_pos)/sizeof(*si_r_pos);
	double si_z_length[] = {27.0, 27.0, 54.0, 54.0}; //Simplified Toy Model
	double si_thick_bar = 0.55/100.*9.37;
	


	for (int ilayer = 0; ilayer < nTrckLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("BARR", ilayer);
		cyl->set_string_param("material" , "G4_Si"            );
		cyl->set_double_param("radius"   , si_r_pos[ilayer]   );
		cyl->set_double_param("thickness", si_thick_bar       );
		cyl->set_double_param("place_z"  , 0                  );
		cyl->set_double_param("length"   , si_z_length[ilayer]);
		cyl->SetActive();
		cyl->SuperDetector("SIMPLEVST");
		cyl->set_color(0,0.5,1);
		g4Reco->registerSubsystem(cyl);	
	}
	#endif	
	
	//---------------------------
	// Black hole to suck loopers out of their misery
	#ifdef _BLACKHOLE_
	//Full Wrap single Cylinder
	wrap_with_cylindrical_blackhole(g4Reco,250,-250,350,true);
	#endif

	#ifdef _RICH_
	//EicFRichSubsystem *RICH = new EicFRichSubsystem("RICH");
	dRICH_Subsystem *RICH = new dRICH_Subsystem("RICH");
	g4Reco->registerSubsystem(RICH);
	#endif
	
	#ifdef _MRICH_
	PHG4mRICHSubsystem *mRICH = new PHG4mRICHSubsystem("mRICH", 0);
	g4Reco->registerSubsystem(mRICH);
	#endif

	//beampipe as implemented by Rey with slight modification
	#ifdef _BEAMPIPE_
	PipeInit();
	double pipe_radius=0;	
	pipe_radius = Pipe(g4Reco, pipe_radius , true);	
	#endif


	// ------------
        #ifdef _MPGD_
        double gap_betweenCZ = 1.5, Gap_betweenlayer = 1.5;
        //double thickness = 0.355199;                                                                                                                
        double thicknessMPGD = 0.36499;
 	int nCZlayer = 2;
        bool use_2Dreadout = true;
        if (use_2Dreadout) {
                gap_betweenCZ = 0;
                nCZlayer = 1;
        }
  //double BMT_r[6] = {20., 20.+nCZlayer*thicknessMPGD+gap_betweenCZ+Gap_betweenlayer, 50-nCZlayer*thicknessMPGD-gap_betweenCZ-Gap_betweenlayer/2, 50\
+Gap_betweenlayer/2, 80-(nCZlayer*thicknessMPGD+gap_betweenCZ)*2-Gap_betweenlayer, 80-nCZlayer*thicknessMPGD-gap_betweenCZ};                          
  	double BMT_r[6] = {    47.7153, 49.5718, 71.8958, 73.7523, 75.6088, 77.4653  };
  	
	for (int iScale = 0; iScale < 6; iScale++)
	{
		BMT_r[iScale] = BMT_r[iScale]*1.2;
	}

	//double BMT_r[4] = {   71.8958, 73.7523, 75.6088, 77.4653  };


        PHG4CylinderStripSubsystem *example01;
        //double bmt_length = (1-exp(-2*prapidity))/exp(-prapidity)*80;
        
	const double prapidity =1.1;
        //Original projective
	double bmt_length_inner = (1-exp(-2*prapidity))/exp(-prapidity)*BMT_r[0] ;
        double bmt_length_outer = (1-exp(-2*prapidity))/exp(-prapidity)*BMT_r[2];
        
	//scale up
	//double bmt_length_inner = (1-exp(-2*prapidity))/exp(-prapidity)*50*1.226 - 10;
        //double bmt_length_outer = (1-exp(-2*prapidity))/exp(-prapidity)*78*1.226 - 17;
        double bmt_length;
        
	//const double prapidity =1.0;
        //double bmt_length_inner, bmt_length_outer;
	//double bmt_length = bmt_length_inner = bmt_length_outer = (1-exp(-2*prapidity))/exp(-prapidity)*80;
	//double bmt_length = 250;
        
	for (int ilayer = 0; ilayer< 6; ilayer++){
                example01 = new PHG4CylinderStripSubsystem("BMT",ilayer);
                example01->set_double_param("radius", BMT_r[ilayer]);
                example01->set_string_param("gas", "myMMGas");
                //example01->set_double_param("steplimits", 300e-4);                                                                                  
                example01->set_double_param("phi0", 15*ilayer);
                example01->set_double_param("gap", gap_betweenCZ);
                example01->SetActive();
                example01->SuperDetector("BMT");
                example01->set_int_param("lengthviarapidity",0);
                
		if (ilayer < 2) bmt_length = bmt_length_inner;
		else bmt_length = bmt_length_outer;
		//bmt_length = (1-exp(-2*prapidity))/exp(-prapidity)*BMT_r[ilayer];
		
		example01->set_double_param("length", bmt_length);
                example01->set_double_param("deadzone", 0.2);
                example01->set_int_param("nhit", 2);
                example01->OverlapCheck(true);
                example01->set_int_param("use_2Dreadout",use_2Dreadout);
                g4Reco->registerSubsystem(example01);
                //example01->Print();                                                                                                                 
        }
        #endif

#ifdef _INNERGEMS_
        // Forward GEM tracker module(s);  
        auto fgt = new EicRootGemSubsystem("INNERGEM");
        {
	 	 fgt->SetActive(true);
         	{
	 	//fgt->CheckOverlap();                                                                                                                     \
          	//fgt->SetTGeoGeometryCheckPrecision(0.000001 * etm::um);                                                                                  \
            	// See other GemModule class data in GemGeoParData.h;                                                                                    \
            	// Compose sectors; parameters are:                                                                                                      \
            	//   - layer description (obviously can mix different geometries);                                                                       \
            	//   - azimuthal segmentation;                                                                                                            
            	//   - gas volume center radius;                                                                                                         \
            	//   - Z offset from 0.0 (default);                                                                                                      \
            	//   - azimuthal rotation from 0.0 (default);                                                                                            \
		
		//FullGEMParameters() will calculate the parameters to define the geometry of the GEM disk based on the Z location, minimum eta coverage, inner radius clearance, and number of modules
	        //Array definition: Params[] = {ActiveHeight, CenterRadius, TopWidth, BottomWidth, Z, NModules};
	   
		//New Design for projective Central tracker with eta = 1.1
	    	double GEM_Z[3] = {40, 85, 130};
	    	//double GEM_Z[4] = {40, 70, 100, 130};
		
		array<double,6> Params = FullGEMParameters(10*GEM_Z[0], 1.15, 10*(RadiusFromZEta(GEM_Z[0], 2.5)+2), 12);
		MakeGEM(Params, fgt);
		Params = FullGEMParameters(10*GEM_Z[1], 1.15, 10*(RadiusFromZEta(GEM_Z[1], 2.5)+2), 12);
		MakeGEM(Params, fgt);
		Params = FullGEMParameters(10*GEM_Z[2], 1.15, 10*(RadiusFromZEta(GEM_Z[2], 2.5)+2), 12);
		MakeGEM(Params, fgt);
		//Params = FullGEMParameters(10*GEM_Z[3], 1.15, 10*(RadiusFromZEta(GEM_Z[3], 2.5)+2), 12);
		//MakeGEM(Params, fgt);
	    

			
		Params = FullGEMParameters(-10*GEM_Z[0], 1.15, 10*(RadiusFromZEta(GEM_Z[0], 2.5)+2), 12);
		MakeGEM(Params, fgt);
		Params = FullGEMParameters(-10*GEM_Z[1], 1.15, 10*(RadiusFromZEta(GEM_Z[1], 2.5)+2), 12);
		MakeGEM(Params, fgt);
		Params = FullGEMParameters(-10*GEM_Z[2], 1.15, 10*(RadiusFromZEta(GEM_Z[2], 2.5)+2), 12);
		MakeGEM(Params, fgt);
		//Params = FullGEMParameters(-10*GEM_Z[3], 1.15, 10*(RadiusFromZEta(GEM_Z[3], 2.5)+2), 12);
		//MakeGEM(Params, fgt);
		
	
		}
 	        g4Reco->registerSubsystem(fgt);
	}
#endif

#ifdef _OUTERGEMS_
        // Forward GEM tracker module(s);  
        auto fgt2 = new EicRootGemSubsystem("OUTERGEM");
        {
	 	fgt2->SetActive(true);
         	{
		//FullGEMParameters() will calculate the parameters to define the geometry of the GEM disk based on the Z location, minimum eta coverage, inner radius clearance, and number of modules
	        //Array definition: Params[] = {ActiveHeight, CenterRadius, TopWidth, BottomWidth, Z, NModules};

	    	//Far Hadron Side GEM disk
	    	array<double,6> Params = FullGEMParameters(3450, 1.2, 210, 12);
	    	MakeGEM(Params, fgt2);
	    	Params[4]=Params[4]+50; //Copying previous parameters but shifting in Z
	    	//MakeGEM(Params, fgt2); 

	    	//Far Electron Side GEM disk
	    	Params = FullGEMParameters(-1300, 1.2, 110, 12);
	    	//MakeGEM(Params, fgt2);
			
			cout << "Top Width: " << Params[2] << endl;
       			cout << "Length: " << Params[0] + 210 << endl; 
		}
          	g4Reco->registerSubsystem(fgt2);
	}

#endif


  // Detailed vertex Si tracker from EicToyModel
  // EicRoot vertex tracker; be aware: "VST" will also become a SuperDetector name;
  #ifdef _EICTOYVST_
  auto vst = new EicRootVstSubsystem("VST");
  {
    vst->SetGeometryType(EicGeoParData::NoStructure);
    vst->SetActive(true);

    // Barrel layers; hits belonging to these layers will be labeled internally
    // according to the sequence of these calls;
    {
      auto ibcell = new MapsMimosaAssembly();
      // See other MapsMimosaAssembly class POD entries in MapsMimosaAssembly.h;
      ibcell->SetDoubleVariable("mAssemblyBaseWidth", 17.5 * etm::mm);
  
      // Compose barrel layers; parameters are:
      //  - cell assembly type;
      //  - number of staves in this layer;
      //  - number of chips in a stave;
      //  - chip center installation radius;
      //  - additional stave slope around beam line direction; [degree];
      //  - layer rotation around beam axis "as a whole"; [degree];
      vst->AddBarrelLayer(ibcell, 1*3*12,  1*9, 1*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 2*3*12,  1*9, 2*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 3*3*12,  2*9, 3*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 4*3*12,  2*9, 4*3*23.4 * etm::mm, 12.0, 0.0);
    }

    g4Reco->registerSubsystem(vst);
  }
#endif







	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);


  	G4HitNtuple *hits = new G4HitNtuple("Hits");
  	hits->AddNode("SVTX",0);
  	hits->AddNode("BMT",1);
  	hits->AddNode("CZBMT", 2);
  	se->registerSubsystem(hits);


	//---------------------------
	// fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//---------------------------
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("BARR");
	kalman->set_trackmap_out_name("SvtxTrackMap");


	#ifdef _SIVTX_
	// add Vertexing Layers
	kalman->add_phg4hits(
			"G4HIT_SVTX",				// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,					// radial-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),		// azimuthal-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),		// z-resolution [cm]
			1,					// efficiency,
			0					// noise hits
			);
	#endif

	#ifdef _SIBARR_
	// add Barrel Layers
	kalman->add_phg4hits(
			"G4HIT_BARR",                   	// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,                           	// radial-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// azimuthal-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// z-resolution [cm]
			1,                              	// efficiency,
			0                               	// noise hits
			);
	#endif
	
	#ifdef _SIMPLEVST_
	// add Barrel Layers
	kalman->add_phg4hits(
			"G4HIT_SIMPLEVST",                   	// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,                           	// radial-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// azimuthal-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// z-resolution [cm]
			1,                              	// efficiency,
			0                               	// noise hits
			);
	#endif

	#ifdef _SIDISKS_
	//  add Disk Layers
	kalman->add_phg4hits(
			"G4HIT_FBST",				// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Vertical_Plane,
			pix_size_dis/10000./sqrt(12.),		// radial-resolution [cm]
			pix_size_dis/10000./sqrt(12.),		// azimuthal-resolution [cm]
			999.,                       		// z-resolution [cm]
			1,                          		// efficiency,
			0                           		// noise hits
			);	
	#endif

	#ifdef _MPGD_
	//2D Readout 
    	if(use_2Dreadout)
	{
	    kalman->add_phg4hits(
	        "G4HIT_BMT",                //      const std::string& phg4hitsNames,
	        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
	        2.5/2/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type
	        150e-4,                       //        azimuthal-resolution [cm]
	        150e-4,                           //      z-resolution [cm]
	        1,                           //      efficiency,
	        0                            //      noise hits
		);
	}
	else
	{
	kalman->add_phg4hits(
        	"G4HIT_CZBMT",                //      const std::string& phg4hitsNames,
        	PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
        	2.5/2/sqrt(12),                      //       radial-resolution [cm], only used for Vertical Plane Detector Type
        	150e-4,                       //        azimuthal-resolution [cm]
        	150e-4,                           //      z-resolution [cm]
        	1,                           //      efficiency,
        	0                            //      noise hits
    		);	
	}
	#endif


	#ifdef _INNERGEMS_
        // GEM tracker hits;
        kalman->add_phg4hits(fgt->GetG4HitName(),
                             PHG4TrackFastSim::Vertical_Plane,
                             250e-4, //999. // radial-resolution [cm] (this number is not used in cylindrical geometry)                             
                             50e-4,        // azimuthal (arc-length) resolution [cm]    
                             999., //70e-4       // longitudinal (z) resolution [cm]                                                                   
                             1,// efficiency (fraction)                                                                        
                             0);// hit noise   
	#endif
	#ifdef _OUTERGEMS_
        // GEM tracker hits;
        kalman->add_phg4hits(fgt2->GetG4HitName(),
                             PHG4TrackFastSim::Vertical_Plane,
                             250e-4, //999. // radial-resolution [cm] (this number is not used in cylindrical geometry)                             
                             50e-4,        // azimuthal (arc-length) resolution [cm]    
                             999., //70e-4       // longitudinal (z) resolution [cm]                                                                   
                             1,// efficiency (fraction)                                                                        
                             0);// hit noise   
	#endif
	
	#ifdef _TPC_
        // GridPix tracker hits;
        kalman->add_phg4hits(
			     "G4HIT_GAS",
                             PHG4TrackFastSim::Cylinder,
                             999., // radial-resolution [cm] (this number is not used in cylindrical geometry)                             
                             100e-4,        // azimuthal (arc-length) resolution [cm]    
                             200e-4, //70e-4       // longitudinal (z) resolution [cm]                                                                   
                             0.9,// efficiency (fraction)                                                                        
                             0);// hit noise   
	#endif

	#ifdef _EICTOYVST_
    // Silicon tracker hits;
    kalman->add_phg4hits(vst->GetG4HitName(),		// const std::string& phg4hitsNames
			 PHG4TrackFastSim::Cylinder,	// const DETECTOR_TYPE phg4dettype
			 999.,				// radial-resolution [cm] (this number is not used in cylindrical geometry)
			 // 20e-4/sqrt(12) cm = 5.8e-4 cm, to emulate 20x20 um pixels;
			 //5.8e-4,			// azimuthal (arc-length) resolution [cm]
			 pix_size_vtx/10000./sqrt(12.),		// azimuthal-resolution [cm]
			 //5.8e-4,			// longitudinal (z) resolution [cm]
			 pix_size_vtx/10000./sqrt(12.),		// longitudinal (z) resolution [cm]
			 1,				// efficiency (fraction)
			 0);				// hit noise
	#endif
	
	if(do_projections){
		//NOTE: The output is in cm
		kalman->add_cylinder_state(projname1, projradius1);     // projection on cylinder (DIRC)
		kalman->add_zplane_state  (projname2, projzpos2  );     // projection on vertical planes
		kalman->add_zplane_state  (projname3, projzpos3  );     // projection on vertical planes
		kalman->add_zplane_state  (projname4, projzpos4  );     // projection on vertical planes
	}
	
	//kalman->Verbosity(10);
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_vertex_xy_resolution(0);
	kalman->set_vertex_z_resolution(0);
	kalman->enable_vertexing(false); // this is false by default
	kalman->set_vertex_min_ndf(2);

	se->registerSubsystem(kalman);

	std::string outputFile = "./Output/TEMP/"+(std::string)(out_name)+std::string(B_label)+"_FastSimEval.root";

	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(outputFile);
	if(do_projections){
		fast_sim_eval->AddProjection(projname1);
		fast_sim_eval->AddProjection(projname2);
		fast_sim_eval->AddProjection(projname3);
		fast_sim_eval->AddProjection(projname4);
	}
	se->registerSubsystem(fast_sim_eval);
        //se->registerSubsystem(new TrackFastSimEval());

	// ======================================================================================================
	// IOManagers...
	const std::string dst_name = "DSTFolder/"+std::string(out_name)+std::string(B_label)+"_G4GEM.root";	
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
	out->Verbosity(0);
	se->registerOutputManager(out);

	Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);

	if (nEvents <= 0) return;

	se->run(nEvents);
	se->End();
	delete se;

	gSystem->Exit(0);
}

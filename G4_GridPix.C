#ifndef MACRO_G4GRIDPIX_C
#define MACRO_G4GRIDPIX_C

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <cmath>

R__LOAD_LIBRARY(libg4detectors.so)



namespace G4GRIDPIX
{
  double cage_X0 = 0.01; //X0 = 1%
  double cathod_X0 = 0.005; //X0 = 0.5%
  double kapton_X0 = 28.57; //cm
  double cage_thickness = cage_X0*kapton_X0;
  double cathod_thickness = cathod_X0*kapton_X0;

  double beampipe_clearance = 3.64; 
 
  int n_gas_layers = 125; //14 cm drift
  //int gas_layers = 180; //20 cm drift
  double inner_radius = 23.5;
  double outer_radius = 37.5;
  double cage_length = 100.0;
  string tpcgas = "P10";
  
}  // namespace G4GridPix


void GridPixSetup(PHG4Reco *g4Reco)
{

	
  PHG4CylinderSubsystem* cyl;
  int n_layer = 0;
  double radius = G4GRIDPIX::inner_radius - G4GRIDPIX::cage_thickness;
// inner field cage
  cyl = new PHG4CylinderSubsystem("INCAGE", 0);
  cyl->set_double_param("radius", radius);
  cyl->set_double_param("length", G4GRIDPIX::cage_length);
  cyl->set_string_param("material", "G4_KAPTON");
  cyl->set_double_param("thickness",G4GRIDPIX::cage_thickness);
  cyl->SuperDetector("INCAGE");
  cyl->Verbosity(0);
  g4Reco->registerSubsystem(cyl);
  cyl->OverlapCheck();

  radius += G4GRIDPIX::cage_thickness;

//gas layers
  double tpc_layer_thickness = (G4GRIDPIX::outer_radius - G4GRIDPIX::inner_radius)/G4GRIDPIX::n_gas_layers; //thickness per layer
  for(int ilayer=0; ilayer < G4GRIDPIX::n_gas_layers; ilayer++)
    {
      cout << "gas layer: " << ilayer << " radius: " << radius << endl;
      cyl = new PHG4CylinderSubsystem("GAS",n_layer);
      cyl->set_double_param("radius", radius);
      cyl->set_double_param("length", G4GRIDPIX::cage_length);
      cyl->set_string_param("material", G4GRIDPIX::tpcgas);
      cyl->set_double_param("thickness",tpc_layer_thickness);
      cyl->SetActive();
      cyl->SuperDetector("GAS");
      g4Reco->registerSubsystem(cyl);
      
      radius += tpc_layer_thickness;
      n_layer++;
    }
  cout << "outer gas radius: " << radius << endl;
// outer field cage
  cyl = new PHG4CylinderSubsystem("OUTCAGE", 0);
  cyl->set_double_param("radius", radius);
  cyl->set_double_param("length", G4GRIDPIX::cage_length);
  cyl->set_string_param("material", "G4_KAPTON");
  cyl->set_double_param("thickness",G4GRIDPIX::cage_thickness);
  cyl->SuperDetector("OUTCAGE");
  cyl->Verbosity(0);
  g4Reco->registerSubsystem(cyl);
  cyl->OverlapCheck();

  radius += G4GRIDPIX::cage_thickness;
  cout << "Outer TPC radius: " << radius << endl; 

  //readout
  double readout_Si_thickness = 0.03; //300 um Si
  double readout_FR4_thickness = 0.65; //6.5 mm Si
  double readout_length = readout_Si_thickness;

  cyl = new PHG4CylinderSubsystem("READOUT", 0);
  cyl->set_double_param("radius", G4GRIDPIX::beampipe_clearance);
  cyl->set_string_param("material", "G4_Si");
  cyl->set_double_param("thickness",radius- G4GRIDPIX::beampipe_clearance);
  cyl->set_double_param("length",readout_length);
  cyl->set_double_param("place_z",G4GRIDPIX::cage_length/2.0 + readout_Si_thickness/2.0);
  cyl->SuperDetector("READOUT");
  cyl->Verbosity(0);
  g4Reco->registerSubsystem(cyl);
  cyl->OverlapCheck();

  readout_length +=readout_FR4_thickness;
  cyl = new PHG4CylinderSubsystem("READOUT", 1);
  cyl->set_double_param("radius", G4GRIDPIX::beampipe_clearance);
  cyl->set_string_param("material", "FR4");
  cyl->set_double_param("thickness",radius-G4GRIDPIX::beampipe_clearance);
  cyl->set_double_param("length",readout_length);
  //cyl->set_double_param("place_z",G4GRIDPIX::cage_length/2.0 + readout_Si_thickness + readout_FR4_thickness);
  cyl->set_double_param("place_z",G4GRIDPIX::cage_length/2.0 + readout_Si_thickness +readout_FR4_thickness);
  cyl->SuperDetector("READOUT");
  cyl->Verbosity(0);
  g4Reco->registerSubsystem(cyl);
  cyl->OverlapCheck();

  cyl = new PHG4CylinderSubsystem("CATHOD", 0);
  cyl->set_double_param("radius", G4GRIDPIX::beampipe_clearance);
  cyl->set_string_param("material", "G4_KAPTON");
  cyl->set_double_param("thickness",radius - G4GRIDPIX::beampipe_clearance);
  cyl->set_double_param("length",G4GRIDPIX::cathod_thickness);
  cyl->set_double_param("place_z",-G4GRIDPIX::cage_length/2.0 - G4GRIDPIX::cage_thickness);
  cyl->SuperDetector("CATHOD");
  cyl->Verbosity(0);
  g4Reco->registerSubsystem(cyl);
  cyl->OverlapCheck();

  return;	
	
}
#endif

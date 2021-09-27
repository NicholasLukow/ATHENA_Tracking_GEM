#include "Al_support_Detector.h"
#include <phparameter/PHParameters.h>
#include <g4main/PHG4Detector.h>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;
using namespace std;
//____________________________________________________________________________..
Al_support_Detector::Al_support_Detector(PHG4Subsystem *subsys,
		PHCompositeNode *Node,
		PHParameters *parameters,
		const std::string &dnam)
	: PHG4Detector(subsys, Node, dnam)
	  , m_Params(parameters)
{
}
//_______________________________________________________________
int Al_support_Detector::IsInDetector(G4VPhysicalVolume *volume) const
{
	set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if (iter != m_PhysicalVolumesSet.end())
	{
		return 1;
	}
	return 0;
}
//_______________________________________________________________
void Al_support_Detector::ConstructMe(G4LogicalVolume *logicWorld)
{
	//begin implement your own here://
	// Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
	addDetectorSection( logicWorld , "Al_pos_z" ,  1. );
	addDetectorSection( logicWorld , "Al_neg_z" , -1. );
	//end implement your own here://
	return;
}
// ======================================================================================================
void Al_support_Detector::Print(const std::string &what) const
{
	cout << "Al_support_ Detector:" << endl;
	if (what == "ALL" || what == "VOLUME")
	{
		cout << "Version 0.1" << endl;
		cout << "Parameters:" << endl;
		m_Params->Print();
	}
	return;
}
// ======================================================================================================
void  Al_support_Detector::addDetectorSection( G4LogicalVolume *logicWorld , std::string name , double sign){

	//double z_det[3] = {20.,46.75, 121.};
        //double rin  [3] = {7.13 * 20 / 22.0, 19.0,19.0};
	//double z_det[3] = {16.8,25.375, 121.};
        //double rin  [3] = {16.47 * 16.8 / 22.0, 19.0,19.0};
	
	//Projective 2.5 cone (expanded to outermost GEM)
	//double z_det[3] = {22, 55, 143};
        //double rin  [3] = {3.63, 9.10,23.64};
	
	//Projective 2.5 cone (expanded to outermost GEM) mRICH at -137
	//double z_det_temp[4] = {22, 153, -22, -135};
        //double rin_temp  [4] = {3.63, 25.29, 3.63, 22.31};

	//Symmetric 145, B2P2N2, projective at 2
	//double z_det_temp[4] = {85, 145, -85, -145};
        //double rin_temp  [4] = {23.44, 39.98, 23.44, 39.98};
	
	//Asymmetric, B2P2N2, projective at 2
	double z_det_temp[4] = {85, 170, -85, -150};
        double rin_temp  [4] = {23.44, 46.87, 23.44, 41.36};
	
	//clyinder instead of projective cone
	//double z_det_temp[4] = {85, 170, -85, -150};
        //double rin_temp  [4] = {41.36, 41.36, 41.36, 41.36};

	
	//mRICH at -200
	//double z_det_temp[4] = {22, 153, -22, -198};
        //double rin_temp  [4] = {3.63, 25.29, 3.63, 32.73};
	
	double z_det[2] = {0, 0};//153, -22, -135};
        double rin  [2] = {0,0};//3.63, 25.29, 3.63, 22.31};
	
	//double z_det[3] = {22, 55, 121.};
        //double rin  [3] = {3.63, 9.10,20.0};
	
	double rout [2] = {0};
	const int nzplanes = sizeof(z_det)/sizeof(*z_det);
	//for(int i = 0 ; i < nzplanes ; i++){
	for(int i = 0 ; i < 2 ; i++){
		int index = i + 1-sign;
		rout[i] = rin_temp[index] + 0.5;
		z_det[i] = z_det_temp[index]*10.;
		rin  [i] = rin_temp[index]*10.;
		rout [i] = rout[i]*10.;
	}

	G4Material * G4_mat = G4Material::GetMaterial("G4_Al");

	G4RotationMatrix *rotm = new G4RotationMatrix();
	rotm->rotateX(0);
	rotm->rotateY(0);
	rotm->rotateZ(0);

	G4Color G4_color = G4Color(G4Colour::Yellow());

	G4VSolid *G4_polycone = new G4Polycone(name,0,360*degree,nzplanes,z_det,rin,rout);
	G4LogicalVolume *logical = new G4LogicalVolume(G4_polycone,G4_mat, "Al_support_Logical");
	G4VisAttributes *vis_1 = new G4VisAttributes(G4_color);
	vis_1->SetForceSolid(true);
	logical->SetVisAttributes(vis_1);

	G4VPhysicalVolume *phy_1 = new G4PVPlacement(rotm, G4ThreeVector(0,0,0), logical , "Al_support_", logicWorld, 0, false, OverlapCheck());

	// add it to the list of placed volumes so the IsInDetector method picks them up
	m_PhysicalVolumesSet.insert(phy_1);
}

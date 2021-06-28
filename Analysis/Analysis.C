#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//---------------------------------------------------------------------------------------------------------------------------------------  
//============================================== Change Inputs Here =====================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  
//Here are the only parameters that need changing for most analyses (others are changed in the arguments given... See shell script)
//The following parts correspond to the different simulations you have run and want to process data for
//Give the number of detector configurations as well as the name used to identify it in the rootfile name and a longer description for labels
#define _NDet_ 4
//std::string DetVers[_NDet_] = {"Nominal", "NoGems" , "NoEndcapGem", "OneEndcapGem", "TwoEndcapGem", "NoFarGem"};
//const char *DetectorFullName[_NDet_] = {"Nominal Hybrid", "No GEM Disks", "No Endcap GEM Disks", "One Endcap GEM Disk", "Two Endcap GEM Disks", "Only Endcap GEM Disks"};
std::string DetVers[_NDet_] = { "Nominal3Endcap", "Nominal1Endcap","NoFar3Endcap", "NoFar1Endcap"};
const char *DetectorFullName[_NDet_] = {"3 Endcap GEM Disk", "1 Endcap GEM Disk", "3 Endcap GEM Disks, No Far GEM Disks", "1 Endcap GEM Disk, No Far GEM Disks"};

//Give the number of field maps as well as the name used to identify it in the rootfile name and a longer description for labels
#define _NBField_ 1
//std::string BField[_NBField_] = {"ATHENA", "NewBeAST", "Beast", "B_3.0T"};
//const char *FieldMapName[_NBField_] = {"05-28 BeAST Field Map Update", "05-07 BeAST Field Map Update", "Original BeAST Field Map", "Uniform 3.0T"};
std::string BField[_NBField_] = { "ATHENA"};
const char *FieldMapName[_NBField_] = {"ATHENA"};

//The maximum value for the y-axis. This will automatically be increased to be larger than the largest y-value if it is set too low
double ymax =  3.5;

//---------------------------------------------------------------------------------------------------------------------------------------  
//=======================================================================================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  


void Analysis( int _NEta_ = 1, double EtaMin = -1, double EtaMax = 1, int _NP_ = 29, double PMin = 1, double PMax =30, int _PT_ = 0, int _NY_ = 0, int _PROJ_ = 0, int _PVAR_ = 1, int _DRAWFITS_ = 0, int _DRAWPWGReq_ = 0 )
{

	//defining some variables/arrays
	double AngleValues[_NEta_ + 1];
	double MomentumValues[_NP_ + 1];
	int DimensionArray[2] = { _NP_, _NEta_};

	int _NDIM_;
	double _XMIN_;
	double _XMAX_;
	int _NPLOTS_;
	double effpt[_NDet_][_NBField_][_NEta_][_NP_];
	double effDenompt[_NDet_][_NBField_][_NEta_][_NP_];	
	double eff[_NDet_][_NBField_][_NEta_][_NP_];
	double effDenom[_NDet_][_NBField_][_NEta_][_NP_];	

	if (_NEta_ == 1 && _NP_ > 1) _NDIM_ = 0;
	if (_NP_ == 1 && _NEta_ > 1) _NDIM_ = 1;
	_XMIN_ = (1-_NDIM_)*(PMin-1)+(_NDIM_*(EtaMin-0.5));
	_XMAX_ = (1-_NDIM_)*(PMax+1)+(_NDIM_*(EtaMax+0.5));
	_NPLOTS_ = _NDet_*_NBField_;
	double EtaWidth = (EtaMax - EtaMin)/_NEta_;
	double PWidth = (PMax - PMin)/_NP_;
	
	for (int i = 0; i < _NEta_ + 1; i++)
	{
		AngleValues[i] = EtaMin + i*EtaWidth;
	}
	for (int i = 0; i < _NP_ + 1; i++)
	{
		MomentumValues[i] = PMin + i*PWidth;
	}
	//Read the output tree from the simulation and create histograms binned in eta and p/pt and saves them in an output histogram.root file
	//Make the histograms

	//Defining Histograms and variables
	TH1D *h_momRes[_NDet_][_NBField_][_NEta_][_NP_], *h_ptRes[_NDet_][_NBField_][_NEta_][_NP_];
	TH2D *h_nHits_momBin[_NDet_][_NBField_][_NEta_][_NP_], *h_nHits_ptBin[_NDet_][_NBField_][_NEta_][_NP_];

	TH1D *h_dca2d[_NDet_][_NBField_][_NEta_][_NP_], *h_dca2dpt[_NDet_][_NBField_][_NEta_][_NP_]; 
	
	TH1D *h_DIRCPointRes[_NDet_][_NBField_][_NEta_][_NP_], *h_DIRCPointRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	TH1D *h_DIRCThetaRes[_NDet_][_NBField_][_NEta_][_NP_], *h_DIRCThetaRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	TH1D *h_DIRCPhiRes[_NDet_][_NBField_][_NEta_][_NP_], *h_DIRCPhiRespt[_NDet_][_NBField_][_NEta_][_NP_];	

	TH1D *h_FORPointRes[_NDet_][_NBField_][_NEta_][_NP_], *h_FORPointRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	TH1D *h_FORThetaRes[_NDet_][_NBField_][_NEta_][_NP_], *h_FORThetaRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	TH1D *h_FORPhiRes[_NDet_][_NBField_][_NEta_][_NP_], *h_FORPhiRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	
	TH1D *h_BACKPointRes[_NDet_][_NBField_][_NEta_][_NP_], *h_BACKPointRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	TH1D *h_BACKThetaRes[_NDet_][_NBField_][_NEta_][_NP_], *h_BACKThetaRespt[_NDet_][_NBField_][_NEta_][_NP_];	
	TH1D *h_BACKPhiRes[_NDet_][_NBField_][_NEta_][_NP_], *h_BACKPhiRespt[_NDet_][_NBField_][_NEta_][_NP_];

	for (int iDet = 0; iDet < _NDet_; iDet++)
	{
		for (int iB = 0; iB < _NBField_; iB++)
		{
			for (int iEta = 0; iEta < _NEta_; iEta++)
			{
				for (int iP = 0; iP < _NP_; iP++)
				{
					//Creating the histograms with unique names, these names will be used later to grab the proper histogram for analysis
					h_momRes[iDet][iB][iEta][iP] = new TH1D(Form("h_momRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("dp/p for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 400, -0.15, 0.15);
					h_ptRes[iDet][iB][iEta][iP] = new TH1D(Form("h_ptRes_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("dpt/pt for %0.2lf < pt < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 400, -0.15, 0.15);
					h_nHits_momBin[iDet][iB][iEta][iP] = new TH2D(Form("h_nHits_MomentumBin_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf", DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Hits in each detector for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 5, 0, 5, 10, 0, 10);
					h_nHits_ptBin[iDet][iB][iEta][iP] = new TH2D(Form("h_nHits_PtBin_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf", DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Hits in each detector for %0.2lf < pt < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 5, 0, 5, 10, 0, 10);

					h_dca2d[iDet][iB][iEta][iP] = new TH1D(Form("h_dca2d_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("DCA2D for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -0.1, 0.1);
					h_dca2dpt[iDet][iB][iEta][iP] = new TH1D(Form("h_dca2dpt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("DCA2D for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -0.1, 0.1);
					
					double PHIRANGE = 0.01;
					double THETARANGE = 0.01;
					double POINTRANGE = 0.8;

					h_DIRCPointRes[iDet][iB][iEta][iP] = new TH1D(Form("h_DIRCPointRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Pointing Resolution at DIRC for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*POINTRANGE, POINTRANGE);
					h_DIRCPointRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_DIRCPointRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Pointing Resolution at DIRC for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*POINTRANGE, POINTRANGE);
					h_DIRCThetaRes[iDet][iB][iEta][iP] = new TH1D(Form("h_DIRCThetaRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Theta Resolution at DIRC for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*THETARANGE, THETARANGE);
					h_DIRCThetaRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_DIRCThetaRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Theta Resolution at DIRC for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*THETARANGE, THETARANGE);
					h_DIRCPhiRes[iDet][iB][iEta][iP] = new TH1D(Form("h_DIRCPhiRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Phi Resolution at DIRC for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*PHIRANGE, PHIRANGE);
					h_DIRCPhiRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_DIRCPhiRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Phi Resolution at DIRC for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*PHIRANGE, PHIRANGE);
					
					h_FORPointRes[iDet][iB][iEta][iP] = new TH1D(  Form("h_FORPointRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),    Form("Pointing Resolution at Hadron Calorimeter for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*POINTRANGE, POINTRANGE);
					h_FORPointRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_FORPointRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Pointing Resolution at Hadron Calorimeter for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*POINTRANGE, POINTRANGE);
					h_FORThetaRes[iDet][iB][iEta][iP] = new TH1D(  Form("h_FORThetaRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),       Form("Theta Resolution at Hadron Calorimeter for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*THETARANGE, THETARANGE);
					h_FORThetaRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_FORThetaRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),    Form("Theta Resolution at Hadron Calorimeter for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*THETARANGE, THETARANGE);
					h_FORPhiRes[iDet][iB][iEta][iP] = new TH1D(    Form("h_FORPhiRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),           Form("Phi Resolution at Hadron Calorimeter for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*PHIRANGE, PHIRANGE);
					h_FORPhiRespt[iDet][iB][iEta][iP] = new TH1D(  Form("h_FORPhiRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),        Form("Phi Resolution at Hadron Calorimeter for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*PHIRANGE, PHIRANGE);
					
					h_BACKPointRes[iDet][iB][iEta][iP] = new TH1D(  Form("h_BACKPointRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),    Form("Pointing Resolution at EM Calorimeter for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*POINTRANGE, POINTRANGE);
					h_BACKPointRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_BACKPointRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Pointing Resolution at EM Calorimeter for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*POINTRANGE, POINTRANGE);
					h_BACKThetaRes[iDet][iB][iEta][iP] = new TH1D(  Form("h_BACKThetaRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),       Form("Theta Resolution at EM Calorimeter for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*THETARANGE, THETARANGE);
					h_BACKThetaRespt[iDet][iB][iEta][iP] = new TH1D(Form("h_BACKThetaRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),    Form("Theta Resolution at EM Calorimeter for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*THETARANGE, THETARANGE);
					h_BACKPhiRes[iDet][iB][iEta][iP] = new TH1D(    Form("h_BACKPhiRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),           Form("Phi Resolution at EM Calorimeter for %0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*PHIRANGE, PHIRANGE);
					h_BACKPhiRespt[iDet][iB][iEta][iP] = new TH1D(  Form("h_BACKPhiRespt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]),        Form("Phi Resolution at EM Calorimeter for %0.2lf < p_t < %0.2lf and %0.2lf < #eta < %0.2lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1*PHIRANGE, PHIRANGE);
					
					
					//setting names for the hit histogram
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(1,"Si Vtx");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(2,"Si Barrel");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(3,"Barrel MPGD Tracker");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(4,"GEM Disks");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(5,"Si Disks");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetYaxis()->SetTitle("Hits");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetTitle("Detector");

					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(1,"Si Vtx");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(2,"Si Barrel");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(3,"Barrel MPGD Tracker");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(4,"GEM Disks");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(5,"Si Disks");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetYaxis()->SetTitle("Hits");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetTitle("Detector");


					effpt[iDet][iB][iEta][iP] = 0;
					effDenompt[iDet][iB][iEta][iP] = 0;	
					eff[iDet][iB][iEta][iP] = 0;
					effDenom[iDet][iB][iEta][iP] = 0;

				}
			}
		}
	}
	


	//declaring variables to be extracted from tree
   
  	Int_t           event;
   	Int_t           gtrackID;
   	Int_t           gflavor;
   	Float_t         gpx;
   	Float_t         gpy;
   	Float_t         gpz;
   	Float_t         gvx;
   	Float_t         gvy;
   	Float_t         gvz;
   	Float_t         gvt;
   	Int_t           trackID;
   	Int_t           charge;
   	Int_t           nhits;
   	Float_t         px;
   	Float_t         py;
   	Float_t         pz;
   	Float_t         pcax;
   	Float_t         pcay;
   	Float_t         pcaz;
   	Float_t         dca2d;
   	Int_t	   nHits_FGT;
   	Int_t	   nHits_SVTX;
   	Int_t 	   nHits_FBST;
   	Int_t	   nHits_BARR;
   	Int_t 	   nHits_BMT;
	Float_t	   DIRC_proj_x;
	Float_t	   DIRC_proj_y;
	Float_t	   DIRC_proj_z;
	Float_t	   DIRC_proj_lt;
	Float_t	   DIRC_proj_px;
	Float_t	   DIRC_proj_py;
	Float_t	   DIRC_proj_pz;
	Float_t	   DIRC_x;
	Float_t	   DIRC_y;
	Float_t	   DIRC_z;
	Float_t	   DIRC_lt;
	Float_t	   DIRC_px;
	Float_t	   DIRC_py;
	Float_t	   DIRC_pz;
	Float_t	   FOR_proj_x;
	Float_t	   FOR_proj_y;
	Float_t	   FOR_proj_z;
	Float_t	   FOR_proj_lt;
	Float_t	   FOR_proj_px;
	Float_t	   FOR_proj_py;
	Float_t	   FOR_proj_pz;
	Float_t	   FOR_x;
	Float_t	   FOR_y;
	Float_t	   FOR_z;
	Float_t	   FOR_lt;
	Float_t	   FOR_px;
	Float_t	   FOR_py;
	Float_t	   FOR_pz;
	Float_t	   BACK_proj_x;
	Float_t	   BACK_proj_y;
	Float_t	   BACK_proj_z;
	Float_t	   BACK_proj_lt;
	Float_t	   BACK_proj_px;
	Float_t	   BACK_proj_py;
	Float_t	   BACK_proj_pz;
	Float_t	   BACK_x;
	Float_t	   BACK_y;
	Float_t	   BACK_z;
	Float_t	   BACK_lt;
	Float_t	   BACK_px;
	Float_t	   BACK_py;
	Float_t	   BACK_pz;




//Starting the loop to extract information from the root trees (one root file per detecotr/field combination)
	for (int iDet = 0; iDet < _NDet_; iDet++)
	{
		for (int iB = 0; iB < _NBField_; iB++)
		{

			TString infile = Form("../Output/%s_%s_FastSimEval.root", DetVers[iDet].c_str(), BField[iB].c_str());

			//check if file exists
			if(gSystem->AccessPathName(infile)) {std::cout << "file " << infile << " does not exist" << std::endl; return;} 
			else std::cout << "Found File: " << infile << std::endl;
			//get the file
			TFile myFile(infile.Data());
			//get the tree
			TTree *tree;
			myFile.GetObject("tracks", tree);

			//get number of entries
			   Int_t nentries = (Int_t)tree->GetEntries();

			//set the branch addresses
			   tree->SetBranchAddress("event", &event);
			   tree->SetBranchAddress("gtrackID", &gtrackID);
			   tree->SetBranchAddress("gflavor", &gflavor);
			   tree->SetBranchAddress("gpx", &gpx);
			   tree->SetBranchAddress("gpy", &gpy);
			   tree->SetBranchAddress("gpz", &gpz);
			   tree->SetBranchAddress("gvx", &gvx);
			   tree->SetBranchAddress("gvy", &gvy);
			   tree->SetBranchAddress("gvz", &gvz);
			   tree->SetBranchAddress("gvt", &gvt);
			   tree->SetBranchAddress("trackID", &trackID);
			   tree->SetBranchAddress("charge", &charge);
			   tree->SetBranchAddress("nhits", &nhits);
			   tree->SetBranchAddress("px", &px);
			   tree->SetBranchAddress("py", &py);
			   tree->SetBranchAddress("pz", &pz);
			   tree->SetBranchAddress("pcax", &pcax);
			   tree->SetBranchAddress("pcay", &pcay);
			   tree->SetBranchAddress("pcaz", &pcaz);
			   tree->SetBranchAddress("dca2d", &dca2d);
			   tree->SetBranchAddress("nHit_G4HIT_FGT", &nHits_FGT);
			   tree->SetBranchAddress("nHit_G4HIT_SVTX", &nHits_SVTX);
			   tree->SetBranchAddress("nHit_G4HIT_FBST", &nHits_FBST);
			   tree->SetBranchAddress("nHit_G4HIT_BARR", &nHits_BARR);
			   tree->SetBranchAddress("nHit_G4HIT_BMT", &nHits_BMT);
			   
			   tree->SetBranchAddress("DIRC_proj_x", &DIRC_proj_x);
			   tree->SetBranchAddress("DIRC_proj_y", &DIRC_proj_y);
			   tree->SetBranchAddress("DIRC_proj_z", &DIRC_proj_z);
			   tree->SetBranchAddress("DIRC_proj_t", &DIRC_proj_lt);
			   tree->SetBranchAddress("DIRC_proj_px", &DIRC_proj_px);
			   tree->SetBranchAddress("DIRC_proj_py", &DIRC_proj_py);
			   tree->SetBranchAddress("DIRC_proj_pz", &DIRC_proj_pz);
			   tree->SetBranchAddress("DIRC_x", &DIRC_x);
			   tree->SetBranchAddress("DIRC_y", &DIRC_y);
			   tree->SetBranchAddress("DIRC_z", &DIRC_z);
			   tree->SetBranchAddress("DIRC_t", &DIRC_lt);
			   tree->SetBranchAddress("DIRC_px", &DIRC_px);
			   tree->SetBranchAddress("DIRC_py", &DIRC_py);
			   tree->SetBranchAddress("DIRC_pz", &DIRC_pz);
			   tree->SetBranchAddress("FOR_proj_x", &FOR_proj_x);
			   tree->SetBranchAddress("FOR_proj_y", &FOR_proj_y);
			   tree->SetBranchAddress("FOR_proj_z", &FOR_proj_z);
			   tree->SetBranchAddress("FOR_proj_t", &FOR_proj_lt);
			   tree->SetBranchAddress("FOR_proj_px", &FOR_proj_px);
			   tree->SetBranchAddress("FOR_proj_py", &FOR_proj_py);
			   tree->SetBranchAddress("FOR_proj_pz", &FOR_proj_pz);
			   tree->SetBranchAddress("FOR_x", &FOR_x);
			   tree->SetBranchAddress("FOR_y", &FOR_y);
			   tree->SetBranchAddress("FOR_z", &FOR_z);
			   tree->SetBranchAddress("FOR_t", &FOR_lt);
			   tree->SetBranchAddress("FOR_px", &FOR_px);
			   tree->SetBranchAddress("FOR_py", &FOR_py);
			   tree->SetBranchAddress("FOR_pz", &FOR_pz);
			   tree->SetBranchAddress("BACK_proj_x", &BACK_proj_x);
			   tree->SetBranchAddress("BACK_proj_y", &BACK_proj_y);
			   tree->SetBranchAddress("BACK_proj_z", &BACK_proj_z);
			   tree->SetBranchAddress("BACK_proj_t", &BACK_proj_lt);
			   tree->SetBranchAddress("BACK_proj_px", &BACK_proj_px);
			   tree->SetBranchAddress("BACK_proj_py", &BACK_proj_py);
			   tree->SetBranchAddress("BACK_proj_pz", &BACK_proj_pz);
			   tree->SetBranchAddress("BACK_x", &BACK_x);
			   tree->SetBranchAddress("BACK_y", &BACK_y);
			   tree->SetBranchAddress("BACK_z", &BACK_z);
			   tree->SetBranchAddress("BACK_t", &BACK_lt);
			   tree->SetBranchAddress("BACK_px", &BACK_px);
			   tree->SetBranchAddress("BACK_py", &BACK_py);
			   tree->SetBranchAddress("BACK_pz", &BACK_pz);


			// run the "analysis" simply getting tracks and filling the residual histograms
			for (int i = 0; i < nentries; i++)
			{
				tree->GetEntry(i);
				
				//retrieving the momenta
				TVector3 truthP(gpx, gpy, gpz);
				TVector3 recoP(px, py, pz);

				//hits and projections of tracks 
				TVector3 DIRCPos(DIRC_x, DIRC_y, DIRC_z);
				TVector3 DIRCProjPos(DIRC_proj_x, DIRC_proj_y, DIRC_proj_z);
				TVector3 FORPos(FOR_x, FOR_y, FOR_z);
				TVector3 FORProjPos(FOR_proj_x, FOR_proj_y, FOR_proj_z);
				TVector3 BACKPos(BACK_x, BACK_y, BACK_z);
				TVector3 BACKProjPos(BACK_proj_x, BACK_proj_y, BACK_proj_z);
				//hits and projections of tracks 
				TVector3 DIRCMom(DIRC_px, DIRC_py, DIRC_pz);
				TVector3 DIRCProjMom(DIRC_proj_px, DIRC_proj_py, DIRC_proj_pz);
				TVector3 FORMom(FOR_px, FOR_py, FOR_pz);
				TVector3 FORProjMom(FOR_proj_px, FOR_proj_py, FOR_proj_pz);
				TVector3 BACKMom(BACK_px, BACK_py, BACK_pz);
				TVector3 BACKProjMom(BACK_proj_px, BACK_proj_py, BACK_proj_pz);
				//Loop over eta as to store information in the properly indexed histogram
				for (int iEta = 0; iEta < _NEta_; iEta++)
				{
					//theres probably a better way than loops and if-statements to determine the proper indexed histogram for a given track, but this works
					if (truthP.Eta() >= AngleValues[iEta] && truthP.Eta() <= AngleValues[iEta+1]) 
					{
						//Loop over momentum as to store information in the properly indexed histogram
						for (int iP = 0; iP < _NP_; iP++)
						{				
							//filling histograms which are binned in momentum
							if (truthP.Mag() >= MomentumValues[iP] && truthP.Mag() <= MomentumValues[iP+1])
							{
								effDenom[iDet][iB][iEta][iP] +=1;
								//making sure there is a reconstructed track
								if (trackID != -9999)
								{
									//Calculations for the resolutions at DIRC radius
									double dPhi = DIRCPos.Phi() - DIRCProjPos.Phi();
									double dZ = DIRC_z - DIRC_proj_z;
									double Radius = DIRCPos.Mag()*TMath::Sin(DIRCPos.Theta());
									double CylS = TMath::Sqrt( (Radius*dPhi)*(Radius*dPhi) + dZ*dZ ); //Calculating the distance between the points along the cylindrical surface. Angle are so small that it is effectively equivalent to the magnitude of the difference of the position vectors
									//default values for hit locations are -9999, no hit results in a position magnitude of ~17318.8cm
									if (DIRCPos.Mag() < 17318)
									{
										h_DIRCPointRes[iDet][iB][iEta][iP]->Fill(CylS);				
										h_DIRCThetaRes[iDet][iB][iEta][iP]->Fill(DIRCMom.Theta() - DIRCProjMom.Theta());
										h_DIRCPhiRes[iDet][iB][iEta][iP]->Fill(DIRCMom.Phi() - DIRCProjMom.Phi());
									}
		
									if (FORPos.Mag() < 17318)
									{
										h_FORPointRes[iDet][iB][iEta][iP]->Fill( (FORPos-FORProjPos).Mag() );				
										h_FORThetaRes[iDet][iB][iEta][iP]->Fill(FORMom.Theta() - FORProjMom.Theta());
										h_FORPhiRes[iDet][iB][iEta][iP]->Fill(FORMom.Phi() - FORProjMom.Phi());
									}
									
									if (BACKPos.Mag() < 17318)
									{
										h_BACKPointRes[iDet][iB][iEta][iP]->Fill( (BACKPos-BACKProjPos).Mag() );
										h_BACKThetaRes[iDet][iB][iEta][iP]->Fill(BACKMom.Theta() - BACKProjMom.Theta());
										h_BACKPhiRes[iDet][iB][iEta][iP]->Fill(BACKMom.Phi() - BACKProjMom.Phi());
									}	
									

									h_momRes[iDet][iB][iEta][iP]->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(0., nHits_SVTX);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(1., nHits_BARR);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(2., nHits_BMT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(3., nHits_FGT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(4., nHits_FBST);
								
									eff[iDet][iB][iEta][iP] += 1;	
									h_dca2d[iDet][iB][iEta][iP]->Fill(dca2d);
								
								}
							}

							//filling histograms which are binned in pt
							if (truthP.Pt() >= MomentumValues[iP] && truthP.Pt() <= MomentumValues[iP+1])
							{
								effDenompt[iDet][iB][iEta][iP] +=1;
								//making sure there is a reconstructed track
								if (trackID != -9999)
								{
									//Calculations for the resolutions at DIRC radius
									
									double dPhi = DIRCPos.Phi() - DIRCProjPos.Phi();
									double dZ = DIRC_z - DIRC_proj_z;
									double Radius = DIRCPos.Mag()*TMath::Sin(DIRCPos.Theta());
									double CylS = TMath::Sqrt( (Radius*dPhi)*(Radius*dPhi) + dZ*dZ );
									//default values for hit locations are -9999, no hit results in a position magnitude of ~17318.8cm
									if (DIRCPos.Mag() < 17318)
									{
										h_DIRCPointRespt[iDet][iB][iEta][iP]->Fill(CylS);
										h_DIRCThetaRespt[iDet][iB][iEta][iP]->Fill(DIRCPos.Theta() - DIRCProjPos.Theta());
										h_DIRCPhiRespt[iDet][iB][iEta][iP]->Fill(dPhi);	
									}
									
									if (FORPos.Mag() < 17318)
									{
										h_FORPointRespt[iDet][iB][iEta][iP]->Fill( (DIRCPos-DIRCProjPos).Mag() );	
										h_FORThetaRespt[iDet][iB][iEta][iP]->Fill(DIRCPos.Theta() - DIRCProjPos.Theta());
										h_FORPhiRespt[iDet][iB][iEta][iP]->Fill(dPhi);
									}

									if (BACKPos.Mag() < 17318)
									{
										h_BACKPointRespt[iDet][iB][iEta][iP]->Fill( (DIRCPos-DIRCProjPos).Mag() );	
										h_BACKThetaRespt[iDet][iB][iEta][iP]->Fill(DIRCPos.Theta() - DIRCProjPos.Theta());
										h_BACKPhiRespt[iDet][iB][iEta][iP]->Fill(dPhi);				
									}

									h_ptRes[iDet][iB][iEta][iP]->Fill( (recoP.Pt() - truthP.Pt())/truthP.Pt());  
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(0., nHits_SVTX);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(1., nHits_BARR);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(2., nHits_BMT);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(3., nHits_FGT);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(4., nHits_FBST);
									
									effpt[iDet][iB][iEta][iP] += 1;	
									h_dca2dpt[iDet][iB][iEta][iP]->Fill(dca2d);
								}

							}
						}//end loop over momentum values

					} //end Eta if
				} //end loop over eta values
			} //end loop over events


			//Open an output histogram file and write the histograms
			TFile outfile("./Output/Histogram.root", "UPDATE");
			for (int iEta = 0; iEta < _NEta_; iEta++)
			{
				for (int iP = 0; iP < _NP_; iP++)
				{
					h_momRes[iDet][iB][iEta][iP]->Write();
					h_ptRes[iDet][iB][iEta][iP]->Write();
					h_nHits_ptBin[iDet][iB][iEta][iP]->Write();
					h_nHits_momBin[iDet][iB][iEta][iP]->Write();

					h_dca2d[iDet][iB][iEta][iP]->Write();
					h_dca2dpt[iDet][iB][iEta][iP]->Write();
					
					h_DIRCPointRes[iDet][iB][iEta][iP]->Write();
					h_DIRCPointRespt[iDet][iB][iEta][iP]->Write();
					h_DIRCThetaRes[iDet][iB][iEta][iP]->Write();
					h_DIRCThetaRespt[iDet][iB][iEta][iP]->Write();
					h_DIRCPhiRes[iDet][iB][iEta][iP]->Write();
					h_DIRCPhiRespt[iDet][iB][iEta][iP]->Write();

					h_FORPointRes[iDet][iB][iEta][iP]->Write();
					h_FORPointRespt[iDet][iB][iEta][iP]->Write();
					h_FORThetaRes[iDet][iB][iEta][iP]->Write();
					h_FORThetaRespt[iDet][iB][iEta][iP]->Write();
					h_FORPhiRes[iDet][iB][iEta][iP]->Write();
					h_FORPhiRespt[iDet][iB][iEta][iP]->Write();

					h_BACKPointRes[iDet][iB][iEta][iP]->Write();
					h_BACKPointRespt[iDet][iB][iEta][iP]->Write();
					h_BACKThetaRes[iDet][iB][iEta][iP]->Write();
					h_BACKThetaRespt[iDet][iB][iEta][iP]->Write();
					h_BACKPhiRes[iDet][iB][iEta][iP]->Write();
					h_BACKPhiRespt[iDet][iB][iEta][iP]->Write();
					

				}
			}
			//close the output histogram file
			outfile.Close();
			//close the current file to prepare for new file to be opened
			myFile.Close();
		}//end loop over B Field options
	}//end loop over the detector configurations

//-------END MAKE HISTOGRAMS-------

//Reads the histogram file and fits the distributions, stores the widths in arrays, and uses the arrays to create plots
//Make the plots

	//Strings used for getting histograms and labeling axes. Should be no need to change these unless a change is made to the Fun4All macro.
	std::string ProjectionLocation[3] = {"DIRC", "FOR" , "BACK"};
	std::string LocationFullName[3] = {"DIRC", "RICH" , "EMCal"};
	std::string ProjectionVariable[3] = {"Point", "Theta" , "Phi"};
	std::string VariableFullName[3] = {"Pointing", "Theta" , "Phi"};

	//Calculating the bin centers for plotting
	double MomentumBinCenter[_NP_];
  	for (int i = 0; i < _NP_; i++)
  	{
		MomentumBinCenter[i] = 0.5*(MomentumValues[i]+MomentumValues[i+1]);
  	}
	double AngleBinCenter[_NEta_];
  	for (int i = 0; i < _NEta_; i++)
  	{
		AngleBinCenter[i] = 0.5*(AngleValues[i]+AngleValues[i+1]);
  	}

	//getting the histograms created in the earlier step
	TFile *histoFile = new TFile("./Output/Histogram.root","read");
	

	//Get values from histograms and store in arrays
	double Resolution[_NEta_][_NP_][_NDet_][_NBField_]; 
	double ResErr[_NEta_][_NP_][_NDet_][_NBField_];
	double Efficiency[_NEta_][_NP_][_NDet_][_NBField_]; 



	for (int iD = 0; iD < _NDet_; iD++)
	{
		for (int iA = 0; iA < _NEta_; iA++)
		{
		  	for (int iP = 0; iP < _NP_; iP++)
			{
		        	for (int iB = 0; iB < _NBField_; iB++) 
		        	{
		        		//Construct Histo Name
						TString HistName;

						if(_NY_ == 0)
						{
							if (_PT_) HistName = Form("h_ptRes_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
							else HistName = Form("h_momRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
						}
						else if (_NY_ == 1)
						{
							if (_PT_) HistName = Form("h_dca2dpt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
							else HistName = Form("h_dca2d_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
						}  
						else if (_NY_ == 3)
						{
							if (_PT_) HistName = Form("h_%s%sRes_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",ProjectionLocation[_PROJ_].c_str(),ProjectionVariable[_PVAR_].c_str(), DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
							else HistName = Form("h_%s%sRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",ProjectionLocation[_PROJ_].c_str(),ProjectionVariable[_PVAR_].c_str(), DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
						} 	
						//Get Histogram
		        		//cout << "Getting: " << HistName << endl;
		        		TH1D *Htmp = (TH1D*)histoFile->Get(HistName);
					if(!Htmp && _NY_ != 2){ cout << "Histogram not found: " << HistName << endl; break;}
					else if (Htmp)
					{

						//Extract the value from the histogram and store in array	
			        		TF1 *gausFit = new TF1("gausFit", "gaus");
			        		gausFit->SetParameter(1, 0);
	       		    			Htmp->Fit(gausFit, "Q");
			        		double sigma = gausFit->GetParameter(2);
			        		double sigmaErr = gausFit->GetParError(2);
			        		if (sigma == 0) cout << "ERROR! - Bad Value for : " << HistName << endl;

				    		Resolution[iA][iP][iD][iB] = sigma;
	            				ResErr[iA][iP][iD][iB] = sigmaErr;

						TFile outfile("./Output/Fits.root", "UPDATE");
						Htmp->Write();
						outfile.Close();
			        	}

					//Calculating the efficiency
					if(_PT_ && effDenompt[iD][iB][iA][iP] != 0 ) Efficiency[iA][iP][iD][iB]= effpt[iD][iB][iA][iP]/effDenompt[iD][iB][iA][iP];
 					else if (!_PT_ && effDenom[iD][iB][iA][iP] != 0 ) Efficiency[iA][iP][iD][iB]= eff[iD][iB][iA][iP]/effDenom[iD][iB][iA][iP];
					else {cout << "Division by 0! Setting Efficiency to 0" << endl; Efficiency[iA][iP][iD][iB] = 0;}


			        }//End B Loop
			} // End Momentum Loop
		} // End Angle Loop
	} // End Detector Loop



	//setting options
    gStyle->SetTextSize(0.02);
    gStyle->SetLabelSize(0.04,"xy");
    gStyle->SetFrameFillColor(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    //gStyle->SetTitleSize(0.05,"x");
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.10);
    //gStyle->SetMarkerColor(kBlack);
    //gStyle->SetMarkerStyle(25);
    //gStyle->SetMarkerSize(1.0);  

    gStyle->SetStatBorderSize(0);
    gStyle->SetStatColor(kWhite);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetStatFont(52);
    gStyle->SetStatW(.13);
    gStyle->SetFitFormat("2.1e [%]");
    //gStyle->SetFitFormat("2.1e");


	//create canvas
	TString OutputName = Form("cPlot_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf", PMin, PMax, EtaMin, EtaMax);
	TString OutputTitle = Form("%0.2lf < p < %0.2lf and %0.2lf < #eta < %0.2lf", PMin, PMax, EtaMin, EtaMax);
    TCanvas *c1 = new TCanvas(OutputName, OutputTitle , 0, 0, 1200, 900);
    c1->UseCurrentStyle();
    c1->SetBorderMode(0);
    c1->SetFrameBorderMode(0);
    c1->SetFrameLineColor(kWhite);
    c1->SetFillColor(kWhite);
    c1->cd();


    //PWG Resolution Requirements from the Yellow Report (Table 11.2)
    TF1 *PWGReq[5];
    PWGReq[0] = new TF1("PWGReq_-3.5_-2.5", "TMath::Sqrt((0.10*x)^2+0.5^2)", 1, 30);
    PWGReq[1] = new TF1("PWGReq_-2.5_-1.0", "TMath::Sqrt((0.05*x)^2+0.5^2)", 1, 30);
    PWGReq[2] = new TF1("PWGReq_-1.0_1.0" , "TMath::Sqrt((0.05*x)^2+0.5^2)", 1, 30);
    PWGReq[3] = new TF1("PWGReq_1.0_2.5"  , "TMath::Sqrt((0.05*x)^2+1.0^2)", 1, 30);
    PWGReq[4] = new TF1("PWGReq_2.5_3.5"  , "TMath::Sqrt((0.10*x)^2+2.0^2)", 1, 30);


    //PWG Transverse Pointing Resolution from the Yellow Report (Table 11.2)
    TF1 *PWGReq_DCA2D[6];
    PWGReq_DCA2D[0] = new TF1("PWGReq_DCA2D_-3.0_-2.5", "TMath::Sqrt((30/x)^2+40^2)", 1, 30);
    PWGReq_DCA2D[1] = new TF1("PWGReq_DCA2D_-2.5_-1", "TMath::Sqrt((30/x)^2+20^2)", 1, 30);
    PWGReq_DCA2D[2] = new TF1("PWGReq_DCA2D_-1_1", "TMath::Sqrt((20/x)^2+5^2)", 1, 30);
    PWGReq_DCA2D[3] = new TF1("PWGReq_DCA2D_1_2.5", "TMath::Sqrt((30/x)^2+20^2)", 1, 30);
    PWGReq_DCA2D[4] = new TF1("PWGReq_DCA2D_2.5_3", "TMath::Sqrt((30/x)^2+40^2)", 1, 30);
    PWGReq_DCA2D[5] = new TF1("PWGReq_DCA2D_3_3.5", "TMath::Sqrt((30/x)^2+60^2)", 1, 30);



	double res[DimensionArray[_NDIM_]];
	double reserr[DimensionArray[_NDIM_]];
	int countPlots = 0;

	//double minPt = 99;
	TLine *EfficiencyMinimum = new TLine(_XMIN_, 0.90, _XMAX_, 0.90);
	EfficiencyMinimum->SetLineColor(kRed);
	EfficiencyMinimum->SetLineStyle(9);

	TF1 *ReqFit[_NPLOTS_];
	TGraphErrors *gResPlot[_NPLOTS_];

	//can compare up to 5 different detector configurations before this will give an error
	EColor ColorArray[10] = {kGreen, kBlue, kMagenta, kRed, kCyan, kOrange, kYellow, kPink, kViolet, kBlack};


	for (int iD = 0; iD < _NDet_; iD++)
	{
		for (int iB = 0; iB < _NBField_; iB++) 
		{
			//Index is reset for the new plot
			int index=0;
			if (_NP_*_NEta_ > _NP_+_NEta_-1) {cout << "Incorrect number of bins. Make sure either Eta or Momentum is set to 1 bin" << endl; break;} //sanity check
		 	for (int iP = 0; iP < _NP_; iP++)
		      	{
			  	for (int iA = 0; iA < _NEta_; iA++)
		        {
					if (_NY_ == 0)
					{
						//convert to %
						res[index] = 100*Resolution[iA][iP][iD][iB];
	   	  	  	  		reserr[index] = 100*ResErr[iA][iP][iD][iB];
					}
					else if (_NY_ == 1)
					{
						//convert from cm to um
						res[index] = 10000*Resolution[iA][iP][iD][iB];
	   	  	  	  		reserr[index] = 10000*ResErr[iA][iP][iD][iB];
					}
					else if (_NY_ == 2)
					{
						res[index] = Efficiency[iA][iP][iD][iB];
						reserr[index] = 0;

						//if (res[index] > 0.9 && minPt == 99 ) minPt = MomentumBinCenter[iP];
						//else if (res[index] < 0.9 && minPt < 99) minPt = 99; //incase the previous value was just a fluctuation reset the minPt value
					}
					else if (_NY_ ==3 && _PVAR_ == 0)
					{
						//convert from cm to um
						res[index] = 10000*Resolution[iA][iP][iD][iB];
	   	  	  	  		reserr[index] = 10000*ResErr[iA][iP][iD][iB];
					}
					else if (_NY_ ==3 && (_PVAR_ == 1 || _PVAR_ == 2))
					{
						//convert from rad to mrad
						res[index] = 1000*Resolution[iA][iP][iD][iB];
	   	  	  	  		reserr[index] = 1000*ResErr[iA][iP][iD][iB];
					}
					else {cout << "No valid plotting variable chosen" << endl; break;}
					if ( res[index] > ymax) ymax = res[index];
				  	index++;
				}//End Angle Loop
			}//End Momentum Loop

			//cout << "Minimum Pt: " << minPt << " +/- " << (MomentumBinCenter[0] + MomentumBinCenter[1])/2 << endl;     		

			if(_NY_ == 0 && !_PT_) ReqFit[countPlots] = new TF1(Form("ReqFit_%d", countPlots), " TMath::Sqrt(([0]*x)^2+[1]^2)", 1 , 31);
			else if (_NY_ == 1 && _PT_) ReqFit[countPlots] = new TF1(Form("ReqFit_%d", countPlots), " TMath::Sqrt(([0]/x)^2+[1]^2)", 1 , 31); 

			//the plot is constructed from the calculated bin centers (from above) and the resolution values and errors extracted from the full multidimensional array 	
			if (_NDIM_ == 0) gResPlot[countPlots] = new TGraphErrors(DimensionArray[_NDIM_], MomentumBinCenter , res, 0, reserr );
			else if (_NDIM_ == 1) gResPlot[countPlots] = new TGraphErrors(DimensionArray[_NDIM_], AngleBinCenter , res, 0, reserr  );
	 		else { cout << "Incorrect number of bins. Make sure either Eta or Momentum is set to 1 bin" << endl; break;}

			//setting some plot options
			gResPlot[countPlots]->SetTitle("");
	 		gResPlot[countPlots]->SetMarkerSize(1.3);
	 		gResPlot[countPlots]->SetMarkerStyle(21);
	 		gResPlot[countPlots]->SetMarkerColor(ColorArray[countPlots]);
		    gResPlot[countPlots]->SetLineColor(ColorArray[countPlots]);
		    gResPlot[countPlots]->SetLineWidth(2);
		    gResPlot[countPlots]->SetLineStyle(7);

			if ( (_NY_ == 0 && !_PT_ )|| ( _NY_ == 1 && _PT_))
			{ 	
				//ReqFit[countPlots]->SetParameter(0,0.04);
		     	ReqFit[countPlots]->SetParName(0,"A");
		     	//ReqFit[countPlots]->SetParameter(1,1.0);
		     	ReqFit[countPlots]->SetParName(1,"B");
		     	ReqFit[countPlots]->SetLineColor(ColorArray[countPlots]+1);
			}
			//increment plot index
			countPlots++;

		}//End B Field Loop
	} //End Detector Loop


    //Creating a dummy histogram which will be formatted
	
    TH1D *hdum  = new TH1D(OutputName, "", 15, _XMIN_, _XMAX_);
    hdum->GetXaxis()->SetLabelFont(52);
    hdum->GetYaxis()->SetLabelFont(52);
    hdum->GetXaxis()->SetTitleFont(52);
    hdum->GetYaxis()->SetTitleFont(52);
    if (_NDIM_ == 0 && !_PT_) hdum->GetXaxis()->SetTitle("p [GeV]");
    else if (_NDIM_ == 0 && _PT_) hdum->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (_NDIM_ == 1) hdum->GetXaxis()->SetTitle("#eta");
    hdum->GetXaxis()->SetTitleSize(0.06);
    hdum->GetXaxis()->SetLabelSize(0.05);
    if(_NY_ == 0)
    {
    	if (_PT_) hdum->GetYaxis()->SetTitle("p_{T} resolution   #sigma_{p_{T}} /p_{T}  [%]");
    	else hdum->GetYaxis()->SetTitle("p resolution   #sigma_{p} /p  [%]");
    }
    else if (_NY_ == 1)
    {
    	hdum->GetYaxis()->SetTitle("#sigma_{DCA_{2D}} [#mum]");
    }
    else if (_NY_ == 2)
    {
    	hdum->GetYaxis()->SetTitle("Efficiency");
    }
	else if (_NY_ ==3 && _PVAR_ == 0)
	{
		hdum->GetYaxis()->SetTitle("Pointing Resolution [#mum]");
	}
	else if (_NY_ ==3 && (_PVAR_ == 1 || _PVAR_ == 2))
	{
		hdum->GetYaxis()->SetTitle(Form( "%s Resolution [mRad]" ,VariableFullName[_PVAR_].c_str()));
	}

    hdum->GetYaxis()->SetTitleSize(0.06);
    hdum->GetYaxis()->SetLabelSize(0.05);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(0.75);
    hdum->SetMinimum( 0 );
    hdum->SetMaximum( ymax+0.25*ymax );

    //Title will be automatically set (since either p/eta is integrated over, it will be noted in the title) 
    if (_NY_ ==0 )
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "Resolution for tracks with %0.2lf < #eta < %0.2lf", EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "Resolution for tracks with %0.2lf < p < %0.2lf", PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "Resolution for tracks with %0.2lf < p_T < %0.2lf", PMin, PMax));
    }
    else if (_NY_ == 1 )
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "DCA2D for tracks with %0.2lf < #eta < %0.2lf", EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "DCA2D for tracks with %0.2lf < p < %0.2lf", PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "DCA2D for tracks with %0.2lf < p_T < %0.2lf", PMin, PMax));
    }
    else if (_NY_ == 2 )
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "Efficiency for tracks with %0.2lf < #eta < %0.2lf", EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "Efficiency for tracks with %0.2lf < p < %0.2lf", PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "Efficiency for tracks with %0.2lf < p_T < %0.2lf", PMin, PMax));
    }
	else if (_NY_ == 3)
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "%s Resolution for %s for tracks with %0.2lf < #eta < %0.2lf",VariableFullName[_PVAR_].c_str(), LocationFullName[_PROJ_].c_str(), EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "%s Resolution for %s for tracks with %0.2lf < p < %0.2lf", VariableFullName[_PVAR_].c_str(), LocationFullName[_PROJ_].c_str(), PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "%s Resoltuoin for %s for tracks with %0.2lf < p_T < %0.2lf", VariableFullName[_PVAR_].c_str(), LocationFullName[_PROJ_].c_str(), PMin, PMax));
    }


    //hdum->GetXaxis()->SetNdivisions(408);
    //hdum->GetYaxis()->SetNdivisions(804);
    hdum->GetXaxis()->SetNdivisions(408);
    hdum->GetYaxis()->SetNdivisions(808);
    gPad->SetGrid();
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFrameBorderMode(0);
    gPad->SetFrameLineColor(kWhite);
    hdum->Draw();



	TLegend *legend = new TLegend(.15, .7, .45, .92);
	for (int iPlot = 0; iPlot < _NPLOTS_; iPlot++)
	{ 
	     if(_DRAWFITS_ && ( (_NY_ == 0 && !_PT_ )|| ( _NY_ == 1 && _PT_)) ) gResPlot[iPlot]->Fit(ReqFit[iPlot], "R");
	     gResPlot[iPlot]->Draw("PLsame"); //replace L with C for smooth line

	     //legend->AddEntry(gResPlot[iPlot], Form("Hybrid Detector Resolution for %0.2lf < #eta < %0.2lf", AngleValues[0], AngleValues[1]  ), "P");
	
	     if (_NDet_ > 1 && _NBField_ == 1 ) legend->AddEntry(gResPlot[iPlot], Form("Detector Version: %s", DetectorFullName[iPlot]), "P");
	     else if (_NBField_ > 1 && _NDet_ == 1) legend->AddEntry(gResPlot[iPlot], Form("FieldMap: %s", FieldMapName[iPlot]), "P");
	     else legend->AddEntry(gResPlot[iPlot], "Resolution", "P");
	}
	
	if(_NY_ == 2) EfficiencyMinimum->Draw("same");


    if (_DRAWPWGReq_)
    {
		int WhichPlot = 99;
		if (_NY_ == 0 && !_PT_  )
		{
			if (AngleValues[0] == -3.5 && AngleValues[1] == -2.5) WhichPlot = 0;
			else if (AngleValues[0] == -2.5 && AngleValues[1] == -1.0) WhichPlot = 1;
			else if (AngleValues[0] == -1.0 && AngleValues[1] == 1.0) WhichPlot = 2;
			else if (AngleValues[0] == 1.0 && AngleValues[1] == 2.5) WhichPlot = 3;
			else if (AngleValues[0] == 2.5 && AngleValues[1] == 3.5) WhichPlot = 4;
			else {cout << "No Valid Angle Range for PWG Requirements" << endl;}
			if (WhichPlot != 99)
			{
				PWGReq[WhichPlot]->SetLineColor(kBlack);
				PWGReq[WhichPlot]->Draw("same");
				legend->AddEntry(PWGReq[WhichPlot], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , AngleValues[0], AngleValues[1]));
    	 	}
		}
		if (_NY_ == 1 && _PT_ )
		{
			if (AngleValues[0] == -3. && AngleValues[1] == -2.5) WhichPlot = 0;
			else if (AngleValues[0] == -2.5 && AngleValues[1] == -1.0) WhichPlot = 1;
			else if (AngleValues[0] == -1.0 && AngleValues[1] == 1.0) WhichPlot = 2;
			else if (AngleValues[0] == 1.0 && AngleValues[1] == 2.5) WhichPlot = 3;
			else if (AngleValues[0] == 2.5 && AngleValues[1] == 3.0) WhichPlot = 4;
			else if (AngleValues[0] == 3.0 && AngleValues[1] == 3.5) WhichPlot = 5;
			else {cout << "No Valid Angle Range for PWG Requirements" << endl;}
			if (WhichPlot != 99)
			{
				PWGReq_DCA2D[WhichPlot]->SetLineColor(kBlack);
				PWGReq_DCA2D[WhichPlot]->Draw("same");
				legend->AddEntry(PWGReq_DCA2D[WhichPlot], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , AngleValues[0], AngleValues[1]));
			}	
		}
	}

	legend->SetTextFont(52);
	legend->SetTextSize(0.04);
	legend->SetTextColor(kBlack);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	legend->Draw();

	//Add saving of plots to root file and as pdf to folder
	TFile Plotoutfile("./Plots/Plots.root", "UPDATE");
	c1->Write();
	Plotoutfile.Close();
	c1->SaveAs(Form("./Plots/cPlot_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf.pdf", PMin, PMax, EtaMin, EtaMax));

	


	//Saving the individual fits to a pdf so they can be easily investigated in case of issues

	TFile *FitInFile = new TFile("./Output/Fits.root", "read");
	TCanvas *cPS = new TCanvas("cPS", "", 100, 100, 600, 800);


	for (int iD = 0; iD < _NDet_; iD++)
	{
		for (int iB = 0; iB < _NBField_; iB++) 
		{
			TString SetupName = Form("%s_%s", DetVers[iD].c_str(), BField[iB].c_str());
			TPDF *ps = new TPDF(Form("./Plots/FittedHistograms_%s.pdf", SetupName.Data()));
			cPS->Clear();
			cPS->Divide(2,3);
			int padIndex = 1;
			
		  	for (int iP = 0; iP < _NP_; iP++)
			{
				for (int iA = 0; iA < _NEta_; iA++)
		        	{
					if (padIndex > 6)
					{
						cPS->Clear();
						cPS->Divide(2,3);
						padIndex = 1;
					}
		        		//Construct Histo Name
					TString HistName;
					if(_NY_ == 0)
					{
						if (_PT_) HistName = Form("h_ptRes_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
						else HistName = Form("h_momRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
					}
					else if (_NY_ == 1)
					{
						if (_PT_) HistName = Form("h_dca2dpt_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
						else HistName = Form("h_dca2d_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
					}  
					else if (_NY_ == 3)
					{
						if (_PT_) HistName = Form("h_%s%sRes_%s_%s_Pt_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",ProjectionLocation[_PROJ_].c_str(),ProjectionVariable[_PVAR_].c_str(), DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
						else HistName = Form("h_%s%sRes_%s_%s_P_%0.2lf_%0.2lf_Eta_%0.2lf_%0.2lf",ProjectionLocation[_PROJ_].c_str(),ProjectionVariable[_PVAR_].c_str(), DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
					} 	
					//Get Histogram
		        		TH1D *Htmp = (TH1D*)FitInFile->Get(HistName);
					cPS->cd(padIndex);
					Htmp->Draw();
					cPS->Update();
					padIndex++;
								
				
				}
			}
			ps->Close();
		}
	}




	cout << "DONE!" << endl;
	return;
}

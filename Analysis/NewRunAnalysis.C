#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#define _TSIZE_ 0.06
#define _LSIZE_ 0.05
//---------------------------------------------------------------------------------------------------------------------------------------  
//============================================== Change Inputs Here =====================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  
//Here are the only parameters that need changing for most analyses

//Draws the PWG Requirement line on the plot
#define _DRAWPWGReq_ 0
//Select which requirement:    0: Far-Backward  1: Backward  2: Central  3: Forward  4: Far-Forward
#define _PWGPlot_ 4

//Option to draw functional fit to the plots of the form: f(p) = SQRT( (A[%]*p)^2 + (B[%])^2 )
#define _DRAWFITS_ 0

//perform some rebinning of the residual histograms based on rapidity before fitting to ensure a good fit 
//(see where it is used for more details, may need to tweak the precise behavior yourself depending on the distributions for your histograms) 
#define _REBIN_ 0
//number for rebinning
#define _BINS_ 2

//defines the variable to be plotted on the Y-axis:    0: resolution    1: dca2d    2: efficiency
#define _NY_ 0

// 0: Momentum 1: Pseudo-Rapidity 
//define which variable will be plotted on x axis
#define _NDIM_ 0

//One of the following MUST have only 1 bin (the one not chosen above as the x axis variable)

//Number of Eta Bins, and the minimum and maximum values
#define _NEta_ 1
double EtaMin = -2.5;
double EtaMax = -1.;


//Number of Momentum Bins, and the minimum and maximum values
#define _NP_ 29
double PMin = 1;
double PMax =30;
// Change from momentum to pt ( 0=p  1=pt)
#define _PT_ 1


//The following parts correspond to the different simulations you have run and want to process data for

//Give the number of detector configurations as well as the name used to identify it in the rootfile name and a longer description for labels
#define _NDet_ 1
std::string DetVers[_NDet_] = {"PtBSi"};
const char *DetectorFullName[_NDet_] = {"Berkeley Si Hybrid"};
//std::string DetVers[_NDet_] = {"BSiONLY", "BSiNoGEM","BSiNoBMT", "BerkeleySi"};
//const char *DetectorFullName[_NDet_] = {"Berkeley Silicon Tracker Only", "Berkeley Silicon Hybrid without GEM Disks", "Berkeley Silicon Hybrid without BMT", "Berkeley Silicon Hybrid"};
//std::string DetVers[_NDet_] = {"Nominal", "NewMaterialFwdCombo", "WideGEM5Si"};
//const char *DetectorFullName[_NDet_] = {"Nominal","5 Si Disks + Wide GEMs with RICH+AuLayer", "5 Si Disks + Wide GEMs"};

//Give the number of field maps as well as the name used to identify it in the rootfile name and a longer description for labels
#define _NBField_ 4
std::string BField[_NBField_] = {"ATHENA", "NewBeAST", "Beast", "B_3.0T"};
const char *FieldMapName[_NBField_] = {"05-28 BeAST Field Map Update", "05-07 BeAST Field Map Update", "Original BeAST Field Map", "Uniform 3.0T"};
//std::string BField[_NBField_] = {"Beast", "ATHENA", "B_3.0T"};
//const char *FieldMapName[_NBField_] = {"Beast", "ATHENA", "Uniform 3.0T"};

//Set the number of plots you will compare (Currently it is only written to compare different BFields/Detectors on the same plot)
#define _NPLOTS_ _NBField_

//The maximum value for the y-axis. This will automatically be increased to be larger than the largest y-value if it is set too low
double ymax =0;// 1.6;

//---------------------------------------------------------------------------------------------------------------------------------------  
//=======================================================================================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  



//defining some global variables/arrays
double AngleValues[_NEta_ + 1];
double MomentumValues[_NP_ + 1];
int DimensionArray[2] = { _NP_, _NEta_};
#define _XMIN_ (1-_NDIM_)*(PMin-1)+(_NDIM_*(EtaMin-0.5))
#define _XMAX_ (1-_NDIM_)*(PMax+1)+(_NDIM_*(EtaMax+0.5))


	int effpt[_NDet_][_NBField_][_NEta_][_NP_] = {0};
	int effDenompt[_NDet_][_NBField_][_NEta_][_NP_] = {0};	
	int eff[_NDet_][_NBField_][_NEta_][_NP_] = {0};
	int effDenom[_NDet_][_NBField_][_NEta_][_NP_] = {0};	

void MakeHistogram()
{

	//Defining Histograms and variables
	TH1D *h_momRes[_NDet_][_NBField_][_NEta_][_NP_], *h_ptRes[_NDet_][_NBField_][_NEta_][_NP_];
	TH2D *h_nHits_momBin[_NDet_][_NBField_][_NEta_][_NP_], *h_nHits_ptBin[_NDet_][_NBField_][_NEta_][_NP_];

	TH1D *h_dca2d[_NDet_][_NBField_][_NEta_][_NP_]; 
	//pt binned versions
	TH1D *h_dca2dpt[_NDet_][_NBField_][_NEta_][_NP_]; 
	

	//TEMP
	TH1D *hETA[_NDet_][_NBField_][_NEta_][_NP_];
	TH1D *hrETA[_NDet_][_NBField_][_NEta_][_NP_];

	for (int iDet = 0; iDet < _NDet_; iDet++)
	{
		for (int iB = 0; iB < _NBField_; iB++)
		{
			for (int iEta = 0; iEta < _NEta_; iEta++)
			{
				for (int iP = 0; iP < _NP_; iP++)
				{
					//Creating the histograms with unique names, these names will be used later to grab the proper histogram for analysis
					h_momRes[iDet][iB][iEta][iP] = new TH1D(Form("h_momRes_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("dp/p for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1, 1);
					h_ptRes[iDet][iB][iEta][iP] = new TH1D(Form("h_ptRes_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("dpt/pt for %0.1lf < pt < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1, 1);
					h_nHits_momBin[iDet][iB][iEta][iP] = new TH2D(Form("h_nHits_MomentumBin_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf", DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Hits in each detector for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 5, 0, 5, 10, 0, 10);
					h_nHits_ptBin[iDet][iB][iEta][iP] = new TH2D(Form("h_nHits_PtBin_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf", DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Hits in each detector for %0.1lf < pt < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 5, 0, 5, 10, 0, 10);


					hETA[iDet][iB][iEta][iP] = new TH1D(Form("hETA_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("ETA for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1, 1);
					hrETA[iDet][iB][iEta][iP] = new TH1D(Form("hRECOETA_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("RECO ETA for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1, 1);

					h_dca2d[iDet][iB][iEta][iP] = new TH1D(Form("h_dca2d_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("DCA2D for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -0.1, 0.1);
					h_dca2dpt[iDet][iB][iEta][iP] = new TH1D(Form("h_dca2dpt_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("DCA2D for %0.1lf < p_t < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -0.1, 0.1);
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


			// run the "analysis" simply getting tracks and filling the residual histograms
			for (int i = 0; i < nentries; i++)
			{
				tree->GetEntry(i);
				
				//retrieving the momenta
				TVector3 truthP(gpx, gpy, gpz);
				TVector3 recoP(px, py, pz);

				//Loop over eta as to store information in the properly indexed histogram
				for (int iEta = 0; iEta < _NEta_; iEta++)
				{
					//Loop over momentum as to store information in the properly indexed histogram
					for (int iP = 0; iP < _NP_; iP++)
					{
						//theres probably a better way than loops and if-statements to determine the proper indexed histogram for a given track, but this works
						if (truthP.Eta() >= AngleValues[iEta] && truthP.Eta() <= AngleValues[iEta+1]) 
						{
							//filling histograms which are binned in momentum
							if (truthP.Mag() >= MomentumValues[iP] && truthP.Mag() <= MomentumValues[iP+1])
							{
								effDenom[iDet][iB][iEta][iP] +=1;
								//making sure there is a reconstructed track
								if (trackID != -9999)
								{
									h_momRes[iDet][iB][iEta][iP]->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(0., nHits_SVTX);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(1., nHits_BARR);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(2., nHits_BMT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(3., nHits_FGT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(4., nHits_FBST);
								
									eff[iDet][iB][iEta][iP] += 1;	
									h_dca2d[iDet][iB][iEta][iP]->Fill(dca2d);
								
									hETA[iDet][iB][iEta][iP]->Fill( truthP.Eta() );  
									hrETA[iDet][iB][iEta][iP]->Fill( recoP.Eta() );  
								}
							}
							//filling histograms which are binned in pt
							if (truthP.Pt() >= MomentumValues[iP] && truthP.Pt() <= MomentumValues[iP+1])
							{
								effDenompt[iDet][iB][iEta][iP] +=1;
								//making sure there is a reconstructed track
								if (trackID != -9999)
								{
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
						}//end Eta if

					} //end loop over momentum values
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
					
					hETA[iDet][iB][iEta][iP]->Write();
					hrETA[iDet][iB][iEta][iP]->Write();
					


				}
			}
			//close the output histogram file
			outfile.Close();

			//close the current file to prepare for new file to be opened
			myFile.Close();
		}//end loop over B Field options
	}//end loop over the detector configurations

return;

}

void MakePlot()
{

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
	double DCA2D[_NEta_][_NP_][_NDet_][_NBField_]; 
	double DCA2DErr[_NEta_][_NP_][_NDet_][_NBField_]; 



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
					if (_PT_) HistName = Form("h_ptRes_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
					else HistName = Form("h_momRes_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);

					//Construct DCA Histo Name
					TString DCAHistName;
					if (_PT_) DCAHistName = Form("h_dca2dpt_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
					else DCAHistName = Form("h_dca2d_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);
		        	    	
					//Get Histogram
		        	    	cout << "Getting: " << HistName << endl;
		        	    	TH1D *Htmp = (TH1D*)histoFile->Get(HistName);
					if(!Htmp){ cout << "Histogram not found: " << HistName << endl; break;}
		        	    	
					cout << "Getting: " << DCAHistName << endl;
		        	    	TH1D *Htmp2 = (TH1D*)histoFile->Get(DCAHistName);
					if(!Htmp2){ cout << "Histogram not found: " << DCAHistName << endl; break;}


				    	//Tweaking the histograms to attain a good fit (resolutions degrade at high eta, and the binning becomes too fine)
				    	//This part is related to the art of proper histogram binning. One size does not fit all. You may need to tweak some of this yourself to ensure proper fitting
					if (_REBIN_)
				    	{
				    		if (TMath::Abs(AngleValues[iA]) < 2.5) Htmp->GetXaxis()->SetRangeUser(-0.15, 0.15);
				    		else
				    		{
							Htmp->GetXaxis()->SetRangeUser(-0.3,0.3);
							Htmp->Rebin(_BINS_);
						}
					}
					else  Htmp->GetXaxis()->SetRangeUser(-0.15, 0.15); //for most tracks, the residuals don't exceed 15%
		
					//Extract the value from the histogram and store in array	
			        	TF1 *gausFit = new TF1("gausFit", "gaus");
			        	gausFit->SetParameter(1, 0);
	       		     		Htmp->Fit(gausFit, "Q");
			        	double sigma = gausFit->GetParameter(2);
			        	double sigmaErr = gausFit->GetParError(2);
			        	if (sigma == 0) cout << "ERROR! - Bad Value for : " << HistName << endl;

				    	Resolution[iA][iP][iD][iB] = sigma*100; //make it a percent
	            			ResErr[iA][iP][iD][iB] = sigmaErr*100; //make it a percent

			        	TF1 *gausFit2 = new TF1("gausFit2", "gaus");
			        	gausFit2->SetParameter(1, 0);
	       		     		Htmp2->Fit(gausFit2, "Q");
					double DCAsigma = gausFit2->GetParameter(2);
					double DCAerr = gausFit2->GetParError(2);
										
					DCA2D[iA][iP][iD][iB] = DCAsigma*10000; //convert to um
					DCA2DErr[iA][iP][iD][iB] = DCAerr*10000; // convert to um
					
				//	cout << "=============================================================" << endl;
				//	cout << "EFFICIENCY CALCULATION" << endl;					
				//	cout << "Reco Tracks: " << effpt[iA][iP][iD][iB] << endl;
				//	cout << "Generated Tracks: " << effDenompt[iA][iP][iD][iB] << endl;
				//	cout << "=============================================================" << endl;

					if(_PT_ && effDenompt[iA][iP][iD][iB] != 0 ) Efficiency[iA][iP][iD][iB]= effpt[iA][iP][iD][iB]/effDenompt[iA][iP][iD][iB];
 					else if (!_PT_ && effDenom[iA][iP][iD][iB] != 0 ) Efficiency[iA][iP][iD][iB]= eff[iA][iP][iD][iB]/effDenom[iA][iP][iD][iB];
					else {cout << "Division by 0! Setting Efficiency to 0" << endl; Efficiency[iA][iP][iD][iB] = 0;}
	
					TFile outfile("./Output/Fits.root", "UPDATE");
					Htmp->Write();
					Htmp2->Write();
					outfile.Close();

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
    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 500);
    c1->UseCurrentStyle();
    c1->SetBorderMode(0);
    c1->SetFrameBorderMode(0);
    c1->SetFrameLineColor(kWhite);
    c1->SetFillColor(kWhite);
    c1->cd();


    //PWG Requirements from the Yellow Report (Table 11.25)
    TF1 *PWGReq[5];
    PWGReq[0] = new TF1("PWGReq_-3.5_-2.5", "TMath::Sqrt((0.10*x)^2+0.5^2)", 1, 30);
    PWGReq[1] = new TF1("PWGReq_-2.5_-1.0", "TMath::Sqrt((0.05*x)^2+0.5^2)", 1, 30);
    PWGReq[2] = new TF1("PWGReq_-1.0_1.0" , "TMath::Sqrt((0.05*x)^2+0.5^2)", 1, 30);
    PWGReq[3] = new TF1("PWGReq_1.0_2.5"  , "TMath::Sqrt((0.05*x)^2+1.0^2)", 1, 30);
    PWGReq[4] = new TF1("PWGReq_2.5_3.5"  , "TMath::Sqrt((0.10*x)^2+2.0^2)", 1, 30);



double res[DimensionArray[_NDIM_]];
double reserr[DimensionArray[_NDIM_]];
int countPlots = 0;
TF1 *ReqFit[_NPLOTS_];
TGraphErrors *gResPlot[_NPLOTS_];

//can compare up to 5 different detector configurations before this will give an error
EColor ColorArray[5] = {kGreen, kBlue, kMagenta, kRed, kCyan};


for (int iD = 0; iD < _NDet_; iD++)
{
	for (int iB = 0; iB < _NBField_; iB++) 
	{
		//Index is reset for the new plot
		int index=0;
	 	for (int iP = 0; iP < _NP_; iP++)
	      	{
		  	for (int iA = 0; iA < _NEta_; iA++)
	          	{

			  
				if (_NY_ == 0)
   	  	  	  	{
					res[index] = Resolution[iA][iP][iD][iB];
   	  	  	  		reserr[index] = ResErr[iA][iP][iD][iB];
			  	}
				else if(_NY_ == 1)
				{
					res[index] = DCA2D[iA][iP][iD][iB];
					reserr[index] = DCA2DErr[iA][iP][iD][iB];
				}
				else if (_NY_ == 2)
				{
					res[index] = Efficiency[iA][iP][iD][iB];
					reserr[index] = 0;
				}
				else {cout << "No Valid Choice of Y variable. Exiting..." << endl; break;}
				if ( res[index] > ymax) ymax = res[index];
			  	index++;
			}//End Angle Loop
		}//End Momentum Loop
     		
		ReqFit[countPlots] = new TF1(Form("ReqFit_%d", countPlots), " TMath::Sqrt(([0]*x)^2+[1]^2)", 1 , 31);
		
	
		//the plot is constructed from the calculated bin centers (from above) and the resolution values and errors extracted from the full multidimensional array 	
		if (_NDIM_ == 0 || _NDIM_ == 2) gResPlot[countPlots] = new TGraphErrors(DimensionArray[_NDIM_], MomentumBinCenter , res, 0, reserr );
		else if (_NDIM_ == 1) gResPlot[countPlots] = new TGraphErrors(DimensionArray[_NDIM_], AngleBinCenter , res, 0, reserr  );
 		else { cout << "No plotting variable set. Exiting" << endl; break;}

		//setting some plot options
		gResPlot[countPlots]->SetTitle("");
 		gResPlot[countPlots]->SetMarkerSize(1.3);
 		gResPlot[countPlots]->SetMarkerStyle(21);
 		gResPlot[countPlots]->SetMarkerColor(ColorArray[countPlots]);
	     	gResPlot[countPlots]->SetLineColor(ColorArray[countPlots]);
	     	gResPlot[countPlots]->SetLineWidth(2);
	     	gResPlot[countPlots]->SetLineStyle(7);
	     	ReqFit[countPlots]->SetParameter(0,0.04);
	     	ReqFit[countPlots]->SetParName(0,"A");
	     	ReqFit[countPlots]->SetParameter(1,1.0);
	     	ReqFit[countPlots]->SetParName(1,"B");
	     	ReqFit[countPlots]->SetLineColor(ColorArray[countPlots]+1);

		//increment plot index
		countPlots++;

	}//End B Field Loop
} //End Detector Loop


    //Creating a dummy histogram which will be formatted
    TH1D *hdum  = new TH1D("hdum",   "", 15, _XMIN_, _XMAX_);


    hdum->GetXaxis()->SetLabelFont(52);
    hdum->GetYaxis()->SetLabelFont(52);
    hdum->GetXaxis()->SetTitleFont(52);
    hdum->GetYaxis()->SetTitleFont(52);

    if (_NDIM_ == 0 && !_PT_) hdum->GetXaxis()->SetTitle("p [GeV]");
    else if (_NDIM_ == 0 && _PT_) hdum->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (_NDIM_ == 1) hdum->GetXaxis()->SetTitle("#eta");

    hdum->GetXaxis()->SetTitleSize(_TSIZE_);
    hdum->GetXaxis()->SetLabelSize(_LSIZE_);
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
    
    hdum->GetYaxis()->SetTitleSize(_TSIZE_);
    hdum->GetYaxis()->SetLabelSize(_LSIZE_);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(0.75);
    hdum->SetMinimum( 0 );
    hdum->SetMaximum( ymax+0.25*ymax );

 
    //Title will be automatically set (since either p/eta is integrated over, it will be noted in the title) 
    if (_NY_ ==0 )
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "Resolution for tracks with %0.1lf < #eta < %0.1lf", EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "Resolution for tracks with %0.1lf < p < %0.1lf", PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "Resolution for tracks with %0.1lf < p_T < %0.1lf", PMin, PMax));
    }
    else if (_NY_ == 1 )
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "DCA2D for tracks with %0.1lf < #eta < %0.1lf", EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "DCA2D for tracks with %0.1lf < p < %0.1lf", PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "DCA2D for tracks with %0.1lf < p_T < %0.1lf", PMin, PMax));
    }
    else if (_NY_ == 2 )
    {
	if (_NDIM_ == 0) hdum->SetTitle(Form( "Efficiency for tracks with %0.1lf < #eta < %0.1lf", EtaMin, EtaMax));
	if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "Efficiency for tracks with %0.1lf < p < %0.1lf", PMin, PMax));
	if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "Efficiency for tracks with %0.1lf < p_T < %0.1lf", PMin, PMax));
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
     

     if(_DRAWFITS_) gResPlot[iPlot]->Fit(ReqFit[iPlot], "R");
     gResPlot[iPlot]->Draw("PLsame"); //replace L with C for smooth line

     //legend->AddEntry(gResPlot[iPlot], Form("Hybrid Detector Resolution for %0.1lf < #eta < %0.1lf", AngleValues[0], AngleValues[1]  ), "P");
 
     if (_NDet_ > 1 && _NBField_ == 1 ) legend->AddEntry(gResPlot[iPlot], Form("Detector Version: %s", DetectorFullName[iPlot]), "P");
     else if (_NBField_ > 1 && _NDet_ == 1) legend->AddEntry(gResPlot[iPlot], Form("FieldMap: %s", FieldMapName[iPlot]), "P");
     else legend->AddEntry(gResPlot[iPlot], "Resolution", "P");
}

     if (_DRAWPWGReq_)
     {
	PWGReq[_PWGPlot_]->SetLineColor(kBlack);
	PWGReq[_PWGPlot_]->Draw("same");
	//This legend entry assumes that you have set the eta range to the same range specified in the YR table 11.25
	legend->AddEntry(PWGReq[_PWGPlot_], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , AngleValues[0], AngleValues[1]));
     }

     legend->SetTextFont(52);
     legend->SetTextSize(0.04);
     legend->SetTextColor(kBlack);
     legend->SetFillColor(kWhite);
     legend->SetLineColor(kWhite);
     legend->Draw();


  return ; 

	
}

void NewRunAnalysis()
{

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
	MakeHistogram();
	//Reads the histograms from the histogram.root file and fits them with a gaussian to extract the resolution and plots them
	MakePlot();
	return;

}

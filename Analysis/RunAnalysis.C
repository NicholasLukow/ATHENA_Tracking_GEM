#define Analyze_cxx
#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


//void RunAnalysis(const char* infile = "./InTest_highP.root")
void RunAnalysis( double pmin, double pmax, double etamin, double etamax, const char* detversion, const char* BField)
{



//extract information from the filename
//Alternatively, take the information as input, and construct the filename from it.

TString infile = Form("../Output/P_%0.0lf_%0.0lf_Eta_%0.1lf_%0.1lf_%s_%s_FastSimEval.root", pmin, pmax, etamin, etamax,detversion, BField);

//check if file exists
if(gSystem->AccessPathName(infile))
{
        std::cout << "file " << infile << " does not exist" << std::endl;
	return;
} else 
{
        std::cout << "Found File: " << infile << std::endl;
}

//get the file
TFile myFile(infile.Data());

   
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
//Run the Analysis
	
//Defining Histograms and variables
	TH1D *h_momRes = new TH1D("h_momRes", "(reco_p - truth_p)/truth_p", 1000, -1, 1);
	TH1D *h_ptRes = new TH1D("h_ptRes", "(reco_pt - truth_pt)/truth_pt", 1000, -1, 1);
	
	


for (int i = 0; i < nentries; i++)
{
	tree->GetEntry(i);
	//fill in analysis here
	TVector3 truthP(gpx, gpy, gpz);
	TVector3 recoP(px, py, pz);
	if (trackID != -9999)
	{
		h_momRes->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  	
		h_ptRes->Fill( (recoP.Pt() - truthP.Pt())/truthP.Pt());  	
		
	}	

}
	TFile outfile("./Output/Histogram.root", "UPDATE");

	//Here I should insert lines Getting the 2D histograms I want to update

	
	//add line to fit 1D histogram of delta p to get sigmap
	TF1 *gausfitmom = new TF1("gausfitmom", "gaus");
	TF1 *gausfitpt = new TF1("gausfitpt", "gaus");
	h_momRes->Fit(gausfitmom, "Q");
	h_ptRes->Fit(gausfitpt, "Q");
	double sigmaP = gausfitmom->GetParameter(2);
	double sigmaPt = gausfitpt->GetParameter(2);

	h_momRes->SetNameTitle(Form("h_momRes_P_%0.1f_%0.1lf_Eta_%0.1lf_%0.1lf_%s_%s", pmin, pmax, etamin, etamax, detversion,BField), Form("(reco_p - truth_p)/truth_p for tracks with pt between %0.1f - %0.1lf GeV/c and %0.1lf < eta < %0.1lf - Detector: %s", pmin, pmax, etamin, etamax,detversion) );
	
	h_ptRes->SetNameTitle(Form("h_ptRes_P_%0.1f_%0.1lf_Eta_%0.1lf_%0.1lf_%s_%s",pmin, pmax, etamin, etamax, detversion,BField), Form("(reco_pt - truth_pt)/truth_pt for tracks with pt between %0.1f - %0.1lf GeV/c and %0.1lf < eta < %0.1lf - Detector: %s", pmin, pmax, etamin, etamax,detversion) );

	h_momRes->Write();

	h_ptRes->Write();
	
	outfile.Close();
	//infile.Close();

	

//Close the input file
myFile.Close();
}

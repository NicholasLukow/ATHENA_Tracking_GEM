#define Analyze_cxx

#include <iostream>
#include <fstream>
using namespace std;
#include <TVector3.h>
#include "Analyze.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"

void Analyze::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Analyze.C
//      root> Analyze t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

	//TFile infile(Form("./%s.root", filename));

	ofstream outtextfile;
	outtextfile.open("Out.txt", std::ofstream::app);

	//Setup code here
	TH1D *h_momRes = new TH1D("h_momRes", "Momentum Resolution", 100, -0.1, 0.1);
	TH2D *h2_momRes_vs_mom = new TH2D("h2_momRes_vs_mom", "Momentum Resolution vs Momentum", 41, -0.5, 40.5, 1000, -1, 1); 
	TH2D *h2_momRes_vs_eta = new TH2D("h2_momRes_vs_eta", "Momentum Resolution vs Eta", 70, -3.5, 3.5, 1000, -1, 1); 
	
	TH1D *h_ptRes = new TH1D("h_ptRes", "Pt Resolution", 100, -0.1, 0.1);
	TH2D *h2_ptRes_vs_pt = new TH2D("h2_ptRes_vs_mom", "Pt Resolution vs Pt", 41, -0.5, 40.5, 1000, -1, 1); 
	TH2D *h2_ptRes_vs_eta = new TH2D("h2_ptRes_vs_eta", "Pt Resolution vs Eta", 70, -3.5, 3.5, 1000, -1, 1); 

  	float meanPt = 0;	
	float meanP = 0;
	int trackCount = 0;
	float sigmaP = 0;
	float sigmaPt = 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	TVector3 truthP(gpx, gpy, gpz);
	TVector3 recoP(px, py, pz);
	if (trackID != -9999)
	{
		h_momRes->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  	
		h2_momRes_vs_mom->Fill(truthP.Mag(), (recoP.Mag() - truthP.Mag())/truthP.Mag());  	
		h2_momRes_vs_eta->Fill(truthP.Eta(), (recoP.Mag() - truthP.Mag())/truthP.Mag());  	
		
		h_ptRes->Fill( (recoP.Pt() - truthP.Pt())/truthP.Pt());  	
		h2_ptRes_vs_pt->Fill(truthP.Pt(), (recoP.Pt() - truthP.Pt())/truthP.Pt());  	
		h2_ptRes_vs_eta->Fill(truthP.Eta(), (recoP.Pt() - truthP.Pt())/truthP.Pt());  	
		
		trackCount++;
		meanP += truthP.Mag();
		meanPt += truthP.Pt();
	}
	//Analysis code
	

    }

	//wrap-up code
	TFile outfile("Histogram.root", "UPDATE");
	
	//add line to fit 1D histogram of delta p to get sigmap
	TF1 *gausfitmom = new TF1("gausfitmom", "gaus");
	TF1 *gausfitpt = new TF1("gausfitpt", "gaus");
	h_momRes->Fit(gausfitmom, "NQ");
	h_ptRes->Fit(gausfitpt, "NQ");
	sigmaP = gausfitmom->GetParameter(2);
	sigmaPt = gausfitpt->GetParameter(2);
	meanP = meanP/trackCount;	
	meanPt = meanPt/trackCount;	


	h_momRes->SetNameTitle(Form("h_momRes_%0.1fGeV",meanP), Form("Momentum Resolution for tracks with mean momentum of %0.1f", meanP) );
	h2_momRes_vs_mom->SetNameTitle(Form("h2_momRes_vs_mom_%0.1fGeV",meanP), Form("Momentum Resolution vs Momentum for tracks with mean momentum of %0.1f", meanP) );
	h2_momRes_vs_eta->SetNameTitle(Form("h2_momRes_vs_eta_%0.1fGeV",meanP), Form("Momentum Resolution vs Eta for tracks with mean momentum of %0.1f", meanP) );
	
	h_ptRes->SetNameTitle(Form("h_ptRes_%0.1fGeV",meanPt), Form("Pt Resolution for tracks with mean Pt of %0.1f", meanPt) );
	h2_ptRes_vs_pt->SetNameTitle(Form("h2_ptRes_vs_pt_%0.1fGeV",meanPt), Form("Pt Resolution vs Pt for tracks with mean pt of %0.1f", meanPt) );
	h2_ptRes_vs_eta->SetNameTitle(Form("h2_ptRes_vs_eta_%0.1fGeV",meanPt), Form("Pt Resolution vs Eta for tracks with mean pt of %0.1f", meanPt) );

	h_momRes->Write();
	h2_momRes_vs_mom->Write();
	h2_momRes_vs_eta->Write();

	h_ptRes->Write();
	h2_ptRes_vs_pt->Write();
	h2_ptRes_vs_eta->Write();
	
	outfile.Close();
	//infile.Close();

	
	cout << "TESTING" << endl;
	outtextfile << meanP << ", " << sigmaP << " , " << meanPt << ", " << sigmaPt <<  endl;
	outtextfile.close();
}

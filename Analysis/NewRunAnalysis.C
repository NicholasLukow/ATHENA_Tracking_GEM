#define Analyze_cxx
#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#define _TSIZE_ 0.06
#define _LSIZE_ 0.05
//---------------------------------------------------------------------------------------------------------------------------------------  
//============================================== Change Inputs Here =====================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  

#define _DRAWPWGReq_ 0
#define _PWGPlot_ 3

#define _DRAWFITS_ 1
#define _REBIN_ 1


#define _NBField_ 1
#define _NDet_ 2
#define _NEta_ 1
#define _NP_ 14

//double AngleValues[_NEta_+1] = {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
double AngleValues[_NEta_+1] = {1.0, 3.5};
double MomentumValues[_NP_+1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20};


std::string DetVers[_NDet_] = {"2ffg", "2ffgFull"};
const char *DetectorFullName[_NDet_] = { "Nominal", "NewTestSim"};
std::string BField[_NBField_] = {"ATHENA"};
const char *FieldMapName[_NBField_] = {"ATHENA"};
//std::string BField[_NBField_] = {"Beast", "ATHENA", "B_3.0T"};
//const char *FieldMapName[_NBField_] = {"Beast", "ATHENA", "Uniform 3.0T"};


//For Plot
#define _NDIM_ _NP_
#define _NPLOTS_ _NDet_
#define _YMIN_ 0
#define _YMAX_ 2.5
#define _XMIN_ 0
#define _XMAX_ 21
//---------------------------------------------------------------------------------------------------------------------------------------  
//=======================================================================================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  

void MakeHistogram()
{

	//Defining Histograms and variables

	//Create many histograms based on the bins from the input at the beginning
	TH1D *h_momRes[_NDet_][_NBField_][_NEta_][_NP_], *h_ptRes[_NDet_][_NBField_][_NEta_][_NP_];
	TH2D *h_nHits_momBin[_NDet_][_NBField_][_NEta_][_NP_], *h_nHits_ptBin[_NDet_][_NBField_][_NEta_][_NP_];

	for (int iDet = 0; iDet < _NDet_; iDet++)
	{
		for (int iB = 0; iB < _NBField_; iB++)
		{
			for (int iEta = 0; iEta < _NEta_; iEta++)
			{
				for (int iP = 0; iP < _NP_; iP++)
				{
					h_momRes[iDet][iB][iEta][iP] = new TH1D(Form("h_momRes_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("dp/p for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1, 1);
					h_ptRes[iDet][iB][iEta][iP] = new TH1D(Form("h_ptRes_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("dpt/pt for %0.1lf < pt < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 1000, -1, 1);
					h_nHits_momBin[iDet][iB][iEta][iP] = new TH2D(Form("h_nHits_MomentumBin_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf", DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Hits in each detector for %0.1lf < p < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 4, 0, 4, 10, 0, 10);
					h_nHits_ptBin[iDet][iB][iEta][iP] = new TH2D(Form("h_nHits_PtBin_%s_%s_Pt_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf", DetVers[iDet].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1]), Form("Hits in each detector for %0.1lf < pt < %0.1lf and %0.1lf < #eta < %0.1lf - Detector: %s  Field Map: %s", MomentumValues[iP], MomentumValues[iP+1], AngleValues[iEta], AngleValues[iEta+1], DetectorFullName[iDet], FieldMapName[iB]), 4, 0, 4, 10, 0, 10);

					//setting names for the hit histogram
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(1,"Si Vtx and Barrel");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(2,"Barrel MPGD Tracker");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(3,"GEM Disks");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(4,"Si Disks");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetYaxis()->SetTitle("Hits");
					h_nHits_momBin[iDet][iB][iEta][iP]->GetXaxis()->SetTitle("Detector");

					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(1,"Si Vtx and Barrel");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(2,"Barrel MPGD Tracker");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(3,"GEM Disks");
					h_nHits_ptBin[iDet][iB][iEta][iP]->GetXaxis()->SetBinLabel(4,"Si Disks");
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
   	Int_t	   nHits_FGT;
   	Int_t 	   nHits_FBST;
   	Int_t	   nHits_BARR;
   	Int_t 	   nHits_BMT;




	//StartDetector/Field Loop here! I will have new file for each Field/Detector configuration

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
			   tree->SetBranchAddress("nHit_G4HIT_FGT", &nHits_FGT);
			   tree->SetBranchAddress("nHit_G4HIT_FBST", &nHits_FBST);
			   tree->SetBranchAddress("nHit_G4HIT_BARR", &nHits_BARR);
			   tree->SetBranchAddress("nHit_G4HIT_BMT", &nHits_BMT);
			//Run the Analysis



			for (int i = 0; i < nentries; i++)
			{
				tree->GetEntry(i);
				//fill in analysis here
				TVector3 truthP(gpx, gpy, gpz);
				TVector3 recoP(px, py, pz);
				if (trackID != -9999)
				{
					//Determine the P and Eta index based on values!
					for (int iEta = 0; iEta < _NEta_; iEta++)
					{
						for (int iP = 0; iP < _NP_; iP++)
						{
							//theres almsot assuredly a better way to identify which momentum and eta bin the track belongs too, but this is an expedient solution
							if (truthP.Eta() >= AngleValues[iEta] && truthP.Eta() <= AngleValues[iEta+1]) 
							{
								if (truthP.Mag() >= MomentumValues[iP] && truthP.Mag() <= MomentumValues[iP+1])
								{
									h_momRes[iDet][iB][iEta][iP]->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(0., nHits_BARR);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(1., nHits_BMT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(2., nHits_FGT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(3., nHits_FBST);
								}
								if (truthP.Pt() >= MomentumValues[iP] && truthP.Pt() <= MomentumValues[iP+1])
								{
									h_ptRes[iDet][iB][iEta][iP]->Fill( (recoP.Pt() - truthP.Pt())/truthP.Pt());  
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(0., nHits_BARR);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(1., nHits_BMT);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(2., nHits_FGT);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(3., nHits_FBST);
								}
							} // check if track truth eta is within the values

						} // loop over momentum values
					} // loop over eta values

				} 

			} // end loop over events


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

				}
			}
			//close the output histogram file
			outfile.Close();

			//close the current file to prepare for new file to be opened
			myFile.Close();
		}//End loop over B Field options
	}// End loop over the detector configurations

return;

}

void MakePlot()
{

	double MomentumBinCenter[_NP_];
  	for (int i = 0; i < _NP_; i++)
  	{
		MomentumBinCenter[i] = 0.5*(MomentumValues[i]+MomentumValues[i+1]);
  	}

	//getting the created histograms
	TFile *histoFile = new TFile("./Output/Histogram.root","read");
	

	//Get values from histograms and store in arrays
	double Resolution[_NEta_][_NP_][_NDet_][_NBField_]; //[eta][momentum][detectors]
	double ResErr[_NEta_][_NP_][_NDet_][_NBField_]; //[eta][momentum][detectors]

	for (int iD = 0; iD < _NDet_; iD++)
	{
		for (int iA = 0; iA < _NEta_; iA++)
		{
		  	for (int iP = 0; iP < _NP_; iP++)
			{
		        for (int iB = 0; iB < _NBField_; iB++) 
		        {
					TString HistName = Form("h_momRes_%s_%s_P_%0.1lf_%0.1lf_Eta_%0.1lf_%0.1lf",DetVers[iD].c_str(), BField[iB].c_str(), MomentumValues[iP], MomentumValues[iP+1], AngleValues[iA], AngleValues[iA+1]);

		            //Construct Histo Name
		            //std:string HistName = "hMomentum_"+(std::string)Momentum[iP]+"GeV_"+(std::string)Angle[iA]+"Deg_"+(std::string)BField[iB]+"T";	
			    	//std::string HistName;
		 	    	//if (_PT_) HistName = "h_ptRes_P_"+(std::string)Momentum[iP]+"_"+(std::string)Momentum[iP+1]+"_Eta_"+(std::string)Angle[iA]+"_"+(std::string)Angle[iA+1]+"_"+(std::string)DetVers[iD]+"_"+(std::string)BField[iB];	
		 	    	//else HistName = "h_momRes_P_"+(std::string)Momentum[iP]+"_"+(std::string)Momentum[iP+1]+"_Eta_"+(std::string)Angle[iA]+"_"+(std::string)Angle[iA+1]+"_"+(std::string)DetVers[iD]+"_"+(std::string)BField[iB];	

		            //Get Histogram
		            cout << "Getting: " << HistName << endl;
		            TH1D *Htmp = (TH1D*)histoFile->Get(HistName);
					if(!Htmp) break;


			    //Tweaking the histograms to attain a good fit (resolutions degrade at high eta, and the binning becomes too fine)
			    if (_REBIN_)
			    {
			    	if (TMath::Abs(AngleValues[iA]) < 2.5) Htmp->GetXaxis()->SetRangeUser(-0.1, 0.1);
			    	else if(TMath::Abs(AngleValues[iA]) == 2.5 ) Htmp->GetXaxis()->SetRangeUser(-0.2,0.2);
			    	else
			    	{
					Htmp->GetXaxis()->SetRangeUser(-0.3,0.3);
					Htmp->Rebin(10);
				    }
			    }
			    else  Htmp->GetXaxis()->SetRangeUser(-0.15, 0.15);
	
				//Extract the value from the histogram and store in array	
		        TF1 *gausFit = new TF1("gausFit", "gaus");
		        gausFit->SetParameter(1, 0);
	            Htmp->Fit(gausFit, "Q");
		        double sigma = gausFit->GetParameter(2);
		        double sigmaErr = gausFit->GetParError(2);
		        if (sigma == 0) cout << "ERROR! - Bad Value for : " << HistName << endl;

			    Resolution[iA][iP][iD][iB] = sigma*100; //make it a percent
	            ResErr[iA][iP][iD][iB] = sigmaErr*100; //make it a percent
	
				TFile outfile("./Output/Fits.root", "UPDATE");
				Htmp->Write();
				outfile.Close();

		        }//End B Loop
		    } // End Momentum Loop
		} // End Angle Loop
	} // End Detector Loop


	//Take the array and make the plot

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


    //Old PWG Requirement curve (may need to update)
    //TF1 *PWGReq3p0 = new TF1("PWGReq3p0", " TMath::Sqrt((0.02*x)^2+0.5^2)", 1 , 30);

    TF1 *PWGReq[5];
    PWGReq[0] = new TF1("PWGReq_-3.5_-2.5", "TMath::Sqrt((0.10*x)^2+0.5^2)", 1, 20);
    PWGReq[1] = new TF1("PWGReq_-2.5_-1.0", "TMath::Sqrt((0.05*x)^2+0.5^2)", 1, 20);
    PWGReq[2] = new TF1("PWGReq_-1.0_1.0" , "TMath::Sqrt((0.05*x)^2+0.5^2)", 1, 20);
    PWGReq[3] = new TF1("PWGReq_1.0_2.5"  , "TMath::Sqrt((0.05*x)^2+1.0^2)", 1, 20);
    PWGReq[4] = new TF1("PWGReq_2.5_3.5"  , "TMath::Sqrt((0.10*x)^2+2.0^2)", 1, 20);

	// loop to extract values from array and store in arrays to be turned into plots
//_NDIM_ defines the number of points for the x axis


double res[_NDIM_];
double reserr[_NDIM_];
int countPlots = 0;
TF1 *ReqFit[_NPLOTS_];
TGraphErrors *gResPlot[_NPLOTS_];
EColor ColorArray[5] = {kGreen, kBlue, kMagenta, kRed, kCyan};


for (int iD = 0; iD < _NDet_; iD++)
{
	for (int iB = 0; iB < _NBField_; iB++) 
	{
		int index=0;
	 	for (int iP = 0; iP < _NP_; iP++)
	      	{
		  	for (int iA = 0; iA < _NEta_; iA++)
	          	{
   	  	  	  res[index] = Resolution[iA][iP][iD][iB];
   	  	  	  reserr[index] = ResErr[iA][iP][iD][iB];
			  index++;
			}//End Angle Loop
		}//End Momentum Loop
     		
		ReqFit[countPlots] = new TF1(Form("ReqFit_%d", countPlots), " TMath::Sqrt(([0]*x)^2+[1]^2)", 1 , 31);
		
		//Must edit this if changing the x axis variable of plot
		gResPlot[countPlots] = new TGraphErrors(_NDIM_, MomentumBinCenter , res, 0, reserr );
		//gResPlot[countPlots] = new TGraph(_NDIM_, AngleBinCenter , res  );
 		
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
		countPlots++;

	}//End B Field Loop
} //End Detector Loop
// end loop here

 

    TH1D *hdum  = new TH1D("hdum",   "", 15, _XMIN_, _XMAX_);


    hdum->GetXaxis()->SetLabelFont(52);
    hdum->GetYaxis()->SetLabelFont(52);
    hdum->GetXaxis()->SetTitleFont(52);
    hdum->GetYaxis()->SetTitleFont(52);

    //if (_PT_) hdum->GetXaxis()->SetTitle("p_{T} [GeV]");
    //else hdum->GetXaxis()->SetTitle("p [GeV]");
    hdum->GetXaxis()->SetTitle("p [GeV]");
//      hdum->GetXaxis()->SetTitle("#eta");

    //hdum->GetXaxis()->SetTitle("#eta");
    hdum->GetXaxis()->SetTitleSize(_TSIZE_);
    hdum->GetXaxis()->SetLabelSize(_LSIZE_);
    //if (_PT_) hdum->GetYaxis()->SetTitle("p_{T} resolution   #sigma_{p_{T}} /p_{T}  [%]");
    //else hdum->GetYaxis()->SetTitle("p resolution   #sigma_{p} /p  [%]");
    hdum->GetYaxis()->SetTitle("p resolution   #sigma_{p} /p  [%]");
    
    hdum->GetYaxis()->SetTitleSize(_TSIZE_);
    hdum->GetYaxis()->SetLabelSize(_LSIZE_);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(0.75);
    hdum->SetMinimum( _YMIN_);
    hdum->SetMaximum( _YMAX_);

    //hdum->GetXaxis()->SetNdivisions(408);
    //hdum->GetYaxis()->SetNdivisions(804);
    hdum->GetXaxis()->SetNdivisions(408);
    hdum->GetYaxis()->SetNdivisions(808);

    //c1->cd(1);
    gPad->SetGrid();
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFrameBorderMode(0);
    gPad->SetFrameLineColor(kWhite);


    hdum->Draw();


     //#if _NOW_

TLegend *legend = new TLegend(.15, .70, .45, .92);
for (int iPlot = 0; iPlot < _NPLOTS_; iPlot++)
{ 
     
     if(_DRAWFITS_) gResPlot[iPlot]->Fit(ReqFit[iPlot], "R");
     gResPlot[iPlot]->Draw("PLsame"); //replace L with C for smooth line

     legend->AddEntry(gResPlot[iPlot], Form("Hybrid Detector Resolution for %0.1lf < #eta < %0.1lf", AngleValues[0], AngleValues[1]  ), "P");
     //legend->AddEntry(gResPlot[iPlot], Form("Detector Version: %s", DetectorFullName[iPlot]), "P");
     //legend->AddEntry(gResPlot[iPlot], Form("FieldMap: %s", FieldMapName[iPlot]), "P");
}

     if (_DRAWPWGReq_)
     {
	PWGReq[_PWGPlot_]->SetLineColor(kBlack);
	PWGReq[_PWGPlot_]->Draw("same");
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
	//Read the output tree from the simulation and create histograms binned in eta and p/pt and saves them in an output histogram.root file
	MakeHistogram();
	//Reads the histograms from the histogram.root file and fits them with a gaussian to extract the resolution and plots them
	MakePlot();
	return;
}

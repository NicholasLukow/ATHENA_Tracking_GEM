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
//Here are the only parameters that need changing for most analyses

//Draws the PWG Requirement line on the plot
#define _DRAWPWGReq_ 0
//Select which requirement:    0: Far-Backward  1: Backward  2: Central  3: Forward  4: Far-Forward
#define _PWGPlot_ 3

//Option to draw functional fit to the plots of the form: f(p) = SQRT( (A[%]*p)^2 + (B[%])^2 )
#define _DRAWFITS_ 0

//perform some rebinning of the residual histograms based on rapidity before fitting to ensure a good fit 
//(see where it is used for more details, may need to tweak the precise behavior yourself depending on the distributions for your histograms) 
#define _REBIN_ 1



// 0: Momentum 1: Pseudo-Rapidity 
//define which variable will be plotted on x axis
#define _NDIM_ 0

//One of the following MUST have only 1 bin (the one not chosen above as the x axis variable)

//Number of Eta Bins, and the minimum and maximum values
#define _NEta_ 1
double EtaMin = -3.5;
double EtaMax = 3.5;


//Number of Momentum Bins, and the minimum and maximum values
#define _NP_ 30
double PMin = 1;
double PMax =30;
// Change from momentum to pt ( 0=p  1=pt)
#define _PT_ 0


//The following parts correspond to the different simulations you have run and want to process data for

//Give the number of detector configurations as well as the name used to identify it in the rootfile name and a longer description for labels
#define _NDet_ 2
std::string DetVers[_NDet_] = {"nominal", "WideGEM"};
const char *DetectorFullName[_NDet_] = { "Nominal", "Wide GEM"};

//Give the number of field maps as well as the name used to identify it in the rootfile name and a longer description for labels
#define _NBField_ 1
std::string BField[_NBField_] = {"ATHENA"};
const char *FieldMapName[_NBField_] = {"ATHENA"};
//std::string BField[_NBField_] = {"Beast", "ATHENA", "B_3.0T"};
//const char *FieldMapName[_NBField_] = {"Beast", "ATHENA", "Uniform 3.0T"};

//Set the number of plots you will compare (Currently it is only written to compare different BFields/Detectors on the same plot)
#define _NPLOTS_ _NDet_


//---------------------------------------------------------------------------------------------------------------------------------------  
//=======================================================================================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  



//defining some global variables/arrays
double ymax = 0;
double AngleValues[_NEta_ + 1];
double MomentumValues[_NP_ + 1];
int DimensionArray[2] = { _NP_, _NEta_};
#define _XMIN_ (1-_NDIM_)*(PMin-1)+(_NDIM_*(EtaMin-0.5))
#define _XMAX_ (1-_NDIM_)*(PMax+1)+(_NDIM_*(EtaMax+0.5))



void MakeHistogram()
{

	//Defining Histograms and variables
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
					//Creating the histograms with unique names, these names will be used later to grab the proper histogram for analysis
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
			   tree->SetBranchAddress("nHit_G4HIT_FGT", &nHits_FGT);
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

				//making sure there is a reconstructed track
				if (trackID != -9999)
				{
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
									h_momRes[iDet][iB][iEta][iP]->Fill( (recoP.Mag() - truthP.Mag())/truthP.Mag());  
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(0., nHits_BARR);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(1., nHits_BMT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(2., nHits_FGT);
									h_nHits_momBin[iDet][iB][iEta][iP]->Fill(3., nHits_FBST);
								}
								//filling histograms which are binned in pt
								if (truthP.Pt() >= MomentumValues[iP] && truthP.Pt() <= MomentumValues[iP+1])
								{
									h_ptRes[iDet][iB][iEta][iP]->Fill( (recoP.Pt() - truthP.Pt())/truthP.Pt());  
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(0., nHits_BARR);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(1., nHits_BMT);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(2., nHits_FGT);
									h_nHits_ptBin[iDet][iB][iEta][iP]->Fill(3., nHits_FBST);
								}
							}

						} //end loop over momentum values
					} //end loop over eta values

				} 

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


		        	    	//Get Histogram
		        	    	cout << "Getting: " << HistName << endl;
		        	    	TH1D *Htmp = (TH1D*)histoFile->Get(HistName);
					if(!Htmp){ cout << "Histogram not found: " << HistName << endl; break;}


				    	//Tweaking the histograms to attain a good fit (resolutions degrade at high eta, and the binning becomes too fine)
				    	//This part is related to the art of proper histogram binning. One size does not fit all. You may need to tweak some of this yourself to ensure proper fitting
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
					else  Htmp->GetXaxis()->SetRangeUser(-0.15, 0.15); //for most central rapidity tracks, the residuals don't exceed 15%
		
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
   	  	  	  res[index] = Resolution[iA][iP][iD][iB];
   	  	  	  reserr[index] = ResErr[iA][iP][iD][iB];
			  index++;
			  if ( Resolution[iA][iP][iD][iB] > ymax) ymax = Resolution[iA][iP][iD][iB];
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
    if (_PT_) hdum->GetYaxis()->SetTitle("p_{T} resolution   #sigma_{p_{T}} /p_{T}  [%]");
    else hdum->GetYaxis()->SetTitle("p resolution   #sigma_{p} /p  [%]");
    
    
    hdum->GetYaxis()->SetTitleSize(_TSIZE_);
    hdum->GetYaxis()->SetLabelSize(_LSIZE_);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(0.75);
    hdum->SetMinimum( 0 );
    hdum->SetMaximum( ymax+0.5 );
   
    //Title will be automatically set (since either p/eta is integrated over, it will be noted in the title) 
    if (_NDIM_ == 0) hdum->SetTitle(Form( "Resolution for tracks with %0.1lf < #eta < %0.1lf", EtaMin, EtaMax));
    if (_NDIM_ == 1 && !_PT_) hdum->SetTitle(Form( "Resolution for tracks with %0.1lf < p < %0.1lf", PMin, PMax));
    if (_NDIM_ == 1 && _PT_) hdum->SetTitle(Form( "Resolution for tracks with %0.1lf < p_T < %0.1lf", PMin, PMax));


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



TLegend *legend = new TLegend(.35, .70, .65, .92);
for (int iPlot = 0; iPlot < _NPLOTS_; iPlot++)
{ 
     

     if(_DRAWFITS_) gResPlot[iPlot]->Fit(ReqFit[iPlot], "R");
     gResPlot[iPlot]->Draw("PLsame"); //replace L with C for smooth line

     //legend->AddEntry(gResPlot[iPlot], Form("Hybrid Detector Resolution for %0.1lf < #eta < %0.1lf", AngleValues[0], AngleValues[1]  ), "P");
 
     if (_NDet_ > 1 ) legend->AddEntry(gResPlot[iPlot], Form("Detector Version: %s", DetectorFullName[iPlot]), "P");
     else if (_NBField_ > 1) legend->AddEntry(gResPlot[iPlot], Form("FieldMap: %s", FieldMapName[iPlot]), "P");
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

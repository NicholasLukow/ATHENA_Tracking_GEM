#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTree.h"
#include "TLine.h"


#define _TSIZE_ 0.06
#define _LSIZE_ 0.05


void Plotter(TString infile="./Output/Histogram.root")
{
  gStyle->SetOptStat(1); //1-shows statistics box (entries, mean, std)   0 - hides it
  gStyle->SetOptFit(1); // 1- shows the fit parameters box   0 - hides it

  //gStyle->SetPadGridX(0);
  //gStyle->SetPadGridY(0);




  //gStyle->SetStatW(0.2); gStyle->SetStatH(0.2);
  
  TFile *rootfile = new TFile(infile,"read");
  
  TH1D *Htmp;


//---------------------------------------------------------------------------------------------------------------------------------------  
//============================================== Change Inputs Here =====================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  
#define _NBField_ 2
#define _NDet_ 1
#define _NEta_ 14
double AngleValues[_NEta_+1] = {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
#define _NP_ 14
double MomentumValues[_NP_+1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20};

//std::string DetVers[_NDet_] = {"noGem" ,"1ffg", "2ffg"};
//const char *DetectorFullName[_NDet_] = {"No GEM Disks", "3Endcap/1Post-RICH GEM Disks", "3Endcap/2Post-RICH GEM Disks"};
std::string DetVers[_NDet_] = {"2ffg"};
const char *DetectorFullName[_NDet_] = { "3Endcap/2Post-RICH GEM Disks"};
std::string BField[_NBField_] = {"Beast", "ATHENA"};
const char *FieldMapName[_NBField_] = {"Beast", "ATHENA"};
//---------------------------------------------------------------------------------------------------------------------------------------  
//=======================================================================================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  



  //creating two arrays, one of type double for use in plotting, one of type string for name searching histograms
  std::string Momentum[_NP_+1];
  double MomentumBinCenter[_NP_];
  for (int i = 0; i < _NP_+1; i++)
  {
	Momentum[i] = Form("%0.1lf",MomentumValues[i]);
  	if (i<_NP_) MomentumBinCenter[i] = 0.5*(MomentumValues[i]+MomentumValues[i+1]);
  }
  
  //creating two arrays, one of type double for use in plotting, one of type string for name searching histograms
  std::string Angle[_NEta_+1];
  double AngleBinCenter[_NEta_];
  for (int i = 0; i < _NEta_+1; i++)
  {
	Angle[i] = Form("%0.1lf",AngleValues[i]);
  	if (i<_NEta_) AngleBinCenter[i] = 0.5*(AngleValues[i]+AngleValues[i+1]);
  }




  double Resolution[_NEta_][_NP_][_NDet_][_NBField_]; //[eta][momentum][detectors]

  for (int iD = 0; iD < _NDet_; iD++)
  {
	for (int iA = 0; iA < _NEta_; iA++)
	{
	  for (int iP = 0; iP < _NP_; iP++)
	      {
	          for (int iB = 0; iB < _NBField_; iB++) 
	          {

	            //Construct Histo Name
	            //std:string HistName = "hMomentum_"+(std::string)Momentum[iP]+"GeV_"+(std::string)Angle[iA]+"Deg_"+(std::string)BField[iB]+"T";	
	 	    std:string HistName = "h_ptRes_P_"+(std::string)Momentum[iP]+"_"+(std::string)Momentum[iP+1]+"_Eta_"+(std::string)Angle[iA]+"_"+(std::string)Angle[iA+1]+"_"+(std::string)DetVers[iD]+"_"+(std::string)BField[iB];	

	            //Get Histogram
	            cout << "Getting: " << HistName << endl;
	            Htmp = (TH1D*)rootfile->Get(HistName.c_str());
		    
		    
		    //Tweaking the histograms to attain a good fit (resolutions degrade at high eta, and the binning becomes too fine)
		    if (TMath::Abs(AngleValues[iA]) < 2.5) Htmp->GetXaxis()->SetRangeUser(-0.1, 0.1);
		    else if(TMath::Abs(AngleValues[iA]) == 2.5 ) Htmp->GetXaxis()->SetRangeUser(-0.2,0.2);
		    else
		    {
			Htmp->GetXaxis()->SetRangeUser(-0.3,0.3);
			Htmp->Rebin(10);
		    }
		    
                        //Extract the value from the histogram and store in array	
	                TF1 *gausFit = new TF1("gausFit", "gaus");
	                gausFit->SetParameter(1, 0);
	                Htmp->Fit(gausFit, "Q");
	                double sigma = gausFit->GetParameter(2);
	                if (sigma == 0) cout << "ERROR! - Bad Value for : " << HistName << endl;
			
    		        Resolution[iA][iP][iD][iB] = sigma*100; //make it a percent
  		
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


     //Old PWG Requirement curve (may need to update)
     TF1 *PWGReq3p0 = new TF1("PWGReq3p0", " TMath::Sqrt((0.02*x)^2+0.5^2)", 1 , 30);




// loop to extract values from array and store in arrays to be turned into plots
//_NDIM_ defines the number of points for the x axis
#define _NDIM_ _NP_
#define _NPLOTS_ _NBField_

#define _DRAWFITS_ 1

double res[_NDIM_];
int countPlots = 0;
TF1 *ReqFit[_NPLOTS_];
TGraph *gResPlot[_NPLOTS_];
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
			  index++;
			}//End Angle Loop
		}//End Momentum Loop
     		
		ReqFit[countPlots] = new TF1(Form("ReqFit_%d", countPlots), " TMath::Sqrt(([0]*x)^2+[1]^2)", 1 , 31);
		
		//Must edit this if changing the x axis variable of plot
		gResPlot[countPlots] = new TGraph(_NDIM_, MomentumBinCenter , res  );
 		
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

 

    TH1D *hdum  = new TH1D("hdum",   "", 15, 0.8, 3.5);


    hdum->GetXaxis()->SetLabelFont(52);
    hdum->GetYaxis()->SetLabelFont(52);
    hdum->GetXaxis()->SetTitleFont(52);
    hdum->GetYaxis()->SetTitleFont(52);

    hdum->GetXaxis()->SetTitle("p_{T} [GeV]");
    //hdum->GetXaxis()->SetTitle("#eta");
    hdum->GetXaxis()->SetTitleSize(_TSIZE_);
    hdum->GetXaxis()->SetLabelSize(_LSIZE_);
    hdum->GetYaxis()->SetTitle("p_{T} resolution   #sigma_{p_{T}} /p_{T}  [%]");
    hdum->GetYaxis()->SetTitleSize(_TSIZE_);
    hdum->GetYaxis()->SetLabelSize(_LSIZE_);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(0.75);
    hdum->SetMinimum( 0.00);
    hdum->SetMaximum( 10.0);

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

     //legend->AddEntry(gResPlot[iPlot], Form("Detector Version: %s", DetectorFullName[iPlot]), "P");
     legend->AddEntry(gResPlot[iPlot], Form("FieldMap: %s", FieldMapName[iPlot]), "P");
}
     legend->SetTextFont(52);
     legend->SetTextSize(0.04);
     legend->SetTextColor(kBlack);
     legend->SetFillColor(kWhite);
     legend->SetLineColor(kWhite);
     legend->Draw();


  return ; 
}

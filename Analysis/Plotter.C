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
#define _ETA_ 0 // 0: 1-1.5 1: 1.5-2 2:2-2.5 3: 2.5-3 4: 3-3.5


void Plotter(TString infile="./Histogram.root")
{
  gStyle->SetOptStat(1); //1-shows statistics box (entries, mean, std)   0 - hides it
  gStyle->SetOptFit(1); // 1- shows the fit parameters box   0 - hides it

  //gStyle->SetPadGridX(0);
  //gStyle->SetPadGridY(0);




  //gStyle->SetStatW(0.2); gStyle->SetStatH(0.2);
  
  TFile *rootfile = new TFile(infile,"read");
  
  TH1D *Htmp;
  
 

  //Instead of using this to make the plots, I could use this to extract the resolutions which I will then plug into the code Alexander gave me

  std::string BField[1] = {"Beast"};
  std::string Angle[6] = {"1.0", "1.5", "2.0", "2.5", "3.0", "3.5"};
  std::string Momentum[15] = {"1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0", "10.0", "12.0", "14.0", "16.0", "18.0", "20.0"};


  double Resolution[5][14];


for (int iA = 0; iA < 5; iA++)
{
  for (int iP = 0; iP < 14; iP++)
      {

          for (int iB = 0; iB <1; iB++) 
          {

            //Construct Histo Name
            //std:string HistName = "hMomentum_"+(std::string)Momentum[iP]+"GeV_"+(std::string)Angle[iA]+"Deg_"+(std::string)BField[iB]+"T";	
 	    std:string HistName = "h_ptRes_P_"+(std::string)Momentum[iP]+"_"+(std::string)Momentum[iP+1]+"_Eta_"+(std::string)Angle[iA]+"_"+(std::string)Angle[iA+1];	

            //Get Histogram
            //cout << "Getting: " << HistName << endl;
            Htmp = (TH1D*)rootfile->Get(HistName.c_str());

	    if (iA < 3) Htmp->GetXaxis()->SetRangeUser(-0.1, 0.1);
	    else if(iA ==3) Htmp->GetXaxis()->SetRangeUser(-0.2,0.2);
	    else
		{
		Htmp->GetXaxis()->SetRangeUser(-0.3,0.3);
		Htmp->Rebin(10);
		}

            //Extract Value
            //Htmp = (TH1D*)rootfile->Get(Form("Histograms/Location_%s/%s", Location[LOCATION], HistName));
            TF1 *gausFit = new TF1("gausFit", "gaus");
            gausFit->SetParameter(1, 0);
	    if (iA ==0) Htmp->Fit(gausFit, "");
            else Htmp->Fit(gausFit, "Q");
            double sigma = gausFit->GetParameter(2);
            if (sigma == 0) cout << "ERROR! - Bad Value for : " << HistName << endl;
		Resolution[iA][iP] = sigma*100; //make it a percent
		if (iA == 0) cout << "Resolution of tracks with transverse momentum " << Momentum[iP] << "-" << Momentum[iP+1] << " and Eta " << Angle[iA] << "-" << Angle[iA+1] << " : " << sigma << endl;
  		
		TFile outfile("./Output/Fits.root", "UPDATE");
		Htmp->Write();
		outfile.Close();
	        //Htmp->SetNameTitle( "h_PtFit_P_"+(std::string)Momentum[iP]+"_"+(std::string)Momentum[iP+1]+"_Eta_"+(std::string)Angle[iA]+"_"+(std::string)Angle[iA+1], "(reco_p - truth_p)/truth_p for tracks with pt between "+(std::string)Momentum[iP]+"-"+(std::string)Momentum[iP+1]+" and Eta "+(std::string)Angle[iA]+"-"+(std::string)Angle[iA+1] );


            //Print value to output file for plugging into other code
          }//End B Loop
    } // End Momentum Loop
    
} // End Angle Loop

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

    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 500);
    c1->UseCurrentStyle();
    c1->SetBorderMode(0);
    c1->SetFrameBorderMode(0);
    c1->SetFrameLineColor(kWhite);
    c1->SetFillColor(kWhite);
    c1->cd();


     //Old PWG Requirement curve (may need to update)
     TF1 *PWGReq3p0 = new TF1("PWGReq3p0", " TMath::Sqrt((0.02*x)^2+0.5^2)", 1 , 30);


     TF1 *ReqFit = new TF1("ReqFit", " TMath::Sqrt(([0]*x)^2+[1]^2)", 1 , 31);


     //get values from 2D array to use for plotting
     double res[14];
     double meanpt[14] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19};
     for (int i = 0; i < 14; i++)
     {
   	res[i] = Resolution[_ETA_][i];
     }     

     TGraph *gResPlot;
     gResPlot = new TGraph(14, meanpt , res  );
     gResPlot->SetTitle("Hybrid Detector With GEM Disks");
     gResPlot->SetMarkerSize(1.3);
     gResPlot->SetMarkerStyle(21);
     gResPlot->SetMarkerColor(kGreen);
 


 

    TH1D *hdum  = new TH1D("hdum",   "", 15, -0.01, 32);

    if (_ETA_ == 0) hdum->SetTitle("Pions at #eta ~0.5");

    hdum->GetXaxis()->SetLabelFont(52);
    hdum->GetYaxis()->SetLabelFont(52);
    hdum->GetXaxis()->SetTitleFont(52);
    hdum->GetYaxis()->SetTitleFont(52);

    hdum->GetXaxis()->SetTitle("p_{T} [GeV]");
    hdum->GetXaxis()->SetTitleSize(_TSIZE_);
    hdum->GetXaxis()->SetLabelSize(_LSIZE_);
    hdum->GetYaxis()->SetTitle("p_{T} resolution   #sigma_{p_{T}} /p_{T}  [%]");
    hdum->GetYaxis()->SetTitleSize(_TSIZE_);
    hdum->GetYaxis()->SetLabelSize(_LSIZE_);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(0.75);
    hdum->SetMinimum( 0.00);
    hdum->SetMaximum( 3.50);

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

     gResPlot->SetLineColor(kGreen);
     gResPlot->SetLineWidth(2);
     gResPlot->SetLineStyle(7);
     ReqFit->SetParameter(0,0.04);
     ReqFit->SetParName(0,"A");
     ReqFit->SetParameter(1,1.0);
     ReqFit->SetParName(1,"B");
     ReqFit->SetLineColor(kGreen+2);
     gResPlot->Fit(ReqFit, "R");
     gResPlot->Draw("PL"); //replace L with C for smooth line


     TLegend *legend = new TLegend(.15, .70, .45, .92);
     legend->AddEntry(gResPlot, "Hybrid Detector With GEM Disks", "P");

     legend->SetTextFont(52);
     legend->SetTextSize(0.04);
     legend->SetTextColor(kBlack);
     legend->SetFillColor(kWhite);
     legend->SetLineColor(kWhite);
     legend->Draw();



  return ; 
}

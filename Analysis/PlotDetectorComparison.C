#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//---------------------------------------------------------------------------------------------------------------------------------------  
//============================================== Change Inputs Here =====================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  

//---------------------------------------------------------------------------------------------------------------------------------------  
//=======================================================================================================================================
//---------------------------------------------------------------------------------------------------------------------------------------  


void PlotDetectorComparison()
{
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

    TF1 *PWGReq[5];
    PWGReq[0] = new TF1("PWGReq_-3.5_-2.5", "TMath::Sqrt((0.10*x)^2+0.5^2)", .1, 30);
    PWGReq[1] = new TF1("PWGReq_-2.5_1.0", "TMath::Sqrt((0.05*x)^2+0.5^2)", .1, 30);
    PWGReq[2] = new TF1("PWGReq_1.0_2.5"  , "TMath::Sqrt((0.05*x)^2+1.0^2)", .1, 30);
    PWGReq[3] = new TF1("PWGReq_2.5_3.5"  , "TMath::Sqrt((0.10*x)^2+2.0^2)", .1, 30);
    //PWGReq[2] = new TF1("PWGReq_-1.0_1.0" , "TMath::Sqrt((0.05*x)^2+0.5^2)", .1, 30);
    
    TF1 *PWGReq_DCA2D[6];
    PWGReq_DCA2D[0] = new TF1("PWGReq_DCA2D_-3.0_-2.5", "TMath::Sqrt((30/x)^2+40^2)", .1, 30);
    PWGReq_DCA2D[1] = new TF1("PWGReq_DCA2D_-2.5_-1", "TMath::Sqrt((30/x)^2+20^2)", .1, 30);
    PWGReq_DCA2D[2] = new TF1("PWGReq_DCA2D_-1_1", "TMath::Sqrt((20/x)^2+5^2)", .1, 30);
    PWGReq_DCA2D[3] = new TF1("PWGReq_DCA2D_1_2.5", "TMath::Sqrt((30/x)^2+20^2)", .1, 30);
    PWGReq_DCA2D[4] = new TF1("PWGReq_DCA2D_2.5_3", "TMath::Sqrt((30/x)^2+40^2)", .1, 30);
    PWGReq_DCA2D[5] = new TF1("PWGReq_DCA2D_3_3.5", "TMath::Sqrt((30/x)^2+60^2)", .1, 30);


TFile *infile = new TFile("./Output/TempPlots.root", "read");

const int NUM = 7;

TGraphErrors *graphs[NUM*2];
//double EtaBinEdges[NUM+1] = {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0};
double EtaBinEdges[NUM+1] = {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0};
//double EtaBinEdges[NUM+1] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
for (int i = 0; i < NUM; i++)
{
	//graphs[i] = (TGraphErrors*)infile->Get(Form("Plot_Complete_Eta_%0.1lf_%0.1lf", EtaBinEdges[i], EtaBinEdges[i+1] ));
	graphs[i] = (TGraphErrors*)infile->Get(Form("Plot_OrigHybrid_Eta_%0.1lf_%0.1lf", EtaBinEdges[i], EtaBinEdges[i+1] ));
	//graphs[i+NUM] = (TGraphErrors*)infile->Get(Form("Plot_OneRearGEMRing_Eta_%0.1lf_%0.1lf", EtaBinEdges[i], EtaBinEdges[i+1] ));
	//graphs[i+NUM] = (TGraphErrors*)infile->Get(Form("Plot_NoOuter_Eta_%0.1lf_%0.1lf", EtaBinEdges[i], EtaBinEdges[i+1] ));
	graphs[i+NUM] = (TGraphErrors*)infile->Get(Form("Plot_Projective_Eta_%0.1lf_%0.1lf", EtaBinEdges[i], EtaBinEdges[i+1] ));
}
        EColor ColorArray[10] = {kBlue, kMagenta, kRed, kGreen, kCyan, kOrange, kYellow, kPink, kViolet, kBlack};

	TCanvas *c1 = new TCanvas("c1", "Plots", 0, 0, 1200, 900);
	c1->cd();

    TH1D *hdum  = new TH1D("Plot", "", 15, 0, 31);
    hdum->GetXaxis()->SetLabelFont(52);
    hdum->GetYaxis()->SetLabelFont(52);
    hdum->GetXaxis()->SetTitleFont(52);
    hdum->GetYaxis()->SetTitleFont(52);
    //hdum->GetXaxis()->SetTitle("p_T [GeV]");
    hdum->GetXaxis()->SetTitle("p [GeV]");
    hdum->GetYaxis()->SetTitle("p resolution   #sigma_{p} /p  [%]");
    //hdum->GetYaxis()->SetTitle("#sigma_{DCA_{z}} [#mum]");
    //hdum->GetYaxis()->SetTitle("#sigma_{DCA_{2D}} [#mum]");
    hdum->GetXaxis()->SetTitleSize(0.06);
    hdum->GetXaxis()->SetLabelSize(0.05);
    hdum->GetYaxis()->SetTitleSize(0.06);
    hdum->GetYaxis()->SetLabelSize(0.05);
    hdum->GetXaxis()->SetTitleOffset(0.90);
    hdum->GetYaxis()->SetTitleOffset(1.2);
    hdum->SetMinimum( 0 );
    hdum->SetMaximum( 5 );
    hdum->GetXaxis()->SetNdivisions(408);
    hdum->GetYaxis()->SetNdivisions(808);
    gPad->SetGrid();
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFrameBorderMode(0);
    gPad->SetFrameLineColor(kWhite);
    gPad->SetLeftMargin(0.15);
    hdum->Draw();
    TLegend *legend = new TLegend(.15, .7, .45, .92);

	for (int j = 0; j<NUM;j++)
	{
	 	graphs[j]->SetMarkerSize(2.0);//1.3
	 	graphs[j]->SetMarkerStyle(21);
		graphs[j]->SetMarkerColor(ColorArray[j]);
		graphs[j]->SetLineColor(ColorArray[j]);
		graphs[j]->Draw("PLsame");
	        //legend->AddEntry(graphs[j], Form("%0.1lf < #eta < %0.1lf Baseline2", EtaBinEdges[j], EtaBinEdges[j+1] ), "P");
	        legend->AddEntry(graphs[j], Form("%0.1lf < #eta < %0.1lf Original Hybrid", EtaBinEdges[j], EtaBinEdges[j+1] ), "P");
	 	graphs[j+NUM]->SetMarkerSize(2.0);//1.3
	 	graphs[j+NUM]->SetMarkerStyle(25);
	 	graphs[j+NUM]->SetLineStyle(1);
		graphs[j+NUM]->SetMarkerColor(ColorArray[j]);
		graphs[j+NUM]->SetLineColor(ColorArray[j]);
		graphs[j+NUM]->Draw("PLsame");
	        //legend->AddEntry(graphs[j+NUM], Form("%0.1lf < #eta < %0.1lf One GEM Ring", EtaBinEdges[j], EtaBinEdges[j+1] ), "P");
	        //legend->AddEntry(graphs[j+NUM], Form("%0.1lf < #eta < %0.1lf No Outer GEM Disk", EtaBinEdges[j], EtaBinEdges[j+1] ), "P");
	        legend->AddEntry(graphs[j+NUM], Form("%0.1lf < #eta < %0.1lf Projective", EtaBinEdges[j], EtaBinEdges[j+1] ), "P");
		

	}
	int PWGReqChoice = 0;
	/*
	double EtaReqBinEdges[7] = {-3.5, -2.5, -1.0, 1.0,  2.5, 3.0, 3.5};
	PWGReq_DCA2D[PWGReqChoice]->SetLineColor(kBlack);
	PWGReq_DCA2D[PWGReqChoice]->Draw("same");
	legend->AddEntry(PWGReq_DCA2D[PWGReqChoice], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , EtaReqBinEdges[PWGReqChoice], EtaReqBinEdges[PWGReqChoice+1]));
	PWGReqChoice = 4;
	PWGReq_DCA2D[PWGReqChoice]->SetLineColor(kBlack);
	PWGReq_DCA2D[PWGReqChoice]->SetLineStyle(9);
	PWGReq_DCA2D[PWGReqChoice]->Draw("same");
	legend->AddEntry(PWGReq_DCA2D[PWGReqChoice], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , EtaReqBinEdges[PWGReqChoice], EtaReqBinEdges[PWGReqChoice+1]));
	PWGReqChoice = 5;
	PWGReq_DCA2D[PWGReqChoice]->SetLineColor(kBlack);
	PWGReq_DCA2D[PWGReqChoice]->SetLineStyle(2);
	PWGReq_DCA2D[PWGReqChoice]->Draw("same");
	legend->AddEntry(PWGReq_DCA2D[PWGReqChoice], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , EtaReqBinEdges[PWGReqChoice], EtaReqBinEdges[PWGReqChoice+1]));
	*/

	double EtaReqBinEdges[5] = {-3.5, -2.5, 0, 2.5, 3.5};
	PWGReq[PWGReqChoice]->SetLineColor(kBlack);
	PWGReq[PWGReqChoice]->Draw("same");
	legend->AddEntry(PWGReq[PWGReqChoice], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , EtaReqBinEdges[PWGReqChoice], EtaReqBinEdges[PWGReqChoice+1]));
	PWGReqChoice = 1;
	PWGReq[PWGReqChoice]->SetLineColor(kBlack);
	PWGReq[PWGReqChoice]->SetLineStyle(9);
	PWGReq[PWGReqChoice]->Draw("same");
	legend->AddEntry(PWGReq[PWGReqChoice], Form( "PWG Requirement for %0.1lf < #eta < %0.1lf " , EtaReqBinEdges[PWGReqChoice], EtaReqBinEdges[PWGReqChoice+1]));

legend->Draw();

}

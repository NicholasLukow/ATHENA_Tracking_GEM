#include <fstream>
#include <vector>

Bool_t DEBUG = false;
Double_t DegtoRad = TMath::Pi()/180.0;

Double_t RotateTheta(Double_t ang);

void Myplot_matscan(TString inFile = "./matscan", TString outFile = "matscan.root")
{
  //TString inFilePath = "/eic/u/mposik/Fun4all/macros/macros/Mymatscan/Geometry/";
  TString inFilePath = "./Comparisons/";
  inFile = "OneRing";
  //inFile = "NoOuter";
  //read in matscan results
  fstream fIN;
  fIN.open(Form("%s%s.dat",inFilePath.Data(),inFile.Data()),ios::in);
  if(!fIN) cout << "Can't open " << inFile << endl;

  vector<Double_t> Theta, Length, X0, Lambda, Phi, Len, Eta, ThetaPrime;

  Double_t ith,iphi,ilen,ix0,ilambda,iEta,ithprime;

  while(fIN >> ith >> iphi >> ilen >> ix0 >> ilambda)
    {
      if(ith!=90.0)
        {
          Theta.push_back(ith); Phi.push_back(iphi); 
          ithprime = RotateTheta(ith); ThetaPrime.push_back(ithprime);
          iEta = -1.0*TMath::Log(TMath::Tan(ithprime*DegtoRad/2.0)); Eta.push_back(iEta);
          Len.push_back(ilen); X0.push_back(ix0*100); Lambda.push_back(ilambda);      
          if(DEBUG) cout << "Theta " << ith << " ThetaPrime: " << ithprime << " Eta: " << iEta << " X0 : " << ix0 << endl;
        }
    }
  fIN.close();
  TFile *fRoot = new TFile(Form("%s%s",inFilePath.Data(),outFile.Data()),"UPDATE");
//make Plot
  //TGraph *gmat = new TGraph(Phi.size(),&Phi[0],&X0[0]);
  TGraph *gmat = new TGraph(Eta.size(),&Eta[0],&X0[0]);
  TCanvas *c1 = new TCanvas("c1","mat");
  c1->cd(1);
  gmat->Draw("al");
  //gmat->GetXaxis()->SetRangeUser(0,360);
  //gmat->GetYaxis()->SetRangeUser(0,0.25);
  gmat->GetXaxis()->SetRangeUser(-4,4);
  gmat->GetXaxis()->SetTitle("Eta");
  //gmat->GetXaxis()->SetTitle("Phi");
  gmat->GetYaxis()->SetRangeUser(0,1000);
  gmat->GetYaxis()->SetTitle("Material Budget, X/X_0 [%]");
  gmat->SetTitle("Material Scan of Projective Hybrid");
  //gmat->SetTitle("Material Scan of Al Support for Si Disks");
  //gmat->SetTitle("Material Scan in Phi at Fixed Eta = 3.0");
  //gmat->SetTitle("Material Scan of B-1.0 P-1.0 N-1.0");
  gmat->Write(inFile);



  const int PHIBINS = 1;
  const int ETABINS = 2000;


  //TCanvas *c2 = new TCanvas("c2", "2Dmat");
  //c2->cd(1);
 
  //TH2D *h2d_mat = new TH2D("h2d_mat", "2D Material Scan", PHIBINS, 0, 2*TMath::Pi(), 50, 1.3, 3.4);
  //TH2D *h2d_mat = new TH2D("h2d_mat", "2D Material Scan", PHIBINS, 0, 2*TMath::Pi(), 50, 0.8, 5.0);
  //TH2D *h2d_mat = new TH2D("h2d_mat", "2D Material Scan", PHIBINS, 0, 360, 50, 0.8, 5.0);



/*
  for(int iEta = 0; iEta < 200; iEta++)
  {
	for(int iPhi = 0; iPhi < PHIBINS; iPhi++)
	{
		if( iEta == 37) cout << "Eta: " << Eta[iEta*PHIBINS] << "  Theta: " << Theta[iEta*PHIBINS] << "   Phi: " << Phi[iPhi] << "   X0: " << X0[PHIBINS*iEta+iPhi] << endl;
		h2d_mat->Fill(Phi[iPhi],Eta[PHIBINS*iEta], X0[iEta*PHIBINS+iPhi]);
	}
  }
  //h2d_mat->Draw("colz pol");
  //h2d_mat->Draw("colz");
  h2d_mat->Draw("lego2 pol");
*/

 
}

Double_t RotateTheta(Double_t ang)
{
  Double_t thprime = 90.0 - ang;
  return thprime;
}

#define Analyze_cxx
#include <TVector3.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void RunAnalyze(const char* infile = "./InTest_highP.root")
{
//load the analysis script library
R__LOAD_LIBRARY(Analyze_C.so)

//get the file
TFile myFile(infile);

//get the tree
TTree *tree;
myFile.GetObject("tracks", tree);

//Run the Analysis
Analyze a(tree);
a.Loop();


//Close the input file
myFile.Close();
}

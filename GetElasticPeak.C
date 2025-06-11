// ******************************************************************************** //
//                     Created by J.McMurtry on June 10th, 2025                     //
//                                                                                  //
// Script to analyze elastic yields by applying a global cut and fitting the data   //
// ******************************************************************************** //

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TPaveText.h"

// **** ========== Useful functions ========== ****  
// returns today's date
std::string GetDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);
  std::string yyyy = to_string(1900 + ltm.tm_year);
  std::string mm = to_string(1 + ltm.tm_mon);
  std::string dd = to_string(ltm.tm_mday);
  std::string date = mm + '/' + dd + '/' + yyyy;
  return date;
}

// Wrap long lines of text whenusing pt->AddText
void AddWrappedText(TPaveText* pt, const TString& text, int maxCharsPerLine = 60) {
  int pos = 0;
  while (pos < text.Length()) {
    TString line = text(pos, maxCharsPerLine);
    int lastSpace = text.Last(' ');
    if (lastSpace != kNPOS && pos + maxCharsPerLine < text.Length())
      line = text(pos, lastSpace), pos += lastSpace + 1;
    else
      pos += maxCharsPerLine;
    pt->AddText(line);
  }
}

// Crudely fitting a gaussian and a polynomial to our data
void FitGausQuad( TH1D *htest, double thresh=0.5 ){
  int binmax = htest->GetNbinsX();
  int binlow = 1, binhigh = binmax;

  //double max = htest->GetBinContent(binmax);

  //while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  //while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh);

  TF1 *fitfunc = new TF1("fitfunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2", -5, 5);
  fitfunc->SetParameters(10, 0, 0.1, 0, 0, 0);
  fitfunc->SetParNames("Amp", "Mean", "Sigma", "Offset", "Slope", "Square");
  
  htest->Fit(fitfunc,"q0S","",xlow, xhigh);
  cout << "xlow = " << xlow << ", xhigh = " << xhigh << endl;
  cout << "binlow = " << binlow << ", binhigh = " << binhigh << endl;
}


// **** ========== Main functions ========== **** 
void GetElasticPeak( const char *configfilename, const char *outfilename="ElasticPeak.root"){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  // TStopwatch *sw2 = new TStopwatch();
  sw->Start(); //sw2->Start();  
  
  TChain *C = new TChain("T");
  cout << "Ok, let's get started!" << endl;

  // Reading in the config file
  ifstream infile(configfilename);
  char runlistfile[1000];
  TString currentline, readline;
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endRunlist") ){
    if( !currentline.BeginsWith("#") ){
      sprintf(runlistfile,"%s",currentline.Data());
      ifstream run_list(runlistfile);
      while( readline.ReadLine( run_list ) && !readline.BeginsWith("endlist") ){
	if( !readline.BeginsWith("#") ){
	  cout << readline << "\n";
          C->Add(readline);
	}
      }
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }
  
  cout << "Let's ball" << endl;
  
  // Time to initialize branches we will use:
  const int MAXHEEP = 100;

  // List of variables to cut on
  double dxECAL[MAXHEEP], dyECAL[MAXHEEP];
  double eprime_eth[MAXHEEP], ecalo[MAXHEEP];

  // Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  // the * applies it to all branches, the 0 disables those branches. to enable would need to make 1
  C->SetBranchStatus("*",0);

  // Initializing branches
  C->SetBranchStatus("heep.*",1);
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("g.runnum",1);
  C->SetBranchStatus("earm.ecal.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("sbs.tr.*",1);

  // Filling Arrays
  C->SetBranchAddress("heep.dxECAL", dxECAL);
  C->SetBranchAddress("heep.dyECAL", dyECAL);
  C->SetBranchAddress("heep.eprime_eth", eprime_eth);
  C->SetBranchAddress("heep.ecalo", ecalo);

  TFile *fout = new TFile(outfilename, "RECREATE");

  // Set up the histograms we want to look at:
  TH1D *hdxECAL = new TH1D("hdxECAL", "heep.dxECAL after global cut; heep.dxECAL (m);", 100, -0.06, 0.06);
  TH1D *hdyECAL = new TH1D("hdyECAL", "heep.dyECAL after global cut; heep.dyECAL (m);", 100, -0.06, 0.06);

  TH2D *hdxECAL_v_dyECAL = new TH2D("hdxECAL_v_dyECAL", "heep.dxECAL vs heep.dyECAL ; heep.dyECAL (m); heep.dxECAL (m)", 100, -0.06, 0.06, 100, -0.06, 0.06);

  TH1D *hEdivP = new TH1D("hEdivP", "E/P after global cut; E/P;", 100, 0.0, 1.3);

  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );

  int treenum=0, currenttreenum=0;
  int ngoodevs = 0;
  while( C->GetEntry( nevent++ ) && nevent){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;

      if (GlobalCut) {
        delete GlobalCut;
        GlobalCut = nullptr;
      }
      GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
      GlobalCut->UpdateFormulaLeaves();
    }
   
    if( nevent % 100000 == 0 ) cout << nevent << endl;
    bool passedcut = GlobalCut->EvalInstance(0) != 0;

    if( passedcut ){
      ngoodevs ++;
      hdxECAL->Fill( dxECAL[0]);
      hdyECAL->Fill( dyECAL[0]);
      hdxECAL_v_dyECAL->Fill( dyECAL[0], dxECAL[0]);
      hEdivP->Fill(ecalo[0]/eprime_eth[0]);
    }
    
  }

  TString outfilepdf = outfilename; // Make a pdf file to save these histograms to
  outfilepdf.ReplaceAll(".root",".pdf");

  FitGausQuad( hdxECAL, 0.5);
  TF1 *fitfuncX = (TF1*) (hdxECAL->GetListOfFunctions()->FindObject("fitfunc"));
  cout << "Fit mean X = " << fitfuncX->GetParameter("Mean") << endl;

  FitGausQuad( hdyECAL, 0.5);
  TF1 *fitfuncY = (TF1*) (hdyECAL->GetListOfFunctions()->FindObject("fitfunc"));
  cout << "Fit mean Y = " << fitfuncY->GetParameter("Mean") << endl;
  
  // Make some plots for us to look at
  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(2,2);
  c1->cd(1);    hdxECAL->Draw();    fitfuncX->Draw("same");
  c1->cd(2);    hdyECAL->Draw();    fitfuncY->Draw("same");
  c1->cd(3);    hdxECAL_v_dyECAL->Draw();
  c1->cd(4);    hEdivP->Draw();
  //c1->Update();
  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page


  /**** Summary Canvas ****/
  TCanvas *cSummary = new TCanvas("cSummary","Summary");
  cSummary->cd();
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form(" Date of creation: %s",GetDate().c_str()));
  sw->Stop();
  TText *t1 = pt->AddText(Form("Macro processing time: CPU %.1fs | Real %.1fs",sw->CpuTime(),sw->RealTime()));
  t1->SetTextColor(kBlue);  
  pt->AddText(Form(" Total # events analyzed: %lld ",nevent));
  pt->AddText(Form("Total # events passed global cuts: %lld", ngoodevs));
  TText *t2 = pt->AddText(" Global cuts: ");
  t2->SetTextColor(kRed);
  AddWrappedText(pt, globalcut.GetTitle());
  pt->Draw();
  cSummary->Print(outfilepdf + ")"); // Close out the pdf
  
  
  //  cout << "Total number of events: " << nevent << endl;
  //cout << "Number of events passing elastic cuts: " << ngoodevs << endl;
  //cout << "Our global cut: " << globalcut << endl;
  fout->Write();
  cout << "That's all folks!" << endl;
}

// ******************************************************************************** //
//                     Created by J.McMurtry on August 4th, 2025                    //
//                                                                                  //
// Script to look at ECal/HCal timing and number of events for different thresholds //
// ******************************************************************************** //

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include <map>
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

// **** ========== Main functions ========== **** 
void ECalScan( const char *configfilename, const char *chargefile="runCharge.txt", const char *outfilename="output/ECalADCtime.root"){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();  

  TChain *C = new TChain("T");

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

  TCut globalcut = "sbs.hcal.nclus>0&&earm.ecal.nclus>0";

  cout << "Let's ball" << endl;

   // Time to initialize branches we will use:
  const int MAXBLK = 100;

  // List of variables to cut on
  double atimeHCal[MAXBLK], atimeECal[MAXBLK];
  double runnum[1];

  C->SetBranchStatus("*",0);

  // Initializing branches
  C->SetBranchStatus("earm.ecal.atimeblk",1);
  C->SetBranchStatus("earm.ecal.nclus",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  C->SetBranchStatus("sbs.hcal.nclus",1);
  C->SetBranchStatus("g.runnum",1);
  
  // Filling Arrays
  C->SetBranchAddress("earm.ecal.atimeblk", atimeECal);
  C->SetBranchAddress("sbs.hcal.atimeblk", atimeHCal);
  C->SetBranchAddress("g.runnum", runnum);

  TFile *fout = new TFile(outfilename, "RECREATE");

  TH1D *hclus_atimediff = new TH1D("hclus_atimediff", "ECAL ADC Clus time - HCAL ADC Clus time (ns);", 300, -150, 150);
  
  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );

  int treenum=0, currenttreenum=0;
  int ngoodevs = 0;
  map<int, int> runmap;
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
      if( runmap.find(runnum[0]) == runmap.end() ){
        runmap.insert({runnum[0], 1});
      }
      else{
	runmap[runnum[0]]++;
      }
      ngoodevs ++;
      hclus_atimediff->Fill( atimeECal[0] - atimeHCal[0]);
    }    
  }

  // Let's do some charge corrections
  ifstream fcharge(chargefile);  // read a text file which contains all the total charge values of each run
  int runrow;
  double chargerow;
  map<int, double> chargemap;
  while(fcharge >> runrow >> chargerow){
    //cout << "Read charge: " << runrow << ", " << chargerow << endl;
    if( chargemap.find(runrow) == chargemap.end() ){
        chargemap.insert({runrow, chargerow});
    }
  }
  fcharge.close();

  // find the total charge from the runs we input
  double totalcharge = 0;
  for( const auto& [run, val] : runmap) {
    cout << "Charge = " << chargemap.at(run) << endl;
    totalcharge += chargemap.at(run);
  }
  cout << "total charge = " << totalcharge << endl;

  // Let's try to estimate the background rate around the peak
  double delt = 5; // 20 ns window
  double t0 = 8; // start of time window
  int lbin = hclus_atimediff->FindBin(t0);
  int rbin = hclus_atimediff->FindBin(t0+delt);

  double nbkg = hclus_atimediff->Integral(lbin, rbin);
  cout << "Estimated background rate = " << nbkg/delt << " counts/ns" << endl;

  // Now let's try to estimate the actual number of signal events
  double tplus = 15; // 20 ns window
  double tmin = t0-tplus; // start of time window
  int binmin = hclus_atimediff->FindBin(tmin);
  //int binmax = hclus_atimediff->FindBin(tmin+ttot);
  double ntot = hclus_atimediff->Integral(binmin, lbin);
  double nsig = ntot - nbkg * tplus / delt;
  cout << "Estimated number signal event count = " << nsig << endl;
  cout << "Adjusted per coulomb = " << nsig/totalcharge << endl;

  TString outfilepdf = outfilename; // Make a pdf file to save these histograms to
  outfilepdf.ReplaceAll(".root",".pdf");
  
  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(1,1);
  c1->cd(1);    hclus_atimediff->Draw();

  double ymax = hclus_atimediff->GetMaximum();
  TLine *lline = new TLine(t0, 0, t0, ymax);
  TLine *rline = new TLine(t0+delt, 0, t0+delt, ymax);
  TLine *sline = new TLine(tmin, 0, tmin, ymax);
  lline->SetLineColor(kRed);
  rline->SetLineColor(kRed);
  sline->SetLineColor(kGreen);
  lline->Draw("same");     rline->Draw("same");     sline->Draw("same");
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
  //pt->AddText(Form("Total # events passed global cuts: %lld", ngoodevs));
  //TText *t2 = pt->AddText(" Global cuts: ");
  //t2->SetTextColor(kRed);
  //AddWrappedText(pt, globalcut.GetTitle());
  TText *t3 = pt->AddText(Form(" Average Charge Normalized Number of Events: %.0f",  nsig/totalcharge ));
  t3->SetTextColor(kRed);
  //pt->AddText(Form(" Signal fraction: %.5f ",fsig));
  //TText *t4 = pt->AddText(Form(" Fit Integrated Elastic Yield Y: %.0f", ElasticYieldY));
  //t4->SetTextColor(kRed);
  pt->Draw();
  cSummary->Print(outfilepdf + ")"); // Close out the pdf


  fout->Write();
  cout << "That's all folks!" << endl;
}

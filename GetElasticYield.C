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

// Exclude specified range in fit
Double_t RejectFunc(Double_t *x, Double_t *par){
  if (x[0] > -0.08 && x[0] < 0.09){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Fit our background to the sides away from central peak
TF1 *FitBkgrSide( TH1D *htest){
  TF1 *Bkgr = new TF1("Bkgr", RejectFunc, -5, 5, 3);
  Bkgr->SetParameters(10, 0, -0.1);
  Bkgr->SetNpx(1000);
  Bkgr->SetLineColor(kBlue);
  htest->Fit(Bkgr, "R");
  return Bkgr;
}

// Fit our signal over the background
TF1 *FitGausQuad( TH1D *htest, vector<double> value ){
  int binmax = htest->GetNbinsX();
  int binlow = 1, binhigh = binmax;

  //double max = htest->GetBinContent(binmax);

  //while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  //while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh);

  TF1 *fitfunc = new TF1("fitfunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2", -5, 5);
  fitfunc->SetParameters(10, 0, 0.1, value[0], value[1], value[2]);
  fitfunc->SetParNames("Amp", "Mean", "Sigma", "Offset", "Slope", "Quad");
  fitfunc->FixParameter(3, value[0]);
  fitfunc->FixParameter(4, value[1]);
  fitfunc->FixParameter(5, value[2]);
  fitfunc->SetNpx(1000);
  fitfunc->SetLineColor(kRed);
  
  htest->Fit(fitfunc,"q0S","",xlow, xhigh);
  cout << "xlow = " << xlow << ", xhigh = " << xhigh << endl;
  cout << "binlow = " << binlow << ", binhigh = " << binhigh << endl;
  return fitfunc;
}

// Seperating out signal
TF1* GausOnly( vector<double> value ){
  TF1 *fitgaus = new TF1("fitgaus", "[0]*exp(-0.5*((x-[1])/[2])^2)", -1, 1);
  fitgaus->SetParameters(value[0], value[1], value[2]);
  fitgaus->SetNpx(1000);
  fitgaus->SetLineColor(kBlue);
  //fitgaus->Draw("same");
  return fitgaus;
}

// Sperating out background
TF1* QuadOnly( vector<double> value ){
  TF1 *fitquad = new TF1("fitquad", "[0]+[1]*x+[2]*x^2", -1, 1);
  fitquad->SetParameters(value[0], value[1], value[2]);
  fitquad->SetNpx(1000);
  fitquad->SetLineColor(kMagenta);
  // fitquad->Draw("same");
  return fitquad;
}

//map<int, double> mapbyrun(const vector<pair<int, double>>& entries) {
//  map<int, double> totals;
//}

int Readfilename( const string& filename){
  string pattern = "gep5_fullreplay_";
  size_t start = filename.find(pattern);
  if( start == string::npos) return -1;

  start += pattern.length();
  size_t end = filename.find("_", start);
  string runStr = filename.substr(start, end - start);
  return stoi(runStr);
}


// **** ========== Main functions ========== **** 
void GetElasticYield( const char *configfilename, const char *chargefile="runCharge.txt", const char *outfilename="ElasticPeak.root"){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  // TStopwatch *sw2 = new TStopwatch();
  sw->Start(); //sw2->Start();  
  
  TChain *C = new TChain("T");
  TChain *R = new TChain("TSsbsgep");
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
	  R->Add(readline);
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
  double L1A[MAXHEEP], TS10[MAXHEEP];
  double runnum[1];

  // Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  // the * applies it to all branches, the 0 disables those branches. to enable would need to make 1
  C->SetBranchStatus("*",0);
  R->SetBranchStatus("*",0);

  // Initializing branches
  C->SetBranchStatus("heep.*",1);
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("g.runnum",1);
  C->SetBranchStatus("earm.ecal.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  R->SetBranchStatus("sbsgep.L1A.scalerRate",1);
  R->SetBranchStatus("sbsgep.TS10_EcalScint.scalerRate",1);

  // Filling Arrays
  C->SetBranchAddress("heep.dxECAL", dxECAL);
  C->SetBranchAddress("heep.dyECAL", dyECAL);
  C->SetBranchAddress("heep.eprime_eth", eprime_eth);
  C->SetBranchAddress("heep.ecalo", ecalo);
  C->SetBranchAddress("g.runnum", runnum);
  R->SetBranchAddress("sbsgep.L1A.scalerRate", L1A);
  R->SetBranchAddress("sbsgep.TS10_EcalScint.scalerRate", TS10);
  

  TFile *fout = new TFile(outfilename, "RECREATE");

  // Set up the histograms we want to look at:
  int bins = 200;   // useful for normalizing later on
  double range = 0.36;
  TH1D *hdxECAL = new TH1D("hdxECAL", "heep.dxECAL after global cut; heep.dxECAL (m);", bins, -range/2, range/2);
  TH1D *hdyECAL = new TH1D("hdyECAL", "heep.dyECAL after global cut; heep.dyECAL (m);", bins, -range/2, range/2);

  TH2D *hdxECAL_v_dyECAL = new TH2D("hdxECAL_v_dyECAL", "heep.dxECAL vs heep.dyECAL ; heep.dyECAL (m); heep.dxECAL (m)", bins/2, -range/2, range/2, bins/2, -range/2, range/2);

  TH1D *hEdivP = new TH1D("hEdivP", "E/P after global cut; E/P;", 100, 0.0, 1.3);

  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );

  int treenum=0, currenttreenum=0;
  int ngoodevs = 0;
  map<int, int> runyield;
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
      if( runyield.find(runnum[0]) == runyield.end() ){
        runyield.insert({runnum[0], 1});
      }
      else{
	runyield[runnum[0]]++;
      }
      ngoodevs ++;
      hdxECAL->Fill( dxECAL[0]);
      hdyECAL->Fill( dyECAL[0]);
      hdxECAL_v_dyECAL->Fill( dyECAL[0], dxECAL[0]);
      hEdivP->Fill(ecalo[0]/eprime_eth[0]);
    }
    
  }
  cout << "Print my map:" << endl;
  for (const auto& element : runyield){
    cout << element.first << " : " << element.second << "\n";
  }
  //cout << runyield << endl;

  // Let's do some charge corrections
  ifstream fcharge(chargefile);
  int runrow;
  double chargerow;
  while(fcharge >> runrow >> chargerow){
    cout << "Read charge: " << runrow << ", " << chargerow << endl;
  }
  fcharge.close();

  // Let's do some livetime corrections
  long nts=0;
  double L1Asum=0, liveT=0;
  int treenumR=0, currenttreenumR=0;
  map<int, double> runlivet;
  map<int, int> runlivetN;
  int currentRunNumber = -1;
  while( R->GetEntry( nts++ ) && nts){
    currenttreenumR = R->GetTreeNumber();
    if( nts == 1 || currenttreenumR != treenumR ){
      treenumR = currenttreenumR;
    }

    TFile* currentfile = R->GetCurrentFile();
    string fname = currentfile->GetName();

    currentRunNumber = Readfilename(fname);

    L1Asum += L1A[0];
    liveT += L1A[0] / (TS10[0]+4);

    if( runlivet.find(currentRunNumber) == runlivet.end() ){
        runlivet.insert({currentRunNumber, liveT});
	runlivetN.insert({currentRunNumber, 1});
      }
      else{
	runlivet[currentRunNumber] += liveT;
	runlivetN[currentRunNumber]++;
      }
    //cout << "filename : " << currentRunNumber << endl;
    
  }
  double L1Aavg = L1Asum/nts;
  double liveTavg = liveT/nts;
  cout << "L1A average = " << L1Aavg << endl;
  cout << "Live time average = " << liveTavg << endl;
  
  TString outfilepdf = outfilename; // Make a pdf file to save these histograms to
  outfilepdf.ReplaceAll(".root",".pdf");

  // Let's fit curves to our histograms
  vector<double> bparX, bparY;
  vector<double> sparX, sparY;

  TF1 *fitBkgrSideX = FitBkgrSide(hdxECAL);  // Fit the background at the sides
  bparX.push_back(fitBkgrSideX->GetParameter(0));
  bparX.push_back(fitBkgrSideX->GetParameter(1));
  bparX.push_back(fitBkgrSideX->GetParameter(2));
  TF1 *fitAllX = FitGausQuad(hdxECAL, bparX); // Fit signal over the background
  sparX.push_back(fitAllX->GetParameter(0));
  sparX.push_back(fitAllX->GetParameter(1));
  sparX.push_back(fitAllX->GetParameter(2));
  TF1 *fitSigX = GausOnly(sparX); // Extract only signal and backgroun
  TF1 *fitBkgrX = QuadOnly(bparX);

  
  TF1 *fitBkgrSideY = FitBkgrSide(hdyECAL);
  bparY.push_back(fitBkgrSideY->GetParameter(0));
  bparY.push_back(fitBkgrSideY->GetParameter(1));
  bparY.push_back(fitBkgrSideY->GetParameter(2));
  TF1 *fitAllY = FitGausQuad(hdyECAL, bparY);
  sparY.push_back(fitAllY->GetParameter(0));
  sparY.push_back(fitAllY->GetParameter(1));
  sparY.push_back(fitAllY->GetParameter(2));
  TF1 *fitSigY = GausOnly(sparY);
  TF1 *fitBkgrY = QuadOnly(bparY);
  
  
  // Make some plots for us to look at
  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(2,2);
  c1->cd(1);    hdxECAL->Draw();    fitAllX->Draw("same");    fitBkgrX->Draw("same");    fitSigX->Draw("same");
  c1->cd(2);    hdyECAL->Draw();    fitAllY->Draw("same");    fitBkgrY->Draw("same");    fitSigY->Draw("same");
  c1->cd(3);    hdxECAL_v_dyECAL->Draw();
  c1->cd(4);    hEdivP->Draw();
  //c1->Update();
  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page

  double ElasticYieldX = (fitSigX->Integral(-1,1))/(range/bins);
  double ElasticYieldY = (fitSigY->Integral(-1,1))/(range/bins);

  cout << "Elastic yeild X = " << ElasticYieldX << endl;
  cout << "Elastic yeild Y = " << ElasticYieldY << endl;


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
  TText *t3 = pt->AddText(Form(" Fit Integrated Elastic Yield X: %.0f", ElasticYieldX));
  t3->SetTextColor(kRed);
  TText *t4 = pt->AddText(Form(" Fit Integrated Elastic Yield Y: %.0f", ElasticYieldY));
  t4->SetTextColor(kRed);
  pt->Draw();
  cSummary->Print(outfilepdf + ")"); // Close out the pdf
  
  
  //  cout << "Total number of events: " << nevent << endl;
  //cout << "Number of events passing elastic cuts: " << ngoodevs << endl;
  //cout << "Our global cut: " << globalcut << endl;
  fout->Write();
  cout << "That's all folks!" << endl;
}

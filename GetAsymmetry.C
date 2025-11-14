// ------------------------------------------------------------------------ //
// Script to extract azimuthal distribution and helicity asymmetry for GEp  //
//                                                                          //
// ---------                                                                //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 08-12-2025                  //
// ---------                                                                //
// ** Do not tamper with this sticker! Log any updates to the script above. //
// ------------------------------------------------------------------------ //

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include <map>
#include <cmath>
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
  if (x[0] > -0.04 && x[0] < 0.06){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Fit our background to the sides away from central peak
TF1 *FitBkgrSide( TH1D *htest){
  int firstbin = htest->FindFirstBinAbove(0);
  int lastbin = htest->FindLastBinAbove(0);
  double xlow = htest->GetBinLowEdge(firstbin);
  double xhigh = htest->GetBinLowEdge(lastbin);
  
  TF1 *Bkgr = new TF1("Bkgr", RejectFunc, xlow, xhigh, 3);
  Bkgr->SetParameters(0, 0, -0.1);
  Bkgr->SetNpx(1000);
  Bkgr->SetLineColor(kBlue);
  htest->Fit(Bkgr, "R");
  return Bkgr;
}

// Fit our signal over the background
TF1 *FitGausQuad( TH1D *htest, vector<double> value ){
  int firstbin = htest->FindFirstBinAbove(0);
  int lastbin = htest->FindLastBinAbove(0);

  //int binmax = htest->GetNbinsX();
  //int binlow = 1, binhigh = binmax;
  //double max = htest->GetBinContent(binmax);
  //while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  //while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(firstbin);
  double xhigh = htest->GetBinLowEdge(lastbin);

  TF1 *fitfunc = new TF1("fitfunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2", xlow, xhigh);
  fitfunc->SetParameters(10, 0, 0.1, value[0], value[1], value[2]);
  fitfunc->SetParNames("Amp", "Mean", "Sigma", "Offset", "Slope", "Quad");
  fitfunc->FixParameter(3, value[0]);
  fitfunc->FixParameter(4, value[1]);
  fitfunc->FixParameter(5, value[2]);
  fitfunc->SetNpx(2000);
  fitfunc->SetLineColor(kRed);
  
  htest->Fit(fitfunc,"q0S","",xlow, xhigh);
  //cout << "xlow = " << xlow << ", xhigh = " << xhigh << endl;
  //cout << "binlow = " << binlow << ", binhigh = " << binhigh << endl;
  return fitfunc;
}

// Seperating out signal
TF1* GausOnly( vector<double> value , TH1D *htest){
  int firstbin = htest->FindFirstBinAbove(0);
  int lastbin = htest->FindLastBinAbove(0);
  double xlow = htest->GetBinLowEdge(firstbin);
  double xhigh = htest->GetBinLowEdge(lastbin);
  
  TF1 *fitgaus = new TF1("fitgaus", "[0]*exp(-0.5*((x-[1])/[2])^2)", xlow, xhigh);
  fitgaus->SetParameters(value[0], value[1], value[2]);
  fitgaus->SetNpx(1000);
  fitgaus->SetLineColor(kBlue);
  //fitgaus->Draw("same");
  return fitgaus;
}

// Sperating out background
TF1* QuadOnly( vector<double> value, TH1D *htest ){
  int firstbin = htest->FindFirstBinAbove(0);
  int lastbin = htest->FindLastBinAbove(0);
  double xlow = htest->GetBinLowEdge(firstbin);
  double xhigh = htest->GetBinLowEdge(lastbin);
  
  TF1 *fitquad = new TF1("fitquad", "[0]+[1]*x+[2]*x^2", xlow, xhigh);
  fitquad->SetParameters(value[0], value[1], value[2]);
  fitquad->SetNpx(1000);
  fitquad->SetLineColor(kMagenta);
  // fitquad->Draw("same");
  return fitquad;
}

// Read run number from file name
int Readfilename( const string& filename){
  string pattern = "gep5_fullreplay_";
  size_t start = filename.find(pattern);
  if( start == string::npos) return -1;

  start += pattern.length();
  size_t end = filename.find("_", start);
  string runStr = filename.substr(start, end - start);
  return stoi(runStr);
}

// Structure for storing our info that is specific to each run
struct RunData{
  int events;
  double livetime;
  double charge;
};

// Function to draw asymmetry with lines at each bin
void DrawAsymmetry(TH1* h, double frac = 0.4, int color = kRed, int linewidth = 2, bool drawAxis = true, bool forceYminZero = false) {
    if (!h) return;

    // Force drawing into the current pad
    if (gPad) gPad->cd();

    h->SetLineColor(color);
    h->SetLineWidth(linewidth);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(color);

    if (drawAxis) h->Draw("P"); // only draw axes the first time
    if (forceYminZero) h->SetMinimum(0); //Set axis min to zero if desired

    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; i++) {
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);
        double binWidth = h->GetBinWidth(i);
        double half = 0.5 * binWidth * frac;
        TLine *line = new TLine(x - half, y, x + half, y);
        line->SetLineColor(color);
        line->SetLineWidth(linewidth);
        line->Draw("same");
    }

    h->Draw("P same");
}


// **** ========== Main functions ========== **** 
void GetAsymmetry( const char *configfilename="configElastic.cfg", const char *chargefile="runCharge.txt", const char *outfilename="output/Asymmetry.root"){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0); // No stats box on histograms

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
	  //R->Add(readline);
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

  TCut selectioncut = "";
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endselection") ){
    if( !currentline.BeginsWith("#") ){
      selectioncut += currentline.Data();
    }
  }

  //while ( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ) {
  //  currentline = currentline.Strip(TString::kBoth); // remove whitespace
  //  if (currentline.Length() == 0 || currentline.BeginsWith("#")) continue; // skip blank or commented lines


  cout << "Global cut = " << globalcut.GetTitle() << endl;
  cout << "Selection cut = " << selectioncut.GetTitle() << endl;
  cout << "Let's ball" << endl;
  
  // Time to initialize branches we will use:
  const int MAXHEEP = 100;

  // List of variables to cut on
  double dxECAL[MAXHEEP], dyECAL[MAXHEEP];
  double eprime_eth[MAXHEEP], ecalo[MAXHEEP];
  double runnum[1];
  double dpp[MAXHEEP];
  double sclose[MAXHEEP], zclose[MAXHEEP];
  double ttheta[MAXHEEP], tphi[MAXHEEP];
  double true_hel[1];

  // Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  // the * applies it to all branches, the 0 disables those branches. to enable would need to make 1
  C->SetBranchStatus("*",0);

  // Initializing branches
  C->SetBranchStatus("heep.*",1);
  C->SetBranchStatus("scalhel.true_hel",1);
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("g.runnum",1);
  C->SetBranchStatus("earm.ecal.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("sbs.hcal.nblk*",1);
  C->SetBranchStatus("sbs.hcal.y*",1);
  C->SetBranchStatus("sbs.hcal.x*",1);

  // Filling Arrays
  C->SetBranchAddress("heep.dxECAL", dxECAL);
  C->SetBranchAddress("heep.dyECAL", dyECAL);
  C->SetBranchAddress("heep.eprime_eth", eprime_eth);
  C->SetBranchAddress("heep.ecalo", ecalo);
  C->SetBranchAddress("heep.dpp", dpp);
  C->SetBranchAddress("g.runnum", runnum);
  C->SetBranchAddress("sbs.gemFPP.track.sclose", sclose);
  C->SetBranchAddress("sbs.gemFPP.track.zclose", zclose);
  C->SetBranchAddress("sbs.gemFPP.track.theta", ttheta);
  C->SetBranchAddress("sbs.gemFPP.track.phi", tphi);
  C->SetBranchAddress("scalhel.true_hel", true_hel);
  

  TFile *fout = new TFile(outfilename, "RECREATE");

  // Set up the histograms we want to look at:
  int bins = 200;   // useful for normalizing later on
  double range = 0.25;
  TH1D *hdxECAL = new TH1D("hdxECAL", "heep.dxECAL after global cut; heep.dxECAL (m);", bins, -range/2, range/2);
  TH1D *hdyECAL = new TH1D("hdyECAL", "heep.dyECAL after global cut; heep.dyECAL (m);", bins, -range/2, range/2);

  TH2D *hdxECAL_v_dyECAL = new TH2D("hdxECAL_v_dyECAL", "heep.dxECAL vs heep.dyECAL ; heep.dyECAL (m); heep.dxECAL (m)", bins/2, -range/2, range/2, bins/2, -range/2, range/2);

  TH1D *hEdivP = new TH1D("hEdivP", "E/P after global cut; E/P;", bins, 0.0, 1.3);
  TH1D *hdpp = new TH1D("hdpp", "dpp after global cut; dpp;", bins, -range/2, range/2);

  TH1D *hsclose = new TH1D("hsclose", "Distance of Closest Approach; DOCA (m)", bins, 0.0, 0.015);
  TH1D *hzclose = new TH1D("hzclose", "z of Closest Approach; z_{close} (m)", bins, 0.0, 3.0);
  TH1D *hzclose0 = new TH1D("hzclose0", "z of Closest Approach; z_{close} (m)", bins, 0.0, 3.0);
  TH1D *hzclose1 = new TH1D("hzclose1", "z of Closest Approach; z_{close} (m)", bins, 0.0, 3.0);
  TH1D *httheta = new TH1D("httheta", "Polar Scattering Angle;#theta_{FPP} (degrees)", bins, 0.0, 30.0);
  TH1D *htphi = new TH1D("htphi", "Azimuthal Scattering Angle; #phi (rad); Events", 18, -M_PI, M_PI);
  TH1D *htphipos = new TH1D("htphipos", "Azimuthal Scattering Angle; #phi (rad); Events", 18, -M_PI, M_PI);
  TH1D *htphineg = new TH1D("htphineg", "Azimuthal Scattering Angle; #phi (rad); Events", 18, -M_PI, M_PI);
  TH2D *hzvtheta = new TH2D("hzvtheta", "z vs #theta ;#theta (deg);z_{close} (m)", bins, 0.0, 30.0, bins, 0.0, 3.0);

  htphi->Sumw2();  htphipos->Sumw2();  htphineg->Sumw2();
  
  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
  TTreeFormula *SelectionCut = new TTreeFormula( "SelectionCut", selectioncut.GetTitle(), C );

  int treenum=0, currenttreenum=0;
  int ngoodevs = 0;
  int ngoodasym = 0;
  //map<int, int> runyield;
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

      if (SelectionCut) {
        delete SelectionCut;
        SelectionCut = nullptr;
      }
      SelectionCut = new TTreeFormula( "SelectionCut", selectioncut.GetTitle(), C );
      SelectionCut->UpdateFormulaLeaves();
    }
   
    if( nevent % 100000 == 0 ) cout << nevent << endl;
    bool passedcut = GlobalCut->EvalInstance(0) != 0;
    bool passedselection = SelectionCut->EvalInstance(0) != 0;

    if( passedcut ){
      //if( runyield.find(runnum[0]) == runyield.end() ){
      //  runyield.insert({runnum[0], 1});
      //}
      //else{
      //runyield[runnum[0]]++;
      //}
      ngoodevs ++;
      hdxECAL->Fill( dxECAL[0]);
      hdyECAL->Fill( dyECAL[0]);
      hdxECAL_v_dyECAL->Fill( dyECAL[0], dxECAL[0]);
      hEdivP->Fill(ecalo[0]/eprime_eth[0]);
      hdpp->Fill( dpp[0]);
      hsclose->Fill( sclose[0]);
      hzclose->Fill( zclose[0]);
      httheta->Fill( ttheta[0] * 180 / M_PI );
      hzvtheta->Fill( ttheta[0] * 180 / M_PI, zclose[0]);
      if (ttheta[0] * 180/ M_PI < 1) {
	      hzclose1->Fill(zclose[0]);
      } else {
	      hzclose0->Fill(zclose[0]);
      }

      if (passedselection) {
        ngoodasym ++;
        htphi->Fill( tphi[0]);
        if (true_hel[0] > 0) {
	        htphipos->Fill( tphi[0]);
        } else {
	        htphineg->Fill( tphi[0]);
        }
      }
    }
    
  }
  
  
  TString outfilepdf = outfilename; // Make a pdf file to save these histograms to
  outfilepdf.ReplaceAll(".root",".pdf");

  
  // Make some plots for us to look at
  TCanvas *c1 = new TCanvas("c1","",1600,1200);

  hzclose->SetLineColor(kBlack);
  hzclose->SetLineWidth(2);
  hzclose0->SetLineColor(kRed);
  hzclose0->SetLineWidth(2);
  hzclose0->SetFillColor(kRed - 9);

  hzclose1->SetLineColor(kBlue);
  hzclose1->SetLineWidth(2);
  hzclose1->SetFillColor(kBlue - 9);
  c1->cd();
  c1->Divide(2,2);
  c1->cd(1);    hsclose->Draw();
  c1->cd(2);    httheta->Draw();
  c1->cd(3);    hzclose->Draw();    hzclose0->Draw("same");    hzclose1->Draw("same");
  // Create legend
  TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);  // x1,y1,x2,y2 in NDC
  leg->AddEntry(hzclose,  "All events", "l");
  leg->AddEntry(hzclose0, "#theta < 1 deg", "fl");  // f = filled
  leg->AddEntry(hzclose1, "#theta > 1 deg", "fl");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);   // transparent
  leg->SetTextSize(0.03);
  leg->Draw();
  c1->cd(4);    hzvtheta->Draw();
  c1->Update();
  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page



  // Let's look at the asymmetry now
  TH1D *htphidif = (TH1D*) htphipos->Clone("htphidif");
  htphidif->Add(htphineg, -1);

  TH1D *htphiasym = (TH1D*) htphidif->Clone("htphiasym");
  htphiasym->Divide(htphi); //Our helicity asymmetry
  htphiasym->SetTitle("Helicity Asymmetry (N_{+}-N_{-})/(N_{+}+N_{-}), fit=s_{1}sin(#phi)+c_{1}cos(#phi)");
  htphiasym->GetXaxis()->SetTitle("#phi (rad)");
  htphiasym->GetYaxis()->SetTitle("(N_{+}-N_{-})/(N_{+}+N_{-})");

  TF1 *fitfunc = new TF1("fitfunc", "[0]*sin(x) + [1]*cos(x)", -M_PI, M_PI);  //Prepare our sin+cos fit for the asymmetry
  fitfunc->SetParameters(0.1, 0.1);
  fitfunc->SetParNames("s_{1}", "c_{1}");
  gStyle->SetOptFit(111);
  htphiasym->SetStats(kTRUE);
  
  TCanvas *c2 = new TCanvas("c2","",1600,1200);
  c2->cd();
  c2->Divide(1,2);
  c2->cd(1);
  gPad->cd();
  DrawAsymmetry(htphi,1,kBlack,2,true,true);    // Plotting the helicity separated phi distributions
  DrawAsymmetry(htphipos,1,kRed,2,false,true);
  DrawAsymmetry(htphineg,1,kBlue,2,false,true);
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);  // x1,y1,x2,y2 in NDC
  leg2->AddEntry(htphi,  "N_{+}+N_{-}", "p");
  leg2->AddEntry(htphipos, "N_{+}", "p");
  leg2->AddEntry(htphineg, "N_{-}", "p");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);   // transparent
  leg2->SetTextSize(0.05);
  leg2->Draw();

  c2->cd(2);
  gPad->cd();
  htphiasym->Fit(fitfunc, "R");   // Fitting the asymmetry
  DrawAsymmetry(htphiasym, 1, kBlack, 2);
  int nbins = htphiasym->GetNbinsX(); // Drawing error bars as vertical lines for cleaner visuals
  for (int i = 1; i <= nbins; ++i) {
    double x = htphiasym->GetBinCenter(i);
    double y = htphiasym->GetBinContent(i);
    double err = htphiasym->GetBinError(i);

    TLine *vline = new TLine(x, y - err, x, y + err);
    vline->SetLineColor(kBlack);
    vline->SetLineWidth(1);
    vline->Draw("SAME");
  }
  c2->Update();
  TPaveStats *st = (TPaveStats*)htphiasym->GetListOfFunctions()->FindObject("stats");   //Adjusting stat box for cleaner visuals
  if (st) {
    st->SetX1NDC(0.75); // left
    st->SetX2NDC(0.95); // right
    st->SetY1NDC(0.75); // bottom
    st->SetY2NDC(0.9); // top
    st->Draw();
  }
  c2->Modified();
  c2->Update();
  c2->Print(outfilepdf + "");
  

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
  pt->AddText(Form("Total # events passed selection cuts: %lld", ngoodasym));
  TText *t2 = pt->AddText(" Global cuts: ");
  t2->SetTextColor(kRed);
  AddWrappedText(pt, globalcut.GetTitle());
  TText *t3 = pt->AddText(" Selection cuts: ");
  t3->SetTextColor(kRed);
  AddWrappedText(pt, selectioncut.GetTitle());
  //  TText *t3 = pt->AddText(Form(" Average Charge Normalized Elastic Yield: %.0f", normyield * fsig / totalcharge));
  //t3->SetTextColor(kRed);
  //pt->AddText(Form(" Signal fraction: %.5f ",fsig));
  //TText *t4 = pt->AddText(Form(" Fit Integrated Elastic Yield Y: %.0f", ElasticYieldY));
  //t4->SetTextColor(kRed);
  pt->Draw();
  cSummary->Print(outfilepdf + ")"); // Close out the pdf

  
  
  //  cout << "Total number of events: " << nevent << endl;
  //cout << "Number of events passing elastic cuts: " << ngoodevs << endl;
  //cout << "Our global cut: " << globalcut << endl;
  fout->Write();
  cout << "That's all folks!" << endl;
}

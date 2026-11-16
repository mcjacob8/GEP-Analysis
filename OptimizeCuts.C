// ------------------------------------------------------------------------ //
// Script to optimize several cuts for selecting elastic events from GEp.   //
//
// The basic quality cuts that are applied are:
//     earm.ecal.nblk>2
//     sbs.hcal.nblk>1
//     (sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)
//     (sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)
//     sbs.gemFPP.track.ngoodhits[0]>1
//     sbs.gemFPP.track.sclose<0.005
//
// The following cuts are optimized in this script:
//     heep.dt_ADC
//     heep.dpp
//     (sbs.tr.vz[0]+0.1)
//     sqrt(pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*6.7-sbs.hcal.y+0.0132)/0.0425,2)+pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*6.7-sbs.hcal.x-0.192)/0.0497,2))<=3.5
//     sqrt(pow((heep.dxECAL-0.01+0.006*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL-(0.0029+0.00781*earm.ecal.x+0.0181*pow(earm.ecal.x,2)+0.00378*pow(earm.ecal.x,3)))/0.0249,2))<=3.5

//                                                                          //
// ---------                                                                //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 11-16-2026                  //
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
#include "TLegend.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPaveStats.h"

// Helper function: Fits the core of the peak in a histogram with a Gaussian
// The fit range is the region where the content is above half the maximum bin
TF1* FitCoreGaus(TH1D *h, const char *fname = "fcore") {
  const int nb = h->GetNbinsX();
  if (h->GetEntries() < 100) return nullptr;

  int    maxbin = h->GetMaximumBin();
  double maxval = h->GetBinContent(maxbin);
  double center = h->GetBinCenter(maxbin);
  double half   = 0.5 * maxval;

  // Walk outward from the peak until the content drops below half max
  int lo = maxbin, hi = maxbin;
  while (lo > 1  && h->GetBinContent(lo - 1) > half) lo--;
  while (hi < nb && h->GetBinContent(hi + 1) > half) hi++;

  // Need enough bins to constrain 3 parameters
  if (hi - lo + 1 < 5) return nullptr;

  double fitlo = h->GetXaxis()->GetBinLowEdge(lo);
  double fithi = h->GetXaxis()->GetBinUpEdge(hi);
  double fwhm  = fithi - fitlo;

  TF1 *f = new TF1(fname, "gaus", fitlo, fithi);
  f->SetParameters(maxval, center, fwhm / 2.355);

  int status = h->Fit(f, "RQ");   // R: use function range, Q: quiet, 0: don't draw
  return (status == 0) ? f : nullptr;
}

// Two-stage fit: core Gaussian first, then gaus + background over a wider range.
// bkgOrder: 0 = constant, 1 = linear. Returns the TF1 or nullptr on failure.
TF1* FitGausBkg(TH1D *h, double nsigRange = 5.0, int bkgOrder = 0,
                const char *fname = "fgb") {
  // Stage 1: core fit for starting values
  TF1 *fc = FitCoreGaus(h, "fcore_tmp");
  if (!fc) return nullptr;
  double amp0 = fc->GetParameter(0);
  double mu0  = fc->GetParameter(1);
  double sig0 = fabs(fc->GetParameter(2));
  delete fc;

  // Stage 2: wider range, clipped to the histogram axis
  double xlo = std::max(h->GetXaxis()->GetXmin(), mu0 - nsigRange * sig0);
  double xhi = std::min(h->GetXaxis()->GetXmax(), mu0 + nsigRange * sig0);

  TString form = (bkgOrder == 0) ? "gaus(0) + pol0(3)" : "gaus(0) + pol1(3)";
  TF1 *f = new TF1(fname, form, xlo, xhi);

  // Background seed: average content of the bins at the edges of the fit range
  int blo = h->FindBin(xlo), bhi = h->FindBin(xhi);
  int nedge = std::max(1, (bhi - blo + 1) / 10);
  double bkg = 0;
  for (int i = 0; i < nedge; i++)
    bkg += h->GetBinContent(blo + i) + h->GetBinContent(bhi - i);
  bkg /= (2.0 * nedge);

  f->SetParameter(0, std::max(amp0 - bkg, 1.0));
  f->SetParameter(1, mu0);
  f->SetParameter(2, sig0);
  f->SetParameter(3, bkg);
  if (bkgOrder == 1) f->SetParameter(4, 0);
  f->SetParLimits(2, 0.2 * sig0, 5 * sig0);   // keep sigma from running away

  int status = h->Fit(f, "RQ0");
  return (status == 0) ? f : nullptr;
}

// **** ========== Main functions ========== **** 
void OptimizeCuts( const char *configfilename="configCuts.cfg", const char *outfilename="output/OptimizedCuts.root"){
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0); // No stats box on histograms

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();

  double HCal_Dist = 6.7; // Distance from target to HCal in meters, default Kin3
  TChain *C = new TChain("T");
  cout << "Ok, let's get started!" << endl;

  // Reading in the config file
  ifstream infile(configfilename);
  if( !infile ) return;

  TString currentline;

  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();

	if( skey == "HCal_Dist" ){
	  HCal_Dist = sval.Atof();
	}
        }
        }
    }

  cout << "Let's ball" << endl;
  
  // Time to initialize branches we will use:
  const int MAXHEEP = 100;

  // List of variables to cut on
  double dxECAL[MAXHEEP], dyECAL[MAXHEEP];
  double ecalx[MAXHEEP], ecaly[MAXHEEP];
  double hcalx[MAXHEEP], hcaly[MAXHEEP];
  double FPP_x[MAXHEEP], FPP_xp[MAXHEEP];
  double FPP_y[MAXHEEP], FPP_yp[MAXHEEP];
  double eprime_eth[MAXHEEP], ecalo[MAXHEEP];
  double dt_ADC[MAXHEEP];
  double tr_vz[1];
  double dpp[MAXHEEP], pp[MAXHEEP];
  double sclose[MAXHEEP], zclose[MAXHEEP];

  // Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  // the * applies it to all branches, the 0 disables those branches. to enable would need to make 1
  C->SetBranchStatus("*",0);

  // Initializing branches
  C->SetBranchStatus("heep.*",1);
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("earm.ecal.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("sbs.hcal.nblk*",1);
  C->SetBranchStatus("sbs.hcal.y*",1);
  C->SetBranchStatus("sbs.hcal.x*",1);

  // Filling Arrays
  C->SetBranchAddress("heep.dxECAL", dxECAL);
  C->SetBranchAddress("heep.dyECAL", dyECAL);
  C->SetBranchAddress("earm.ecal.x", ecalx);
  C->SetBranchAddress("earm.ecal.y", ecaly);
  C->SetBranchAddress("sbs.hcal.x", hcalx);
  C->SetBranchAddress("sbs.hcal.y", hcaly);
  C->SetBranchAddress("sbs.gemFPP.track.x", FPP_x);
  C->SetBranchAddress("sbs.gemFPP.track.xp", FPP_xp);
  C->SetBranchAddress("sbs.gemFPP.track.y", FPP_y);
  C->SetBranchAddress("sbs.gemFPP.track.yp", FPP_yp);
  C->SetBranchAddress("sbs.gemFPP.track.sclose", sclose);
  C->SetBranchAddress("sbs.gemFPP.track.zclose", zclose);
  C->SetBranchAddress("heep.eprime_eth", eprime_eth);
  C->SetBranchAddress("heep.ecalo", ecalo);
  C->SetBranchAddress("heep.dpp", dpp);
  C->SetBranchAddress("heep.pp", pp);
  C->SetBranchAddress("heep.dt_ADC", dt_ADC);
  C->SetBranchAddress("sbs.tr.vz", tr_vz);

  TFile *fout = new TFile(outfilename, "RECREATE");

  // Set up the histograms we want to look at:
  int bins = 100;   // useful for normalizing later on
  double range = 0.25; //hard coding these for Kin3 cuts at the moment
  TH1D *hsclose = new TH1D("hsclose", "Distance of Closest Approach; DOCA (m)", bins, 0.0, 0.006);
  TH1D *htr_vz = new TH1D("htr_vz", "Target Vertex Position; z_{tr} (m)", bins, -0.3, 0.2);
  TH1D *hdt_ADC = new TH1D("hdt_ADC", "heep.dt_ADC; t (ns)", bins, -15, 15);

  TH1D *hdpp = new TH1D("hdpp", "pp/pp_th-1; dpp;", bins, -0.1, 0.13);
  TH1D *hpresid_x = new TH1D("hpresid_x", "Proton Residuals in X; Residual (m)", bins, -0.1, 0.5);
  TH1D *hpresid_y = new TH1D("hpresid_y", "Proton Residuals in Y; Residual (m)", bins, -0.3, 0.25);

  //TH1D *hdxECAL = new TH1D("hdxECAL", "heep.dxECAL after global cut; heep.dxECAL (m);", bins, -range/2, range/2);
  //TH1D *hdyECAL = new TH1D("hdyECAL", "heep.dyECAL after global cut; heep.dyECAL (m);", bins, -range/2, range/2);

  //TH2D *hdxECAL_v_dyECAL = new TH2D("hdxECAL_v_dyECAL", "heep.dxECAL vs heep.dyECAL ; heep.dyECAL (m); heep.dxECAL (m)", bins/2, -range/2, range/2, bins/2, -range/2, range/2);

  //TH1D *hEdivP = new TH1D("hEdivP", "E/P after global cut; E/P;", bins, 0.0, 1.3);


  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );

  int treenum=0, currenttreenum=0;
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
      hsclose->Fill( sclose[0]);
      htr_vz->Fill( tr_vz[0]);
      hdt_ADC->Fill( dt_ADC[0]);
    }
  }

  TCanvas *c1 = new TCanvas("c1","",1600,1200);

  c1->cd();
  c1->Divide(2,2);
  c1->cd(1);    hsclose->Draw();
  c1->cd(2);    htr_vz->Draw();

  c1->cd(3);    hdt_ADC->Draw();
  TF1 *f = FitGausBkg(hdt_ADC, 5.0, 0);
  if (f) {
    double mu = f->GetParameter(1), sig = fabs(f->GetParameter(2));
    double nsigma = 3.5;

    cout << "dt_ADC: mu = " << mu << ", sigma = " << sig << endl;

    globalcut += TCut(Form("abs(heep.dt_ADC-%g)<%g", mu, nsigma * sig));
    cout << "Adding cut: abs(heep.dt_ADC-" << mu << ")<" << nsigma * sig << endl;
    f->Draw("same");
  }
  c1->Update();



  //*********************************************************//
  // Begin loop to optimize hcal track projections           //
  //*********************************************************//
  nevent=0;
  //GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
  treenum=0, currenttreenum=0;
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
      hpresid_x->Fill( FPP_x[0] + FPP_xp[0] * HCal_Dist - hcalx[0]);
      hpresid_y->Fill( FPP_y[0] + FPP_yp[0] * HCal_Dist - hcaly[0]);
    }
  }

  TCanvas *c2 = new TCanvas("c2","",1600,1200);
  c2->cd();
  c2->Divide(2,2);
  c2->cd(1);    hpresid_x->Draw();
  TF1 *fpresid_x = FitGausBkg(hpresid_x, 5.0, 0);
  if (fpresid_x) {
    double mu = fpresid_x->GetParameter(1), sig = fabs(fpresid_x->GetParameter(2));
    cout << "presid_x: mu = " << mu << ", sigma = " << sig << endl;
    fpresid_x->Draw("same");
  }
  c2->cd(2);    hpresid_y->Draw();
  TF1 *fpresid_y = FitGausBkg(hpresid_y, 5.0, 0);
  if (fpresid_y) {
    double mu = fpresid_y->GetParameter(1), sig = fabs(fpresid_y->GetParameter(2));
    cout << "presid_y: mu = " << mu << ", sigma = " << sig << endl;
    fpresid_y->Draw("same");
  }
  c2->Update();
  if (fpresid_x && fpresid_y) {
    double nsigma = 3.5;
    globalcut += TCut(Form("sqrt(pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*%g-sbs.hcal.y-%g)/%g,2)+pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*%g-sbs.hcal.x-%g)/%g,2))<=%g",
                            HCal_Dist, fpresid_y->GetParameter(1), nsigma * fabs(fpresid_y->GetParameter(2)),
                            HCal_Dist, fpresid_x->GetParameter(1), nsigma * fabs(fpresid_x->GetParameter(2)),
                            nsigma));
    cout << "Adding cut: sqrt(pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*" << HCal_Dist << "-sbs.hcal.y-" << fpresid_y->GetParameter(1) << ")/" << nsigma * fabs(fpresid_y->GetParameter(2)) << ",2)+pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*" << HCal_Dist << "-sbs.hcal.x-" << fpresid_x->GetParameter(1) << ")/" << nsigma * fabs(fpresid_x->GetParameter(2)) << ",2))<=" << nsigma << endl;
  }



  //Begin next loop to optimize dpp cut
  nevent=0;
  //GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
  treenum=0, currenttreenum=0;
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
      hdpp->Fill( dpp[0]);
    }
  }

  TCanvas *c3 = new TCanvas("c3","",1600,1200);
  c3->cd();
  c3->Divide(2,2);
  c3->cd(1);    hdpp->Draw();
  TF1 *fdpp = FitGausBkg(hdpp, 5.0, 0);
  if (fdpp) {
    double mu = fdpp->GetParameter(1), sig = fabs(fdpp->GetParameter(2));
    double nsigma = 3.5;

    cout << "dpp: mu = " << mu << ", sigma = " << sig << endl;

    globalcut += TCut(Form("abs(heep.dpp-%g)<%g", mu, nsigma * sig));
    cout << "Adding cut: abs(heep.dpp-" << mu << ")<" << nsigma * sig << endl;
    fdpp->Draw("same");
  }
  c3->Update();


  cout << "That's all folks!" << endl;
}
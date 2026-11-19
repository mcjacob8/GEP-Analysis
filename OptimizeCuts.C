// ------------------------------------------------------------------------------------ //
// Script to optimize several cuts for selecting elastic events from GEp.               //
//                                                                                      //
// The basic quality cuts that are applied are:                                         //
//     earm.ecal.nblk>2                                                                 //
//     sbs.hcal.nblk>1                                                                  //
//     (sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)                     //
//     (sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)                   //
//     sbs.gemFPP.track.ngoodhits[0]>1                                                  //
//     sbs.gemFPP.track.sclose<0.005                                                    //
//     (sbs.tr.vz[0]+0.1)                                                               //
//                                                                                      //
// The following cuts are optimized (fit) in this script (in order):                    //
//     heep.dt_ADC                                                                      //
//     Eliptical cut on pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*6.7-sbs.hcal.y),2)  //
//                    + pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*6.7-sbs.hcal.x),2)  //
//     Eliptical cut on pow((heep.dxECAL-POLY(earm.ecal.x)),2)                          //
//                    + pow((heep.dyECAL-POLY(earm.ecal.x)),2)                          //
//     heep.dpp                                                                         //
// ------------------------------------------------------------------------------------ //
// To run this script:                                                                  //
//   analyzer -l 'OptimizeCuts.C("configOptimizeCuts.cfg","output/OptimizedCuts.root")' //
//                                                                                      //
// The config file should have the form:                                                //
//    /path/to/replayed/rootfiles                                                       //
//    endlist                                                                           //
//    #Basic quality cuts to apply                                                      //
//    endcut                                                                            //
//    HCal_Dist #config based on Kin1 or Kin3                                           //
//    endconfig                                                                         //
//                                                                                      //
// ---------                                                                            //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 11-16-2026                              //
// ---------                                                                            //
// ** Do not tamper with this sticker! Log any updates to the script above.             //
// ------------------------------------------------------------------------------------ //

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
  if (h->GetEntries() < 100) { cout << "  core: only " << h->GetEntries() << " entries" << endl; return nullptr; }

  int    maxbin = h->GetMaximumBin();
  double maxval = h->GetBinContent(maxbin);
  double half   = 0.5 * maxval;

  int lo = maxbin, hi = maxbin;
  while (lo > 1  && h->GetBinContent(lo - 1) > half) lo--;
  while (hi < nb && h->GetBinContent(hi + 1) > half) hi++;
  cout << "  core: max bin " << maxbin << " (x=" << h->GetBinCenter(maxbin) << "), max=" << maxval
       << ", half-max bins " << lo << "-" << hi << " of " << nb << endl;

  if (hi - lo + 1 < 5) { cout << "  core: fewer than 5 bins above half max" << endl; return nullptr; }

  double fitlo = h->GetXaxis()->GetBinLowEdge(lo);
  double fithi = h->GetXaxis()->GetBinUpEdge(hi);
  TF1 *f = new TF1(fname, "gaus", fitlo, fithi);
  f->SetParameters(maxval, h->GetBinCenter(maxbin), (fithi - fitlo) / 2.355);

  TFitResultPtr r = h->Fit(f, "RQ0S");
  cout << "  core: status " << r->Status() << ", valid " << r->IsValid()
       << ", mu=" << f->GetParameter(1) << ", sigma=" << f->GetParameter(2) << endl;
  return r->IsValid() ? f : nullptr;
}

// Helper function: Does a two stage fit: core Gaussian first, then gaus + background over a wider range.
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
  f->SetParNames("Amp", "Mean", "Sigma", "Const");

  int status = h->Fit(f, "RQ");
  return (status == 0) ? f : nullptr;
}

// Helper function: Fits a polynomial to the profile (mean of y in each x bin) of h2.
// The profile is returned via pout for plotting. fitxlo/fitxhi restrict the fit range.
TF1* FitCorrelation(TH2D *h2, int order, const char *tag,
                    TProfile **pout = nullptr,
                    double fitxlo = 0, double fitxhi = 0) {
  if (fitxlo >= fitxhi) {
    fitxlo = h2->GetXaxis()->GetXmin();
    fitxhi = h2->GetXaxis()->GetXmax();
  }

  TProfile *prof = h2->ProfileX(Form("prof_%s", tag));
  if (pout) *pout = prof;

  TF1 *fp = new TF1(Form("fpoly_%s", tag), Form("pol%d", order), fitxlo, fitxhi);
  int status = prof->Fit(fp, "RQ0");
  if (status != 0) cout << tag << ": polynomial fit failed (status " << status << ")" << endl;
  return (status == 0) ? fp : nullptr;
}

// Helper function: Builds a 1D residual histogram (y - f(x)) from the 2D histogram by looping over its bins.
// Each bin is filled at its center with weight = bin content. Elimates an extra loop over the entries in the tree.
TH1D* MakeResidual(TH2D *h2, TF1 *fp, const char *name,
                   int nbins, double lo, double hi) {
  TH1D *hr = new TH1D(name, Form("%s Corrected Residual;residual", h2->GetTitle()), nbins, lo, hi);
  hr->Sumw2();
  TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
  TAxis *ar = hr->GetXaxis();

  for (int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    double fx = fp->Eval(ax->GetBinCenter(ix));
    for (int iy = 1; iy <= h2->GetNbinsY(); iy++) {
      double w = h2->GetBinContent(ix, iy);
      if (w <= 0) continue;
      double a = ay->GetBinLowEdge(iy) - fx;
      double b = ay->GetBinUpEdge(iy)  - fx;

      int ia = std::max(1,     ar->FindBin(a));
      int ib = std::min(nbins, ar->FindBin(b));
      for (int ir = ia; ir <= ib; ir++) {
        double ov = std::min(b, ar->GetBinUpEdge(ir)) - std::max(a, ar->GetBinLowEdge(ir));
        if (ov > 0) hr->Fill(ar->GetBinCenter(ir), w * ov / (b - a));
      }
    }
  }
  return hr;
}

// Helper function: Draws vertical lines at mu +/- nsigma*sigma (and optionally at mu) on the current pad.
// Call after the histogram has been drawn, so the pad's y range is known.
// Set draw = false to skip the lines without touching the calling code.
void DrawCutLines(double mu, double sig, double nsigma = 3.0,
                  bool draw = true, int color = kRed) {
  if (!draw) return;
  gPad->Update();                       // creates the stat box if it doesn't exist yet
  double ymin = gPad->GetUymin(), ymax = gPad->GetUymax();
  double xs[2] = { mu - nsigma * sig, mu + nsigma * sig };

  for (int i = 0; i < 2; i++) {
    TLine *l = new TLine(xs[i], ymin, xs[i], ymax);
    l->SetLineColor(color);
    l->SetLineWidth(2);
    l->SetLineStyle(7);
    l->Draw("same");
  }

  // Redraw the stat box last so it sits in front of the lines
  TPaveStats *st = (TPaveStats*) gPad->GetPrimitive("stats");
  if (st) st->Draw();
  gPad->Modified();
}


// **** -------------------------------------------------------------------------- **** //
// **** ----------------------------- Main function ------------------------------ **** //
// **** -------------------------------------------------------------------------- **** //
void OptimizeCuts( const char *configfilename="configCuts.cfg", const char *outfilename="output/OptimizedCuts2.root"){
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  //gStyle->SetOptStat(0); // No stats box on histograms
  gStyle->SetOptFit(1111);
  gStyle->SetStatX(0.93);   // right edge (NDC)
  gStyle->SetStatY(0.93);   // top edge
  gStyle->SetStatW(0.17);   // width
  gStyle->SetStatH(0.15);   // height
  gStyle->SetStatFontSize(0.03);
  bool showCutLines = true;    // one switch for all plots

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
  TH1D *hpresid_x = new TH1D("hpresid_x", "sbs.gemFPP.track.x+sbs.gemFPP.track.xp*6.7-sbs.hcal.x; Residual (m)", bins, -0.1, 0.5);
  TH1D *hpresid_y = new TH1D("hpresid_y", "sbs.gemFPP.track.y+sbs.gemFPP.track.yp*6.7-sbs.hcal.y; Residual (m)", bins, -0.3, 0.25);

  //TH1D *hdxECAL = new TH1D("hdxECAL", "heep.dxECAL after global cut; heep.dxECAL (m);", bins, -range/2, range/2);
  //TH1D *hdyECAL = new TH1D("hdyECAL", "heep.dyECAL after global cut; heep.dyECAL (m);", bins, -range/2, range/2);

  TH2D *hdxECAL_v_ecal_x = new TH2D("hdxECAL_v_ecal_x", "heep.dxECAL vs earm.ecal.x ; earm.ecal.x (m); heep.dxECAL (m)", bins*2, -1.5, 1.5, 40, -0.15, 0.15);
  TH2D *hdyECAL_v_ecal_x = new TH2D("hdyECAL_v_ecal_x", "heep.dyECAL vs earm.ecal.x ; earm.ecal.x (m); heep.dyECAL (m)", bins*2, -1.5, 1.5, 40, -0.2, 0.2);
  TH1D *heresid_x = new TH1D("heresid_x", "Electron Arm Residual in X After Fit; heep.dxECAL (m);", bins, -range/2, range/2);
  TH1D *heresid_y = new TH1D("heresid_y", "Electron Arm Residual in Y After Fit; heep.dyECAL (m);", bins, -range/2, range/2);


  //Histograms for after* all cuts are applied
  TH1D *hpp = new TH1D("hpp", "Proton Momentum; p (GeV);", bins, 4.5, 9);
  TH2D *hdxECAL_v_dyECAL_before = new TH2D("hdxECAL_v_dyECAL", "heep.dxECAL vs heep.dyECAL ; heep.dyECAL (m); heep.dxECAL (m)", bins, -0.1, 0.1, bins, -0.1, 0.1);
  TH2D *hdxECAL_v_dyECAL_after = new TH2D("hdxECAL_v_dyECAL", "heep.dxECAL vs heep.dyECAL ; heep.dyECAL (m); heep.dxECAL (m)", bins, -0.1, 0.1, bins, -0.1, 0.1);
  

  TH1D *hEdivP = new TH1D("hEdivP", "heep.ecalo/heep.eprime_eth; E_{ECal}/E'_{e}(#theta_{e});", bins, 0.0, 2.0);


  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );

  int treenum=0, currenttreenum=0;

  //*********************************************************//
  // Begin first loop to optimize timing cuts                //
  //*********************************************************//
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

  TString outfilepdf = outfilename; // Make a pdf file to save these histograms to
  outfilepdf.ReplaceAll(".root",".pdf");

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

    DrawCutLines(mu, sig, nsigma, showCutLines);
    f->Draw("same");
  }
  c1->Update();
  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page



  //*********************************************************//
  // Begin loop to optimize hcal track projections           //
  //*********************************************************//
  nevent=0;
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

  TCanvas *c2 = new TCanvas("c2","HCal Track Projections",1600,1200);
  c2->cd();
  c2->Divide(2,1);
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
  c2->Print(outfilepdf + "");
  if (fpresid_x && fpresid_y) {
    double nsigma = 3.5;
    globalcut += TCut(Form("sqrt(pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*%g-sbs.hcal.y-%g)/%g,2)+pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*%g-sbs.hcal.x-%g)/%g,2))<=%g",
                            HCal_Dist, fpresid_y->GetParameter(1), fabs(fpresid_y->GetParameter(2)),
                            HCal_Dist, fpresid_x->GetParameter(1), fabs(fpresid_x->GetParameter(2)),
                            nsigma));
    cout << "Adding cut: sqrt(pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*" << HCal_Dist << "-sbs.hcal.y-" << fpresid_y->GetParameter(1) << ")/" << fabs(fpresid_y->GetParameter(2)) << ",2)+pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*" << HCal_Dist << "-sbs.hcal.x-" << fpresid_x->GetParameter(1) << ")/" << fabs(fpresid_x->GetParameter(2)) << ",2))<=" << nsigma << endl;
  }



  //*********************************************************//
  // Begin loop to optimize ecal track projections           //
  //*********************************************************//
  nevent=0;
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
      hdxECAL_v_ecal_x->Fill( ecalx[0], dxECAL[0]);
      hdyECAL_v_ecal_x->Fill( ecalx[0], dyECAL[0]);
      //hpresid_x->Fill( FPP_x[0] + FPP_xp[0] * HCal_Dist - hcalx[0]);
      //hpresid_y->Fill( FPP_y[0] + FPP_yp[0] * HCal_Dist - hcaly[0]);
      hdxECAL_v_dyECAL_before->Fill( dyECAL[0], dxECAL[0]);
    }
  }

  TProfile *gx = nullptr, *gy = nullptr;
  TF1 *fpx = FitCorrelation(hdxECAL_v_ecal_x, 1, "x", &gx);   // linear
  TF1 *fpy = FitCorrelation(hdyECAL_v_ecal_x, 3, "y", &gy);   // cubic

  TCanvas *c3 = new TCanvas("c3","ECal Correlations",1600,1200);
  c3->cd();
  c3->Divide(2,2);

  c3->cd(1);    hdxECAL_v_ecal_x->Draw("colz");
  TH1D *hresx = MakeResidual(hdxECAL_v_ecal_x, fpx, "hresx", 100, -0.3, 0.3);
  TF1 *fx = FitGausBkg(hresx, 5.0, 0);

  if (fpx) { cout << "x: polynomial fit succeeded" << endl; }
  else cout << "x: polynomial fit failed" << endl;

  if (fpx) {
    hdxECAL_v_ecal_x->Draw("colz"); gx->SetMarkerStyle(20); gx->Draw("P same"); fpx->SetLineColor(kRed); fpx->Draw("same");
    c3->cd(2);   hresx->Draw();
  }

  c3->cd(3);    hdyECAL_v_ecal_x->Draw("colz");
  TH1D *hresy = MakeResidual(hdyECAL_v_ecal_x, fpy, "hresy", 100, -0.3, 0.3);
  TF1 *fy = FitGausBkg(hresy, 5.0, 0);
  if (fy) {
    hdyECAL_v_ecal_x->Draw("colz");  gy->SetMarkerStyle(20); gy->Draw("P same"); fpy->SetLineColor(kRed); fpy->Draw("same");
    c3->cd(4);   hresy->Draw();
  }
  c3->Update();
  c3->Print(outfilepdf + "");

  if (fx && fy) {
    double nsigma = 3.5;
    globalcut += TCut(Form("sqrt(pow((heep.dxECAL-(%g+%g*earm.ecal.x))/%g,2)+pow((heep.dyECAL-(%g+%g*earm.ecal.x+%g*pow(earm.ecal.x,2)+%g*pow(earm.ecal.x,3)))/%g,2))<=%g",
                  fpx->GetParameter(0), fpx->GetParameter(1), fy->GetParameter(2),
                  fpy->GetParameter(0), fpy->GetParameter(1), fpy->GetParameter(2), fpy->GetParameter(3), fpy->GetParameter(2),
                  nsigma));                     
    cout << "Adding cut: sqrt(pow((heep.dxECAL-("
     << fpx->GetParameter(0) << "+"
     << fpx->GetParameter(1) << "*earm.ecal.x))/"
     << fx->GetParameter(2) << ",2)+pow((heep.dyECAL-("
     << fpy->GetParameter(0) << "+"
     << fpy->GetParameter(1) << "*earm.ecal.x+"
     << fpy->GetParameter(2) << "*pow(earm.ecal.x,2)+"
     << fpy->GetParameter(3) << "*pow(earm.ecal.x,3)))/"
     << fy->GetParameter(2) << ",2))<=" << nsigma << endl;
  }



  nevent=0;
  treenum=0, currenttreenum=0;
  //*********************************************************//
  // Begin loop to optimize dpp cut                          //
  //*********************************************************//
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

  TCanvas *c4 = new TCanvas("c4","heep.dpp",1600,1200);
  c4->cd();
  c4->Divide(1,1);
  c4->cd(1);    hdpp->Draw();
  TF1 *fdpp = FitGausBkg(hdpp, 5.0, 0);
  if (fdpp) {
    double mu = fdpp->GetParameter(1), sig = fabs(fdpp->GetParameter(2));
    double nsigma = 3.5;

    cout << "dpp: mu = " << mu << ", sigma = " << sig << endl;

    globalcut += TCut(Form("abs(heep.dpp-%g)<%g", mu, nsigma * sig));
    cout << "Adding cut: abs(heep.dpp-" << mu << ")<" << nsigma * sig << endl;
    DrawCutLines(mu, sig, nsigma, showCutLines);
    fdpp->Draw("same");
  }
  c4->Update();
  c4->Print(outfilepdf + "");



  nevent=0;
  treenum=0, currenttreenum=0;
  //*********************************************************//
  // Final loop to apply all optimized cuts                  //
  //*********************************************************//
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
      hpp->Fill( pp[0]);
      hEdivP->Fill( ecalo[0]/eprime_eth[0]);
      hdxECAL_v_dyECAL_after->Fill( dyECAL[0], dxECAL[0]);
    }
  }

  TCanvas *c5 = new TCanvas("c5","heep.pp",1600,1200);
  c5->cd();
  c5->Divide(2,2);
  c5->cd(1);    hpp->Draw();
  c5->cd(2);    hEdivP->Draw();
  c5->cd(3);    hdxECAL_v_dyECAL_before->Draw("colz");
  c5->cd(4);    hdxECAL_v_dyECAL_after->Draw("colz");
  c5->Update();
  c5->Print(outfilepdf + ")");

  cout << "Final global cut: " << globalcut.GetTitle() << endl;
  fout->Write();
  cout << "That's all folks!" << endl;
}
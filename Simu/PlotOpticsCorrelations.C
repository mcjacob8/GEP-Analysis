// ------------------------------------------------------------------------------------------------------- //
// This helper script generates a canvas that allows one to plot different distributions side by side for  //
// purposes of evaluating the behavior of different correlations for different optics models. It starts by //
// performing a fit on the dipole midplane projection residuals as a function of momentum.                 //
//                                                                                                         //
// ---------                                                                                               //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 10-26-2026                                                 //
// ---------                                                                                               //
// ** Do not tamper with this sticker! Log any updates to the script above.                                //
// ------------------------------------------------------------------------------------------------------- //

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"
#include <iostream>

// Usage (from the ROOT prompt):
//   root -l 'PlotOpticsCorrelations.C("optics_tree.root")'
//   root -l 'PlotOpticsCorrelations.C("optics_tree.root", "abs(xtar)<0.05")'   // optional cut
void PlotOpticsCorrelations( const char *filename,
                            const char *cut = "",
                            const char *outpng = "ytar_vs_p.png" ){

  TFile *f = TFile::Open(filename, "READ");
  if( !f || f->IsZombie() ){
    std::cout << "Could not open " << filename << std::endl;
    return;
  }

  TTree *T = (TTree*) f->Get("Tout");
  if( !T ){
    std::cout << "No tree named Tout in " << filename << std::endl;
    return;
  }

  std::cout << "Tout has " << T->GetEntries() << " entries" << std::endl;

  gStyle->SetOptStat(1111);


  // Let's do some plotting for x residuals
  TCanvas *c1 = new TCanvas("c1", "xtar vs p", 900, 700);
  c1->SetRightMargin(0.13);
  //T->Draw("xtar+xptar*(1.59638+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.)):p>>hxtar_p(100, 0.75, 10.25, 100, -0.01, 0.02)", cut, "colz");
  T->Draw("xtar+xptar*(1.618+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.)):p>>hxtar_p(100, 0.75, 10.25, 100, -0.01, 0.02)", cut, "colz");
  TH2D *h1 = (TH2D*) gDirectory->Get("hxtar_p");
  h1->SetTitle("x_{tar} Residual vs p;p (GeV);x_{tar} Residual(m)");
  h1->Draw("colz");

  // mean in each p bin (error = error on the mean) and with "s" (error = standard deviation)
  TProfile *hprof  = h1->ProfileX("hprof");
  TProfile *hprofx = h1->ProfileX("hprofx", 1, -1, "s");

  hprof->SetMarkerStyle(20);
  hprofx->SetLineColor(kRed);
  hprof->Draw("same");
  hprofx->Draw("same");
  TF1 *fmean = new TF1("fmean", "[0] + [1]*pow(x,[2])", 1, 10);
  TF1 *fsig  = new TF1("fsig",  "[0] + [1]*pow(x,[2])", 1, 10);
  hprof->Fit("fmean", "", "same");        // fit the mean vs p, drawn on the 2D plot
  c1->SaveAs("xtar_vs_p.png");

  // standard deviation vs p: put the spread of each bin into a histogram and fit that
  TH1D *hsig = h1->ProjectionX("hsig");
  hsig->Reset();
  for( int i=1; i<=hsig->GetNbinsX(); i++ ){
    double n = hprofx->GetBinEntries(i);
    if( n < 20 ) continue;
    double s = hprofx->GetBinError(i);
    hsig->SetBinContent(i, s);
    hsig->SetBinError(i, s/sqrt(2.0*n));
  }
  new TCanvas("csig", "sigma vs p", 900, 700);
  hsig->SetTitle("standard deviation of x_{tar} Residual vs p;p (GeV);#sigma");
  hsig->SetMarkerStyle(20);
  hsig->Fit("fsig");


  // Let's do some plotting for y residuals
  TCanvas *c2 = new TCanvas("c2", "ytar vs p", 900, 700);
  c2->SetRightMargin(0.13);
  //T->Draw("ytar+yptar*(1.59638+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.)):p>>hytar_p(100, 0.75, 10.25, 100, -0.01, 0.01)", cut, "colz");
  T->Draw("ytar+yptar*(1.618+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.)):p>>hytar_p(100, 0.75, 10.25, 100, -0.01, 0.01)", cut, "colz");
  TH2D *h2 = (TH2D*) gDirectory->Get("hytar_p");
  h2->SetTitle("y_{tar} Residual vs p;p (GeV);y_{tar} Residual(m)");
  h2->Draw("colz");

  TProfile *h2prof  = h2->ProfileX("h2prof");
  TProfile *h2profx = h2->ProfileX("h2profx", 1, -1, "s");

  h2prof->SetMarkerStyle(20);
  h2profx->SetLineColor(kRed);
  h2prof->Draw("same");
  h2profx->Draw("same");
  TF1 *f2mean = new TF1("f2mean", "[0] + [1]*pow(x,[2])", 1, 10);
  TF1 *f2sig  = new TF1("f2sig",  "[0] + [1]*pow(x,[2])", 1, 10);
  h2prof->Fit("f2mean", "", "same");        // fit the mean vs p, drawn on the 2D plot
  c2->SaveAs("ytar_vs_p.png");

  // standard deviation vs p
  TH1D *h2sig = h2->ProjectionX("h2sig");
  h2sig->Reset();
  for( int i=1; i<=h2sig->GetNbinsX(); i++ ){
    double n = h2profx->GetBinEntries(i);
    if( n < 20 ) continue;
    double s = h2profx->GetBinError(i);
    h2sig->SetBinContent(i, s);
    h2sig->SetBinError(i, s/sqrt(2.0*n));
  }
  new TCanvas("csig", "sigma vs p", 900, 700);
  h2sig->SetTitle("standard deviation of y_{tar} Residual vs p;p (GeV);#sigma");
  h2sig->SetMarkerStyle(20);
  h2sig->Fit("f2sig");


  // This last section I will hard coded some plotting for looking at different correlations
  //TFile *f1 = TFile::Open("SBS_optics_fitresult_11.root", "READ");
  //TFile *f2 = TFile::Open("SBS_optics_fitresult_12.root", "READ");
  //TFile *f3 = TFile::Open("SBS_optics_fitresult_13.root", "READ");
  //TFile *f4 = TFile::Open("SBS_optics_fitresult_14.root", "READ");
  //const char *cut1 = "";
  //const char *cut2 = "p>3.0";
  //const char *cut3 = "abs(xtar+xptar*(1.59638+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.))+1.269e-3-0.01907/pow(p,0.74182))<=3.0*(8.6718e-4+0.0081846/pow(p,1.38916))&&abs(ytar+yptar*(1.59638+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.))+3.2995e-5-0.001701/pow(p,2.30806))<=3.0*(8.876e-4+9.04423e-3/pow(p,1.41323))";
  //const char *cut4 = "abs(xtar+xptar*(1.59638+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.))+1.269e-3-0.01907/pow(p,0.74182))<=3.0*(8.6718e-4+0.0081846/pow(p,1.38916))&&abs(ytar+yptar*(1.59638+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.))+3.2995e-5-0.001701/pow(p,2.30806))<=3.0*(8.876e-4+9.04423e-3/pow(p,1.41323))&&p>3.0";

  TFile *f1 = TFile::Open("db_Kin3_fitresults.root", "READ");
  TFile *f2 = TFile::Open("SBS_optics_fitresult_16.root", "READ");
  TFile *f3 = TFile::Open("db_Kin3_fitresults.root", "READ");
  TFile *f4 = TFile::Open("SBS_optics_fitresult_16.root", "READ");
  const char *cut1 = "abs(xtar+xptar*(1.618+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.))+1.41e-2-0.02631/pow(p,0.20046))<=3.0*(1.2532e-4+0.0091/pow(p,1.1492))&&abs(ytar+yptar*(1.618+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.))+6.682e-5-3.532e-3/pow(p,1.73556))<=3.0*(6.597e-4+6.2484e-3/pow(p,1.06111))&&p>5.5&&p<8.0";
  const char *cut2 = "abs(xtar+xptar*(1.618+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.))+1.41e-2-0.02631/pow(p,0.20046))<=3.0*(1.2532e-4+0.0091/pow(p,1.1492))&&abs(ytar+yptar*(1.618+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.))+6.682e-5-3.532e-3/pow(p,1.73556))<=3.0*(6.597e-4+6.2484e-3/pow(p,1.06111))&&p>5.5&&p<8.0";
  const char *cut3 = "abs(xtar+xptar*(1.618+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.))+1.41e-2-0.02631/pow(p,0.20046))<=3.0*(1.2532e-4+0.0091/pow(p,1.1492))&&abs(ytar+yptar*(1.618+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.))+6.682e-5-3.532e-3/pow(p,1.73556))<=3.0*(6.597e-4+6.2484e-3/pow(p,1.06111))&&p>4.5&&p<9.5";  //9.5
  const char *cut4 = "abs(xtar+xptar*(1.618+1.22/2.)-(xfp-xpfp*(2.085-1.22/2.))+1.41e-2-0.02631/pow(p,0.20046))<=3.0*(1.2532e-4+0.0091/pow(p,1.1492))&&abs(ytar+yptar*(1.618+1.22/2.)-(yfp-ypfp*(2.085-1.22/2.))+6.682e-5-3.532e-3/pow(p,1.73556))<=3.0*(6.597e-4+6.2484e-3/pow(p,1.06111))&&p>4.5&&p<9.5";  //9.0


  TTree *T1 = (TTree*) f1->Get("Tout");
  TTree *T2 = (TTree*) f2->Get("Tout");
  TTree *T3 = (TTree*) f3->Get("Tout");
  TTree *T4 = (TTree*) f4->Get("Tout");

  TCanvas *c3 = new TCanvas("c3", "ytarrecon-ytar vs ytar", 1500, 900);
  c3->Divide(2, 2);
  c3->SetRightMargin(0.13);
  c3->cd(1);
  T1->Draw("(ytarrecon-ytar):ytar>>hytar1(100, -0.135, 0.06, 100, -0.02, 0.02)", cut1, "colz");
  TH2D *h10 = (TH2D*) gDirectory->Get("hytar1");
  h10->SetTitle("ytarrecon-ytar:ytar;ytar (m);ytarrecon-ytar(m)");
  h10->Draw("colz");
  c3->cd(2);
  T2->Draw("(ytarrecon-ytar):ytar>>hytar2(100, -0.135, 0.06, 100, -0.02, 0.02)", cut2, "colz");
  TH2D *h20 = (TH2D*) gDirectory->Get("hytar2");
  h20->SetTitle("ytarrecon-ytar:ytar;ytar (m);ytarrecon-ytar(m)");
  h20->Draw("colz");
  c3->cd(3);
  T3->Draw("(ytarrecon-ytar):ytar>>hytar3(100, -0.135, 0.06, 100, -0.02, 0.02)", cut3, "colz");
  TH2D *h30 = (TH2D*) gDirectory->Get("hytar3");
  h30->SetTitle("ytarrecon-ytar:ytar;ytar (m);ytarrecon-ytar(m)");
  h30->Draw("colz");
  c3->cd(4);
  T4->Draw("(ytarrecon-ytar):ytar>>hytar4(100, -0.135, 0.06, 100, -0.02, 0.02)", cut4, "colz");
  TH2D *h40 = (TH2D*) gDirectory->Get("hytar4");
  h40->SetTitle("ytarrecon-ytar:ytar;ytar (m);ytarrecon-ytar(m)");
  h40->Draw("colz");
  c3->Update();
  c3->SaveAs("ytar_vs_ytar.png");

  TCanvas *c4 = new TCanvas("c4", "yptarrecon-yptar vs ytar", 1500, 900);
  c4->Divide(2, 2);
  c4->SetRightMargin(0.13);
  c4->cd(1);
  T1->Draw("(yptarrecon-yptar):ytar>>hyptar1(100, -0.135, 0.06, 100, -0.01, 0.01)", cut1, "colz");
  TH2D *h11 = (TH2D*) gDirectory->Get("hyptar1");
  h11->SetTitle("yptarrecon-yptar:ytar;ytar (m);yptarrecon-yptar(rad)");
  h11->Draw("colz");
  c4->cd(2);
  T2->Draw("(yptarrecon-yptar):ytar>>hyptar2(100, -0.135, 0.06, 100, -0.01, 0.01)", cut2, "colz");
  TH2D *h21 = (TH2D*) gDirectory->Get("hyptar2");
  h21->SetTitle("yptarrecon-yptar:ytar;ytar (m);yptarrecon-yptar(rad)");
  h21->Draw("colz");
  c4->cd(3);
  T3->Draw("(yptarrecon-yptar):ytar>>hyptar3(100, -0.135, 0.06, 100, -0.01, 0.01)", cut3, "colz");
  TH2D *h31 = (TH2D*) gDirectory->Get("hyptar3");
  h31->SetTitle("yptarrecon-yptar:ytar;ytar (m);yptarrecon-yptar(rad)");
  h31->Draw("colz");
  c4->cd(4);
  T4->Draw("(yptarrecon-yptar):ytar>>hyptar4(100, -0.135, 0.06, 100, -0.01, 0.01)", cut4, "colz");
  TH2D *h41 = (TH2D*) gDirectory->Get("hyptar4");
  h41->SetTitle("yptarrecon-yptar:ytar;ytar (m);yptarrecon-yptar(rad)");
  h41->Draw("colz");
  c4->Update();
  c4->SaveAs("yptar_vs_ytar.png");

  TCanvas *c5 = new TCanvas("c5", "ytarrecon-ytar vs p", 1500, 900);
  c5->Divide(2, 2);
  c5->SetRightMargin(0.13);
  c5->cd(1);
  T1->Draw("(ytarrecon-ytar):p>>hytarp1(100, 0.5, 10.5, 100, -0.02, 0.02)", cut1, "colz");
  TH2D *h12 = (TH2D*) gDirectory->Get("hytarp1");
  h12->SetTitle("ytarrecon-ytar:p;p (GeV);ytarrecon-ytar(m)");
  h12->Draw("colz");
  c5->cd(2);
  T2->Draw("(ytarrecon-ytar):p>>hytarp2(100, 0.5, 10.5, 100, -0.02, 0.02)", cut2, "colz");
  TH2D *h22 = (TH2D*) gDirectory->Get("hytarp2");
  h22->SetTitle("ytarrecon-ytar:p;p (GeV);ytarrecon-ytar(m)");
  h22->Draw("colz");
  c5->cd(3);
  T3->Draw("(ytarrecon-ytar):p>>hytarp3(100, 0.5, 10.5, 100, -0.02, 0.02)", cut3, "colz");
  TH2D *h32 = (TH2D*) gDirectory->Get("hytarp3");
  h32->SetTitle("ytarrecon-ytar:p;p (GeV);ytarrecon-ytar(m)");
  h32->Draw("colz");
  c5->cd(4);
  T4->Draw("(ytarrecon-ytar):p>>hytarp4(100, 0.5, 10.5, 100, -0.02, 0.02)", cut4, "colz");
  TH2D *h42 = (TH2D*) gDirectory->Get("hytarp4");
  h42->SetTitle("ytarrecon-ytar:p;p (GeV);ytarrecon-ytar(m)");
  h42->Draw("colz");
  c5->Update();
  c5->SaveAs("ytar_vs_p.png");

  TCanvas *c6 = new TCanvas("c6", "precon/p - 1 vs p", 1500, 900);
  c6->Divide(2, 2);
  c6->SetRightMargin(0.13);
  c6->cd(1);
  T1->Draw("(precon/p - 1):p>>hprecon1(100, 0.5, 10.5, 100, -0.1, 0.1)", cut1, "colz");
  TH2D *h13 = (TH2D*) gDirectory->Get("hprecon1");
  h13->SetTitle("precon/p - 1:p;p (GeV);precon/p - 1");
  h13->Draw("colz");
  c6->cd(2);
  T2->Draw("(precon/p - 1):p>>hprecon2(100, 0.5, 10.5, 100, -0.1, 0.1)", cut2, "colz");
  TH2D *h23 = (TH2D*) gDirectory->Get("hprecon2");
  h23->SetTitle("precon/p - 1:p;p (GeV);precon/p - 1");
  h23->Draw("colz");
  c6->cd(3);
  T3->Draw("(precon/p - 1):p>>hprecon3(100, 0.5, 10.5, 100, -0.1, 0.1)", cut3, "colz");
  TH2D *h33 = (TH2D*) gDirectory->Get("hprecon3");
  h33->SetTitle("precon/p - 1:p;p (GeV);precon/p - 1");
  h33->Draw("colz");
  c6->cd(4);
  T4->Draw("(precon/p - 1):p>>hprecon4(100, 0.5, 10.5, 100, -0.1, 0.1)", cut4, "colz");
  TH2D *h43 = (TH2D*) gDirectory->Get("hprecon4");
  h43->SetTitle("precon/p - 1:p;p (GeV);precon/p - 1");
  h43->Draw("colz");
  c6->Update();
  c6->SaveAs("precon_vs_p.png");


  TCanvas *c7 = new TCanvas("c7", "precon/p - 1 vs p", 1500, 900);
  c7->Divide(2, 2);
  c7->SetRightMargin(0.13);
  c7->cd(1);
  T1->Draw("ytarrecon-ytar>>hy1(100,-0.01,0.01)", cut1);
  TH1D *h14 = (TH1D*) gDirectory->Get("hy1");
  h14->SetTitle("ytarrecon-ytar;ytarrecon-ytar (m);counts");
  h14->Draw();
  c7->cd(2);
  T2->Draw("ytarrecon-ytar>>hy2(100,-0.01,0.01)", cut2);
  TH1D *h24 = (TH1D*) gDirectory->Get("hy2");
  h24->SetTitle("ytarrecon-ytar;ytarrecon-ytar (m);counts");
  h24->Draw();
  c7->cd(3);
  T3->Draw("ytarrecon-ytar>>hy3(100,-0.01,0.01)", cut3);
  TH1D *h34 = (TH1D*) gDirectory->Get("hy3");
  h34->SetTitle("ytarrecon-ytar;ytarrecon-ytar (m);counts");
  h34->Draw();
  c7->cd(4);
  T4->Draw("ytarrecon-ytar>>hy4(100,-0.01,0.01)", cut4);
  TH1D *h44 = (TH1D*) gDirectory->Get("hy4");
  h44->SetTitle("ytarrecon-ytar;ytarrecon-ytar (m);counts");
  h44->Draw();
  c7->Update();


}
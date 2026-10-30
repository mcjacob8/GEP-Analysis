// ------------------------------------------------------------------------------------------------------- //
// This is a simplehelper script for comparing fits for dipole midplane projection residuals for optics    //
// modeling.                                                                                               //
//                                                                                                         //
// ---------                                                                                               //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 8-29-2026                                                  //
// ---------                                                                                               //
// ** Do not tamper with this sticker! Log any updates to the script above.                                //
// ------------------------------------------------------------------------------------------------------- //


#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include <algorithm>

// Draws the four fit functions f(p) = a + b*p^c on a 2x2 canvas.
// A second set of parameters (from your other data) is overlaid on the same pads.
//
// Usage:  root -l draw_fits.C
//
// Fill in the parameters below. Each row is {a, b, c} for one fit.

void CompareFits(){

  // ---- Edit these -----------------------------------------------------------
  const int nfits = 4;

  // Pad titles, one per fit (these are placeholders, rename them)
  const char *titles[nfits] = {
    "Fit 1;p (GeV);value",
    "Fit 2;p (GeV);value",
    "Fit 3;p (GeV);value",
    "Fit 4;p (GeV);value"
  };

  // First data set: {a, b, c} for each fit (placeholder numbers, replace with yours)
  // Original from energy loss cut, wide histogram
  double par3[nfits][3] = {
    { -0.12459, 0.133844, -0.0232155 },
    { 0.00257429, 0.0160591, -3.92505 },
    { 0.0278524, -0.0268968, 0.0193966 },
    { 0.00106795, 0.0107266, -1.34561 }
  };

  // Second data set, drawn over the first. Set draw_second = false to only draw the first.
  // Original cuts I've been using
  bool draw_second = true;
  double par2[nfits][3] = {
    { -1.41e-2, 0.02631, -0.20046 },
    { 1.2532e-4, 0.0091, -1.1492 },
    { -6.682e-5, 3.532e-3, -1.73556 },
    { 6.597e-4, 6.2484e-3, -1.06111 }
  };

  // Cuts from original distributions, wide histograms
  double par5[nfits][3] = {
    { -0.0935463, 0.102709, -0.0300443 },
    { 0.00248319, 0.0159179, -2.54584 },
    { 0.0360191, -0.0352607, 0.012204 },
    { 0.00106652, 0.0161007, -1.52479 }
  };

  // Cuts from original distributions, narrow histograms
  double par4[nfits][3] = {
    { -0.11243, 0.1216, -0.0252001 },
    { 0.00231494, 0.0115821, -2.86041 },
    { 0.0154565, -0.0147923, 0.0256719 },
    { 0.000322836, 0.00785046, -0.978336 }
  };

  // Cuts from new distributions, narrow histograms
  double par1[nfits][3] = {
    { -0.0131841, 0.0254668, -0.211042 },
    { 0.000574698, 0.00997693, -1.52906 },
    { -7.12483e-05, 0.00387567, -1.75892 },
    { -3.19278e-05, 0.00761138, -0.927624 }
  };

  const char *label1 = "New Fits";
  const char *label2 = "Old Fits";
  const char *label3 = "data set 3";
  double pmin = 1, pmax = 10;   // range of p to draw
  // -----------------------------------------------------------------------------

  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("cfits", "fits", 1200, 900);
  c->Divide(2, 2);

  for( int i=0; i<nfits; i++ ){
    c->cd(i+1);

    TF1 *f1 = new TF1( Form("f1_%d", i), "[0] + [1]*pow(x,[2])", pmin, pmax );
    f1->SetParameters( par1[i][0], par1[i][1], par1[i][2] );
    f1->SetLineColor(kBlack);
    f1->SetLineWidth(2);
    f1->SetTitle( titles[i] );

    TF1 *f2 = 0;
    if( draw_second ){
      f2 = new TF1( Form("f2_%d", i), "[0] + [1]*pow(x,[2])", pmin, pmax );
      f2->SetParameters( par2[i][0], par2[i][1], par2[i][2] );
      f2->SetLineColor(kRed);
      f2->SetLineStyle(2);
      f2->SetLineWidth(2);

      // Use a common vertical range so neither curve is clipped
      double lo = std::min( f1->GetMinimum(pmin, pmax), f2->GetMinimum(pmin, pmax) );
      double hi = std::max( f1->GetMaximum(pmin, pmax), f2->GetMaximum(pmin, pmax) );
      double pad = 0.1*(hi - lo);
      f1->SetMinimum( lo - pad );
      f1->SetMaximum( hi + pad );
    }

    f1->Draw();
    if( f2 ) f2->Draw("same");

    TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->AddEntry( f1, label1, "l" );
    if( f2 ) leg->AddEntry( f2, label2, "l" );
    leg->Draw();
  }

  c->SaveAs("fits_comparison.png");
}
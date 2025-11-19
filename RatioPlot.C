// ------------------------------------------------------------------------ //
// Script to proton electromagnetic form factor ratios for final            //
// presentation. It includes data points from various sources.              //
//                                                                          //
// This script should eventually be modified to accept input data files for //
// calculated values and errors from GEp-V analysis.                        //
//                                                                          //
// ---------                                                                //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 11-19-2025                  //
// ---------                                                                //
// ** Do not tamper with this sticker! Log any updates to the script above. //
// ------------------------------------------------------------------------ //
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TText.h"
#include "TLine.h"


void RatioPlot() {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // GEp-V Results
    const int nGEP5 = 2;
    double GEP5R[nGEP5]  = {0, 0};       // Ratio values
    double GEP5Q[nGEP5]  = {5.732, 11.109};  // Q2 values
    double GEP5Rstat[nGEP5] = {0, 0};   // statistical errors
    double GEP5Rsyst[nGEP5] = {0, 0};   // systematic errors
    double GEP5err[nGEP5]; // Combine errors in quadrature
    for(int i = 0; i < nGEP5; i++) {
        GEP5err[i] = sqrt(GEP5Rstat[i]*GEP5Rstat[i] +  GEP5Rsyst[i]*GEP5Rsyst[i]);
    }
    TGraphErrors *GEP5 = new TGraphErrors(nGEP5, GEP5Q, GEP5R, 0, GEP5err);
    GEP5->SetTitle("Form Factor Ratio Results;Q^{2} (GeV^{2});#mu_{p}G_{E}^{P}/G_{M}^{P}");
    GEP5->GetXaxis()->SetLimits(0, 12.5);
    GEP5->GetYaxis()->SetRangeUser(-0.5, 1.25);
    GEP5->SetMarkerStyle(20);
    GEP5->SetMarkerSize(1.2);
    GEP5->SetMarkerColor(kGreen);
    GEP5->SetLineColor(kGreen);

    // GEp-I Results (Punjabi et al, PRC 2005)
    const int nGEP1 = 10;
    double GEP1R[nGEP1]  = {0.979, 0.951, 0.883, 0.798, 0.789, 0.777, 0.747, 0.703, 0.615, 0.606};       // Ratio values
    double GEP1Q[nGEP1]  = {0.49, 0.79, 1.18, 1.48, 1.77, 1.88, 2.13, 2.47, 2.97, 3.47};   // Q2 values
    double GEP1Rstat[nGEP1] = {0.016, 0.012, 0.013, 0.029, 0.024, 0.024, 0.032, 0.023, 0.029, 0.042};   // statistical errors
    double GEP1Rsyst[nGEP1] = {0.006, 0.010, 0.018, 0.026, 0.035, 0.033, 0.034, 0.033, 0.021, 0.014};   // systematic errors
    double GEP1err[nGEP1]; // Combine errors in quadrature
    for(int i = 0; i < nGEP1; i++) {
        GEP1err[i] = sqrt(GEP1Rstat[i]*GEP1Rstat[i] +  GEP1Rsyst[i]*GEP1Rsyst[i]);
    }
    TGraphErrors *GEP1 = new TGraphErrors(nGEP1, GEP1Q, GEP1R, 0, GEP1err);
    GEP1->SetMarkerStyle(20);
    GEP1->SetMarkerSize(1.2);
    GEP1->SetMarkerColor(kBlue);
    GEP1->SetLineColor(kBlue);

     // GEp-II Results (Puckett et al, PRC 2012)
    const int nGEP2 = 4;
    double GEP2R[nGEP2]  = {0.571, 0.517, 0.450, 0.354};       // Ratio values
    double GEP2Q[nGEP2]  = {3.5, 4.0, 4.8, 5.6};  // Q2 values
    double GEP2Rstat[nGEP2] = {0.072, 0.055, 0.052, 0.085};   // statistical errors
    double GEP2Rsyst[nGEP2] = {0.007, 0.009, 0.012, 0.019};   // systematic errors
    double GEP2err[nGEP2]; // Combine errors in quadrature
    for(int i = 0; i < nGEP2; i++) {
        GEP2err[i] = sqrt(GEP2Rstat[i]*GEP2Rstat[i] +  GEP2Rsyst[i]*GEP2Rsyst[i]);
    }
    TGraphErrors *GEP2 = new TGraphErrors(nGEP2, GEP2Q, GEP2R, 0, GEP2err);
    GEP2->SetMarkerStyle(21);
    GEP2->SetMarkerSize(1.2);
    GEP2->SetMarkerColor(kRed);
    GEP2->SetLineColor(kRed);


    // GEp-III Results (Puckett al, PRC 2017)
    const int nGEP3 = 3;
    double GEP3R[nGEP3]  = {0.448, 0.348, 0.145};       // Ratio values
    double GEP3Q[nGEP3]  = {5.2, 6.8, 8.537};  // Q2 values
    double GEP3Rstat[nGEP3] = {0.06, 0.105, 0.175};   // statistical errors
    double GEP3Rsyst[nGEP3] = {0.006, 0.01, 0.024};   // systematic errors
    double GEP3err[nGEP3]; // Combine errors in quadrature
    for(int i = 0; i < nGEP3; i++) {
        GEP3err[i] = sqrt(GEP3Rstat[i]*GEP3Rstat[i] +  GEP3Rsyst[i]*GEP3Rsyst[i]);
    }
    TGraphErrors *GEP3 = new TGraphErrors(nGEP3, GEP3Q, GEP3R, 0, GEP3err);
    GEP3->SetMarkerStyle(22);
    GEP3->SetMarkerSize(1.2);
    GEP3->SetMarkerColor(kBlack);
    GEP3->SetLineColor(kBlack);

    // Create TLatex object
    
    TLatex *latex = new TLatex();
    latex->SetNDC();                 // coordinates from 0–1
    latex->SetTextSize(0.175);         // big font
    latex->SetTextFont(62);           // bold Helvetica
    latex->SetTextAlign(22);          // center text at the coordinates
    latex->SetTextAngle(45);          // rotate 45 degrees
    latex->SetTextColorAlpha(kBlue, 0.3);



    // Draw everything
    TCanvas *c = new TCanvas("c", "Scatter with Errors", 800, 600);
    c->SetLeftMargin(0.15);    // Adjust spacing for axis labels
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    GEP5->GetXaxis()->SetTitleSize(0.05);    // Adjust title and label sizes for better visibility
    GEP5->GetYaxis()->SetTitleSize(0.05);
    GEP5->GetXaxis()->SetLabelSize(0.04);    // size of the tick labels
    GEP5->GetYaxis()->SetLabelSize(0.04);
    GEP5->GetYaxis()->SetNdivisions(505);    // Adjust tick divisions for cleaner look
    GEP5->GetXaxis()->SetNdivisions(505);
    GEP5->GetXaxis()->CenterTitle();         // Center the axis title
    GEP5->GetYaxis()->CenterTitle();

    GEP5->Draw("AP");   // Draw everything to start

    TLine *xaxis = new TLine(0, 0, 12.5, 0);  // Add an x-axis at y=0 for reference
    xaxis->SetLineColor(kBlack);
    xaxis->SetLineWidth(1);
    xaxis->Draw();

    // Begin Plotting different datasets
    TLegend *leg = new TLegend(0.15, 0.15, 0.4, 0.35);  // x1,y1,x2,y2 in NDC
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    GEP1->Draw("P SAME");
    leg->AddEntry(GEP1, "GEp-I", "p");
    GEP2->Draw("P SAME");
    leg->AddEntry(GEP2, "GEp-II", "p");
    GEP3->Draw("P SAME");
    leg->AddEntry(GEP3, "GEp-III", "p");
    GEP5->Draw("P SAME");
    leg->AddEntry(GEP5, "GEp-V", "p");
    leg->Draw();


    latex->DrawLatex(0.55, 0.525, "Preliminary");  // Optional "Preliminary" stamp
    c->SaveAs("output/RatioResults.pdf");  // Optional save as PDF
}
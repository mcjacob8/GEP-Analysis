// ------------------------------------------------------------------------ //
// This script is for comparing replayed simultation results with either    //
// other simulation or with real data. Orginally designed to look at timing //
//                                                                          //
// ---------                                                                //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 6-23-2026                   //
// ---------                                                                //
// ** Do not tamper with this sticker! Log any updates to the script above. //
// ------------------------------------------------------------------------ //

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TKey.h"
#include "TClass.h"
#include "TLegend.h"
#include "TString.h"
#include "TDirectory.h"

using namespace std;



// === Main comparison macro ===
void SimuCompare(const char *root1 = "/volatile/halla/sbs/mcjacob/GEP/Kin3/Simu/May27/replayed-GEP3mod_preinit_all.root", const char *root2 = "/volatile/halla/sbs/mcjacob/GEP/Kin3/Simu/Jun22/replayed_GEP3mod_preinit_all.root", 
                    TString outpdf = "../output/SimuCompare.pdf") {

    TFile *f1 = TFile::Open(root1);
    TFile *f2 = TFile::Open(root2);

    if(!f1 || !f2){
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    // Histograms to compare
    std::vector<TString> histNames = {
        "hsbs_gemFT_clust_Utime",
        "hsbs_gemFT_clust_Vtime"
    };


    for(int m = 0; m < 14; m++) {
        histNames.push_back(Form("hsbs_gemFT_m%d_TimeU_all",  m));
        histNames.push_back(Form("hsbs_gemFT_m%d_TimeU_good", m));
        histNames.push_back(Form("hsbs_gemFT_m%d_TimeV_all",  m));
        histNames.push_back(Form("hsbs_gemFT_m%d_TimeV_good", m));
    }

    TCanvas *c = new TCanvas("c","c",1600,1200);
    gStyle->SetOptStat(0);
    //c->Print(outpdf + "("); //Open PDF

    c->Clear();
    c->Divide(2,1);

    {
        c->cd(1);
        TH1 *h1 = (TH1*)f1->Get("hsbs_gemFT_clust_Utime");
        TH1 *h2 = (TH1*)f2->Get("hsbs_gemFT_clust_Utime");

        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);

        if(h1->Integral()>0) h1->Scale(1./h1->Integral());
        if(h2->Integral()>0) h2->Scale(1./h2->Integral());

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
        leg->SetTextSize(0.03);
        leg->AddEntry(h1, root1, "l");
        leg->AddEntry(h2, root2, "l");
        leg->Draw();

        c->cd(2);
        h1 = (TH1*)f1->Get("hsbs_gemFT_clust_Vtime");
        h2 = (TH1*)f2->Get("hsbs_gemFT_clust_Vtime");

        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);

        if(h1->Integral()>0) h1->Scale(1./h1->Integral());
        if(h2->Integral()>0) h2->Scale(1./h2->Integral());

        h1->Draw("hist");
        h2->Draw("hist same");

        TLegend *leg2 = new TLegend(0.55, 0.70, 0.88, 0.88);
        leg2->SetTextSize(0.03);
        leg2->AddEntry(h1, root1, "l");
        leg2->AddEntry(h2, root2, "l");
        leg2->Draw();
    }

    //c->Print(outpdf);
    c->Print(outpdf + "(");

    std::vector<TString> categories = {
        "TimeU_all",
        "TimeU_good",
        "TimeV_all",
        "TimeV_good"
    };

    for(const auto& cat : categories){

        c->Clear();
        c->Divide(4,4);

        for(int m=0; m<14; m++){

            c->cd(m+1);

            TString hname =
                Form("hsbs_gemFT_m%d_%s", m, cat.Data());

            TH1 *h1 = (TH1*)f1->Get(hname);
            TH1 *h2 = (TH1*)f2->Get(hname);

            if(!h1 || !h2) continue;

            h1->SetLineColor(kBlue);
            h2->SetLineColor(kRed);

            h1->SetLineWidth(2);
            h2->SetLineWidth(2);

            if(h1->Integral()>0) h1->Scale(1./h1->Integral());
            if(h2->Integral()>0) h2->Scale(1./h2->Integral());

            double maxy =
                std::max(h1->GetMaximum(),
                         h2->GetMaximum());

            h1->SetMaximum(1.1*maxy);

            h1->SetTitle(Form("Module %d: %s", m, cat.Data()));

            h1->Draw("hist");
            h2->Draw("hist same");
        }

        c->Print(outpdf);
    }

    c->Clear();
    c->Print(outpdf + ")");
    f1->Close();
    f2->Close();
}


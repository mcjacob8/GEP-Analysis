// ------------------------------------------------------------------------ //
// This script is for comparing the results of the GetTrackingCutsFast      //
// script, which produced many DB paramaters for the GEMs for ADC threshold //
// and timing cuts. It is useful for comparing how these paramaters might   //
// change between different run conditions (e.g. high voltages, latency...) //
//                                                                          //
// ---------                                                                //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 11-11-2025                  //
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

// Function to parse a .dat file
map<string, vector<double>> ReadParamFile(const string &filename) {
    map<string, vector<double>> data;
    ifstream fin(filename);
    if (!fin) {
        cerr << "Error opening " << filename << endl;
        return data;
    }

    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue; // skip comments

        string key, eq;
        stringstream ss(line);
        ss >> key >> eq;
        if (eq != "=") continue;

        vector<double> vals;
        double val;
        while (ss >> val) vals.push_back(val);

        data[key] = vals;
    }
    return data;
}

// Helper: collect all TH1 histograms recursively from a ROOT file
void CollectHistograms(TDirectory *dir, map<string, TH1*> &histMap, const string &prefix = "") {
    TIter nextkey(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TH1")) {
            string name = prefix + obj->GetName();
            histMap[name] = (TH1*)obj;
        } else if (obj->InheritsFrom("TDirectory")) {
            CollectHistograms((TDirectory*)obj, histMap, prefix + obj->GetName() + "/");
        }
    }
}

// Fit a Gaussian to the histogram using FWHM method (taken from GetTrackingCutsFast_GEp_AJRP.C)
void FitGaus_FWHM( TH1 *htest, double thresh=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh+1);

  htest->Fit("gaus","qS","",xlow, xhigh);
}

// === Main comparison macro ===
void TrackingCompare(const char *root1 = "output/GEPtrackingcuts.root", const char *root2 = "output/GEPtrackingcuts_optics.root", 
                    const char *FT_file1 = "output/GEPtrackingcuts_FT.dat", const char *FT_file2 = "output/GEPtrackingcuts_optics_FT.dat", 
                    const char *FPP_file1 = "output/GEPtrackingcuts_FPP.dat", const char *FPP_file2 = "output/GEPtrackingcuts_optics_FPP.dat",
                    const char *outpdf = "output/TrackingCompare.pdf") {
    // Read in all of our input files   
    TFile *f1 = TFile::Open(root1);
    TFile *f2 = TFile::Open(root2);
    if (!f1 || !f2) {
        cerr << "Error: Cannot open one of the ROOT files." << endl;
        return;
    }
    map<string, TH1*> hists1, hists2;
    CollectHistograms(f1, hists1);
    CollectHistograms(f2, hists2);
    auto data1FT = ReadParamFile(FT_file1);
    auto data2FT = ReadParamFile(FT_file2);
    auto data1FPP = ReadParamFile(FPP_file1);
    auto data2FPP = ReadParamFile(FPP_file2);

    cout << "Loaded " << data1FT.size() << " parameters from " << FT_file1 << endl;
    cout << "Loaded " << data2FT.size() << " parameters from " << FT_file2 << endl;
    cout << "Loaded " << data1FPP.size() << " parameters from " << FPP_file1 << endl;
    cout << "Loaded " << data2FPP.size() << " parameters from " << FPP_file2 << endl;
    cout << "Loaded " << hists1.size() << " histograms from " << root1 << endl;
    cout << "Loaded " << hists2.size() << " histograms from " << root2 << endl;

    // --- Define histogram names to compare ---
    vector<string> histNames = {
        "hADCmaxsampU", "hADCmaxsampV", "hADCmaxstripU", "hADCmaxstripV",
        "hADCU", "hADCV", "hDeconvADCUmaxstrip", "hDeconvADCVmaxstrip",
        "hDeconvADCUclust", "hDeconvADCVclust", "hADCasym", "hADCasymDeconv",
        "hUtimeMaxStrip", "hVtimeMaxStrip", "hUtime", "hVtime",
        "hUtimeMaxStripDeconv", "hVtimeMaxStripDeconv", "hUtimeDeconv", "hVtimeDeconv",
        "hTavg_corr", "hdeltat", "hdeltatDeconv" // fill all 23
    };

    // --- Define modules to compare ---
    //vector<string> modules = {"_FTmod0","_FTmod1"};
    //vector<string> modules = {"_FTmod0","_FTmod1","_FTmod2","_FTmod3","_FTmod4","_FTmod5","_FTmod6","_FTmod7","_FTmod8","_FTmod9"
    //      ,"_FTmod10","_FTmod11","_FTmod12","_FTmod13"}; // add FT all modules
    //vector<string> modules = {"_FPPmod0","_FPPmod1","_FPPmod2","_FPPmod3","_FPPmod4","_FPPmod5","_FPPmod6","_FPPmod7","_FPPmod8","_FPPmod9"
    //      ,"_FPPmod10","_FPPmod11","_FPPmod12","_FPPmod13","_FPPmod14","_FPPmod15","_FPPmod16","_FPPmod17","_FPPmod18","_FPPmod19"
    //      ,"_FPPmod20","_FPPmod21","_FPPmod22","_FPPmod23","_FPPmod24","_FPPmod25","_FPPmod26","_FPPmod27","_FPPmod28","_FPPmod29"
    //      ,"_FPPmod30","_FPPmod31"}; // add FPP all modules
    vector<string> modules = {"_FTmod0","_FTmod1","_FTmod2","_FTmod3","_FTmod4","_FTmod5","_FTmod6","_FTmod7","_FTmod8","_FTmod9"
            ,"_FTmod10","_FTmod11","_FTmod12","_FTmod13", "_FPPmod0","_FPPmod1","_FPPmod2","_FPPmod3","_FPPmod4","_FPPmod5","_FPPmod6","_FPPmod7","_FPPmod8","_FPPmod9"
            ,"_FPPmod10","_FPPmod11","_FPPmod12","_FPPmod13","_FPPmod14","_FPPmod15","_FPPmod16","_FPPmod17","_FPPmod18","_FPPmod19"
            ,"_FPPmod20","_FPPmod21","_FPPmod22","_FPPmod23","_FPPmod24","_FPPmod25","_FPPmod26","_FPPmod27","_FPPmod28","_FPPmod29"
            ,"_FPPmod30","_FPPmod31"}; // add FT and FPP all modules


    // --- Create canvas for overlay plotting ---
    TCanvas *c = new TCanvas("c", "Histogram Comparison", 800, 600);
    c->Print(Form("%s[", outpdf)); // begin multi-page PDF

    int plotsPerPage = 12; // number of pads per page (to match asthetic output of GetTrackingCutsFast)
    int nHist = histNames.size();

    // Loop over histograms to generate overlaid plots
    for (const string &mod : modules) {
        for (int i = 0; i < nHist; ++i) {
            const string &name = histNames[i];
            string fullName = name + mod;

            auto it1 = hists1.find(fullName);
            auto it2 = hists2.find(fullName);

            if (it1 == hists1.end() || it2 == hists2.end()) {
                cout << "  Warning: " << fullName << " not found in one of the files!" << endl;
                continue;
            }

            TH1 *h1 = (TH1*)it1->second->Clone();
            TH1 *h2 = (TH1*)it2->second->Clone();

            // --- Normalize histogram and record scale factor ---
            double scale1 = (h1->Integral() != 0) ? 1.0 / h1->Integral() : 1.0;
            double scale2 = (h2->Integral() != 0) ? 1.0 / h2->Integral() : 1.0;

            h1->Scale(scale1);
            h2->Scale(scale2);

            h1->SetLineColor(kBlue);
            h2->SetLineColor(kRed);
            h1->SetLineWidth(2);
            h2->SetLineWidth(2);

            double maxY = max(h1->GetMaximum(), h2->GetMaximum());
            h1->SetMaximum(1.1 * maxY);
            h1->SetTitle((fullName + ";X;Normalized entries").c_str());

            // --- Draw in grid ---
            if (i % plotsPerPage == 0) {
                c->Clear();
                c->Divide(4, 3); // 4x3 grid = 12 plots per page
            }

            c->cd((i % plotsPerPage) + 1);

            // Only attempt fit if both histograms have enough events
            const double minEntries = 300.0;
            if (h1->GetEntries() > minEntries && h2->GetEntries() > minEntries) {
                TF1 *fit1 = nullptr;
                TF1 *fit2 = nullptr;

                // Determine if fitting gaussian or landau based on histogram index
                if (i == 0 || i == 1 || i == 6 || i == 7) {
                    h1->Fit("landau","qS","",0,1000.);
                    h2->Fit("landau","qS","",0,1000.);
                    fit1 = h1->GetFunction("landau");
                    fit2 = h2->GetFunction("landau");
                } else if (i == 2 || i == 3 || i == 8 || i == 9) {
                    h1->Fit("landau","qS","",0,3000.);
                    h2->Fit("landau","qS","",0,3000.);
                    fit1 = h1->GetFunction("landau");
                    fit2 = h2->GetFunction("landau");
                } else if (i == 4 || i == 5) {
                    h1->Fit("landau","qS","",0,7500.);
                    h2->Fit("landau","qS","",0,7500.);
                    fit1 = h1->GetFunction("landau");
                    fit2 = h2->GetFunction("landau");
                } else{
                    FitGaus_FWHM(h1, 0.3);
                    FitGaus_FWHM(h2, 0.3);
                    fit1 = h1->GetFunction("gaus");
                    fit2 = h2->GetFunction("gaus");
                }

                // Customize fits
                if (fit1 && fit2) {
                    fit1->SetLineColor(kBlack + 2);
                    fit1->SetLineStyle(2);
                    fit1->SetLineWidth(2);
                    fit2->SetLineColor(kBlack);
                    fit2->SetLineStyle(2);
                    fit2->SetLineWidth(2);
                }

                h1->Draw("HIST");
                h2->Draw("HIST SAME");
    
                if (fit1) fit1->Draw("SAME");
                if (fit2) fit2->Draw("SAME");

                TLegend *leg = new TLegend(0.75, 0.550, 0.98, 0.75);
                leg->AddEntry(h1, "LH2", "l");
                leg->AddEntry(h2, "Optics", "l");
                if (fit1) leg->AddEntry(fit1, "LH2 fit", "l");
                if (fit2) leg->AddEntry(fit2, "Optics fit", "l");
                leg->SetTextSize(0.04);
                leg->Draw();

            } else {
                // If not enough entries, just draw histograms without fits
                h1->Draw("HIST");
                h2->Draw("HIST SAME");

                TLegend *leg = new TLegend(0.75, 0.550, 0.98, 0.75);
                leg->AddEntry(h1, "LH2", "l");
                leg->AddEntry(h2, "Optics", "l");
                leg->SetTextSize(0.04);
                leg->Draw();
            }

            if ((i + 1) % plotsPerPage == 0 || i == nHist - 1) {
                c->Print(Form("%s", outpdf)); // save the full page
            }
        }
    }   

    //c->Print("TrackingCompare.pdf]"); // OPTIONAL end multi-page PDF (OR add in later plots)

    // --- Now compare parameter .dat files ---
    // --- Define difference histograms ---
    map<string, TH1D*> histMap;
    histMap["maxstrip_t0"]              = new TH1D("hDiff_maxstrip_t0", "maxstrip_t0;#Delta value;Count", 100, -10, 10);
    histMap["maxstrip_tsigma"]          = new TH1D("hDiff_maxstrip_tsigma", "maxstrip_tsigma;#Delta value;Count", 100, -6, 6);
    histMap["maxstrip_t0_deconv"]       = new TH1D("hDiff_maxstrip_t0_deconv", "maxstrip_t0_deconv;#Delta value;Count", 100, -10, 10);
    histMap["maxstrip_tsigma_deconv"]   = new TH1D("hDiff_maxstrip_tsigma_deconv", "maxstrip_tsigma_deconv;#Delta value;Count", 100, -10, 10);
    histMap["HitTimeMean"]              = new TH1D("hDiff_HitTimeMean", "HitTimeMean;#Delta value;Count", 50, -7, 5);
    histMap["HitTimeSigma"]             = new TH1D("hDiff_HitTimeSigma", "HitTimeSigma;#Delta value;Count", 100, -10, 10);
    histMap["HitTimeMeanDeconv"]        = new TH1D("hDiff_HitTimeMeanDeconv", "HitTimeMeanDeconv;#Delta value;Count", 100, -10, 10);
    histMap["HitTimeSigmaDeconv"]       = new TH1D("hDiff_HitTimeSigmaDeconv", "HitTimeSigmaDeconv;#Delta value;Count", 100, -10, 10);
    histMap["deltat_sigma"]             = new TH1D("hDiff_deltat_sigma", "deltat_sigma;#Delta value;Count", 100, -10, 4);
    histMap["deltat_sigma_deconv"]      = new TH1D("hDiff_deltat_sigma_deconv", "deltat_sigma_deconv;#Delta value;Count", 100, -15, 5);
    histMap["sigma_tcorr"]              = new TH1D("hDiff_sigma_tcorr", "sigma_tcorr;#Delta value;Count", 50, -3, 3);
    histMap["threshold_sample"]         = new TH1D("hDiff_threshold_sample", "threshold_sample;#Delta value;Count", 100, -100, 100);
    histMap["threshold_stripsum"]       = new TH1D("hDiff_threshold_stripsum", "threshold_stripsum;#Delta value;Count", 100, -100, 100);
    histMap["threshold_clustersum"]     = new TH1D("hDiff_threshold_clustersum", "threshold_clustersum;#Delta value;Count", 100, -100, 100);
    histMap["threshold_maxcombo_deconv"]= new TH1D("hDiff_threshold_maxcombo_deconv", "threshold_maxcombo_deconv;#Delta value;Count", 100, -10, 10);
    histMap["ADCasym_sigma"]            = new TH1D("hDiff_ADCasym_sigma", "ADCasym_sigma;#Delta value;Count", 50, -0.5, 0.5);

    int nCompared = 0;
    map<string, vector<double>> data1 = data1FT;       // Make a copy of the maps and append FT and FPP results together
    map<string, vector<double>> data2 = data2FT;
    data1.insert(data1FPP.begin(), data1FPP.end());    // Can comment these lines out and edit previous ones to look at only FT or FPP
    data2.insert(data2FPP.begin(), data2FPP.end());


    // Loop through parameters
    for (const auto &entry : data1) {
        nCompared++;
        const string &key = entry.first;
        const vector<double> &vals1 = entry.second;

        auto it2 = data2.find(key);
        if (it2 == data2.end()) continue;

        const vector<double> &vals2 = it2->second;

        // Extract base parameter name (everything after last '.')
        string paramName = key.substr(key.find_last_of('.') + 1);
        // cout << "Comparing parameter: " << paramName << endl;

        // Skip parameters we don’t have a histogram for
        auto histIt = histMap.find(paramName);
        if (histIt == histMap.end()) continue;
        TH1D *hist = histIt->second;

        // Fill with differences
        for (size_t i = 0; i < vals1.size() && i < vals2.size(); i++) {
            double diff = vals2[i] - vals1[i];
            hist->Fill(diff);
        }
    }

    // Draw results on new canvas
    TCanvas *c2 = new TCanvas("cDiffs", "Parameter Differences", 1200, 800);
    c2->Divide(4,4); // Total of 16 parameters being plotted

    int pad = 1;
    for (auto &pair : histMap) {
        c2->cd(pad++);
        pair.second->Draw();
        if (pad > 16) {
            c2->Print(Form("%s", outpdf));  // save these plots to the pdf as well
            break;
        }
    }

    c2->Print(Form("%s]", outpdf)); // end multi-page PDF

    
    cout << "\nCompared " << nCompared << " common parameters." << endl;

    // Clean up
    // c->Close();
    // c2->Close();
    f1->Close();
    f2->Close();
}
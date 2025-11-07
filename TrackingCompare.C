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
#include "TKey.h"
#include "TClass.h"
#include "TLegend.h"
#include "TString.h"
#include "TDirectory.h"

using namespace std;


// Helper: parse a line like "sbs.gemFT.m11.HitTimeSigmaDeconv = 21.0469 19.7812"
bool ParseLine(const string &line, string &key, vector<double> &values) {
    if (line.empty() || line[0] == '#') return false;

    size_t eqPos = line.find('=');
    if (eqPos == string::npos) return false;

    key = line.substr(0, eqPos);
    // trim whitespace from key
    key.erase(0, key.find_first_not_of(" \t"));
    key.erase(key.find_last_not_of(" \t") + 1);

    string valuePart = line.substr(eqPos + 1);
    stringstream ss(valuePart);
    double val;
    values.clear();
    while (ss >> val) values.push_back(val);

    return !values.empty();
}

// Load a .dat file into a map<string, vector<double>>
map<string, vector<double>> LoadDatFile(const string &filename) {
    map<string, vector<double>> data;
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return data;
    }

    string line, key;
    vector<double> values;
    while (getline(infile, line)) {
        if (ParseLine(line, key, values)) {
            data[key] = values;
        }
    }
    infile.close();
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

// === Main comparison macro ===
void TrackingCompare(const char *root1 = "output/GEPtrackingcuts.root", const char *root2 = "output/GEPtrackingcuts_optics.root", const char *FT_file1 = "output/GEPtrackingcuts_FT.dat", const char *FT_file2 = "output/GEPtrackingcuts_FT.dat", const char *FPP_file1 = "output/GEPtrackingcuts_FPP.dat", const char *FPP_file2 = "output/GEPtrackingcuts_FPP.dat") {
    TFile *f1 = TFile::Open(root1);
    TFile *f2 = TFile::Open(root2);
    if (!f1 || !f2) {
        cerr << "Error: Cannot open one of the ROOT files." << endl;
        return;
    }
    map<string, TH1*> hists1, hists2;
    CollectHistograms(f1, hists1);
    CollectHistograms(f2, hists2);

    auto data1 = LoadDatFile(FT_file1);
    auto data2 = LoadDatFile(FT_file2);
    auto data3 = LoadDatFile(FPP_file1);
    auto data4 = LoadDatFile(FPP_file2);

    cout << "Loaded " << data1.size() << " parameters from " << FT_file1 << endl;
    cout << "Loaded " << data2.size() << " parameters from " << FT_file2 << endl;
    cout << "Loaded " << data3.size() << " parameters from " << FPP_file1 << endl;
    cout << "Loaded " << data4.size() << " parameters from " << FPP_file2 << endl;
    cout << "Loaded " << hists1.size() << " histograms from " << root1 << endl;
    cout << "Loaded " << hists2.size() << " histograms from " << root2 << endl;

    vector<string> histNames = {
        "hADCmaxsampU", "hADCmaxsampV", "hADCmaxstripU", "hADCmaxstripV",
        "hADCU", "hADCV", "hDeconvADCUmaxstrip", "hDeconvADCVmaxstrip",
        "hDeconvADCUclust", "hDeconvADCVclust", "hADCasym", "hADCasymDeconv",
        "hUtimeMaxStrip", "hVtimeMaxStrip", "hUtime", "hVtime",
        "hUtimeMaxStripDeconv", "hVtimeMaxStripDeconv", "hUtimeDeconv", "hVtimeDeconv",
        "hTavg_corr", "hdeltat", "hdeltatDeconv" // fill all 23
    };

    // --- Define modules
    vector<string> modules = {"_FTmod0","_FTmod1"}; // add all modules


    TCanvas *c = new TCanvas("c", "Histogram Comparison", 800, 600);
    c->Print("TrackingCompare.pdf["); // begin multi-page PDF

    int nCompared = 0;

    int plotsPerPage = 12; // number of pads per page
    int nHist = histNames.size();

    
    for (const string &mod : modules) {
        //const string &name = histNames[i];
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

            if (h1->Integral() != 0) h1->Scale(1.0 / h1->Integral());
            if (h2->Integral() != 0) h2->Scale(1.0 / h2->Integral());

            h1->SetLineColor(kBlue);
            h2->SetLineColor(kRed);
            h1->SetLineWidth(2);
            h2->SetLineWidth(2);

            double maxY = max(h1->GetMaximum(), h2->GetMaximum());
            h1->SetMaximum(1.1 * maxY);
            h1->SetTitle((fullName + ";X;Normalized entries").c_str());

            // --- Draw in grid
            if (i % plotsPerPage == 0) {
                c->Clear();
                c->Divide(4, 3); // 4x3 grid = 12 plots per page
            }

            c->cd((i % plotsPerPage) + 1);
            h1->Draw("HIST");
            h2->Draw("HIST SAME");

            // --- Overlay fit functions if they exist ---
            TF1 *fit1 = nullptr;
            TF1 *fit2 = nullptr;

            if (h1->GetListOfFunctions() && h1->GetListOfFunctions()->GetSize() > 0)
                fit1 = (TF1*)h1->GetListOfFunctions()->First();

            if (h2->GetListOfFunctions() && h2->GetListOfFunctions()->GetSize() > 0)
                fit2 = (TF1*)h2->GetListOfFunctions()->First();

            if (fit1) {
                fit1->SetLineColor(kBlue + 2);
                fit1->SetLineStyle(2);
                fit1->SetLineWidth(100);
                fit1->Draw("SAME");
            }
            if (fit2) {
                fit2->SetLineColor(kRed + 2);
                fit2->SetLineStyle(2);
                fit2->SetLineWidth(100);
                fit2->Draw("SAME");
            }

            TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
            leg->AddEntry(h1, root1, "l");
            leg->AddEntry(h2, root2, "l");
            if (fit1) leg->AddEntry(fit1, (string(root1) + " fit").c_str(), "l");
            if (fit2) leg->AddEntry(fit2, (string(root2) + " fit").c_str(), "l");
            leg->SetTextSize(0.04);
            leg->Draw();

            if ((i + 1) % plotsPerPage == 0 || i == nHist - 1) {
                c->Print("TrackingCompare.pdf"); // save the full page
            }
        }
    }   

    c->Print("TrackingCompare.pdf]"); // end multi-page PDF

    cout << "Compared and plotted " << nCompared << " matching histograms." << endl;

    f1->Close();
    f2->Close();

    // Histograms for value differences (optional)
    TH1D *hDiff = new TH1D("hDiff", "Parameter Value Differences;#Delta value;Count", 100, -10, 10);

    nCompared = 0;

    // Loop over all parameters in file1
    for (const auto &entry : data1) {
        const string &param = entry.first;
        const vector<double> &vals1 = entry.second;

        auto it2 = data2.find(param);
        if (it2 == data2.end()) {
            cout << "Warning: " << param << " not found in " << FT_file2 << endl;
            continue;
        }

        const vector<double> &vals2 = it2->second;

        cout << param << endl;
        for (size_t i = 0; i < vals1.size() && i < vals2.size(); i++) {
            double diff = vals2[i] - vals1[i];
            cout << "   value[" << i << "]: " << vals1[i] << " -> " << vals2[i]
                 << " (diff = " << diff << ")" << endl;
            hDiff->Fill(diff);
        }
        nCompared++;
    }

    cout << "\nCompared " << nCompared << " common parameters." << endl;

    // Plot difference histogram
    //TCanvas *c = new TCanvas("c", "Tracking Parameter Comparison", 800, 600);
    //hDiff->Draw();
}
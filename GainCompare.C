// A simple script to compare gain coefficents between different calibrations
// Input 1 = new, Input 2 = old
// Usage: GainCompare("output/GEMFPP_gainmatch_temp.dat", "output/GEMFPP_gainmatch_pass1.dat");
void GainCompare(const char* gainfile1 = "output/GEMFPP_gainmatch_temp.dat", const char* gainfile2 = "output/GEMFPP_gainmatch_pass1.dat") {
  
  // Helper function to read gain file into a map
  auto readGainFile = [](const char* filename) {
    std::map<TString, std::vector<double>> gains;
    std::ifstream fin(filename);
    if (!fin.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return gains;
    }

    TString line;
    while (line.ReadLine(fin)) {
      if (line.BeginsWith("#") || line.IsWhitespace()) continue; // skip comments or blanks
      TObjArray* tokens = line.Tokenize(" ");
      if (tokens->GetEntries() < 3) continue;

      TString label = ((TObjString*)tokens->At(0))->GetString(); // e.g. sbs.gemFPP.m0.vgain
      TString eq = ((TObjString*)tokens->At(1))->GetString();    // should be "="

      std::vector<double> values;
      for (int i = 2; i < tokens->GetEntries(); i++) {
        values.push_back(((TObjString*)tokens->At(i))->GetString().Atof());
      }
      gains[label] = values;
      delete tokens;
    }
    return gains;
  };

  auto gain1 = readGainFile(gainfile1);
  auto gain2 = readGainFile(gainfile2);

  // Plot comparison
  gStyle->SetOptStat(1110);
  TCanvas* c = new TCanvas("c", "Gain Comparison", 1000, 800);
  c->Divide(2, 1);

  TGraph* g1 = new TGraph;
  TGraph* g2 = new TGraph;

  TH1D * hGain = new TH1D("hGain", "New Gain/Old Gain; nG/oG;", 100, 0.5, 1.7);

  int point = 0;
  for (auto& [label, vals1] : gain1) {
    // Check if the second file has the same key
    if (gain2.find(label) == gain2.end()) {
      std::cout << "Warning: " << label << " not found in " << gainfile2 << std::endl;
      continue;
    }
    auto vals2 = gain2[label];

    int n = vals1.size();
    if ((int)vals2.size() != n) {
      std::cout << "Warning: size mismatch for " << label << std::endl;
      n = std::min(n, (int)vals2.size());
    }

    for (int i = 0; i < n; i++) {
      g1->SetPoint(point, point+1, vals1[i]);
      g2->SetPoint(point, point+1, vals2[i]);
      hGain->Fill( vals1[i] / vals2[i]);
      point++;
    }

  }

  c->cd(1);
  g1->SetTitle("Gain Coefficient; APV; Coefficient");
  g1->SetMarkerStyle(20);   // filled circle
  g1->SetMarkerSize(0.5); 
  g1->SetMarkerColor(kBlue);
  g1->SetLineStyle(0);      // disable line

  g2->SetMarkerStyle(21);   // filled square
  g2->SetMarkerSize(0.5); 
  g2->SetMarkerColor(kRed);
  g2->SetLineStyle(0);      // disable line
  
  g1->Draw("AP");
  g2->Draw("P SAME");

  auto leg = new TLegend(0.65, 0.8, 0.9, 0.9);
  leg->AddEntry(g1, "New Gain", "p");
  leg->AddEntry(g2, "Old Gain", "p");
  leg->Draw();

  c->cd(2);
  hGain->SetStats(kTRUE);
  hGain->Draw("HIST");
  gPad->Update();
  c->Update();
}

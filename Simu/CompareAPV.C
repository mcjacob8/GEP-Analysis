// ------------------------------------------------------------------------------------------------------- //
// Script for comparing spectra on an APV basis, namely for looking at varaitions in timing distributions. //
// Used to generate offsets in hit timing for each APV to improve simulation agreement.                    //  
//                                                                                                         //
// Usage:                                                                                                  //
//   root -l                                                                                               //
//   .L CompareAPV.C+                                                                                      //
//   CompareAPV("configElastic.cfg", "output/APV_time.root");                                              //
//                                                                                                         //
// Arguments:                                                                                              //
//   configfilename : Configuration file containing list of ROOT files to process and cuts to use          //
//   outfilename    : Output file for storing comparison results                                           //
//                                                                                                         //
// ---------                                                                                               //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 6-13-2026                                                  //
// ---------                                                                                               //
//  Jacob McMurtry, rby2vw@virginia.edu EDITED 6-13-2026 - CM overlay from new TTree variables             //
//  Jacob McMurtry, rby2vw@virginia.edu EDITED 6-17-2026 - Overhauled script to produce timing of all APVs //
// ---------                                                                                               //
// ** Do not tamper with this sticker! Log any updates to the script above.                                //
// ------------------------------------------------------------------------------------------------------- //

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLine.h>
#include <TStyle.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>


//Structure to hold unique APV information: module id, vtp, fiber, adcID
struct apvInfoGEM {
	int gemid, vtp, fiber, adcID;

	void print() const {
	std::cout 
	<< "GEMId: " << gemid 
	<< " VTP: " << vtp
	<< " Fiber: " << fiber
	<< " ADC ID: " << adcID
	<< std::endl;
	};
	bool operator<(const apvInfoGEM& other) const {
	if (gemid != other.gemid) return gemid < other.gemid;
	if (vtp != other.vtp) return vtp < other.vtp;
	if (fiber != other.fiber) return fiber < other.fiber;
	return adcID < other.adcID;
	}
};

// Begin Main function
void CompareAPV( const char *configfilename="configElastic.cfg", const char *outfilename="../output/APV_time.root"){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  //gStyle->SetOptStat(0); // No stats box on histograms
  
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

  const int MAXHEEP = 100;

  // List of variables that will be plotted
  double UTime[MAXHEEP], VTime[MAXHEEP];
  double ADC_ID_U[MAXHEEP], ADC_ID_V[MAXHEEP];
  double mpd_U[MAXHEEP], mpd_V[MAXHEEP];
  double crate_U[MAXHEEP], crate_V[MAXHEEP];
  double module[MAXHEEP];
  double trackindexFT[MAXHEEP];
  double nstripuFT[MAXHEEP], nstripvFT[MAXHEEP];

  // Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  // the * applies it to all branches, the 0 disables those branches. to enable would need to make 1
  C->SetBranchStatus("*",0);

  // Initializing branches
  C->SetBranchStatus("g.runnum",1);
  C->SetBranchStatus("heep.*",1);
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("sbs.gemFT.hit.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("earm.ecal.*",1);
  C->SetBranchStatus("sbs.hcal.*",1);


  // Filling Arrays
  C->SetBranchAddress("sbs.gemFT.hit.mpd_V", mpd_V);
  C->SetBranchAddress("sbs.gemFT.hit.crate_V", crate_V);
  C->SetBranchAddress("sbs.gemFT.hit.adc_id_V", ADC_ID_V);
  C->SetBranchAddress("sbs.gemFT.hit.Vtime", VTime);
  C->SetBranchAddress("sbs.gemFT.hit.mpd_U", mpd_U);
  C->SetBranchAddress("sbs.gemFT.hit.crate_U", crate_U);
  C->SetBranchAddress("sbs.gemFT.hit.adc_id_U", ADC_ID_U);
  C->SetBranchAddress("sbs.gemFT.hit.Utime", UTime);
  C->SetBranchAddress("sbs.gemFT.hit.module", module);
  C->SetBranchAddress("sbs.gemFT.hit.trackindex", trackindexFT);
  C->SetBranchAddress("sbs.gemFT.hit.nstripu",nstripuFT);
  C->SetBranchAddress("sbs.gemFT.hit.nstripv",nstripvFT);

  TFile *fout = new TFile(outfilename, "RECREATE");

  // Set up the histograms we want to look at:
  int bins = 100;   // useful for normalizing later on
  double range = 0.25;
  //Histograms to look at/compare APV timing distributions with
  TH2D *hFTM2_APV_v_Vtime = new TH2D("hFTM2_APV_v_Vtime", "FT M2 APV vs hit.Vtime ; hit.adc_id_V; hit.Vtime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);
  TH2D *hFTM2_APV_v_Utime = new TH2D("hFTM2_APV_v_Utime", "FT M2 APV vs hit.Utime ; hit.adc_id_U; hit.Utime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);

  TH2D *hFTM3_APV_v_Vtime = new TH2D("hFTM3_APV_v_Vtime", "FT M3 APV vs hit.Vtime ; hit.adc_id_V; hit.Vtime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);
  TH2D *hFTM3_APV_v_Utime = new TH2D("hFTM3_APV_v_Utime", "FT M3 APV vs hit.Utime ; hit.adc_id_U; hit.Utime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);
  
  TH2D *hFTM4_APV_v_Vtime = new TH2D("hFTM4_APV_v_Vtime", "FT M4 APV vs hit.Vtime ; hit.adc_id_V; hit.Vtime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);
  TH2D *hFTM4_APV_v_Utime = new TH2D("hFTM4_APV_v_Utime", "FT M4 APV vs hit.Utime ; hit.adc_id_U; hit.Utime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);
  
  TH2D *hFTM5_APV_v_Vtime = new TH2D("hFTM5_APV_v_Vtime", "FT M5 APV vs hit.Vtime ; hit.adc_id_V; hit.Vtime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);
  TH2D *hFTM5_APV_v_Utime = new TH2D("hFTM5_APV_v_Utime", "FT M5 APV vs hit.Utime ; hit.adc_id_U; hit.Utime (ns)", 30, -0.5, 29.5, bins/2, 0, 125);



  /**************************************************/
  /*********** Begin Event Loop and Setup ***********/
  /**************************************************/
  long nevent=0;
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
  int treenum=0, currenttreenum=0;
  int ngoodevs = 0;

  std::map<apvInfoGEM, TH1D*> hTime;

  while( C->GetEntry( nevent++ ) && nevent){
    currenttreenum = C->GetTreeNumber();

    // Setting up our cuts
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;
      if (GlobalCut) {
        delete GlobalCut;
        GlobalCut = nullptr;
      }
      GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
      GlobalCut->UpdateFormulaLeaves();
    }
    if( nevent % 10000 == 0 ) cout << nevent << endl;
    bool passedcut = GlobalCut->EvalInstance(0) != 0;

    if( passedcut ){        //Fill histograms for events passing a looser cut first to check quality of later cuts
      for (int i = 0; i < MAXHEEP; i++) {

        // Create histogram if first time seeing this APV
        if( int(trackindexFT[i]) == 0 && nstripuFT[i]>1&&nstripvFT[i]>1 ){
          apvInfoGEM key;
          key.gemid = module[i];
          key.vtp   = crate_U[i];
          key.fiber = mpd_U[i];
          key.adcID = ADC_ID_U[i];

          if (hTime.find(key) == hTime.end()) {
            TString hname = Form(
              "hTime_M%d_VTP%d_F%d_ADC%d",
              key.gemid,
              key.vtp,
              key.fiber,
              key.adcID
            );
            TString htitle = Form(
              "Module %d VTP %d Fiber %d ADC %d;Time;Counts",
              key.gemid,
              key.vtp,
              key.fiber,
              key.adcID
            );
            hTime[key] = new TH1D(
              hname,
              htitle,
              200,   // bins
              -50,   // xmin
              150    // xmax
            );
          }
          hTime[key]->Fill(UTime[i]);
          
          // Repeat for V side
          key.vtp   = crate_V[i];
          key.fiber = mpd_V[i];
          key.adcID = ADC_ID_V[i];

          if (hTime.find(key) == hTime.end()) {
            TString hname = Form(
              "hTime_M%d_VTP%d_F%d_ADC%d",
              key.gemid,
              key.vtp,
              key.fiber,
              key.adcID
            );
            TString htitle = Form(
              "Module %d VTP %d Fiber %d ADC %d;Time;Counts",
              key.gemid,
              key.vtp,
              key.fiber,
              key.adcID
            );
            hTime[key] = new TH1D(
              hname,
              htitle,
              200,   // bins
              -50,   // xmin
              150    // xmax
            );
          }
          hTime[key]->Fill(VTime[i]);
        }
      }
    }
  }

  // Now we fit these histograms and store the mean time values
  std::map<apvInfoGEM, double> APVMean;

  for (auto const& entry : hTime) {
    const apvInfoGEM& key = entry.first;
    TH1D* h = entry.second;

    if (h->GetEntries() < 1000) continue;

    double peak = h->GetBinCenter(h->GetMaximumBin());
    TF1 gaus("gaus","gaus", peak - 13, peak + 13);
    h->Fit(&gaus,"RQ");

    double mean = gaus.GetParameter(1);

    APVMean[key] = mean-80;
  }


  //Create a root file to save these histograms to
  for (auto const& histpair : hTime) {
    histpair.second->Write();
  }
  fout->Close();

  TString outfiletxt = outfilename;
  outfiletxt.ReplaceAll(".root", ".txt");
  std::ofstream outfile(outfiletxt.Data());

  outfile << "# Module VTP Fiber APV MeanTime\n";

  for (auto const& entry : APVMean) {

    const apvInfoGEM& key = entry.first;
    double mean = entry.second;

    outfile
        << key.gemid << " "
        << key.vtp   << " "
        << key.fiber << " "
        << key.adcID << " "
        << mean
        << "\n";
  }

  outfile.close();

}
  
    
    
// Old components of this script that may be useful for future reference:
      //M2 V strip if fiber = 8 -> Fill
      //M2 U strip if fiber = 10 -> Fill
      //M3 V strip if fiber = 0 -> Fill
      //M3 U strip if fiber = 2 -> Fill

      //Some sort of psuedo code:
      //Loop over the hit.Utime vector
        //if nmodule[j] == 2 and mpd[j] == 8
          //Fill (adcid vs Utime)
          //Fill (adcid vs Utime)
          //etc.

      /*for (int i = 0; i < MAXHEEP; i++) {

        int mod = module[i];
        if( int(trackindexFT[i]) == 0 && nstripuFT[i]>1&&nstripvFT[i]>1 ){
          switch(mod) {
            case 2: {
              int offsetV = (mpd_V[i] == 9) ? 15 : 0;
              hFTM2_APV_v_Vtime->Fill(offsetV + ADC_ID_V[i], VTime[i]);
              int offsetU = (mpd_U[i] == 11) ? 15 : 0;
              hFTM2_APV_v_Utime->Fill(offsetU + ADC_ID_U[i], UTime[i]);
              break;
            }
            case 3: {
              int offsetV = (mpd_V[i] == 1) ? 15 : 0;
              hFTM3_APV_v_Vtime->Fill(offsetV + ADC_ID_V[i], VTime[i]);
              int offsetU = (mpd_U[i] == 3) ? 15 : 0;
              hFTM3_APV_v_Utime->Fill(offsetU + ADC_ID_U[i], UTime[i]);
              break;
            }
            case 4: {
              int offsetV = (mpd_V[i] == 5) ? 15 : 0;
              hFTM4_APV_v_Vtime->Fill(offsetV + ADC_ID_V[i], VTime[i]);
              int offsetU = (mpd_U[i] == 8) ? 15 : 0;
              hFTM4_APV_v_Utime->Fill(offsetU + ADC_ID_U[i], UTime[i]);
              break;
            }
            case 5: {
              int offsetV = (mpd_V[i] == 10) ? 15 : 0;
              hFTM5_APV_v_Vtime->Fill(offsetV + ADC_ID_V[i], VTime[i]);
              int offsetU = (mpd_U[i] == 0) ? 15 : 0;
              hFTM5_APV_v_Utime->Fill(offsetU + ADC_ID_U[i], UTime[i]);
              break;
            }
            default:
              continue;
          }
        }*/


        /*
        // Module 2
        if (module[i]==2 && mpd_V[i] == 8){
          hFTM2_APV_v_Vtime->Fill(ADC_ID_V[i], VTime[i]);
        } else if(module[i]==2 && mpd_V[i] == 9){
          hFTM2_APV_v_Vtime->Fill(15+ADC_ID_V[i], VTime[i]);
        } else if (module[i]==2 && mpd_U[i] == 10){
          hFTM2_APV_v_Utime->Fill(ADC_ID_U[i], UTime[i]);
        } else if(module[i]==2 && mpd_U[i] == 11){
          hFTM2_APV_v_Utime->Fill(15+ADC_ID_U[i], UTime[i]);
        // Module 3
        } else if (module[i]==3 && mpd_V[i] == 0){
          hFTM3_APV_v_Vtime->Fill(ADC_ID_V[i], VTime[i]);
        } else if(module[i]==3 && mpd_V[i] == 1){
          hFTM3_APV_v_Vtime->Fill(15+ADC_ID_V[i], VTime[i]);
        } else if (module[i]==3 && mpd_U[i] == 2){
          hFTM3_APV_v_Utime->Fill(ADC_ID_U[i], UTime[i]);
        } else if(module[i]==3 && mpd_U[i] == 3){
          hFTM3_APV_v_Utime->Fill(15+ADC_ID_U[i], UTime[i]);
        // Module 4
        } else if (module[i]==4 && mpd_V[i] == 4){
          hFTM4_APV_v_Vtime->Fill(ADC_ID_V[i], VTime[i]);
        } else if(module[i]==4 && mpd_V[i] == 5){
          hFTM4_APV_v_Vtime->Fill(15+ADC_ID_V[i], VTime[i]);
        } else if (module[i]==4 && mpd_U[i] == 6){
          hFTM4_APV_v_Utime->Fill(ADC_ID_U[i], UTime[i]);
        } else if(module[i]==4 && mpd_U[i] == 8){
          hFTM4_APV_v_Utime->Fill(15+ADC_ID_U[i], UTime[i]);
        // Module 5
        } else if (module[i]==5 && mpd_V[i] == 9){
          hFTM5_APV_v_Vtime->Fill(ADC_ID_V[i], VTime[i]);
        } else if(module[i]==5 && mpd_V[i] == 10){
          hFTM5_APV_v_Vtime->Fill(15+ADC_ID_V[i], VTime[i]);
        } else if (module[i]==5 && mpd_U[i] == 11){
          hFTM5_APV_v_Utime->Fill(ADC_ID_U[i], UTime[i]);
        } else if(module[i]==5 && mpd_U[i] == 0){
          hFTM5_APV_v_Utime->Fill(15+ADC_ID_U[i], UTime[i]);

        } else {
          continue; // Skip if not one of the specified conditions
        }*/
/*      }
    }
  }

  //Let's do some plotting and saving to a file
  TString outfilepdf = outfilename; // Make a pdf file to save these histograms to
  outfilepdf.ReplaceAll(".root",".pdf");

  TProfile *profM2V = hFTM2_APV_v_Vtime->ProfileX("profM2V", 1, -1, "s");
  TProfile *profM2U = hFTM2_APV_v_Utime->ProfileX("profM2U", 1, -1, "s");
  TProfile *profM3V = hFTM3_APV_v_Vtime->ProfileX("profM3V", 1, -1, "s");
  TProfile *profM3U = hFTM3_APV_v_Utime->ProfileX("profM3U", 1, -1, "s");
  TProfile *profM4V = hFTM4_APV_v_Vtime->ProfileX("profM4V", 1, -1, "s");
  TProfile *profM4U = hFTM4_APV_v_Utime->ProfileX("profM4U", 1, -1, "s");
  TProfile *profM5V = hFTM5_APV_v_Vtime->ProfileX("profM5V", 1, -1, "s");
  TProfile *profM5U = hFTM5_APV_v_Utime->ProfileX("profM5U", 1, -1, "s");
  profM2V->SetLineColor(kRed);
  profM2U->SetLineColor(kRed);
  profM3V->SetLineColor(kRed);
  profM3U->SetLineColor(kRed);
  profM4V->SetLineColor(kRed);
  profM4U->SetLineColor(kRed);
  profM5V->SetLineColor(kRed);
  profM5U->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(2,2);
  c1->cd(1);    hFTM2_APV_v_Vtime->Draw("COLZ");    profM2V->Draw("sames");
  gPad->SetLogz(1);
  c1->cd(2);    hFTM3_APV_v_Vtime->Draw("COLZ");    profM3V->Draw("sames");
  gPad->SetLogz(1);
  c1->cd(3);    hFTM4_APV_v_Vtime->Draw("COLZ");    profM4V->Draw("sames");
  gPad->SetLogz(1);
  c1->cd(4);    hFTM5_APV_v_Vtime->Draw("COLZ");    profM5V->Draw("sames");
  gPad->SetLogz(1);
  c1->Update();
  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page
  //c1->Print(outfilepdf);

  TCanvas *c2 = new TCanvas("c2","",1600,1200);
  c2->Divide(2,2);
  c2->cd(1);    hFTM2_APV_v_Utime->Draw("COLZ");    profM2U->Draw("sames");
  gPad->SetLogz(1);
  c2->cd(2);    hFTM3_APV_v_Utime->Draw("COLZ");    profM3U->Draw("sames");
  gPad->SetLogz(1);
  c2->cd(3);    hFTM4_APV_v_Utime->Draw("COLZ");    profM4U->Draw("sames");
  gPad->SetLogz(1);
  c2->cd(4);    hFTM5_APV_v_Utime->Draw("COLZ");    profM5U->Draw("sames");
  gPad->SetLogz(1);
  c2->Update();
  c2->Print(outfilepdf + ")");  //Close the pdf

  fout->Write();
  cout << "That's all folks!" << endl;
}*/
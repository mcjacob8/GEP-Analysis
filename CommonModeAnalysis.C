// ------------------------------------------------------------------------------------------------------- //
// Script for GEM common mode studies: Intended to create plots of APV raw                                 //
// frames with common mode values overlaid for diagnostic studies.                                         //
//                                                                                                         //
// Usage:                                                                                                  //
//   root -l                                                                                               //
//   .L plot_gem_apv_commonmode.C+                                                                         //
//   plot_gem_apv_commonmode("replay.root", 0, 4, 4, 4);                                                   //
//                                                                                                         //
// Arguments:                                                                                              //
//   filename : SBS-Replay ROOT output file                                                                //
//   event    : event number                                                                               //
//   module   : GEM module index (m0, m1, ...)                                                             //
//   mpd      : MPD (fiber) number (module specific)                                                       //
//   apv      : APV starting number                                                                        //
//                                                                                                         //
// For this script to be useful, full readout events with no zero suppression are needed (every 1/100      //
// events, or data taken in pedestal mode, or events that fail ONLINE zero suppression).                   //
// Data needs to be replayed in a specific manner for this to work properly:                               //
//   replay_gep.C                                                       //                                 //
//     In the argument of the main replay script, requiretrack = 0, nontrackingmode = 1, and dogems = 1.   //
//   replay_gep.cdef                                                                                       //
//     Make sure there are no cuts that require tracks (e.g. #GoodFrontTrack)                              //     
//   db_sbs.gemF*.dat                                                                                      //
//     Ensure that sbs.gemF*.zerosuppress = 0                                                              //
//   replay_F*GEM_gep.odef                                                                                 //
//     Uncomment block for strip variables of the desired module                                           //  
//   SBSGEMModule.cxx/h                                                                                    //
//     Ensure that common-mode variables are stored in the output ROOT tree                                //  
//     (i.e. use proper version that includes sbs.gemF*.m*.strip.CMsamples etc.)                           //   
//                                                                                                         //
// ---------                                                                                               //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 12-15-2025                                                 //
// ---------                                                                                               //
// ** Do not tamper with this sticker! Log any updates to the script above.                                //
//  Jacob McMurtry, rby2vw@virginia.edu EDITED 1-8-2026 - CM overlay from new TTree variables              //
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

//void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEM_strip/GEP3/LH2/rootfiles/gep5_fullreplay_6085_stream0_0_seg70_70.root", Long64_t event =41, int module = 4,int mpd = 4, int APV = 2)
// Run 5514 - Pedestal Run
//void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/mcjacob/GEP/5514/rootfiles/gep5_replayed_5514_stream0_0_seg1_1_firstevent0_nevent1000.root", Long64_t event = 2, int module = 4,int mpd = 4, int APV = 4)
// Run 3175 - 1 uA FRO run, done using zerosuppress = 0
//void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/mcjacob/GEP/3175/rootfiles/gep5_replayed_3175_stream2_2_seg1_1_firstevent0_nevent1000.root", Long64_t event = 3, int module = 4,int mpd = 4, int APV = 4)
// Run 6075 - 22 uA run
//void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/mcjacob/GEP/rootfiles/gep5_replayed_6075_stream2_2_seg0_0_firstevent0_nevent5000.root", Long64_t event = 121, int module = 4,int mpd = 4, int APV = 4)
// Run 6075 "non tracking mode"
//void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/mcjacob/GEP/Jan6/test4/rootfiles/gep5_replayed_3175_stream1_1_seg0_0_firstevent0_nevent5000.root", Long64_t event = 10, int module = 4,int mpd = 4, int APV = 4)
//void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/mcjacob/GEP/Jan6/test4/rootfiles/gep5_replayed_6075_stream1_1_seg0_0_firstevent0_nevent5000.root", Long64_t event = 168, int module = 4,int mpd = 4, int APV = 4)
void CommonModeAnalysis(const char *filename="/volatile/halla/sbs/mcjacob/GEP/Jan7/rootfiles/gep5_replayed_6075_stream1_1_seg0_0_firstevent0_nevent5000.root", Long64_t event = 67, int module = 4,int mpd = 4, int APV = 4)
{

  gStyle->SetOptStat(0);

  const int NSAMP = 6;
  const int NSTRIP = 128;   // Number of strips on APV


  cout << "Open seseame?" << std::endl;

  // ---------------------------------------------------------
  // Open file and tree
  // ---------------------------------------------------------
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file " << filename << std::endl;
    return;
  }
  //TChain *T = new TChain("T");
  //T->Add(filename);
  TTree *T = (TTree*)f->Get("T");
  if (!T) {
    std::cerr << "Tree T not found" << std::endl;
    return;
  }

  // ---------------------------------------------------------
  // Branches (GEp-V naming)
  // ---------------------------------------------------------

  const int MAXSTRIP = 6*128;

  int MAXSTRIPS = 100000;
  int test = 3840*2 +2;

  Double_t   strip_mpd[MAXSTRIPS];
  Double_t   adc_id[test];
  Double_t   strip_istrip[MAXSTRIPS];
  Double_t   enable_cm[MAXSTRIPS];
  Double_t   cm_good[MAXSTRIPS];
  Double_t   BUILD_ALL_SAMPLES[test];

  Double_t adc_samples[6*3840*2 +6 ];  // 1D packed
  Double_t rawadc_samples[6*3840*2 +6];
  Double_t CMsamples[6*3840*2 +6];
  Double_t   nstrips;
  Double_t   IsU[test];
  Double_t   IsV[test];

  cout << "Huh?" << std::endl;

  TString prefix = Form("sbs.gemFT.m%d.strip", module);

  T->SetBranchStatus("*",0);
  //T->SetBranchStatus((prefix + ".*").Data(), 1);
  T->SetBranchStatus(prefix + ".*", 1);

  T->SetBranchAddress(prefix + ".nstripsfired", &nstrips);

  T->SetBranchAddress(prefix + ".mpd",        strip_mpd);
  T->SetBranchAddress(prefix + ".adc_id",     adc_id);
  T->SetBranchAddress(prefix + ".istrip",     strip_istrip);
  T->SetBranchAddress(prefix + ".IsU",        IsU);
  T->SetBranchAddress(prefix + ".IsV",        IsV);
  T->SetBranchAddress(prefix + ".ENABLE_CM",  enable_cm);
  T->SetBranchAddress(prefix + ".CM_GOOD",    cm_good);
  T->SetBranchAddress(prefix + ".BUILD_ALL_SAMPLES",    BUILD_ALL_SAMPLES);
  T->SetBranchAddress(prefix + ".ADCsamples", adc_samples);
  T->SetBranchAddress(prefix + ".rawADCsamples", rawadc_samples);
  T->SetBranchAddress(prefix + ".CMsamples", CMsamples);

  cout << "Analyzing file: " << filename << std::endl;

  // ---------------------------------------------------------
  // Load event
  // ---------------------------------------------------------
  T->GetEntry(event);

  cout << "seg fault here?" << std::endl;

  // ---------------------------------------------------------
  // Canvas
  // ---------------------------------------------------------
  TCanvas *c = new TCanvas("c",
    Form("GEM m%d MPD %d Event %lld", module, mpd, event),
    1200, 800);
  c->Divide(3,2);
  //c->Divide(2,1);

    cout << "seg fault here 2?" << std::endl;

  // ---------------------------------------------------------
  // Loop over samples
  // ---------------------------------------------------------

  
  //TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", "sbs.gemFT.m4.strip.ENABLE_CM == 0", T );
  /*TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", "sbs.gemFT.m4.strip.nstripsfired > 6000", T );
  int nevent = 0;
  while( T->GetEntry( nevent++ ) && nevent){
    //GlobalCut = new TTreeFormula( "GlobalCut", "sbs.gemFT.m4.strip.ENABLE_CM == 0", T );
    GlobalCut = new TTreeFormula( "GlobalCut", "sbs.gemFT.m4.strip.BUILD_ALL_SAMPLES > 0 && sbs.gemFT.m4.strip.CM_GOOD == 1", T );
    bool passedcut = GlobalCut->EvalInstance(0) != 0;
    if( passedcut ){
      cout << "success" << std::endl;
      cout << nevent << std::endl;
      //break;
    }

  }*/

  // Strip variables are stored in a 1D array where the index = isamp + NSAMP*istrip
  // Where isamp is the time sample index (i.e. 0,1,2,3,4,5), NSAMP is the number of time samples (6), and istrip is the index of firedstrips (not to be confused with strip number strip.istrip[istrip]))
  for (int apvs = 0; apvs < 6; apvs++) {
    cout << "nstrips " << nstrips << std::endl;
    TH1D *h_APVFrame = new TH1D("h_APVframe", Form("Module %d MPD %d APV %d;Strip;ADC", module, mpd, apvs+APV), NSTRIP*NSAMP, 0, NSTRIP*NSAMP);
    TH1D *h_APVFrameCM = new TH1D("h_APVframeCM", Form("Module %d MPD %d APV %d;Strip;ADC", module, mpd, apvs+APV), NSTRIP*NSAMP, 0, NSTRIP*NSAMP);
    for (int isamp = 0; isamp < NSAMP; isamp++) {
    
      for (int i = 0; i < nstrips; i++) {   // Loop over number of strips fired
        //for (int i = 0; i < 3000; i++) {
        int strip = (int)strip_istrip[i];
        double adc = adc_samples[i*NSAMP + isamp];
        double rawadc = rawadc_samples[i*NSAMP + isamp];
        double cm = CMsamples[i*NSAMP + isamp];

        if (strip_mpd[i] == mpd && adc_id[i] == apvs+APV) { // Check if strip is in the correct MPD and APV
          //cout << strip << std::endl;
          int local = (strip - 1) % 128 + 1;
          //cout << "local " << local << std::endl;
          h_APVFrame->SetBinContent(local+isamp*NSTRIP, rawadc);
          h_APVFrameCM->SetBinContent(local+isamp*NSTRIP, cm);
        }
      }
    }
    double ymin = -100.0;
    double ymax = h_APVFrame->GetMaximum() + 100.0;

    c->cd(apvs+1);
    h_APVFrame->SetMinimum(ymin);
    h_APVFrame->Draw("LP");
    h_APVFrameCM->SetLineColor(kRed);
    h_APVFrameCM->SetLineWidth(2);
    h_APVFrameCM->Draw("L same");

    TLine *l1 = new TLine(128, ymin, 128, ymax);
    l1->SetLineStyle(2);   l1->SetLineWidth(1);
    l1->Draw("same");
    TLine *l2 = new TLine(128*2, ymin, 128*2, ymax);
    l2->SetLineStyle(2);   l2->SetLineWidth(1);
    l2->Draw("same");
    TLine *l3 = new TLine(128*3, ymin, 128*3, ymax);
    l3->SetLineStyle(2);   l3->SetLineWidth(1);
    l3->Draw("same");
    TLine *l4 = new TLine(128*4, ymin, 128*4, ymax);
    l4->SetLineStyle(2);   l4->SetLineWidth(1);
    l4->Draw("same");
    TLine *l5 = new TLine(128*5, ymin, 128*5, ymax);
    l5->SetLineStyle(2);   l5->SetLineWidth(1);
    l5->Draw("same");
  }
  /*TH1D *h_APVFrame = new TH1D("h_APVframe", Form("Module %d MPD %d APV %d;Strip;ADC", module, mpd, APV), NSTRIP*NSAMP, 0, NSTRIP*NSAMP);
  for (int isamp = 0; isamp < NSAMP; isamp++) {
    
    for (int i = 0; i < nstrips; i++) {   // Loop over number of strips fired

      int strip = (int)strip_istrip[i];
      double adc = adc_samples[i*NSAMP + isamp];

      if (strip_mpd[i] == mpd && adc_id[i] == APV) { // Check if strip is in the correct MPD and APV
        cout << strip << std::endl;
        int local = (strip - 1) % 128 + 1;
        cout << "local " << local << std::endl;
        h_APVFrame->SetBinContent(local+isamp*NSTRIP, adc);

      }

    }

  }
  c->cd(1);
  h_APVFrame->SetMinimum(-100.0);
  h_APVFrame->Draw("LP");*/
  c->Update();
}


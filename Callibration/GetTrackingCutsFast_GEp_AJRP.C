#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/QuantFuncMathCore.h"

using namespace ROOT;
using namespace Math;

void DrawHist(TH1D *hist, vector<double> value, int flag = 0){

  //flag = 0 for calorimeter plots
  //flag = 1 for ADC plots
  //flag = 2 for time plots

  hist->Draw();

  TPaveText *pt;
  if(flag == 0){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    pt->AddText(Form("Center = %g",value[0]));
    pt->AddText(Form("Width = %g",value[1]));
  }
  else if(flag == 1){
    pt = new TPaveText(0.55,0.75,0.75,0.88,"ndc");
    pt->AddText(Form("Threshold = %g",value[0]));
  }
  else if(flag == 2){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    if(value.size() == 2){
      pt->AddText(Form("Mean = %g",value[0]));
      pt->AddText(Form("#sigma = %g",value[1]));
    }
    if(value.size() == 1){
      pt->AddText(Form("Cut = %g",value[0]));
    }
 }
  pt->SetFillColor(0);
  pt->Draw("same");


}

void FitGaus_FWHM( TH1D *htest, double thresh=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh+1);

  htest->Fit("gaus","qS","",xlow, xhigh);
}

void GetTrackingCutsFast_GEp_AJRP( const char *configfilename, const char *outfilename="output/GEPtrackingcuts_10uA.root", int nmodules=14, int nmodules_back=32, double thresh=0.005, double nsig_tstrip=4.5, double nsig_dt=5.0){

 
  
  ifstream infile(configfilename);
  
  TChain *C = new TChain("T");

  TString currentline;
  while( infile >> currentline && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( infile >> currentline && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  //Initialize branches:
  
  //what branches do we need? Tracking stuff:
  //int ntracks=0;
  const int MAXNTRACKS=10;
  //int nhits=0;
  const int MAXNHITS = 100;
  const int MAXNCLUST = 1;
  const int MAXNTRIG = 10;

  C->SetBranchStatus("*",0);

  //Branches to enable:
  
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("sbs.gemFT.hit.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("sbs.gemFPP.hit.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("heep.*",1);

  //enable these top-level ECAL branches for cuts only:
  C->SetBranchStatus("earm.ecal.e",1);
  C->SetBranchStatus("earm.ecal.nblk",1);
  C->SetBranchStatus("earm.ecal.x",1);
  C->SetBranchStatus("earm.ecal.y",1);
  C->SetBranchStatus("earm.ecal.eblk",1);
  C->SetBranchStatus("earm.ecal.atimeblk",1);

  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.nblk",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.eblk",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  
  //Comment irrelevant bigbite branches
  // double EPS, ESH;
  // C->SetBranchStatus("bb.etot_over_p",1);
  // C->SetBranchStatus("bb.ps.e",1);
  // C->SetBranchStatus("bb.sh.e",1);
  // C->SetBranchAddress("bb.ps.e",&EPS);
  // C->SetBranchAddress("bb.sh.e",&ESH);

  // C->SetBranchStatus("bb.grinch_tdc.clus.*",1);
  
  //Track variables for FT:
  double ntracks;
  double tracknhits[MAXNTRACKS];
  double trackngoodhits[MAXNTRACKS];
  double trackchi2ndf[MAXNTRACKS];
  double trackchi2ndf_hitquality[MAXNTRACKS];
  double trackt0[MAXNTRACKS];
  double xfp[MAXNTRACKS], yfp[MAXNTRACKS], thfp[MAXNTRACKS], phfp[MAXNTRACKS];
  double xtar[MAXNTRACKS], ytar[MAXNTRACKS], thtar[MAXNTRACKS], phtar[MAXNTRACKS];
  double rxfp[MAXNTRACKS], ryfp[MAXNTRACKS], rthfp[MAXNTRACKS], rphfp[MAXNTRACKS];
  double dxfp[MAXNTRACKS], dyfp[MAXNTRACKS], dthfp[MAXNTRACKS], dphfp[MAXNTRACKS];
  double p[MAXNTRACKS],px[MAXNTRACKS],py[MAXNTRACKS],pz[MAXNTRACKS],vz[MAXNTRACKS]; 
  //double xfcp[MAXNCLUST],yfcp[MAXNCLUST],zfcp[MAXNCLUST],xbcp[MAXNCLUST],ybcp[MAXNCLUST],zbcp[MAXNCLUST];
  //add variables for theta, phi, sclose, zclose in case tracker is a polarimeter tracker:

  //Track variables for back tracker:
  double ntracks_back;
  double tracknhitsback[MAXNTRACKS],trackngoodhitsback[MAXNTRACKS],trackchi2ndfback[MAXNTRACKS];
  double trackchi2ndf_hitqualityback[MAXNTRACKS],trackt0back[MAXNTRACKS];
  double xback[MAXNTRACKS],yback[MAXNTRACKS],xpback[MAXNTRACKS],ypback[MAXNTRACKS];
  double theta[MAXNTRACKS],phi[MAXNTRACKS],sclose[MAXNTRACKS],zclose[MAXNTRACKS];

  C->SetBranchAddress("sbs.gemFT.track.ntrack", &ntracks );
 
  
  C->SetBranchAddress("sbs.gemFPP.track.ntrack", &ntracks_back );
  C->SetBranchAddress("sbs.gemFPP.track.sclose", sclose );
  
  //Hit variables for FT: 
  double ngoodhits_FT;
  //Track hit variables:
  double modFT[MAXNHITS];
  double nstripuFT[MAXNHITS],nstripvFT[MAXNHITS];
  double trackindexFT[MAXNHITS];
  double ADCmaxsampU_FT[MAXNHITS], ADCmaxsampV_FT[MAXNHITS];
  double ADCmaxstripU_FT[MAXNHITS], ADCmaxstripV_FT[MAXNHITS];
  double ADCU_FT[MAXNHITS], ADCV_FT[MAXNHITS], ADCavg_FT[MAXNHITS], ADCasym_FT[MAXNHITS];
  //what's relevant for deconvolution is the max two-sample combo and the cluster-summed max two-sample combo:
  double DeconvADCU_FT[MAXNHITS], DeconvADCV_FT[MAXNHITS], DeconvADCmaxcomboU_FT[MAXNHITS], DeconvADCmaxcomboV_FT[MAXNHITS];
  double ADCasymDeconv_FT[MAXNHITS];
  double isampmaxUstrip_FT[MAXNHITS], isampmaxVstrip_FT[MAXNHITS];
  double icombomaxUstripDeconv_FT[MAXNHITS], icombomaxVstripDeconv_FT[MAXNHITS];

  double UtimeMaxStrip_FT[MAXNHITS], VtimeMaxStrip_FT[MAXNHITS];
  double UtimeMaxStripDeconv_FT[MAXNHITS], VtimeMaxStripDeconv_FT[MAXNHITS];
  double Utime_FT[MAXNHITS], Vtime_FT[MAXNHITS];
  double UtimeDeconv_FT[MAXNHITS], VtimeDeconv_FT[MAXNHITS];

  double HitTavg_corr_FT[MAXNHITS];
  double deltat_FT[MAXNHITS], deltatDeconv_FT[MAXNHITS];

  C->SetBranchAddress("sbs.gemFT.hit.ngoodhits",&ngoodhits_FT);
  C->SetBranchAddress("sbs.gemFT.hit.module",modFT);
  C->SetBranchAddress("sbs.gemFT.hit.nstripu",nstripuFT);
  C->SetBranchAddress("sbs.gemFT.hit.nstripv",nstripvFT);
  C->SetBranchAddress("sbs.gemFT.hit.trackindex",trackindexFT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampU",ADCmaxsampU_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampV",ADCmaxsampV_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCmaxstripU",ADCmaxstripU_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCmaxstripV",ADCmaxstripV_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCU",ADCU_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCV",ADCV_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCavg",ADCavg_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCasym",ADCasym_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCmaxcomboUclust_deconv",DeconvADCU_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCmaxcomboVclust_deconv",DeconvADCV_FT);
  C->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxcomboU",DeconvADCmaxcomboU_FT);
  C->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxcomboV",DeconvADCmaxcomboV_FT);
  C->SetBranchAddress("sbs.gemFT.hit.ADCasym_deconv",ADCasymDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.isampmaxUstrip",isampmaxUstrip_FT);
  C->SetBranchAddress("sbs.gemFT.hit.isampmaxVstrip",isampmaxVstrip_FT);
  C->SetBranchAddress("sbs.gemFT.hit.icombomaxUstripDeconv",icombomaxUstripDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.icombomaxVstripDeconv",icombomaxVstripDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStrip",UtimeMaxStrip_FT);
  C->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStrip",VtimeMaxStrip_FT);
  C->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStripDeconv",UtimeMaxStripDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStripDeconv",VtimeMaxStripDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.Utime",Utime_FT);
  C->SetBranchAddress("sbs.gemFT.hit.Vtime",Vtime_FT);
  C->SetBranchAddress("sbs.gemFT.hit.UtimeDeconv",UtimeDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.VtimeDeconv",VtimeDeconv_FT);
  C->SetBranchAddress("sbs.gemFT.hit.Tavg_corr",HitTavg_corr_FT);
  C->SetBranchAddress("sbs.gemFT.hit.deltat",deltat_FT);
  C->SetBranchAddress("sbs.gemFT.hit.deltat_deconv",deltatDeconv_FT);
  
  double ngoodhits_FPP;

  double modFPP[MAXNHITS];
  double nstripuFPP[MAXNHITS],nstripvFPP[MAXNHITS];
  double trackindexFPP[MAXNHITS];
  double ADCmaxsampU_FPP[MAXNHITS], ADCmaxsampV_FPP[MAXNHITS];
  double ADCmaxstripU_FPP[MAXNHITS], ADCmaxstripV_FPP[MAXNHITS];
  double ADCU_FPP[MAXNHITS], ADCV_FPP[MAXNHITS], ADCavg_FPP[MAXNHITS], ADCasym_FPP[MAXNHITS];
  //what's relevant for deconvolution is the max two-sample combo and the cluster-summed max two-sample combo:
  double DeconvADCU_FPP[MAXNHITS], DeconvADCV_FPP[MAXNHITS], DeconvADCmaxcomboU_FPP[MAXNHITS], DeconvADCmaxcomboV_FPP[MAXNHITS];

  double ADCasymDeconv_FPP[MAXNHITS];
  double isampmaxUstrip_FPP[MAXNHITS], isampmaxVstrip_FPP[MAXNHITS];
  double icombomaxUstripDeconv_FPP[MAXNHITS], icombomaxVstripDeconv_FPP[MAXNHITS];

  double UtimeMaxStrip_FPP[MAXNHITS], VtimeMaxStrip_FPP[MAXNHITS];
  double UtimeMaxStripDeconv_FPP[MAXNHITS], VtimeMaxStripDeconv_FPP[MAXNHITS];
  double Utime_FPP[MAXNHITS], Vtime_FPP[MAXNHITS];
  double UtimeDeconv_FPP[MAXNHITS], VtimeDeconv_FPP[MAXNHITS];

  double HitTavg_corr_FPP[MAXNHITS];
  double deltat_FPP[MAXNHITS], deltatDeconv_FPP[MAXNHITS];


  C->SetBranchAddress("sbs.gemFPP.hit.ngoodhits",&ngoodhits_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.module",modFPP);
  C->SetBranchAddress("sbs.gemFPP.hit.nstripu",nstripuFPP);
  C->SetBranchAddress("sbs.gemFPP.hit.nstripv",nstripvFPP);
  C->SetBranchAddress("sbs.gemFPP.hit.trackindex",trackindexFPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampU",ADCmaxsampU_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampV",ADCmaxsampV_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCmaxstripU",ADCmaxstripU_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCmaxstripV",ADCmaxstripV_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCU",ADCU_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCV",ADCV_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCavg",ADCavg_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCasym",ADCasym_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCmaxcomboUclust_deconv",DeconvADCU_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCmaxcomboVclust_deconv",DeconvADCV_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxcomboU",DeconvADCmaxcomboU_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxcomboV",DeconvADCmaxcomboV_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.ADCasym_deconv",ADCasymDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.isampmaxUstrip",isampmaxUstrip_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.isampmaxVstrip",isampmaxVstrip_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.icombomaxUstripDeconv",icombomaxUstripDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.icombomaxVstripDeconv",icombomaxVstripDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStrip",UtimeMaxStrip_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStrip",VtimeMaxStrip_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStripDeconv",UtimeMaxStripDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStripDeconv",VtimeMaxStripDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.Utime",Utime_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.Vtime",Vtime_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.UtimeDeconv",UtimeDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.VtimeDeconv",VtimeDeconv_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.Tavg_corr",HitTavg_corr_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.deltat",deltat_FPP);
  C->SetBranchAddress("sbs.gemFPP.hit.deltat_deconv",deltatDeconv_FPP);

  TFile *fout = new TFile(outfilename,"RECREATE");


  //Front Tracker histos:
  TH2D *hADCmaxsampU_vs_modFT = new TH2D("hADCmaxsampU_vs_modFT", "FT U strips; module; Max U strip max sample ADC", nmodules, -0.5, nmodules-0.5, 400,0,4000);
  TH2D *hADCmaxsampV_vs_modFT = new TH2D("hADCmaxsampV_vs_modFT", "FT V strips; module; Max V strip max sample ADC", nmodules, -0.5, nmodules-0.5, 400,0,4000);
  TH2D *hADCmaxstripU_vs_modFT = new TH2D("hADCmaxstripU_vs_modFT", "FT U strips; module; Max U strip ADC sum", nmodules, -0.5, nmodules-0.5, 1000,0,10000);
  TH2D *hADCmaxstripV_vs_modFT = new TH2D("hADCmaxstripV_vs_modFT", "FT V strips; module; Max V strip ADC sum", nmodules, -0.5, nmodules-0.5, 1000,0,10000);
  TH2D *hADCU_vs_modFT = new TH2D("hADCU_vs_modFT", "FT U cluster sum; module; cluster ADC sum", nmodules, -0.5, nmodules-0.5, 1000,0,20000);
  TH2D *hADCV_vs_modFT = new TH2D("hADCV_vs_modFT", "FT V cluster sum; module; cluster ADC sum", nmodules, -0.5, nmodules-0.5, 1000,0,20000);
  

  TH2D *hADCavg_vs_modFT = new TH2D("hADCavg_vs_modFT", "FT hit cluster sums; module; (ADCU+ADCV)/2", nmodules,-0.5,nmodules-0.5, 1000,0,20000);
  TH2D *hADCasym_vs_modFT = new TH2D("hADCasym_vs_modFT", "FT hit cluster sum asymmetry; module; (ADCU-ADCV)/(ADCU+ADCV)",nmodules,-0.5,nmodules-0.5,250,-1.0,1.0);
 
  TH2D *hDeconvADCmaxcomboUstrip_vs_modFT = new TH2D("hDeconvADCmaxcomboUstrip_vs_modFT","FT max U strip; module; Max Deconv. two-sample combo",nmodules,-0.5,nmodules-0.5,400,0,4000);
  TH2D *hDeconvADCmaxcomboVstrip_vs_modFT = new TH2D("hDeconvADCmaxcomboVstrip_vs_modFT","FT max V strip; module; Max Deconv. two-sample combo",nmodules,-0.5,nmodules-0.5,400,0,4000);

  TH2D *hDeconvADCmaxcomboUclust_vs_modFT = new TH2D("hDeconvADCmaxcomboUclust_vs_modFT","FT U cluster sum; module; Max Deconv. two-sample combo",nmodules,-0.5,nmodules-0.5,1000,0,10000);
  TH2D *hDeconvADCmaxcomboVclust_vs_modFT = new TH2D("hDeconvADCmaxcomboVclust_vs_modFT","FT V cluster sum; module; Max Deconv. two-sample combo",nmodules,-0.5,nmodules-0.5,1000,0,10000);

  TH2D *hADCasymDeconv_vs_modFT = new TH2D("hADCasymDeconv_vs_modFT","FT Deconv ADC asym;module;(ADCU-ADCV)/(ADCU+ADCV)",nmodules,-0.5,nmodules-0.5,250,-1,1);

  TH2D *hisampmaxUstrip_vs_modFT = new TH2D("hisampmaxUstrip_vs_modFT","FT max U strip; module; max time sample",nmodules,-0.5,nmodules-0.5,6,-0.5,5.5);
  TH2D *hisampmaxVstrip_vs_modFT = new TH2D("hisampmaxVstrip_vs_modFT","FT max V strip; module; max time sample",nmodules,-0.5,nmodules-0.5,6,-0.5,5.5);

  TH2D *hicombomaxUstrip_vs_modFT = new TH2D("hicombomaxUstrip_vs_modFT","FT max U strip; module; max Deconv combo (isamp = 2nd sample of combo)",nmodules,-0.5,nmodules-0.5,7,-0.5,6.5);
  TH2D *hicombomaxVstrip_vs_modFT = new TH2D("hicombomaxVstrip_vs_modFT","FT max V strip; module; max Deconv combo (isamp = 2nd sample of combo)",nmodules,-0.5,nmodules-0.5,7,-0.5,6.5);

  TH2D *hUtimeMaxStrip_vs_modFT = new TH2D("hUtimeMaxStrip_vs_modFT","FT max U strip; module; Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,0,150);
  TH2D *hVtimeMaxStrip_vs_modFT = new TH2D("hVtimeMaxStrip_vs_modFT","FT max V strip; module; Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,0,150);

  TH2D *hUtimeMaxStripDeconv_vs_modFT = new TH2D("hUtimeMaxStripDeconv_vs_modFT","FT max U strip; module; Deconv Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,-150,150);
  TH2D *hVtimeMaxStripDeconv_vs_modFT = new TH2D("hVtimeMaxStripDeconv_vs_modFT","FT max V strip; module; Deconv Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,-150,150);

  TH2D *hUtime_vs_modFT = new TH2D("hUtime_vs_modFT","FT U cluster mean; module; Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,0,150);
  TH2D *hVtime_vs_modFT = new TH2D("hVtime_vs_modFT","FT V cluster mean; module; Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,0,150);

  TH2D *hUtimeDeconv_vs_modFT = new TH2D("hUtimeDeconv_vs_modFT","FT U cluster mean; module; Deconv Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,-150,150);
  TH2D *hVtimeDeconv_vs_modFT = new TH2D("hVtimeDeconv_vs_modFT","FT V cluster mean; module; Deconv Hit time (ns)",nmodules,-0.5,nmodules-0.5,300,-150,150);

  TH2D *hTavg_corr_vs_modFT = new TH2D("hTavg_corr_vs_modFT", "FT hits, average corrected time (ns); module; T_{corr} (ns)",nmodules,-0.5,nmodules-0.5,300,-150,150);

  TH2D *hdeltat_vs_modFT = new TH2D("hdeltat_vs_modFT", "FT hits; module; #Deltat (ns)", nmodules,-0.5,nmodules-0.5,300,-75,75);
  TH2D *hdeltatDeconv_vs_modFT = new TH2D("hdeltatDeconv_vs_modFT", "FT hits; module; Deconv #Deltat (ns)", nmodules,-0.5,nmodules-0.5,300,-150,150);


  //Back tracker histos:
 
  TH2D *hADCmaxsampU_vs_modFPP = new TH2D("hADCmaxsampU_vs_modFPP", "FPP U strips; module; Max U strip max sample ADC", nmodules_back, -0.5, nmodules_back-0.5, 400,0,4000);
  TH2D *hADCmaxsampV_vs_modFPP = new TH2D("hADCmaxsampV_vs_modFPP", "FPP V strips; module; Max V strip max sample ADC", nmodules_back, -0.5, nmodules_back-0.5, 400,0,4000);
  TH2D *hADCmaxstripU_vs_modFPP = new TH2D("hADCmaxstripU_vs_modFPP", "FPP U strips; module; Max U strip ADC sum", nmodules_back, -0.5, nmodules_back-0.5, 1000,0,10000);
  TH2D *hADCmaxstripV_vs_modFPP = new TH2D("hADCmaxstripV_vs_modFPP", "FPP V strips; module; Max V strip ADC sum", nmodules_back, -0.5, nmodules_back-0.5, 1000,0,10000);
  TH2D *hADCU_vs_modFPP = new TH2D("hADCU_vs_modFPP", "FPP U cluster sum; module; Cluster ADC sum", nmodules_back, -0.5, nmodules_back-0.5, 1000,0,20000);
  TH2D *hADCV_vs_modFPP = new TH2D("hADCV_vs_modFPP", "FPP V cluster sum; module; Cluster ADC sum", nmodules_back, -0.5, nmodules_back-0.5, 1000,0,20000);
  

  TH2D *hADCavg_vs_modFPP = new TH2D("hADCavg_vs_modFPP", "FPP hit cluster sums; module; (ADCU+ADCV)/2", nmodules_back,-0.5,nmodules_back-0.5, 1000,0,20000);
  TH2D *hADCasym_vs_modFPP = new TH2D("hADCasym_vs_modFPP", "FPP hit cluster sum asymmetry; module; (ADCU-ADCV)/(ADCU+ADCV)",nmodules_back,-0.5,nmodules_back-0.5,250,-1.0,1.0);
 
  TH2D *hDeconvADCmaxcomboUstrip_vs_modFPP = new TH2D("hDeconvADCmaxcomboUstrip_vs_modFPP","FPP max U strip; module; Max Deconv. two-sample combo",nmodules_back,-0.5,nmodules_back-0.5,400,0,4000);
  TH2D *hDeconvADCmaxcomboVstrip_vs_modFPP = new TH2D("hDeconvADCmaxcomboVstrip_vs_modFPP","FPP max V strip; module; Max Deconv. two-sample combo",nmodules_back,-0.5,nmodules_back-0.5,400,0,4000);

  TH2D *hDeconvADCmaxcomboUclust_vs_modFPP = new TH2D("hDeconvADCmaxcomboUclust_vs_modFPP","FPP U cluster sum; module; Max Deconv. two-sample combo",nmodules_back,-0.5,nmodules_back-0.5,1000,0,10000);
  TH2D *hDeconvADCmaxcomboVclust_vs_modFPP = new TH2D("hDeconvADCmaxcomboVclust_vs_modFPP","FPP V cluster sum; module; Max Deconv. two-sample combo",nmodules_back,-0.5,nmodules_back-0.5,1000,0,10000);

  TH2D *hADCasymDeconv_vs_modFPP = new TH2D("hADCasymDeconv_vs_modFPP","FPP Deconv ADC asym;module;(ADCU-ADCV)/(ADCU+ADCV)",nmodules_back,-0.5,nmodules_back-0.5,250,-1,1);

  TH2D *hisampmaxUstrip_vs_modFPP = new TH2D("hisampmaxUstrip_vs_modFPP","FPP max U strip; module; max time sample",nmodules_back,-0.5,nmodules_back-0.5,6,-0.5,5.5);
  TH2D *hisampmaxVstrip_vs_modFPP = new TH2D("hisampmaxVstrip_vs_modFPP","FPP max V strip; module; max time sample",nmodules_back,-0.5,nmodules_back-0.5,6,-0.5,5.5);

  TH2D *hicombomaxUstrip_vs_modFPP = new TH2D("hicombomaxUstrip_vs_modFPP","FPP max U strip; module; max Deconv combo (isamp = 2nd sample of combo)",nmodules_back,-0.5,nmodules_back-0.5,7,-0.5,6.5);
  TH2D *hicombomaxVstrip_vs_modFPP = new TH2D("hicombomaxVstrip_vs_modFPP","FPP max V strip; module; max Deconv combo (isamp = 2nd sample of combo)",nmodules_back,-0.5,nmodules_back-0.5,7,-0.5,6.5);

  TH2D *hUtimeMaxStrip_vs_modFPP = new TH2D("hUtimeMaxStrip_vs_modFPP","FPP max U strip; module; Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,0,150);
  TH2D *hVtimeMaxStrip_vs_modFPP = new TH2D("hVtimeMaxStrip_vs_modFPP","FPP max V strip; module; Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,0,150);

  TH2D *hUtimeMaxStripDeconv_vs_modFPP = new TH2D("hUtimeMaxStripDeconv_vs_modFPP","FPP max U strip; module; Deconv Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,-150,150);
  TH2D *hVtimeMaxStripDeconv_vs_modFPP = new TH2D("hVtimeMaxStripDeconv_vs_modFPP","FPP max V strip; module; Deconv Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,-150,150);

  TH2D *hUtime_vs_modFPP = new TH2D("hUtime_vs_modFPP","FPP U cluster mean; module; Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,0,150);
  TH2D *hVtime_vs_modFPP = new TH2D("hVtime_vs_modFPP","FPP V cluster mean; module; Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,0,150);

  TH2D *hUtimeDeconv_vs_modFPP = new TH2D("hUtimeDeconv_vs_modFPP","FPP U cluster mean; module; Deconv Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,-150,150);
  TH2D *hVtimeDeconv_vs_modFPP = new TH2D("hVtimeDeconv_vs_modFPP","FPP V cluster mean; module; Deconv Hit time (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,-150,150);

  TH2D *hTavg_corr_vs_modFPP = new TH2D("hTavg_corr_vs_modFPP", "FPP hits, average corrected time (ns); module; T_{corr} (ns)",nmodules_back,-0.5,nmodules_back-0.5,300,-150,150);

  TH2D *hdeltat_vs_modFPP = new TH2D("hdeltat_vs_modFPP", "FPP hits; module; #Deltat (ns)", nmodules_back,-0.5,nmodules_back-0.5,300,-75,75);
  TH2D *hdeltatDeconv_vs_modFPP = new TH2D("hdeltatDeconv_vs_modFPP", "FPP hits; module; Deconv #Deltat (ns)", nmodules_back,-0.5,nmodules_back-0.5,300,-150,150);
  


  long nevent=0;

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  int treenum=0, currenttreenum=0;

  while( C->GetEntry( nevent++ ) && nevent){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    if( nevent % 1000 == 0 ) cout << nevent << endl;
    
    bool passedcut = GlobalCut->EvalInstance(0) != 0;

    // cout << "passed cut, ntracks, ngoodhits = " << passedcut << ", "
    // 	 << ntracks << ", " << ngoodhits << endl;


    
    if( passedcut && int(ntracks) >= 1){

      //Grab trigger TDC time: 
      // double ttrig = 0.0;
      // for( int itrig=0; itrig<Ntrig; itrig++ ){
      // 	if( int(trig_elemID[itrig]) == 5 ){
      // 	  ttrig = trig_tdc[itrig];
      // 	}
      // }
      
      for( int ihit=0; ihit<int(ngoodhits_FT); ihit++ ){
	if( int(trackindexFT[ihit]) == 0 && nstripuFT[ihit]>1&&nstripvFT[ihit]>1 ){
	  
	  hADCmaxsampU_vs_modFT->Fill(modFT[ihit],ADCmaxsampU_FT[ihit]);
	  hADCmaxsampV_vs_modFT->Fill(modFT[ihit],ADCmaxsampV_FT[ihit]);
	  hADCmaxstripU_vs_modFT->Fill(modFT[ihit],ADCmaxstripU_FT[ihit]);
	  hADCmaxstripV_vs_modFT->Fill(modFT[ihit],ADCmaxstripV_FT[ihit]);
	  hADCU_vs_modFT->Fill( modFT[ihit], ADCU_FT[ihit] );
	  hADCV_vs_modFT->Fill( modFT[ihit], ADCV_FT[ihit] );

	  hADCavg_vs_modFT->Fill( modFT[ihit], ADCavg_FT[ihit]);
	  hADCasym_vs_modFT->Fill( modFT[ihit], ADCasym_FT[ihit]);

	  hDeconvADCmaxcomboUstrip_vs_modFT->Fill( modFT[ihit], DeconvADCmaxcomboU_FT[ihit] );
	  hDeconvADCmaxcomboVstrip_vs_modFT->Fill( modFT[ihit], DeconvADCmaxcomboV_FT[ihit] );

	  hDeconvADCmaxcomboUclust_vs_modFT->Fill( modFT[ihit], DeconvADCU_FT[ihit] );
	  hDeconvADCmaxcomboVclust_vs_modFT->Fill( modFT[ihit], DeconvADCV_FT[ihit] );

	  hADCasymDeconv_vs_modFT->Fill( modFT[ihit], ADCasymDeconv_FT[ihit] );

	  hisampmaxUstrip_vs_modFT->Fill( modFT[ihit], isampmaxUstrip_FT[ihit] );
	  hisampmaxVstrip_vs_modFT->Fill( modFT[ihit], isampmaxVstrip_FT[ihit] );
	  
	  hicombomaxUstrip_vs_modFT->Fill( modFT[ihit], icombomaxUstripDeconv_FT[ihit] );
	  hicombomaxVstrip_vs_modFT->Fill( modFT[ihit], icombomaxVstripDeconv_FT[ihit] );


	  hUtimeMaxStrip_vs_modFT->Fill( modFT[ihit], UtimeMaxStrip_FT[ihit] );
	  hVtimeMaxStrip_vs_modFT->Fill( modFT[ihit], VtimeMaxStrip_FT[ihit] );

	  hUtimeMaxStripDeconv_vs_modFT->Fill( modFT[ihit], UtimeMaxStripDeconv_FT[ihit] );
	  hVtimeMaxStripDeconv_vs_modFT->Fill( modFT[ihit], VtimeMaxStripDeconv_FT[ihit] );

	  hUtime_vs_modFT->Fill( modFT[ihit], Utime_FT[ihit] );
	  hVtime_vs_modFT->Fill( modFT[ihit], Vtime_FT[ihit] );

	  hUtimeDeconv_vs_modFT->Fill( modFT[ihit], UtimeDeconv_FT[ihit] );
	  hVtimeDeconv_vs_modFT->Fill( modFT[ihit], VtimeDeconv_FT[ihit] );
	  
	  hTavg_corr_vs_modFT->Fill( modFT[ihit], HitTavg_corr_FT[ihit] );
	  hdeltat_vs_modFT->Fill( modFT[ihit], deltat_FT[ihit] );
	  hdeltatDeconv_vs_modFT->Fill( modFT[ihit], deltatDeconv_FT[ihit] );
	  
	}
      }

      if( ntracks_back > 0 && sclose[0] < 0.015 ){
      
	for( int ihit=0; ihit<int(ngoodhits_FPP); ihit++ ){
	  if( int(trackindexFPP[ihit]) == 0 && nstripuFPP[ihit]>1&&nstripvFPP[ihit]>1 ){
	  
	    hADCmaxsampU_vs_modFPP->Fill(modFPP[ihit],ADCmaxsampU_FPP[ihit]);
	    hADCmaxsampV_vs_modFPP->Fill(modFPP[ihit],ADCmaxsampV_FPP[ihit]);
	    hADCmaxstripU_vs_modFPP->Fill(modFPP[ihit],ADCmaxstripU_FPP[ihit]);
	    hADCmaxstripV_vs_modFPP->Fill(modFPP[ihit],ADCmaxstripV_FPP[ihit]);
	    hADCU_vs_modFPP->Fill( modFPP[ihit], ADCU_FPP[ihit] );
	    hADCV_vs_modFPP->Fill( modFPP[ihit], ADCV_FPP[ihit] );

	    hADCavg_vs_modFPP->Fill( modFPP[ihit], ADCavg_FPP[ihit]);
	    hADCasym_vs_modFPP->Fill( modFPP[ihit], ADCasym_FPP[ihit]);

	    hDeconvADCmaxcomboUstrip_vs_modFPP->Fill( modFPP[ihit], DeconvADCmaxcomboU_FPP[ihit] );
	    hDeconvADCmaxcomboVstrip_vs_modFPP->Fill( modFPP[ihit], DeconvADCmaxcomboV_FPP[ihit] );

	    hDeconvADCmaxcomboUclust_vs_modFPP->Fill( modFPP[ihit], DeconvADCU_FPP[ihit] );
	    hDeconvADCmaxcomboVclust_vs_modFPP->Fill( modFPP[ihit], DeconvADCV_FPP[ihit] );

	    hADCasymDeconv_vs_modFPP->Fill( modFPP[ihit], ADCasymDeconv_FPP[ihit] );

	    hisampmaxUstrip_vs_modFPP->Fill( modFPP[ihit], isampmaxUstrip_FPP[ihit] );
	    hisampmaxVstrip_vs_modFPP->Fill( modFPP[ihit], isampmaxVstrip_FPP[ihit] );
	  
	    hicombomaxUstrip_vs_modFPP->Fill( modFPP[ihit], icombomaxUstripDeconv_FPP[ihit] );
	    hicombomaxVstrip_vs_modFPP->Fill( modFPP[ihit], icombomaxVstripDeconv_FPP[ihit] );


	    hUtimeMaxStrip_vs_modFPP->Fill( modFPP[ihit], UtimeMaxStrip_FPP[ihit] );
	    hVtimeMaxStrip_vs_modFPP->Fill( modFPP[ihit], VtimeMaxStrip_FPP[ihit] );

	    hUtimeMaxStripDeconv_vs_modFPP->Fill( modFPP[ihit], UtimeMaxStripDeconv_FPP[ihit] );
	    hVtimeMaxStripDeconv_vs_modFPP->Fill( modFPP[ihit], VtimeMaxStripDeconv_FPP[ihit] );

	    hUtime_vs_modFPP->Fill( modFPP[ihit], Utime_FPP[ihit] );
	    hVtime_vs_modFPP->Fill( modFPP[ihit], Vtime_FPP[ihit] );

	    hUtimeDeconv_vs_modFPP->Fill( modFPP[ihit], UtimeDeconv_FPP[ihit] );
	    hVtimeDeconv_vs_modFPP->Fill( modFPP[ihit], VtimeDeconv_FPP[ihit] );
	  
	    hTavg_corr_vs_modFPP->Fill( modFPP[ihit], HitTavg_corr_FPP[ihit] );
	    hdeltat_vs_modFPP->Fill( modFPP[ihit], deltat_FPP[ihit] );
	    hdeltatDeconv_vs_modFPP->Fill( modFPP[ihit], deltatDeconv_FPP[ihit] );
	  
	  }
	}
      }
    }
  }

  //// Here we create a pdf which will be filled in the following loop ////

  TString outfilepdf = outfilename;
  outfilepdf.ReplaceAll(".root",".pdf");

  TString outfilepdf_time = outfilename;
  outfilepdf_time.ReplaceAll(".root","_timecuts.pdf");

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  // c1->Divide(2,3);
  // c1->cd(1);
  // DrawHist(hdxfcp,{hdxfcp->GetMean(),hdxfcp->GetRMS()*4.5});
  // c1->cd(2);
  // DrawHist(hdyfcp,{hdyfcp->GetMean(),hdyfcp->GetRMS()*4.5});
  // c1->cd(3);
  // DrawHist(hdxbcp,{hdxbcp->GetMean(),hdxbcp->GetRMS()*4.5});
  // c1->cd(4);
  // DrawHist(hdybcp,{hdybcp->GetMean(),hdybcp->GetRMS()*4.5});
  // c1->cd(5);
  // DrawHist(hdthcp,{hdthcp->GetMean(),hdthcp->GetRMS()*4.5});
  // c1->cd(6);
  // DrawHist(hdphcp,{hdphcp->GetMean(),hdphcp->GetRMS()*4.5});

  //  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page

  ///////////////////////////////////////////////////////////////////
  double nsigma = nsig_tstrip;
  
  double maxstrip_t0[nmodules+nmodules_back][2];
  double maxstrip_tsigma[nmodules+nmodules_back][2];
  //double maxstrip_tcut[nmodules+nmodules_back][2];

  double maxstrip_t0_deconv[nmodules+nmodules_back][2];
  double maxstrip_tsigma_deconv[nmodules+nmodules_back][2];
  //double maxstrip_tcut_deconv[nmodules+nmodules_back][2];

  //  double maxstrip_t0_fit[nmodules+nmodules_back][2];
  // double maxstrip_tsigma_fit[nmodules+nmodules_back][2];
  
  double tmean[nmodules+nmodules_back][4];
  double tsigma[nmodules+nmodules_back][4];

  double sigma_tcorr[nmodules+nmodules_back];
  double sigma_dt[nmodules+nmodules_back], sigma_dtdeconv[nmodules+nmodules_back];

  double thresh_sample[nmodules+nmodules_back];
  double thresh_strip[nmodules+nmodules_back];
  double thresh_cluster[nmodules+nmodules_back];

  double thresh_maxcombo_deconv[nmodules+nmodules_back];
  double thresh_clustersum_deconv[nmodules+nmodules_back];
  
  double ADCasym_sigma[nmodules+nmodules_back];
  
  TString fname_dbFT = outfilename;
  fname_dbFT.ReplaceAll(".root","_FT.dat");
  
  ofstream dbFT(fname_dbFT.Data());

  TString fname_dbFPP = outfilename;
  fname_dbFPP.ReplaceAll(".root","_FPP.dat");

  ofstream dbFPP(fname_dbFPP.Data());
  
  // ADC (software) thresholds, timing offsets and cuts:

  //  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  //For each module, we want the following plots: 
  // max U and V strip times (2 plots)
  // max U and V strip deconvoluted times (2 plots)
  // mean U and V cluster times (2 plots)
  // mean U and V cluster deconvoluted times (2 plots)
  // average corrected time (1 plot)
  // regular and deconvoluted hit time difference (2 plots)
  // U and V strip ADC max sample (2 plots)
  // U and V strip ADC strip sum (2 plots)
  // U and V strip ADC cluster sums (2 plots)
  // U and V strip deconv ADC max two-sample combo (2 plots)
  // U and V cluster-summed deconv ADC max two-sample combo ( 2 plots)
  // ADC asymmetry (regular and deconv) 2 plots
  // 12 plots for ADC-related stuff:
  // I suppose we could also do correlation coefficient plots, but
  // we don't do that yet:
  // In total then, we have
  //let's take default values for modules without enough statistics:
  for( int imod=0; imod<nmodules; imod++ ){

    for( int i=0; i<2; i++ ){
      
      maxstrip_t0[imod][i] = 85.0;
      maxstrip_tsigma[imod][i] = 10.0;
      maxstrip_t0_deconv[imod][i] = 40.0;
      maxstrip_tsigma_deconv[imod][i] = 20.0;

      tmean[imod][i] = 85.0;
      tsigma[imod][i] = 10.0;
      tmean[imod][i+2] = 40.0;
      tsigma[imod][i+2] = 20.0;
    }

    sigma_tcorr[imod] = 8.0;
    sigma_dt[imod] = 5.0;
    sigma_dtdeconv[imod] = 15.0;

    thresh_sample[imod] = 40.0; //default
    thresh_strip[imod] = 120.0; //default
    thresh_cluster[imod] = 240.0; //default
    thresh_maxcombo_deconv[imod] = 40.0; //default
    thresh_clustersum_deconv[imod] = 60.0; //default
    //Start with timing:

    c1->Clear();
    c1->Divide(4,3,0.001,0.001);

    TString histname;
    
    TH1D *htUtemp, *htVtemp;
    
    htUtemp = hUtimeMaxStrip_vs_modFT->ProjectionY(histname.Format("hUtimeMaxStrip_FTmod%d",imod),imod+1,imod+1);
    htVtemp = hVtimeMaxStrip_vs_modFT->ProjectionY(histname.Format("hVtimeMaxStrip_FTmod%d",imod),imod+1,imod+1);						   
						   
    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(1);
      FitGaus_FWHM( htUtemp, 0.3 );
      //htUtemp->DrawCopy();
      htUtemp->SetTitle(histname.Format("FT mod %d; max U strip time (ns);",imod));
      htUtemp->Draw();
      c1->cd(2);
      FitGaus_FWHM( htVtemp, 0.3 );
      htVtemp->SetTitle(histname.Format("FT mod %d; max V strip time (ns);",imod));
      htVtemp->Draw();

      maxstrip_t0[imod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma[imod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      maxstrip_t0[imod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma[imod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      
      
    }

    htUtemp = hUtime_vs_modFT->ProjectionY(histname.Format("hUtime_FTmod%d",imod),imod+1,imod+1);
    htVtemp = hVtime_vs_modFT->ProjectionY(histname.Format("hVtime_FTmod%d",imod),imod+1,imod+1);

    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(3);
      FitGaus_FWHM( htUtemp,0.3 );
      htUtemp->SetTitle(histname.Format("FT mod %d; U cluster mean time (ns);",imod));
      htUtemp->Draw();
      c1->cd(4);
      FitGaus_FWHM( htVtemp,0.3 );
      htVtemp->SetTitle(histname.Format("FT mod %d; V cluster mean time (ns);",imod));
      htVtemp->Draw();

      tmean[imod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[imod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      tmean[imod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[imod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }

    htUtemp = hUtimeMaxStripDeconv_vs_modFT->ProjectionY(histname.Format("hUtimeMaxStripDeconv_FTmod%d",imod),imod+1,imod+1);
    htVtemp = hVtimeMaxStripDeconv_vs_modFT->ProjectionY(histname.Format("hVtimeMaxStripDeconv_FTmod%d",imod),imod+1,imod+1);
    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(5);
      FitGaus_FWHM( htUtemp,0.3);
      htUtemp->SetTitle(histname.Format("FT mod %d; Max U strip deconv time (ns);",imod));
      htUtemp->Draw();
      c1->cd(6);
      FitGaus_FWHM( htVtemp,0.3);
      htVtemp->SetTitle(histname.Format("FT mod %d; Max V strip deconv time (ns);",imod));
      htVtemp->Draw();

      maxstrip_t0_deconv[imod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma_deconv[imod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      maxstrip_t0_deconv[imod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma_deconv[imod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      
    }

    htUtemp = hUtimeDeconv_vs_modFT->ProjectionY(histname.Format("hUtimeDeconv_FTmod%d",imod),imod+1,imod+1);
    htVtemp = hVtimeDeconv_vs_modFT->ProjectionY(histname.Format("hVtimeDeconv_FTmod%d",imod),imod+1,imod+1);
    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(7);
      FitGaus_FWHM( htUtemp,0.3);
      htUtemp->SetTitle(histname.Format("FT mod %d; Deconv. U cluster time (ns);",imod));
      htUtemp->Draw();
      c1->cd(8);
      FitGaus_FWHM( htVtemp,0.3);
      htVtemp->SetTitle(histname.Format("FT mod %d; Deconv V cluster time (ns);",imod));
      htVtemp->Draw();

      tmean[imod][2] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[imod][2] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      tmean[imod][3] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[imod][3] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }

    TH1D *htemp = hTavg_corr_vs_modFT->ProjectionY(histname.Format("hTavg_corr_FTmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() >= 300 ){
      c1->cd(9);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FT mod %d; Hit average corrected time (ns);",imod));
      htemp->Draw();

      //We don't care about the mean for this variable
      sigma_tcorr[imod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }

    htemp = hdeltat_vs_modFT->ProjectionY(histname.Format("hdeltat_FTmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() >= 300 ){
      c1->cd(10);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FT mod %d; #Deltat(ns);",imod));
      htemp->Draw();

      sigma_dt[imod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      
    }

    htemp = hdeltatDeconv_vs_modFT->ProjectionY(histname.Format("hdeltatDeconv_FTmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() >= 300 ){
      c1->cd(11);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FT mod %d; #Deltat_{Deconv}(ns);",imod));
      htemp->Draw();

      sigma_dtdeconv[imod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }
    
    TString pdffilename = outfilename;
    pdffilename.ReplaceAll(".root",".pdf");

    if( imod == 0 ) pdffilename += "(";

    //if( imod+1 == nmodules ) pdffilename += ")";
    c1->Print(pdffilename,"pdf");

    c1->Clear();

    c1->Divide(4,3,.001,.001);
    
    //Second page: ADC stuff:
    
    TH1D *htempU = hADCmaxsampU_vs_modFT->ProjectionY(histname.Format("hADCmaxsampU_FTmod%d",imod),imod+1,imod+1);
    TH1D *htempV = hADCmaxsampV_vs_modFT->ProjectionY(histname.Format("hADCmaxsampV_FTmod%d",imod),imod+1,imod+1);

    double MPV_U,Sigma_U;
    double MPV_V,Sigma_V;
    
    if( htempU->GetEntries() >= 300 ){
      c1->cd(1)->SetLogy();
      htempU->SetTitle(histname.Format("FT mod %d; Max U strip max ADC sample;",imod));
      htempU->Fit("landau","qS","",0,1000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(2)->SetLogy();
      htempV->SetTitle(histname.Format("FT mod %d; Max V strip max ADC sample;",imod));
      htempV->Fit("landau","qS","",0,1000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_sample[imod] = std::min(threshU,threshV);

      std::cout << "FT mod " << imod << ", (threshU,threshV,thresh_samp)=(" << threshU << ", "
		<< threshV << ", " << thresh_sample[imod] << ")" << std::endl;
    }

    htempU = hADCmaxstripU_vs_modFT->ProjectionY(histname.Format("hADCmaxstripU_FTmod%d",imod),imod+1,imod+1);
    htempV = hADCmaxstripV_vs_modFT->ProjectionY(histname.Format("hADCmaxstripV_FTmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(3)->SetLogy();

      htempU->SetTitle(histname.Format("FT mod %d; Max U strip ADC sum;",imod));
      htempU->Fit("landau","qS","",0,3000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(4)->SetLogy();
      htempV->SetTitle(histname.Format("FT mod %d; Max V strip ADC sum;",imod));
      htempV->Fit("landau","qS","",0,3000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_strip[imod] = std::min(threshU,threshV);

      std::cout << "FT mod " << imod << ", (threshU,threshV,thresh_strip)=(" << threshU << ", "
       		<< threshV << ", " << thresh_strip[imod] << ")" << std::endl;
      
    }

    htempU = hADCU_vs_modFT->ProjectionY(histname.Format("hADCU_FTmod%d",imod),imod+1,imod+1);
    htempV = hADCV_vs_modFT->ProjectionY(histname.Format("hADCV_FTmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(5)->SetLogy();

      htempU->SetTitle(histname.Format("FT mod %d; U cluster ADC sum;",imod));
      htempU->Fit("landau","qS","",0,7500.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(6)->SetLogy();
      htempV->SetTitle(histname.Format("FT mod %d; V cluster ADC sum;",imod));
      htempV->Fit("landau","qS","",0,7500.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_cluster[imod] = std::min(threshU,threshV);

      std::cout << "FT mod " << imod << ", (threshU,threshV,thresh_cluster)=(" << threshU << ", "
       		<< threshV << ", " << thresh_cluster[imod] << ")" << std::endl;
      
    }


    htempU = hDeconvADCmaxcomboUstrip_vs_modFT->ProjectionY(histname.Format("hDeconvADCUmaxstrip_FTmod%d",imod),imod+1,imod+1);
    htempV = hDeconvADCmaxcomboVstrip_vs_modFT->ProjectionY(histname.Format("hDeconvADCVmaxstrip_FTmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(7)->SetLogy();

      htempU->SetTitle(histname.Format("FT mod %d; max U strip max deconv combo;",imod));
      htempU->Fit("landau","qS","",0,1000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(8)->SetLogy();
      htempV->SetTitle(histname.Format("FT mod %d; max V strip max deconv combo;",imod));
      htempV->Fit("landau","qS","",0,1000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_maxcombo_deconv[imod] = std::min(threshU,threshV);

      std::cout << "FT mod " << imod << ", (threshU,threshV,thresh_maxcombo_deconv)=(" << threshU << ", "
       		<< threshV << ", " << thresh_maxcombo_deconv[imod] << ")" << std::endl;
      
    }

    htempU = hDeconvADCmaxcomboUclust_vs_modFT->ProjectionY(histname.Format("hDeconvADCUclust_FTmod%d",imod),imod+1,imod+1);
    htempV = hDeconvADCmaxcomboVclust_vs_modFT->ProjectionY(histname.Format("hDeconvADCVclust_FTmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(9)->SetLogy();

      htempU->SetTitle(histname.Format("FT mod %d; U cluster deconv ADC sum;",imod));
      htempU->Fit("landau","qS","",0,3000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(10)->SetLogy();
      htempV->SetTitle(histname.Format("FT mod %d; V cluster deconv ADC sum;",imod));
      htempV->Fit("landau","qS","",0,3000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_clustersum_deconv[imod] = std::min(threshU,threshV);

      std::cout << "FT mod " << imod << ", (threshU,threshV,thresh_clust_deconv)=(" << threshU << ", "
       		<< threshV << ", " << thresh_clustersum_deconv[imod] << ")" << std::endl;
      
    }

    htemp = hADCasym_vs_modFT->ProjectionY(histname.Format("hADCasym_FTmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() > 300 ){
      c1->cd(11);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FT mod %d; ADC asym;",imod));
      htemp->Draw();

      ADCasym_sigma[imod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }
    
    htemp = hADCasymDeconv_vs_modFT->ProjectionY(histname.Format("hADCasymDeconv_FTmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() > 300 ){
      c1->cd(12);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FT mod %d; Deconv ADC asym;",imod));
      htemp->Draw();
    }
    
    c1->Print(pdffilename,"pdf");
  }


  for( int imod=0; imod<nmodules_back; imod++ ){

    int mod = imod + nmodules; //"global" module number
    for( int i=0; i<2; i++ ){
      
      maxstrip_t0[mod][i] = 85.0;
      maxstrip_tsigma[mod][i] = 10.0;
      maxstrip_t0_deconv[mod][i] = 40.0;
      maxstrip_tsigma[mod][i] = 20.0;

      tmean[mod][i] = 85.0;
      tsigma[mod][i] = 10.0;
      tmean[mod][i+2] = 40.0;
      tsigma[mod][i+2] = 20.0;
    }

    sigma_tcorr[mod] = 8.0;
    sigma_dt[mod] = 5.0;
    sigma_dtdeconv[mod] = 15.0;

    thresh_sample[mod] = 40.0; //default
    thresh_strip[mod] = 120.0; //default
    thresh_cluster[mod] = 240.0; //default
    thresh_maxcombo_deconv[mod] = 40.0; //default
    thresh_clustersum_deconv[mod] = 60.0; //default
    
    //Start with timing:

    
    
    c1->Clear();
    c1->Divide(4,3,0.001,0.001);

    TString histname;
    
    TH1D *htUtemp, *htVtemp;
    
    htUtemp = hUtimeMaxStrip_vs_modFPP->ProjectionY(histname.Format("hUtimeMaxStrip_FPPmod%d",imod),imod+1,imod+1);
    htVtemp = hVtimeMaxStrip_vs_modFPP->ProjectionY(histname.Format("hVtimeMaxStrip_FPPmod%d",imod),imod+1,imod+1);						   
						   
    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(1);
      FitGaus_FWHM( htUtemp, 0.3 );
      //htUtemp->DrawCopy();
      htUtemp->SetTitle(histname.Format("FPP mod %d; max U strip time (ns);",imod));
      htUtemp->Draw();
      c1->cd(2);
      FitGaus_FWHM( htVtemp, 0.3 );
      htVtemp->SetTitle(histname.Format("FPP mod %d; max V strip time (ns);",imod));
      htVtemp->Draw();

      maxstrip_t0[mod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma[mod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      maxstrip_t0[mod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma[mod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      
      
    }

    htUtemp = hUtime_vs_modFPP->ProjectionY(histname.Format("hUtime_FPPmod%d",imod),imod+1,imod+1);
    htVtemp = hVtime_vs_modFPP->ProjectionY(histname.Format("hVtime_FPPmod%d",imod),imod+1,imod+1);

    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(3);
      FitGaus_FWHM( htUtemp,0.3 );
      htUtemp->SetTitle(histname.Format("FPP mod %d; U cluster mean time (ns);",imod));
      htUtemp->Draw();
      c1->cd(4);
      FitGaus_FWHM( htVtemp,0.3 );
      htVtemp->SetTitle(histname.Format("FPP mod %d; V cluster mean time (ns);",imod));
      htVtemp->Draw();

      tmean[mod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[mod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      tmean[mod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[mod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }

    htUtemp = hUtimeMaxStripDeconv_vs_modFPP->ProjectionY(histname.Format("hUtimeMaxStripDeconv_FPPmod%d",imod),imod+1,imod+1);
    htVtemp = hVtimeMaxStripDeconv_vs_modFPP->ProjectionY(histname.Format("hVtimeMaxStripDeconv_FPPmod%d",imod),imod+1,imod+1);
    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(5);
      FitGaus_FWHM( htUtemp,0.3);
      htUtemp->SetTitle(histname.Format("FPP mod %d; Max U strip deconv time (ns);",imod));
      htUtemp->Draw();
      c1->cd(6);
      FitGaus_FWHM( htVtemp,0.3);
      htVtemp->SetTitle(histname.Format("FPP mod %d; Max V strip deconv time (ns);",imod));
      htVtemp->Draw();

      maxstrip_t0_deconv[mod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma_deconv[mod][0] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      maxstrip_t0_deconv[mod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      maxstrip_tsigma_deconv[mod][1] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      
    }

    htUtemp = hUtimeDeconv_vs_modFPP->ProjectionY(histname.Format("hUtimeDeconv_FPPmod%d",imod),imod+1,imod+1);
    htVtemp = hVtimeDeconv_vs_modFPP->ProjectionY(histname.Format("hVtimeDeconv_FPPmod%d",imod),imod+1,imod+1);
    if( htUtemp->GetEntries() >= 300 ){
      c1->cd(7);
      FitGaus_FWHM( htUtemp,0.3);
      htUtemp->SetTitle(histname.Format("FPP mod %d; Deconv. U cluster time (ns);",imod));
      htUtemp->Draw();
      c1->cd(8);
      FitGaus_FWHM( htVtemp,0.3);
      htVtemp->SetTitle(histname.Format("FPP mod %d; Deconv V cluster time (ns);",imod));
      htVtemp->Draw();

      tmean[mod][2] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[mod][2] = ( (TF1*) (htUtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");

      tmean[mod][3] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      tsigma[mod][3] = ( (TF1*) (htVtemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }

    TH1D *htemp = hTavg_corr_vs_modFPP->ProjectionY(histname.Format("hTavg_corr_FPPmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() >= 300 ){
      c1->cd(9);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FPP mod %d; Hit average corrected time (ns);",imod));
      htemp->Draw();

      //We don't care about the mean for this variable
      sigma_tcorr[mod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }

    htemp = hdeltat_vs_modFPP->ProjectionY(histname.Format("hdeltat_FPPmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() >= 300 ){
      c1->cd(10);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FPP mod %d; #Deltat(ns);",imod));
      htemp->Draw();

      sigma_dt[mod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      
    }

    htemp = hdeltatDeconv_vs_modFPP->ProjectionY(histname.Format("hdeltatDeconv_FPPmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() >= 300 ){
      c1->cd(11);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FPP mod %d; #Deltat_{Deconv}(ns);",imod));
      htemp->Draw();

      sigma_dtdeconv[mod] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }
    
    TString pdffilename = outfilename;
    pdffilename.ReplaceAll(".root",".pdf");

    if( imod == 0 ) pdffilename += "(";

    //if( imod+1 == nmodules_back ) pdffilename += ")";
    c1->Print(pdffilename,"pdf");


    c1->Clear();
    
    c1->Divide(4,3,.001,.001);
    
    //Second page: ADC stuff:
    
    TH1D *htempU = hADCmaxsampU_vs_modFPP->ProjectionY(histname.Format("hADCmaxsampU_FPPmod%d",imod),imod+1,imod+1);
    TH1D *htempV = hADCmaxsampV_vs_modFPP->ProjectionY(histname.Format("hADCmaxsampV_FPPmod%d",imod),imod+1,imod+1);

    double MPV_U,Sigma_U;
    double MPV_V,Sigma_V;
    
    if( htempU->GetEntries() >= 300 ){
      c1->cd(1)->SetLogy();
      htempU->SetTitle(histname.Format("FPP mod %d; Max U strip max ADC sample;",imod));
      htempU->Fit("landau","qS","",0,1000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(2)->SetLogy();
      htempV->SetTitle(histname.Format("FPP mod %d; Max V strip max ADC sample;",imod));
      htempV->Fit("landau","qS","",0,1000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_sample[imod+nmodules] = std::min(threshU,threshV);

      std::cout << "FPP mod " << imod << ", (threshU,threshV,thresh_samp)=(" << threshU << ", "
		<< threshV << ", " << thresh_sample[imod+nmodules] << ")" << std::endl;
    }

    htempU = hADCmaxstripU_vs_modFPP->ProjectionY(histname.Format("hADCmaxstripU_FPPmod%d",imod),imod+1,imod+1);
    htempV = hADCmaxstripV_vs_modFPP->ProjectionY(histname.Format("hADCmaxstripV_FPPmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(3)->SetLogy();

      htempU->SetTitle(histname.Format("FPP mod %d; Max U strip ADC sum;",imod));
      htempU->Fit("landau","qS","",0,3000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(4)->SetLogy();
      htempV->SetTitle(histname.Format("FPP mod %d; Max V strip ADC sum;",imod));
      htempV->Fit("landau","qS","",0,3000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_strip[imod+nmodules] = std::min(threshU,threshV);

      std::cout << "FPP mod " << imod << ", (threshU,threshV,thresh_strip)=(" << threshU << ", "
       		<< threshV << ", " << thresh_strip[imod+nmodules] << ")" << std::endl;
      
    }

    htempU = hADCU_vs_modFPP->ProjectionY(histname.Format("hADCU_FPPmod%d",imod),imod+1,imod+1);
    htempV = hADCV_vs_modFPP->ProjectionY(histname.Format("hADCV_FPPmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(5)->SetLogy();

      htempU->SetTitle(histname.Format("FPP mod %d; U cluster ADC sum;",imod));
      htempU->Fit("landau","qS","",0,7500.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(6)->SetLogy();
      htempV->SetTitle(histname.Format("FPP mod %d; V cluster ADC sum;",imod));
      htempV->Fit("landau","qS","",0,7500.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_cluster[imod+nmodules] = std::min(threshU,threshV);

      std::cout << "FPP mod " << imod << ", (threshU,threshV,thresh_cluster)=(" << threshU << ", "
       		<< threshV << ", " << thresh_cluster[imod+nmodules] << ")" << std::endl;
      
    }


    htempU = hDeconvADCmaxcomboUstrip_vs_modFPP->ProjectionY(histname.Format("hDeconvADCUmaxstrip_FPPmod%d",imod),imod+1,imod+1);
    htempV = hDeconvADCmaxcomboVstrip_vs_modFPP->ProjectionY(histname.Format("hDeconvADCVmaxstrip_FPPmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(7)->SetLogy();

      htempU->SetTitle(histname.Format("FPP mod %d; max U strip max deconv combo;",imod));
      htempU->Fit("landau","qS","",0,1000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(8)->SetLogy();
      htempV->SetTitle(histname.Format("FPP mod %d; max V strip max deconv combo;",imod));
      htempV->Fit("landau","qS","",0,1000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_maxcombo_deconv[imod+nmodules] = std::min(threshU,threshV);

      std::cout << "FPP mod " << imod << ", (threshU,threshV,thresh_maxcombo_deconv)=(" << threshU << ", "
       		<< threshV << ", " << thresh_maxcombo_deconv[imod+nmodules] << ")" << std::endl;
      
    }

    htempU = hDeconvADCmaxcomboUclust_vs_modFPP->ProjectionY(histname.Format("hDeconvADCUclust_FPPmod%d",imod),imod+1,imod+1);
    htempV = hDeconvADCmaxcomboVclust_vs_modFPP->ProjectionY(histname.Format("hDeconvADCVclust_FPPmod%d",imod),imod+1,imod+1);
    if( htempU->GetEntries() > 300 ){
      c1->cd(9)->SetLogy();

      htempU->SetTitle(histname.Format("FPP mod %d; U cluster deconv ADC sum;",imod));
      htempU->Fit("landau","qS","",0,3000.);
      htempU->Draw();

      MPV_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_U = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      c1->cd(10)->SetLogy();
      htempV->SetTitle(histname.Format("FPP mod %d; V cluster deconv ADC sum;",imod));
      htempV->Fit("landau","qS","",0,3000.);
      htempV->Draw();

      MPV_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      Sigma_V = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("Sigma");

      double threshU = MPV_U + ROOT::Math::landau_quantile( thresh, Sigma_U );
      double threshV = MPV_V + ROOT::Math::landau_quantile( thresh, Sigma_V );
      
      thresh_clustersum_deconv[imod+nmodules] = std::min(threshU,threshV);

      std::cout << "FPP mod " << imod << ", (threshU,threshV,thresh_clust_deconv)=(" << threshU << ", "
       		<< threshV << ", " << thresh_clustersum_deconv[imod+nmodules] << ")" << std::endl;
      
    }

    htemp = hADCasym_vs_modFPP->ProjectionY(histname.Format("hADCasym_FPPmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() > 300 ){
      c1->cd(11);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FPP mod %d; ADC asym;",imod));
      htemp->Draw();

      ADCasym_sigma[imod+nmodules] = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
    }
    
    htemp = hADCasymDeconv_vs_modFPP->ProjectionY(histname.Format("hADCasymDeconv_FPPmod%d",imod),imod+1,imod+1);
    if( htemp->GetEntries() > 300 ){
      c1->cd(12);
      FitGaus_FWHM( htemp, 0.3 );
      htemp->SetTitle(histname.Format("FPP mod %d; Deconv ADC asym;",imod));
      htemp->Draw();
    }

    if( imod+1 == nmodules_back ) pdffilename += ")";
    
    c1->Print(pdffilename,"pdf");
    
  }
  

  vector<TString> parnames;
  parnames.push_back( "maxstrip_t0" );
  parnames.push_back( "maxstrip_tsigma" );
  parnames.push_back( "maxstrip_t0_deconv" );
  parnames.push_back( "maxstrip_tsigma_deconv" );
  parnames.push_back( "HitTimeMean" );
  parnames.push_back( "HitTimeSigma" );
  parnames.push_back( "HitTimeMeanDeconv" );
  parnames.push_back( "HitTimeSigmaDeconv" );
  //parnames.push_back( "deltat_sigma" );
  //parnames.push_back( "deltat_cut" );
  //parnames.push_back( "deltat_sigma_deconv" );
  //parnames.push_back( "deltat_cut_deconv" );
  //parnames.push_back( "sigma_tcorr" );

  TString dbline;
  
  for( int ipar=0; ipar<parnames.size(); ipar++ ){
    for( int imod=0; imod<nmodules; imod++ ){
      double par0,par1;
      switch(ipar){
      case 0:
	par0 = maxstrip_t0[imod][0];
	par1 = maxstrip_t0[imod][1];
	break;
      case 1:
	par0 = maxstrip_tsigma[imod][0];
	par1 = maxstrip_tsigma[imod][1];
	break;
      case 2:
	par0 = maxstrip_t0_deconv[imod][0];
	par1 = maxstrip_t0_deconv[imod][1];
	break;
      case 3:
	par0 = maxstrip_tsigma_deconv[imod][0];
	par1 = maxstrip_tsigma_deconv[imod][1];
	break;
      case 4:
	par0 = tmean[imod][0];
	par1 = tmean[imod][1];
	break;
      case 5:
	par0 = tsigma[imod][0];
	par1 = tsigma[imod][1];
	break;
      case 6:
	par0 = tmean[imod][2];
	par1 = tmean[imod][3];
	break;
      case 7:
	par0 = tsigma[imod][2];
	par1 = tsigma[imod][3];
	break;
      default:
	par0 = 0.0;
	par1 = 0.0;
	break;
      }
	
      dbFT << dbline.Format( "sbs.gemFT.m%d.%s = %12.6g %12.6g", imod, parnames[ipar].Data(), par0, par1) << endl;
    }

    for( int imod=0; imod<nmodules_back; imod++ ){
      double par0,par1;
      switch(ipar){
      case 0:
	par0 = maxstrip_t0[imod+nmodules][0];
	par1 = maxstrip_t0[imod+nmodules][1];
	break;
      case 1:
	par0 = maxstrip_tsigma[imod+nmodules][0];
	par1 = maxstrip_tsigma[imod+nmodules][1];
	break;
      case 2:
	par0 = maxstrip_t0_deconv[imod+nmodules][0];
	par1 = maxstrip_t0_deconv[imod+nmodules][1];
	break;
      case 3:
	par0 = maxstrip_tsigma_deconv[imod+nmodules][0];
	par1 = maxstrip_tsigma_deconv[imod+nmodules][1];
	break;
      case 4:
	par0 = tmean[imod+nmodules][0];
	par1 = tmean[imod+nmodules][1];
	break;
      case 5:
	par0 = tsigma[imod+nmodules][0];
	par1 = tsigma[imod+nmodules][1];
	break;
      case 6:
	par0 = tmean[imod+nmodules][2];
	par1 = tmean[imod+nmodules][3];
	break;
      case 7:
	par0 = tsigma[imod+nmodules][2];
	par1 = tsigma[imod+nmodules][3];
	break;
      default:
	par0 = 0.0;
	par1 = 0.0;
	break;
      }
      dbFPP << dbline.Format( "sbs.gemFPP.m%d.%s = %12.6g %12.6g", imod, parnames[ipar].Data(), par0, par1) << endl;
    }

    
    dbFT << endl;

    dbFPP << endl;
  }

  parnames.clear();
  parnames.push_back( "deltat_sigma" );
  parnames.push_back( "deltat_sigma_deconv" );
  parnames.push_back( "sigma_tcorr" );
  parnames.push_back( "threshold_sample" );
  parnames.push_back( "threshold_stripsum" );
  parnames.push_back( "threshold_clustersum" );
  parnames.push_back( "threshold_maxcombo_deconv" );
  parnames.push_back( "threshold_clustersum_deconv" );
  parnames.push_back( "ADCasym_sigma" );
  
  for( int ipar=0; ipar<parnames.size(); ipar++ ){
    for( int imod=0; imod<nmodules; imod++ ){
      double par;
      switch(ipar){
      case 2:
	par = sigma_tcorr[imod];
	break;
      case 0:
	par = sigma_dt[imod];
	break;
      case 1:
	par = sigma_dtdeconv[imod];
	break;
      case 3:
	par = thresh_sample[imod];
	break;
      case 4:
	par = thresh_strip[imod];
	break;
      case 5:
	par = thresh_cluster[imod];
	break;
      case 6:
	par = thresh_maxcombo_deconv[imod];
	break;
      case 7:
	par = thresh_clustersum_deconv[imod];
	break;
      case 8:
	par = ADCasym_sigma[imod];
	break;
      default:
	par = 0.0;
	break;
      }

      dbFT << dbline.Format( "sbs.gemFT.m%d.%s = %12.6g", imod, parnames[ipar].Data(), par ) << endl;
      
    }
    dbFT << endl;

    for( int imod=0; imod<nmodules_back; imod++ ){
      double par;
      switch(ipar){
      case 2:
	par = sigma_tcorr[imod+nmodules];
	break;
      case 0:
	par = sigma_dt[imod+nmodules];
	break;
      case 1:
	par = sigma_dtdeconv[imod+nmodules];
	break;
      case 3:
	par = thresh_sample[imod+nmodules];
	break;
      case 4:
	par = thresh_strip[imod+nmodules];
	break;
      case 5:
	par = thresh_cluster[imod+nmodules];
	break;
      case 6:
	par = thresh_maxcombo_deconv[imod+nmodules];
	break;
      case 7:
	par = thresh_clustersum_deconv[imod+nmodules];
	break;
      case 8:
	par = ADCasym_sigma[imod+nmodules];
	break;
      default:
	par = 0.0;
	break;
      }

      dbFPP << dbline.Format( "sbs.gemFPP.m%d.%s = %12.6g", imod, parnames[ipar].Data(), par ) << endl;
      
    }
    dbFPP << endl;
    
  }
    
  

  fout->Write();

  
}

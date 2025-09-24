#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "GEM_cosmic_tracks.C"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TSystem.h"
#include <iostream>
#include <fstream> 
#include "TROOT.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVectorD.h"
#include "TMatrixD.h"
//#include "Math/Functor.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCut.h"
#include "TTreeFormula.h"
// bunch of header files that are needed throughout the script

//What is the purpose of each one of these vector variables?

//Global variables used by the chi2 function:
vector<double> AsymALL, dAsymALL;
vector<int> APVX_asymALL, APVY_asymALL;
vector<vector<double> > weightX, weightY;
vector<int> countX,countY;
//Number of possible APVs in X and Y directions. Looks like in beam transport coordinates
//These variables get overridden by the "strip config" file:
int nAPVmaxX = 12;
int nAPVmaxY = 10;

//I presume a chi2 function which is a helper for the gain match function? Need to understand each input parameter
void chi2_FCN( int &npar, double *gin, double &f, double *par, int flag ){

  //cout << "calculating chi2..." << endl;
  //why do we assume this for the start parameter?
  double chi2 = 0.0;
  for( int i=0; i<AsymALL.size(); i++ ){
    //Define loop variables iy and ix as y and x APV asymmetry respectively?
    // this is the combination of X and Y APV cards for the i'th entry in the asymmetry array
    int iy = APVY_asymALL[i];
    int ix = APVX_asymALL[i];
    //What is this variable for? this variable becomes true if certain quality checks are passed and this particular combination should be included in the fit
    bool useasym = false;
    // Gain starting values? why do we pick 1.0 for start? 1 is default
    double Gx=1.0, Gy=1.0;
    //what is the purpose of each conditional statement?
    if( iy >= 0 && iy < nAPVmaxY && ix >= 0 && ix < nAPVmaxX ){
      //Need to understand par before can understand this variable assignment
      Gx = par[ix + nAPVmaxY];
      Gy = par[iy];
      //why does this variable change? because if we are here then we have a good asymmetry for this combination of X and Y APVs
      useasym = true;
    } else if( iy >= 0 && iy < nAPVmaxY ){ //this is a module-average asymmetry for Y: using weighted average over all X APVs for Gx:
      double sum_Gx = 0.0;
      double sum_Wx = 0.0;
      for( int j = 0; j<nAPVmaxX; j++ ){
	sum_Gx += weightX[iy][j] * par[j+nAPVmaxY];
	sum_Wx += weightX[iy][j];
      }
      //Which definition is actually be used. What is the mechanic here?

      //This is replacing Gx with the statistics-weighted average of all the X APVs:
      Gx = sum_Gx/sum_Wx;
      Gy = par[iy];

      //only use the module-average asymmetry in the fit if there were no X APVs with at least 100 events on this Y APV:
      // 
      // if( countX[iy] == 0 )

      //Let's see what happens if we use the module-average asymmetry here:
      if( countX[iy] == 0 ) useasym = true;
      
    } else if( ix >= 0 && ix < nAPVmaxX ){ //this is a module-average asymmetry for X: using weighted average over all Y APVs for Gy:
      double sum_Gy = 0.0;
      double sum_Wy = 0.0;
      for( int j=0; j<nAPVmaxY; j++ ){
	sum_Gy += weightY[ix][j] * par[j];
	sum_Wy += weightY[ix][j];
      }

      Gx = par[ix+nAPVmaxY];
      Gy = sum_Gy/sum_Wy;

      if( countY[ix] == 0 ) useasym = true;
	
    }
    //The actual asymmetry  
    double Atheory = (Gx - Gy)/(Gx + Gy);
    //it is the sum over all measured asymmetries of the difference between the measured and "expected" asymmetries, which are funcitons of the APV gain coeffs:
    if( useasym ) chi2 += pow( (AsymALL[i] - Atheory)/dAsymALL[i], 2 );
    
  }

  //cout << "chi2 = " << chi2 << endl;
  
  f = chi2;
}


//void GEM_GainMatch(const char *infilename, int nmodules, const char *fname_stripconfig, const char *detname="bb.gem", double chi2cut=10.0, double ADCcut = 1000.0, double target_ADC=4500.0){
void GEM_GainMatch_gep( const char *configfilename, const char *outfname="output/GEM_gainmatch_temp.root" ){

  ifstream configfile(configfilename);

  if( !configfile ) return;
  
  TChain *C = new TChain("T");

  TString currentline;
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );
  
  //convert runnum to a char so that way it can be used in names
  //std::string runnum_temp = std::to_string(runnum);
  //const char *runnum_char = runnum_temp.c_str();
  //cout <<"My number " << runnum_char << endl;

  //Read strip configuration:

  
  
  //setups up an input stream for the strip config file
  //ifstream stripconfigfile(fname_stripconfig);
  // defines an emptry current line for the input stream
  //TString currentline;
  //boolean parameters to help debug reading in the strip config file
  int nmodules=14;

  bool gotxconfig = false;
  bool gotyconfig = false;
  //Places to store the x and y strip information based on the total number of modules. Think the scheme is from target to back tracker. Top to bottom?
  vector<int> nstripx_mod(nmodules);
  vector<int> nstripy_mod(nmodules);

  TString detname = "sbs.gemFT";

  double target_ADC = 4000.;

  double ADCcut = 1000.0;
  double ccor_cut = 0.7;
  double deltat_cut = 25.0;

  double ADCasym_cut = 1.0;
  
  //Parses stripconfigfile and gets x/u strip information
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
    if( currentline.BeginsWith("detname") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
        detname = ( (TObjString*) (*tokens)[1] )->GetString();
      }
    }

    if( currentline.BeginsWith("nmodules") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	nmodules = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	nstripx_mod.resize(nmodules);
	nstripy_mod.resize(nmodules);
      }
    }

    if( currentline.BeginsWith("ADCcut") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	ADCcut = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
      }
    }

    if( currentline.BeginsWith("ADCasym_cut") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	ADCasym_cut = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
      }
    }

    if( currentline.BeginsWith("ccor_cut") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	ccor_cut = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
      }
    }

    if( currentline.BeginsWith("deltat_cut") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	deltat_cut = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
      }
    }
    
    if( currentline.BeginsWith("mod_nstripu") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= nmodules+1 ){
	for( int i=1; i<=nmodules; i++ ){
	  //this line stores the strip info in a vector but it seems like the type is mismatched?
	  nstripx_mod[i-1] = ( (TObjString*) (*tokens)[i] )->GetString().Atoi();
	}
	gotxconfig = true;
      }
    }
    //parses stripconfigfile and gets y/v strip information
    if( currentline.BeginsWith( "mod_nstripv" ) ){
      TObjArray *tokens = currentline.Tokenize( " " );
      if( tokens->GetEntries() >= nmodules+1 ){
	for( int i=1; i<=nmodules; i++ ){
	  //this line stores the strip info in a vector but again type mismatch?
	  nstripy_mod[i-1] = ( (TObjString*) (*tokens)[i] )->GetString().Atoi();
	  
	}
	gotyconfig = true;
      }
    }
  }
 

  //Error message if problem with strip config file
  if( !( gotxconfig && gotyconfig ) ){
    cout << "module strip config not properly defined, quitting..." << endl;
    return;
  }
  // parameters and loop to find the max number of x and y strips respectively
  int nstripxmax = 0, nstripymax=0; 
  for( int i=0; i<nmodules; i++ ){
    nstripxmax = nstripx_mod[i] > nstripxmax ? nstripx_mod[i] : nstripxmax;
    nstripymax = nstripy_mod[i] > nstripymax ? nstripy_mod[i] : nstripymax; 
  }
  
  //something for root?
  gROOT->ProcessLine(".x ~/rootlogon.C");

  //Arbitrarily choose 10000 max hits?
  UInt_t MAXNHITS=1000;
  //Not used?
  //TString fname(infilename);
  //Remind me why TChain is important
  UInt_t MAXNTRACKS=100;
  
  // TChain *C = new TChain("T");

  // C->Add( infilename );
  //Do some file name management. This will of course need to change if you get files from not ewertz volatile

 
  //C->Print();
  //some variables that are useful, for some reason
  TString branchname;
  double ngoodhits;
  double ntracks;
  double besttrack;
  //create a bunch of vector<double>, presumably to store info?
  vector<double> tracknhits(MAXNTRACKS);
  vector<double> trackChi2NDF(MAXNTRACKS);
  vector<double> hit_trackindex(MAXNHITS);
  vector<double> hit_module(MAXNHITS);
  vector<double> hit_layer(MAXNHITS);
  vector<double> hit_nstripu(MAXNHITS);
  vector<double> hit_nstripv(MAXNHITS);
  vector<double> hit_ustripmax(MAXNHITS);
  vector<double> hit_ustriplo(MAXNHITS);
  vector<double> hit_ustriphi(MAXNHITS);
  vector<double> hit_vstripmax(MAXNHITS);
  vector<double> hit_vstriplo(MAXNHITS);
  vector<double> hit_vstriphi(MAXNHITS);
  vector<double> hit_ADCU(MAXNHITS);
  vector<double> hit_ADCV(MAXNHITS);
  vector<double> hit_Ugain(MAXNHITS);
  vector<double> hit_Vgain(MAXNHITS);
  vector<double> hit_ADCavg(MAXNHITS);
  vector<double> hit_ADCasym(MAXNHITS);
  vector<double> hit_ccor_clust(MAXNHITS);
  vector<double> hit_ccor_strip(MAXNHITS);
  vector<double> hit_deltat(MAXNHITS);

  //a map for branch names and a vector for variable names
  vector<double> hit_ADCmaxsampU(MAXNHITS);
  vector<double> hit_ADCmaxsampV(MAXNHITS);
  vector<double> hit_ADCmaxstripU(MAXNHITS);
  vector<double> hit_ADCmaxstripV(MAXNHITS);
  //vector<double> Epreshower(MAXNHITS);
  //vector<double> Eshower(MAXNHITS);
  //double Epreshower,Eshower;
  
  //double xfcp,yfcp,zfcp,xbcp,ybcp,zbcp;
  vector<double> xfcp(MAXNTRACKS);
  vector<double> yfcp(MAXNTRACKS);
  vector<double> zfcp(MAXNTRACKS);
  vector<double> xbcp(MAXNTRACKS);
  vector<double> ybcp(MAXNTRACKS);
  vector<double> zbcp(MAXNTRACKS);
  
  vector<double> track_vz(MAXNTRACKS);
  vector<double> track_p(MAXNTRACKS);


  map<TString,TString> branchnames;
  vector<TString> varnames;
  varnames.push_back("track.ntrack");
  varnames.push_back("track.nhits");
  varnames.push_back("track.besttrack");
  varnames.push_back("track.chi2ndf");
  varnames.push_back("hit.ngoodhits");
  varnames.push_back("hit.trackindex");
  varnames.push_back("hit.module");
  varnames.push_back("hit.layer");
  varnames.push_back("hit.nstripu");
  varnames.push_back("hit.nstripv");
  varnames.push_back("hit.ustripmax");
  varnames.push_back("hit.ustriplo");
  varnames.push_back("hit.ustriphi");
  varnames.push_back("hit.vstripmax");
  varnames.push_back("hit.vstriplo");
  varnames.push_back("hit.vstriphi");
  varnames.push_back("hit.ADCU");
  varnames.push_back("hit.ADCV");
  varnames.push_back("hit.Ugain");
  varnames.push_back("hit.Vgain");
  varnames.push_back("hit.ADCavg");
  varnames.push_back("hit.ADCasym");

  //Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  //cout << "disabling all branches...";
  //the * applies it to all branches, the 0 disables those branches. to enable would need to make 1

  varnames.push_back("hit.ADCmaxsampU");
  varnames.push_back("hit.ADCmaxsampV");
  varnames.push_back("hit.ADCmaxstripU");
  varnames.push_back("hit.ADCmaxstripV");

  varnames.push_back("hit.ccor_clust");
  varnames.push_back("hit.ccor_strip");
  varnames.push_back("hit.deltat");

  cout << "disabling all branches...";
  

  C->SetBranchStatus("*",0);

  //cout << "done." << endl;

  for( int i=0; i<varnames.size(); i++ ){
    //What is actually going on in this for loop. It looks like branch and variable names are getting printed. But also setting branches up in the TChain
    branchnames[varnames[i]] = branchname.Format("%s.%s",detname.Data(),varnames[i].Data());
    // cout << "Branch " << i << " name = " << branchnames[varnames[i]] << endl;
    C->SetBranchStatus( branchnames[varnames[i]].Data(), 1 );
  }
  //Populating data in the TChain branchs?
  // cout << "Setting branch addresses: ";

  //C->SetBranchStatus("bb.tr.*",1);
  //C->SetBranchStatus("bb.tr.vz",1);
  //These are redundant.
  //C->SetBranchStatus("bb.gem.track.chi2ndf",1);
  //C->SetBranchStatus("bb.gem.track.nhits",1);
  //C->SetBranchStatus("bb.etot_over_p",1);
  //C->SetBranchStatus("bb.tr.p",1);
  //C->SetBranchStatus("bb.ps.e",1);
  //C->SetBranchStatus("bb.sh.e",1);
  //C->SetBranchStatus("bb.x_bcp",1);
  //C->SetBranchStatus("bb.y_bcp",1);
  //C->SetBranchStatus("bb.z_bcp",1);
  //C->SetBranchStatus("bb.x_fcp",1);
  //C->SetBranchStatus("bb.y_fcp",1);
  //C->SetBranchStatus("bb.z_fcp",1);


  //For the moment, we're only activating these branches for the purpose of defining cuts. The ones used by the script will be activated below:
  
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.x_bcp",1);
  C->SetBranchStatus("sbs.y_bcp",1);
  C->SetBranchStatus("sbs.z_bcp",1);

  C->SetBranchStatus("sbs.x_fcp",1);
  C->SetBranchStatus("sbs.y_fcp",1);
  C->SetBranchStatus("sbs.z_fcp",1);

  C->SetBranchStatus("heep.dyECAL",1);
  C->SetBranchStatus("heep.dxECAL",1);
  C->SetBranchStatus("heep.dpp",1);
  C->SetBranchStatus("heep.dt_ADC",1);
  C->SetBranchStatus("g.runnum",1);


//  C->SetBranchStatus("sbs.gemFT.track.nhits",1);

  cout << "Setting branch addresses: ";

  C->SetBranchAddress( branchnames["track.ntrack"].Data(), &ntracks );
  C->SetBranchAddress( branchnames["track.besttrack"].Data(), &besttrack );
  C->SetBranchAddress( branchnames["track.nhits"].Data(), &(tracknhits[0]) );
  C->SetBranchAddress( branchnames["track.chi2ndf"].Data(), &(trackChi2NDF[0]) );
  C->SetBranchAddress( branchnames["hit.ngoodhits"].Data(), &ngoodhits );
  C->SetBranchAddress( branchnames["hit.trackindex"].Data(), &(hit_trackindex[0]) );
  C->SetBranchAddress( branchnames["hit.module"].Data(), &(hit_module[0]) );
  C->SetBranchAddress( branchnames["hit.layer"].Data(), &(hit_layer[0]) );
  C->SetBranchAddress( branchnames["hit.nstripu"].Data(), &(hit_nstripu[0]) );
  C->SetBranchAddress( branchnames["hit.nstripv"].Data(), &(hit_nstripv[0]) );
  C->SetBranchAddress( branchnames["hit.ustripmax"].Data(), &(hit_ustripmax[0]) );
  C->SetBranchAddress( branchnames["hit.ustriplo"].Data(), &(hit_ustriplo[0]) );
  C->SetBranchAddress( branchnames["hit.ustriphi"].Data(), &(hit_ustriphi[0]) );
  C->SetBranchAddress( branchnames["hit.vstripmax"].Data(), &(hit_vstripmax[0]) );
  C->SetBranchAddress( branchnames["hit.vstriplo"].Data(), &(hit_vstriplo[0]) );
  C->SetBranchAddress( branchnames["hit.vstriphi"].Data(), &(hit_vstriphi[0]) );
  C->SetBranchAddress( branchnames["hit.ADCU"].Data(), &(hit_ADCU[0]) );
  C->SetBranchAddress( branchnames["hit.ADCV"].Data(), &(hit_ADCV[0]) );
  C->SetBranchAddress( branchnames["hit.Ugain"].Data(), &(hit_Ugain[0]) );
  C->SetBranchAddress( branchnames["hit.Vgain"].Data(), &(hit_Vgain[0]) );
  C->SetBranchAddress( branchnames["hit.ADCavg"].Data(), &(hit_ADCavg[0]) );
  C->SetBranchAddress( branchnames["hit.ADCasym"].Data(), &(hit_ADCasym[0]) );

  // cout << "done." << endl;
  //Setup the output file name 
  //changed output file name to include run number for file organization
  //  TString outfilename = Form("GEM_GainMatch_output/GainRatios_%s_%s.root",detname.Data(),runnum_char);
  
  if (detname == "sbs.gemFPP") { 
    outfname="output/GEMFPP_gainmatch_temp.root";
    C->SetBranchStatus("sbs.gemFT.track.nhits",1);
  } else {
    outfname="output/GEMFT_gainmatch_temp.root";
  }

  TString outfilename(outfname);
  
  C->SetBranchAddress( branchnames["hit.ADCmaxsampU"].Data(), &(hit_ADCmaxsampU[0]) );
  C->SetBranchAddress( branchnames["hit.ADCmaxsampV"].Data(), &(hit_ADCmaxsampV[0]) );
  C->SetBranchAddress( branchnames["hit.ADCmaxstripU"].Data(), &(hit_ADCmaxstripU[0]) );
  C->SetBranchAddress( branchnames["hit.ADCmaxstripV"].Data(), &(hit_ADCmaxstripV[0]) );
  C->SetBranchAddress( branchnames["hit.ccor_clust"].Data(), &(hit_ccor_clust[0]) );
  C->SetBranchAddress( branchnames["hit.ccor_strip"].Data(), &(hit_ccor_strip[0]) );
  C->SetBranchAddress( branchnames["hit.deltat"].Data(), &(hit_deltat[0]) );

  // These variables are important for whatever the global cut setin the config file is
  // Currently there seems to be a segmentation violation issue when switching bb to sbs
  C->SetBranchAddress( "sbs.tr.vz", &(track_vz[0]) );
  C->SetBranchAddress( "sbs.tr.p", &(track_p[0]) );
  //C->SetBranchAddress( "bb.ps.e",&Epreshower);
  //C->SetBranchAddress( "bb.sh.e",&Eshower);
  C->SetBranchAddress( "sbs.x_fcp", &(xfcp[0]) );
  C->SetBranchAddress( "sbs.y_fcp", &(yfcp[0]) );
  C->SetBranchAddress( "sbs.z_fcp", &(zfcp[0]) );
  C->SetBranchAddress( "sbs.x_bcp", &(xbcp[0]) );
  C->SetBranchAddress( "sbs.y_bcp", &(ybcp[0]) );
  C->SetBranchAddress( "sbs.z_bcp", &(zbcp[0]) );
  cout << "done." << endl;
  
  //TString outfilename = Form("GainRatios_%s_temp.root",detname.Data()); 


  

  //  outfilename.Prepend("GainRatios_");

  // cout << "Out file name = " << outfilename << endl;
  //Next couple of lines is to output the info to files
  TFile *fout = new TFile( outfilename.Data(), "RECREATE" );
  
  outfilename.ReplaceAll(".root",".txt");

  ofstream outfile(outfilename.Data());
  //changed output file name to include run number for file organization

  outfilename.ReplaceAll(".txt",".dat");
  
  TString dbfilename(outfilename.Data());
  //dbfilename.Form( "GEM_GainMatch_output/GEM_GainMatchResults_%s_%s.dat",detname.Data(),runnum_char );
  
  ofstream outfile_db(dbfilename.Data());
  //what is going on here with the nAPVmaxX/Y?
  nAPVmaxX = nstripxmax/128;
  if( nstripxmax % 128 > 0 ) nAPVmaxX++;

  nAPVmaxY = nstripymax/128;
  if( nstripymax % 128 > 0 ) nAPVmaxY++;
  //Now defining TClonesArrays?
  TClonesArray *hADCasym_vs_APVXY = new TClonesArray( "TH1D", nAPVmaxX*nAPVmaxY*nmodules );
  TClonesArray *hADCasym_vs_APVX = new TClonesArray( "TH1D", nAPVmaxX*nmodules );
  TClonesArray *hADCasym_vs_APVY = new TClonesArray( "TH1D", nAPVmaxY*nmodules );
  //nested for loops for asymmetry calculations? Need some clarification here
  for( int imodule=0; imodule<nmodules; imodule++ ){
    TString hname;
    for( int iAPVx = 0; iAPVx < nAPVmaxX; iAPVx++ ){
      hname.Form( "hADCasym_vs_APVX%d_mod%d", iAPVx, imodule );

      if( iAPVx*128 < nstripx_mod[imodule] ){
      
	new( (*hADCasym_vs_APVX)[iAPVx+nAPVmaxX*imodule] ) TH1D( hname.Data(), "", 250, -1.01, 1.01 );

      }
      
      for( int iAPVy = 0; iAPVy < nAPVmaxY; iAPVy++ ){

	if( iAPVy*128 < nstripy_mod[imodule] ){
	  if( iAPVx==0 ){
	    hname.Form( "hADCasym_vs_APVY%d_mod%d", iAPVy, imodule );
	    new( (*hADCasym_vs_APVY)[iAPVy+nAPVmaxY*imodule] ) TH1D( hname.Data(), "", 250, -1.01, 1.01 );
	  }
	
	  hname.Form("hADCasym_vs_APV_mod%d_x%d_y%d",imodule,iAPVx,iAPVy);
	  
	  new( (*hADCasym_vs_APVXY)[iAPVy+nAPVmaxY*iAPVx+nAPVmaxX*nAPVmaxY*imodule] ) TH1D( hname.Data(), "", 250,-1.01,1.01);
	}
      }
    }
  }
  
  //GEM_cosmic_tracks *T = new GEM_cosmic_tracks(C);
  
  long nevent=0;
  //make some 1D and 2D histograms
  TH2D *hADCasym_module = new TH2D("hADCasym_module","ADC Asymm vs Module Before Gain;Module;ADC Asymm",nmodules,-0.5,nmodules-0.5,500,-1.01,1.01);
  TH2D *hNstripX_module = new TH2D("hNstripX_module","",nmodules,-0.5,nmodules-0.5,12,0.5,12.5);
  TH2D *hNstripY_module = new TH2D("hNstripY_module","",nmodules,-0.5,nmodules-0.5,12,0.5,12.5);
  TH2D *hADCasym_APV = new TH2D("hADCasym_APV","ADC Asymm vs APV Before Gain;APV;ADC Asymm",nAPVmaxX*nAPVmaxY*nmodules,-0.5,nAPVmaxX*nAPVmaxY*nmodules-0.5,500,-1.01,1.01);
  TH2D *hADCasym_APV_corrected = new TH2D("hADCasym_APV_corrected","ADC Asymm vs APV After Gain;APV;ADC Asymm",nAPVmaxX*nAPVmaxY*nmodules,-0.5,nAPVmaxX*nAPVmaxY*nmodules-0.5,500,-1.01,1.01);

  TH2D *hADCavg_module = new TH2D("hADCavg_module","ADC Average vs Module Before Gain;Module;0.5*(Clust U ADC + Clust V ADC)",nmodules,-0.5,nmodules-0.5,1500,0,30000);
  TH1D *hADCavg_allhits = new TH1D("hADCavg_allhits","",1500,0,30000);

  TH2D *hStripADCsumU_module = new TH2D("hStripADCsumU_module","",nmodules,-0.5,nmodules-0.5,1500,0,15000);
  TH2D *hStripADCmaxU_module = new TH2D("hStripADCmaxU_module","",nmodules,-0.5,nmodules-0.5,1500,0,3000);
  TH2D *hStripADCsumV_module = new TH2D("hStripADCsumV_module","",nmodules,-0.5,nmodules-0.5,1500,0,15000);
  TH2D *hStripADCmaxV_module = new TH2D("hStripADCmaxV_module","",nmodules,-0.5,nmodules-0.5,1500,0,3000);
  
  //int nAPVmax = 

  int treenum=0, currenttreenum=0;
  //cout << "starting event loop:" << endl;
  //Need to understand what this loop is for and what every conditional is for
  while( C->GetEntry( nevent++ ) ){

    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passedcut = GlobalCut->EvalInstance(0) != 0;
    
    if( ngoodhits > MAXNHITS ) continue;
    if( ntracks > MAXNHITS ) continue;
    
    if( nevent % 100000 == 0 ) cout << "event " << nevent << endl;
    //loop over hits:
    int itrack = int(besttrack);
    //cout << "itrack = " << itrack << endl;

    int NTRACKS = int(ntracks);

    //cout << "ntracks = " << NTRACKS << endl;
    if(passedcut && NTRACKS >= 1 && itrack == 0){
      
      int nhits = int(ngoodhits);
      //	cout << "nhits = " << nhits << endl;
      
      for( int ihit=0; ihit<nhits; ihit++ ){
	int tridx = int( hit_trackindex[ihit] );

	//cout << "ihit, tridx = " << ihit << ", " << tridx;

	
	if( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) >= ADCcut && hit_ADCU[ihit] >= ADCcut && hit_ADCV[ihit] >= ADCcut && tridx == itrack && hit_nstripu[ihit]>1 && hit_nstripv[ihit]>1 && hit_ccor_clust[ihit] >= ccor_cut && fabs(hit_deltat[ihit]) <= deltat_cut ){
	  // cout << ", ADCavg[ihit] = " << hit_ADCavg[ihit]
	  // 	 << ", ADCasym[ihit] = " << hit_ADCasym[ihit]
	  // 	 << ", (nstripu,nstripv) = (" << hit_nstripu[ihit] << ", " << hit_nstripv[ihit] << ")"
	  // 	 << ", (ustriplo,ustriphi,ustripmax)=(" << hit_ustriplo[ihit] << ", " << hit_ustriphi[ihit] << ", " << hit_ustripmax[ihit] << ")"
	  // 	 << ", (vstriplo,vstriphi,vstripmax)=(" << hit_vstriplo[ihit] << ", " << hit_vstriphi[ihit] << ", " << hit_vstripmax[ihit] << ")"
	  // 	 << ", module = " << hit_module[ihit] << endl;

	  double hit_ADCasym_nogain = (hit_ADCU[ihit]/hit_Ugain[ihit]-hit_ADCV[ihit]/hit_Vgain[ihit]) / (hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]);
	  hADCavg_allhits->Fill( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) );
	  hADCavg_module->Fill( hit_module[ihit], 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) );
	  hADCasym_module->Fill( hit_module[ihit], hit_ADCasym_nogain );
	  hNstripX_module->Fill( hit_module[ihit], hit_nstripu[ihit] );
	  hNstripY_module->Fill( hit_module[ihit], hit_nstripv[ihit] );

	  //Define a bunch of variables, but what are they for?


	  hStripADCsumU_module->Fill( hit_module[ihit], hit_ADCmaxstripU[ihit] );
	  hStripADCsumV_module->Fill( hit_module[ihit], hit_ADCmaxstripV[ihit] );

	  hStripADCmaxU_module->Fill( hit_module[ihit], hit_ADCmaxsampU[ihit] );
	  hStripADCmaxV_module->Fill( hit_module[ihit], hit_ADCmaxsampV[ihit] );


	  int ixlo = hit_ustriplo[ihit];
	  int ixhi = hit_ustriphi[ihit];
	  int ixmax = hit_ustripmax[ihit];

	  int iylo = hit_vstriplo[ihit];
	  int iyhi = hit_vstriphi[ihit];
	  int iymax = hit_vstripmax[ihit];

	  int xAPVmax = ixmax/128;
	  int yAPVmax = iymax/128;
	  int xAPVlo = ixlo/128;
	  int yAPVlo = iylo/128;

	  int xAPVhi = ixhi/128;
	  int yAPVhi = iyhi/128;

	  int module = int(hit_module[ihit]);
	  //what is the purpose of this conditional and assignment?
	  //To make sure all strips on this hit are on the same APV card, and
	  // to make sure we don't cause a seg fault by accessing an array index that's
	  // out-of-bounds:
	  if( xAPVlo == xAPVmax && xAPVhi == xAPVmax &&
	      yAPVlo == yAPVmax && yAPVhi == yAPVmax &&
	      hit_nstripu[ihit] >= 2 && hit_nstripv[ihit] >= 2 &&
	      xAPVmax*128 < nstripx_mod[module] && yAPVmax*128<nstripy_mod[module] ){

	    //cout << "filling histograms " << endl;

	    //cout << "apvxy hist index = " << yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*module << endl;
	    // cout << "apvx hist index = " << xAPVmax + nAPVmaxX*module << endl;
	    // cout << "apvy hist index = " << yAPVmax + nAPVmaxY*module << endl;
	    //More 1D histograms    
	    ( (TH1D*) (*hADCasym_vs_APVXY)[yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*module] )->Fill( hit_ADCasym_nogain );

	    ( (TH1D*) (*hADCasym_vs_APVX)[xAPVmax + nAPVmaxX*module] )->Fill( hit_ADCasym_nogain );
	    ( (TH1D*) (*hADCasym_vs_APVY)[yAPVmax + nAPVmaxY*module] )->Fill( hit_ADCasym_nogain );

	    hADCasym_APV->Fill(yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*module,hit_ADCasym_nogain);
	    
	    //cout << "done histogram fill" << endl;

	    
	  }
	}
	////
      }
    }
  }
  cout << "Event " << nevent << " done" << endl;
  cout << "Entries " << hADCavg_allhits->GetEntries() << " done" << endl;

  //Fitting landau's to something
  TFitResultPtr fitadcall = hADCavg_allhits->Fit("landau","qS","",2000.0,25000.);

  double MPV_all = ( (TF1*) hADCavg_allhits->GetListOfFunctions()->FindObject("landau") )->GetParameter("MPV");

  cout << "All hits ADC peak position = " << MPV_all << endl;

  
  

  double asympeak[nmodules];
  double R[nmodules];
  
  TH1D *htemp;

  //outfile << "mod_RYX     ";
  //Setup the TCanvas. Should lookup what each one of these functions is
  TCanvas *c1 = new TCanvas("c1","c1",2000,1000);
  c1->Divide(2,1);
  c1->cd(1);

  hADCasym_module->Draw("colz");
  
  gPad->Modified();
  c1->Update();
  
  c1->cd(2);
  //creates vector to store relative gain info
  vector<double> RelativeGainByModule(nmodules);
  TGraphErrors* g = new TGraphErrors();
  int point = 0;
  //what is going on in this for loop?
  for( int i=0; i<nmodules; i++ ){
    TString hnametemp;
    hnametemp.Form("ADCasym_module%d",i);
    htemp = hADCasym_module->ProjectionY(hnametemp.Data(), i+1, i+1 );

    if( htemp->GetEntries() >= 100 ){
      
      double max = htemp->GetMaximum();
      int binmax,binlo,binhi;
      
      binmax = htemp->GetMaximumBin();
      binlo = binmax;
      binhi = binmax;
      
      while( htemp->GetBinContent(binlo) >= 0.5*max && binlo > 1 ){binlo--;}
      while( htemp->GetBinContent(binhi) >= 0.5*max && binhi < 500 ){binhi++; }
      
      
      
      htemp->Fit("gaus","qS","",htemp->GetBinCenter(binlo),htemp->GetBinCenter(binhi) );
      
      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      asympeak[i] = fitfunc->GetParameter(1);
      R[i] = (1.0-asympeak[i])/(1.0+asympeak[i]);

      //         TGraphErrors* g = new TGraphErrors();
      //int point = 0;
    //    for(int i=0; i<nmodules; i++){
    //      if(htemp->GetEntries() >= 100){ // only for bins you fitted
      double x = hADCasym_module->GetXaxis()->GetBinCenter(i+1); // x-position
      double y = asympeak[i];  // mean from Gaussian
      double ex = 0;           // no x error
      double ey = fitfunc->GetParameter(2); // sigma from fit

      g->SetPoint(point, x, y);
      g->SetPointError(point, ex, ey);
      point++;
      
      c1->cd(1);
      hADCasym_module->Draw("colz");
      g->SetMarkerStyle(21);
      g->SetMarkerColor(kRed);
      g->Draw("P SAME"); // overlay points
      
      gPad->Modified();
      c1->Update();
      
      //      gSystem->Sleep(250);
      TString outstring;
      //outfile << outstring.Format(" %12.6g ", R[i] );

      //Within each module, we measure an asymmetry for all possible combinations of X APV (or "U") and Y APV (or "V"):

      //Let us find the gain coefficients for each individual APV card that minimize the chi^2 defined as the
      // sum of squared differences between the modified ADC asymmetries and zero:
      // If APV card X has relative gain Gx and card Y has relative gain Gy, then the observed ADC asymmetry would be:
      // ASYM = ( Gx - Gy )/(Gx + Gy ) = (1 - Ryx_i)/(1 + Ryx_i)
      // Let's choose some reference APVX as having gain = 1. So all gains are measured relative to this common
      // reference. 
      // Ryx = (1-A)/(1+A),
      // For sufficiently small A, Ryx - 1 = (1-A)/(1+A) - 1 = -2A/(1+A) ~= -2A
      //We have asymmetries A_ij 
    }
    //This is where the main gain matching part happens as described by the above. Need to understand what is going on here
    hnametemp.Form("ADCdist_module%d", i );

    htemp = hADCavg_module->ProjectionY( hnametemp.Data(), i+1, i+1 );

    if( htemp->GetEntries() >= 100. ){

      TFitResultPtr ADCfit_module = htemp->Fit("landau","qS","",ADCcut,25000.0);
      double MPV_mod = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");

      // cout << "module " << i << " MPV = " << MPV_mod << endl;
      // cout << "module " << i << " relative gain = " << MPV_mod / target_ADC << endl;
      RelativeGainByModule[i] = MPV_mod/MPV_all;
    } else {
      RelativeGainByModule[i] = 1.0;
    }
  }


  outfile_db << endl << "#Module internal relative gains by APV card: " << endl;

  map<int,vector<double> > Xgain_by_module, Ygain_by_module;
  
  for( int i=0; i<nmodules; i++ ){
    int xAPV_ref = nAPVmaxX/2;
      
    //First, determine all the Y gains relative to the reference APV,
    // A_{ij} = (1-R_{ij})/(1+R_{ij})
    // Ryx_{ij] = (1-A_{ij})/(1+A_{ij}) = Gyi/Gxj
    // 
    //chi^2 = sum_i,j (A_{ij}^(measured) - A_{ij}^{calc})^2/(Delta A_{ij}^{measured})^2
    // Is there a way to make this problem linear?
    // Maybe:
    //For sufficiently small values of R_{yx} - 1, we have:
    
    //A_{ij}^{calc} = (1 - Ryx_{ij})/(1+ Ryx_{ij}) 

    vector<double> Asym, dAsym;
    vector<double> Ryx, dRyx;

    vector<double> Xgain, dXgain, Ygain, dYgain; //by APV within this module:

    double minAsymerr = 1000.0;

    int count=0;

    AsymALL.clear();
    dAsymALL.clear();

    TMinuit *gainfit = new TMinuit( nAPVmaxX+nAPVmaxY );

    TString fitname;
    gainfit->SetName( fitname.Format("gainfit_module%d",i) ); 
    
    gainfit->SetFCN(chi2_FCN);
    
    double arglist[10];
    
    int ierflg=0;
    
    for( int ipar=0; ipar<nAPVmaxY; ipar++ ){
      TString parname;
      parname.Form("Gain_APVY%d",ipar);
      
      gainfit->mnparm( ipar, parname.Data(), 1.0, 0.1,0,0,ierflg ); 
      
    }
    
    for( int ipar=0; ipar<nAPVmaxX; ipar++ ){
      TString parname;
      parname.Form("Gain_APVX%d",ipar);
      
      gainfit->mnparm( ipar+nAPVmaxY, parname.Data(), 1.0, 0.1,0,0,ierflg ); 
      
    }
    AsymALL.clear();
    dAsymALL.clear();

    APVX_asymALL.clear();
    APVY_asymALL.clear();
    
    weightX.clear();
    weightY.clear();

    weightX.resize(nAPVmaxY);
    weightY.resize(nAPVmaxX);
    //Shall we just create a separate TMinuit object for each module? Why TF not?

    countX.resize(nAPVmaxY);
    countY.resize(nAPVmaxX); 

    xAPV_ref = -1;
    double maxnevent = 0.0;
    
    //Card-average asymmetries:
    for( int ix=0; ix<nAPVmaxX; ix++ ){
      weightY[ix].resize(nAPVmaxY);

      countY[ix] = 0;


      if( ix*128 < nstripx_mod[i] ){
      
	TH1D *htemp = ( (TH1D*) (*hADCasym_vs_APVX)[ix+nAPVmaxX*i] );
	if( htemp->GetEntries() >= 100 ){
	  
	  if( xAPV_ref < 0 || htemp->GetEntries() > maxnevent ){
	    xAPV_ref = ix;
	    maxnevent = htemp->GetEntries();
	  }
	  
	  double Amean = htemp->GetMean();
	  double dAmean = htemp->GetMeanError();
	  
	  double Arms = htemp->GetRMS();
	  double dArms = htemp->GetRMSError();
	  
	  htemp->Fit("gaus","SQ","",Amean-Arms, Amean+Arms);
	  
	  double Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	  double dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	  double Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	  double dArms_fit = htemp->GetFunction("gaus")->GetParError(2);
	  
	  //re-run the fit with a tighter range set by the result of the first fit:
	  
	  htemp->Fit( "gaus", "SQ", "", Amean_fit - 2.0*Arms_fit, Amean_fit + 2.0*Arms_fit );
	  
	  Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	  dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	  Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	  dArms_fit = htemp->GetFunction("gaus")->GetParError(2);
	  
	  double Ryxtemp = (1.0-Amean_fit)/(1.0+Amean_fit);
	  
	  double Aup = Amean_fit + dAmean_fit;
	  double Adown = Amean_fit - dAmean_fit;
	  
	  double Ryx_Aup = (1.0 - Aup)/(1.0 + Aup );
	  double Ryx_Adown = (1.0 - Adown)/(1.0 + Adown);
	  
	  AsymALL.push_back( Amean_fit );
	  dAsymALL.push_back( dAmean_fit );
	  
	  APVX_asymALL.push_back( ix );
	  APVY_asymALL.push_back( -1 );
	} else { //If less than 100 entries for this entire APV, fix the gain to 1:
	  gainfit->FixParameter( ix + nAPVmaxY );
	}
      } else {
	gainfit->FixParameter( ix + nAPVmaxY );
      }
    }

    //always fix the gain of the X APV with largest statistics to 1, to fix the absolute normalization
    gainfit->FixParameter( xAPV_ref + nAPVmaxY );
    
    for( int iy=0; iy<nAPVmaxY; iy++ ){
      weightX[iy].resize( nAPVmaxX );

      countX[iy] = 0;

      if( iy * 128 < nstripy_mod[i] ){
	
	TH1D *htemp = ( (TH1D*) (*hADCasym_vs_APVY)[iy+nAPVmaxY*i] );
	if( htemp->GetEntries() >= 100 ){
	  double Amean = htemp->GetMean();
	  double dAmean = htemp->GetMeanError();

	  double Arms = htemp->GetRMS();
	  double dArms = htemp->GetRMSError();

	  htemp->Fit("gaus","SQ","",Amean-Arms, Amean+Arms);
	 
	  double Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	  double dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);

	  double Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	  double dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	  //re-run the fit with a tighter range set by the result of the first fit:

	  htemp->Fit( "gaus", "SQ", "", Amean_fit - 2.0*Arms_fit, Amean_fit + 2.0*Arms_fit );

	  Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	  dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	  Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	  dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	  double Ryxtemp = (1.0-Amean_fit)/(1.0+Amean_fit);

	  double Aup = Amean_fit + dAmean_fit;
	  double Adown = Amean_fit - dAmean_fit;
	  
	  double Ryx_Aup = (1.0 - Aup)/(1.0 + Aup );
	  double Ryx_Adown = (1.0 - Adown)/(1.0 + Adown);

	  AsymALL.push_back( Amean_fit );
	  dAsymALL.push_back( dAmean_fit );

	  APVX_asymALL.push_back( -1 );
	  APVY_asymALL.push_back( iy );
	} else {
	  gainfit->FixParameter( iy );
	}
      } else {
	gainfit->FixParameter( iy );
      }
    }
    
    for( int ix = 0; ix<nAPVmaxX; ix++ ){
      
      for( int iy = 0; iy<nAPVmaxY; iy++ ){

	if( ix * 128 < nstripx_mod[i] && iy*128 < nstripy_mod[i] ){
	
	  TH1D *htemp = ( (TH1D*) (*hADCasym_vs_APVXY)[iy + nAPVmaxY*ix + nAPVmaxY*nAPVmaxX*i] );

	
	
	  weightX[iy][ix] = htemp->GetEntries();
	  weightY[ix][iy] = htemp->GetEntries();
	
	  if( htemp->GetEntries() >= 100 ){

	    countY[ix]++;
	    countX[iy]++;
	  
	    double Amean = htemp->GetMean();
	    double dAmean = htemp->GetMeanError();

	    double Arms = htemp->GetRMS();
	    double dArms = htemp->GetRMSError();

	    htemp->Fit("gaus","SQ","",Amean-Arms, Amean+Arms);
	 
	    double Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	    double dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);

	    double Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	    double dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	    //re-run the fit with a tighter range set by the result of the first fit:

	    htemp->Fit( "gaus", "SQ", "", Amean_fit - 2.0*Arms_fit, Amean_fit + 2.0*Arms_fit );

	    Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	    dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	    Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	    dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	    double Ryxtemp = (1.0-Amean_fit)/(1.0+Amean_fit);

	    double Aup = Amean_fit + dAmean_fit;
	    double Adown = Amean_fit - dAmean_fit;
	  
	    double Ryx_Aup = (1.0 - Aup)/(1.0 + Aup );
	    double Ryx_Adown = (1.0 - Adown)/(1.0 + Adown);

	  
	  
	    Asym.push_back( Amean_fit );
	    dAsym.push_back( dAmean_fit );
	    Ryx.push_back( Ryxtemp );
	    dRyx.push_back( 0.5*fabs(Ryx_Aup - Ryx_Adown) );

	    //ONLY populate the list of asymmetries IF we had a successful fit:
	    AsymALL.push_back( Amean_fit );
	    dAsymALL.push_back( dAmean_fit );

	    APVX_asymALL.push_back( ix );
	    APVY_asymALL.push_back( iy );
	  
	    if( count == 0 || dAmean_fit < minAsymerr ){
	      minAsymerr = dAmean_fit;
	      xAPV_ref = ix;
	      count++;
	    }
	  } else {
	    Asym.push_back( -10.0 );
	    dAsym.push_back( 1.0 );
	    Ryx.push_back( 1.0 );
	    dRyx.push_back( 1.0 );
	  }
	} else {
	  weightX[iy][ix] = 0.0;
	  weightY[ix][iy] = 0.0;
	}
      }
    }

    arglist[0] = 5000;
    arglist[1] = 0.1;
    gainfit->mnexcm("MIGRAD",arglist,2,ierflg);
    
    //First calculate all the Y gains relative to the reference x APV

    double Gmod = RelativeGainByModule[i];
    
    outfile << "mod_Ygain  " << i << "   " << nAPVmaxY << "     ";

    TString varname;
    varname.Form("%s.m%d.vgain = ",detname.Data(),i);

    outfile_db << varname;
    
    for( int iy=0; iy<nAPVmaxY; iy++ ){
      int yidx= iy + nAPVmaxY * xAPV_ref;
      
      // Ygain.push_back( Ryx[yidx] );
      // dYgain.push_back( dRyx[yidx] );

      double Gy, dGy;
      gainfit->GetParameter( iy, Gy, dGy );

      // d(1/x) = -1/x^2 dx = dx / x * (1/x)
      
      Ygain.push_back( 1.0/Gy );
      dYgain.push_back( dGy/(Gy*Gy) );
      
      //cout << "module " << i << ", Y APV " << iy << ", Relative gain = " << Ygain.back() << " +/- " << dYgain.back() << endl;

      outfile << Ygain.back() << "  ";

      if( iy*128 < nstripy_mod[i] ){
	outfile_db << Ygain.back() << "  ";
	Ygain_by_module[i].push_back( Ygain.back() );
      }
    }
    //Output file info again
    outfile << endl;
    outfile_db << endl;
    
    outfile << "mod_Xgain  " << i << "   " << nAPVmaxX << "     ";
    varname.Form("%s.m%d.ugain = ",detname.Data(),i);

    outfile_db << varname;
    
    for( int ix=0; ix<nAPVmaxX; ix++ ){
      //looks like it is a loop just to print out the above gain info
      double Gx, dGx;
      gainfit->GetParameter( ix + nAPVmaxY, Gx, dGx );
      Xgain.push_back( 1.0/Gx );
      dXgain.push_back( dGx/(Gx*Gx) );
      // if( ix != xAPV_ref ){
      // 	double sum_xgain = 0.0;
      // 	double sum_weights = 0.0;
      // 	//to get relative X gains, compute a weighted average of the relative X gains from
      // 	//the relative Y gains: 
      // 	for( int iy=0; iy<nAPVmaxY; iy++ ){
      // 	  int idx = iy + nAPVmaxY*ix;
      // 	  if( fabs( Asym[idx] ) <= 1.0 ){ //then we have successfully fit the asymmetry for this combination:
      // 	    double ygaintemp = Ygain[iy];

      // 	    double dygaintemp = dYgain[iy];

      // 	    double xgaintemp = Ygain[iy]/Ryx[idx];
      // 	    double dxgaintemp = fabs(xgaintemp)*sqrt(pow( dygaintemp/ygaintemp, 2 ) + pow( dRyx[idx]/Ryx[idx],2));

      // 	    sum_xgain += xgaintemp * pow( dxgaintemp, -2 );
      // 	    sum_weights += pow( dxgaintemp, -2 );
      // 	  }
      // 	}

      // 	if( sum_xgain > 0.0 ){
      // 	  Xgain.push_back( sum_xgain/sum_weights );
      // 	  dXgain.push_back( 1.0/sqrt(sum_weights) );
      // 	} else {
      // 	  Xgain.push_back( 1.0 );
      // 	  dXgain.push_back( 1.0 );
      // 	}
      // } else {
      // 	Xgain.push_back( 1.0 );
      // 	dXgain.push_back( 0.0 );
      // }
      //cout << "module " << i << ", X APV " << ix << ", Relative gain = " << Xgain.back() << " +/- " << dXgain.back() << endl;

      outfile << Xgain.back() << "  ";
      if( ix*128 < nstripx_mod[i] ){
	outfile_db << Xgain.back() << "  ";
	Xgain_by_module[i].push_back( Xgain.back() );
      }
    }
    outfile << endl;
    outfile_db << endl;

    //outfile_db << "# Module " << i << " average gain relative to target ADC of " << target_ADC << " = " << Gmod << endl;
  }
  //Define TGraph and Histogram for all gain coefficients. U/X and V/Y gain coefficient vs APV
  TGraph *Gain_APV_all = new TGraph();
  Gain_APV_all->SetName("gGainCoefficients_APV");
  Gain_APV_all->SetTitle("Gain Coefficients vs APV");
  Gain_APV_all->SetMarkerStyle(kFullDotLarge);
  Gain_APV_all->SetLineColor(0);
	
  TGraph *XGain_APV_all = new TGraph();
  XGain_APV_all->SetName("gUX_GainCoefficients_APV");
  XGain_APV_all->SetTitle("U/X Gain Coefficients vs APV");
  XGain_APV_all->SetMarkerStyle(kFullDotLarge);
  XGain_APV_all->SetLineColor(0);
	
  TGraph *YGain_APV_all = new TGraph();
  YGain_APV_all->SetName("gVY_GainCoefficients_APV");
  YGain_APV_all->SetTitle("V/Y Gain Coefficients vs APV");
  YGain_APV_all->SetMarkerStyle(kFullDotLarge);
  YGain_APV_all->SetLineColor(0);
	
  //gStyle->SetOptStat("neMRou");
  TH1D *Gain_histo_all = new TH1D("hGainCoefficient_Histo","",200,0,2);
  TH1D *XGain_histo_all = new TH1D("hUX_GainCoefficient_Histo","",200,0,2);	
  TH1D *YGain_histo_all = new TH1D("hVY_GainCoefficient_Histo","",200,0,2);
	
  int bestcount = 0;
  int mycountx = 0;
  int mycounty = 0;
  for(int mod = 0; mod < nmodules; mod++ ){
    vector<double> Xtemp = Xgain_by_module[mod];
    vector<double> Ytemp = Ygain_by_module[mod];
    int mymodx = Xtemp.size();
    int mymody = Ytemp.size();
    //Define TGraph and Histogram per module for gain coefficients
    TString mynamex;
    mynamex.Form("gUX_GainCoefficients_APV_mod%d",mod);
    TString mynamey;
    mynamey.Form("gVY_GainCoefficients_APV_mod%d",mod);
	
    TString myNameX;
    myNameX.Form("hUX_GainCoefficient_Histo_mod%d",mod);
    TString myNameY;
    myNameY.Form("hVY_GainCoefficient_Histo_mod%d",mod);
	
    TGraph *XGain_APV_mod = new TGraph(mymodx);
    XGain_APV_mod->SetName(mynamex);
    XGain_APV_mod->SetTitle(mynamex);
    XGain_APV_mod->SetMarkerStyle(kFullDotLarge);
    XGain_APV_mod->SetLineColor(0);

    TGraph *YGain_APV_mod = new TGraph(mymody);
    YGain_APV_mod->SetName(mynamey);
    YGain_APV_mod->SetTitle(mynamey);
    YGain_APV_mod->SetMarkerStyle(kFullDotLarge);
    YGain_APV_mod->SetLineColor(0);
	
    TH1D *XGain_histo_mod = new TH1D(myNameX,"",200,0,2);
    TH1D *YGain_histo_mod = new TH1D(myNameY,"",200,0,2);	
    for(int x = 0; x < Xtemp.size(); x++){
      double myXgain = Xtemp.at(x);
      //Fill both histos properly
	
      XGain_histo_all -> Fill(myXgain);
      XGain_APV_mod -> SetPoint(x,x+1,myXgain);
      XGain_histo_mod -> Fill(myXgain);
      Gain_histo_all ->Fill(myXgain);
      Gain_APV_all ->SetPoint(bestcount,bestcount+1,myXgain);
      XGain_APV_all -> SetPoint(mycountx,mycountx+1,myXgain);
		
      mycountx++;
      bestcount++;
    }
    for(int y = 0; y < Ytemp.size(); y++){
      double myYgain = Ytemp.at(y);
      //Fill both histos properly
               
      YGain_histo_all -> Fill(myYgain);
      YGain_APV_mod -> SetPoint(y,y+1,myYgain);
      YGain_histo_mod ->Fill(myYgain);
      Gain_histo_all ->Fill(myYgain);
      Gain_APV_all ->SetPoint(bestcount,bestcount+1,myYgain);
      YGain_APV_all -> SetPoint(mycounty,mycounty+1,myYgain);

      mycounty++;
      bestcount++;
    }
    XGain_APV_mod->Draw("AP");
    XGain_APV_mod->Write();
    YGain_APV_mod->Draw("AP");
    YGain_APV_mod->Write();       
	
  }
  Gain_APV_all->Draw("AP");
  Gain_APV_all->Write();
  XGain_APV_all->Draw("AP");
  XGain_APV_all->Write();
  YGain_APV_all->Draw("AP");
  YGain_APV_all->Write();
	

  




  //We will want a 2nd loop over the data to make plots of corrected ADC spectra:

  //int nAPVmax = 

  TH1D *hADCavg_allhits_corrected = new TH1D("hADCavg_allhits_corrected","",1500,0,30000);
  TH1D *hADCasym_allhits_corrected = new TH1D("hADCasym_allhits_corrected","",250,-1.01,1.01);
  TH2D *hADCavg_module_corrected = new TH2D("hADCavg_module_corrected","ADC Average vs Module After Gain;Module;0.5*(Clust U ADC + Clust V ADC)",nmodules,-0.5,nmodules-0.5,1500,0,30000);
  TH2D *hADCasym_module_corrected = new TH2D("hADCasym_module_corrected","ADC Asymm vs Module After Gain;Module;ADC Asymm",nmodules,-0.5,nmodules-0.5,250,-1.01,1.01);

  TH2D *hADCasym_vs_ADCavg_allhits_corrected = new TH2D("hADCasym_vs_ADCavg_allhits_corrected",";ADC average; ADC asym",1000,0,30000,250,-1.01,1.01);
  
  TH2D *hStripADCsumU_module_corrected = new TH2D("hStripADCsumU_module_corrected","",nmodules,-0.5,nmodules-0.5,1500,0,15000);
  TH2D *hStripADCmaxU_module_corrected = new TH2D("hStripADCmaxU_module_corrected","",nmodules,-0.5,nmodules-0.5,1500,0,3000);
  TH2D *hStripADCsumV_module_corrected = new TH2D("hStripADCsumV_module_corrected","",nmodules,-0.5,nmodules-0.5,1500,0,15000);
  TH2D *hStripADCmaxV_module_corrected = new TH2D("hStripADCmaxV_module_corrected","",nmodules,-0.5,nmodules-0.5,1500,0,3000);

  TH2D *hADC_UV_allhits_corrected = new TH2D("hADC_UV_allhits_corrected","Cluster sum ;ADCU;ADCV",250,0,25000,250,0,25000);
  TH2D *hADC_UVmaxstrip_allhits_corrected = new TH2D("hADC_UVmaxstrip_allhits_corrected","Max Strip sum;ADCU;ADCV",250,0,15000,250,0,15000);
  TH2D *hADC_UVmaxsamp_allhits_corrected = new TH2D("hADC_UVmaxsamp_allhits_corrected","Max Strip max sample;ADCU;ADCV",250,0,3000,250,0,3000);
  
  TClonesArray *hADC_UV_allhits_corrected_mod = new TClonesArray( "TH2D", nmodules);
  TClonesArray *hADC_UVmaxstrip_allhits_corrected_mod = new TClonesArray( "TH2D", nmodules );
  TClonesArray *hADC_UVmaxsamp_allhits_corrected_mod = new TClonesArray( "TH2D", nmodules );
  TClonesArray *hADCasym_vs_ADCavg_allhits_corrected_mod = new TClonesArray( "TH2D", nmodules );
  for(int j =0; j<nmodules;j++){
    TString datname, Datname, datName;
    datname.Form("hADC_UVmaxstrip_allhits_corrected_mod%d",j);
    Datname.Form("hADC_UV_allhits_corrected_mod%d",j);
    datName.Form("hADC_UVmaxsamp_allhits_corrected_mod%d",j);
    new( (*hADC_UVmaxstrip_allhits_corrected_mod)[j] ) TH2D(datname.Data(),"Max Strip sum;ADCU;ADCV",250,0,15000,250,0,15000);
    new( (*hADC_UV_allhits_corrected_mod)[j] ) TH2D(Datname.Data(),"Cluster sum ;ADCU;ADCV",250,0,25000,250,0,25000);
    new( (*hADC_UVmaxsamp_allhits_corrected_mod)[j] ) TH2D(datName.Data(),"Max Strip max sample;ADCU;ADCV",250,0,3000,250,0,3000);	new( (*hADCasym_vs_ADCavg_allhits_corrected_mod)[j] ) TH2D( datname.Format( "hADCasym_vs_ADCavg_corrected_mod%d", j ), ";cluster ADC average; ADC asymmetry",250,0,25000,250,-1.01,1.01);
  } 


  nevent = 0;

  treenum = 0;
  currenttreenum = 0;
  
  vector<double> AppliedFactorsX(nmodules,0.0);
  vector<double> AppliedFactorsY(nmodules,0.0);
  vector<double> HitCounts(nmodules,0.0);

  double AppliedFactorsAllX = 0.0;
  double AppliedFactorsAllY = 0.0;
  double HitCountsAll = 0.0;
  
  //cout << "starting event loop:" << endl;
  while( C->GetEntry( nevent++ ) ){

    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passedcut = GlobalCut->EvalInstance(0) != 0;
    
    if( ngoodhits > MAXNHITS ) continue;
    if( ntracks > MAXNHITS ) continue;
    
    if( nevent % 100000 == 0 ) cout << "event " << nevent << endl;
    //loop over hits:
    int itrack = int(besttrack);
    //cout << "itrack = " << itrack << endl;

    int NTRACKS = int(ntracks);

    //cout << "ntracks = " << NTRACKS << endl;
    
    if( passedcut && NTRACKS > 0 ){

      //      if( trackChi2NDF[itrack] < chi2cut && tracknhits[itrack] > 3 ){
      int nhits = int(ngoodhits);
      //	cout << "nhits = " << nhits << endl;
	
      for( int ihit=0; ihit<nhits; ihit++ ){
	int tridx = int( hit_trackindex[ihit] );

	//cout << "ihit, tridx = " << ihit << ", " << tridx;
	  
	//if( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) >= ADCcut && tridx == itrack ){
	//if( tridx == itrack )
	//if( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) >= ADCcut && hit_ADCU[ihit] >= ADCcut && hit_ADCV[ihit] >= ADCcut && tridx == itrack && hit_nstripu[ihit]>1 && hit_nstripv[ihit]>1 && hit_ccor_clust[ihit] >= ccor_cut && fabs(hit_deltat[ihit]) <= deltat_cut ){
	if( hit_ADCU[ihit] >= ADCcut && hit_ADCV[ihit] >= ADCcut && tridx == itrack && hit_nstripu[ihit]>1 && hit_nstripv[ihit]>1 && hit_ccor_clust[ihit] >= ccor_cut && fabs(hit_deltat[ihit]) <= deltat_cut ){
	  //if( tridx == itrack ){
	  // cout << ", ADCavg[ihit] = " << hit_ADCavg[ihit]
	  // 	 << ", ADCasym[ihit] = " << hit_ADCasym[ihit]
	  // 	 << ", (nstripu,nstripv) = (" << hit_nstripu[ihit] << ", " << hit_nstripv[ihit] << ")"
	  // 	 << ", (ustriplo,ustriphi,ustripmax)=(" << hit_ustriplo[ihit] << ", " << hit_ustriphi[ihit] << ", " << hit_ustripmax[ihit] << ")"
	  // 	 << ", (vstriplo,vstriphi,vstripmax)=(" << hit_vstriplo[ihit] << ", " << hit_vstriphi[ihit] << ", " << hit_vstripmax[ihit] << ")"
	  // 	 << ", module = " << hit_module[ihit] << endl;
	  // hADCavg_allhits->Fill( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) );
	  // hADCavg_module->Fill( hit_module[ihit], 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]+hit_ADCV[ihit]/hit_Vgain[ihit]) );
	  // hADCasym_module->Fill( hit_module[ihit], hit_ADCasym[ihit] );
	  // hNstripX_module->Fill( hit_module[ihit], hit_nstripu[ihit] );
	  // hNstripY_module->Fill( hit_module[ihit], hit_nstripv[ihit] );

	  // hStripADCsumU_module->Fill( hit_module[ihit], hit_ADCmaxstripU[ihit] );
	  // hStripADCsumV_module->Fill( hit_module[ihit], hit_ADCmaxstripV[ihit] );

	  // hStripADCmaxU_module->Fill( hit_module[ihit], hit_ADCmaxsampU[ihit] );
	  // hStripADCmaxV_module->Fill( hit_module[ihit], hit_ADCmaxsampV[ihit] );

	  int ixlo = hit_ustriplo[ihit];
	  int ixhi = hit_ustriphi[ihit];
	  int ixmax = hit_ustripmax[ihit];

	  int iylo = hit_vstriplo[ihit];
	  int iyhi = hit_vstriphi[ihit];
	  int iymax = hit_vstripmax[ihit];

	  int xAPVmax = ixmax/128;
	  int yAPVmax = iymax/128;
	  int xAPVlo = ixlo/128;
	  int yAPVlo = iylo/128;

	  int xAPVhi = ixhi/128;
	  int yAPVhi = iyhi/128;

	  int module = int(hit_module[ihit]);
	  if( xAPVlo == xAPVmax && xAPVhi == xAPVmax &&
	      yAPVlo == yAPVmax && yAPVhi == yAPVmax &&
	      hit_nstripu[ihit] >= 2 && hit_nstripv[ihit] >= 2 &&
	      xAPVmax*128 < nstripx_mod[module] && yAPVmax*128<nstripy_mod[module] ){

	    double Xgaintemp = Xgain_by_module[module][xAPVmax];
	    double Ygaintemp = Ygain_by_module[module][yAPVmax];

	    AppliedFactorsX[module] += Xgaintemp;
	    AppliedFactorsY[module] += Ygaintemp;
	    HitCounts[module] += 1.0;

	    AppliedFactorsAllX += Xgaintemp;
	    AppliedFactorsAllY += Ygaintemp;
	    HitCountsAll += 1.0;
	      
	    hADCavg_allhits_corrected->Fill( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) );
	    hADCavg_module_corrected->Fill( hit_module[ihit], 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) );
	    hADCasym_APV_corrected->Fill(yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*hit_module[ihit],(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp-hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) / (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp));
	    
	    hADCasym_allhits_corrected->Fill( (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp-hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) /
					      (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) );

	    hADCasym_module_corrected->Fill( module,
					     (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp-hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) /
					     (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) );

	    hADCasym_vs_ADCavg_allhits_corrected->Fill( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp),
							(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp-hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) /
							(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) );

	    // hADCasym_module->Fill( hit_module[ihit], hit_ADCasym[ihit] );
	    // hNstripX_module->Fill( hit_module[ihit], hit_nstripu[ihit] );
	    // hNstripY_module->Fill( hit_module[ihit], hit_nstripv[ihit] );

	    hStripADCsumU_module_corrected->Fill( hit_module[ihit], hit_ADCmaxstripU[ihit]*Xgaintemp );
	    hStripADCsumV_module_corrected->Fill( hit_module[ihit], hit_ADCmaxstripV[ihit]*Ygaintemp );

	    hStripADCmaxU_module_corrected->Fill( hit_module[ihit], hit_ADCmaxsampU[ihit]*Xgaintemp );
	    hStripADCmaxV_module_corrected->Fill( hit_module[ihit], hit_ADCmaxsampV[ihit]*Ygaintemp );

	     
	    ( (TH2D*) (*hADC_UVmaxstrip_allhits_corrected_mod)[module] )->Fill( hit_ADCmaxstripU[ihit]*Xgaintemp, hit_ADCmaxstripV[ihit]*Ygaintemp );
	    ( (TH2D*) (*hADC_UV_allhits_corrected_mod)[module] )->Fill(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp, hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp );
	    ( (TH2D*) (*hADC_UVmaxsamp_allhits_corrected_mod)[module] )->Fill(hit_ADCmaxsampU[ihit]*Xgaintemp, hit_ADCmaxsampV[ihit]*Ygaintemp);

	    ( (TH2D*) (*hADCasym_vs_ADCavg_allhits_corrected_mod)[module] )->Fill( 0.5*(hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp), (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp-hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) / (hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp+hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp) );
	    
	    hADC_UV_allhits_corrected->Fill( hit_ADCU[ihit]/hit_Ugain[ihit]*Xgaintemp, hit_ADCV[ihit]/hit_Vgain[ihit]*Ygaintemp );
	    hADC_UVmaxstrip_allhits_corrected->Fill( hit_ADCmaxstripU[ihit]*Xgaintemp, hit_ADCmaxstripV[ihit]*Ygaintemp );
	    hADC_UVmaxsamp_allhits_corrected->Fill( hit_ADCmaxsampU[ihit]*Xgaintemp, hit_ADCmaxsampV[ihit]*Ygaintemp );
	      
	    // if( xAPVlo == xAPVmax && xAPVhi == xAPVmax &&
	    // 	yAPVlo == yAPVmax && yAPVhi == yAPVmax &&
	    // 	hit_nstripu[ihit] >= 2 && hit_nstripv[ihit] >= 2 &&
	    // 	xAPVmax*128 < nstripx_mod[module] && yAPVmax*128<nstripy_mod[module] ){

	    //   // cout << "filling histograms " << endl;

	    //   // cout << "apvxy hist index = " << yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*module << endl;
	    //   // cout << "apvx hist index = " << xAPVmax + nAPVmaxX*module << endl;
	    //   // cout << "apvy hist index = " << yAPVmax + nAPVmaxY*module << endl;
	      
	    //   ( (TH1D*) (*hADCasym_vs_APVXY)[yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*module] )->Fill( hit_ADCasym[ihit] );

	    //   ( (TH1D*) (*hADCasym_vs_APVX)[xAPVmax + nAPVmaxX*module] )->Fill( hit_ADCasym[ihit] );
	    //   ( (TH1D*) (*hADCasym_vs_APVY)[yAPVmax + nAPVmaxY*module] )->Fill( hit_ADCasym[ihit] );
	    //   //cout << "done histogram fill" << endl;
	  }
	  // }
	}
	////
      }
    }
  }
  cout << "Event " << nevent << " done" << endl;


  //Adjust module ABSOLUTE gain AFTER correction for RELATIVE gain:

  outfile_db << endl << "# Module average gains relative to target ADC peak position of " << target_ADC << endl;

  TFitResultPtr ADCfit_all = hADCavg_allhits_corrected->Fit("landau","qS","",2000,25000.0);

  MPV_all = ( (TF1*) hADCavg_allhits_corrected->GetListOfFunctions()->FindObject("landau") )->GetParameter("MPV");

  double thresh_sample = 100.0;
  double thresh_strip = 200.0;
  double thresh_cluster = 500.0;
  
  double avgfactorXall = AppliedFactorsAllX/HitCountsAll;
  double avgfactorYall = AppliedFactorsAllY/HitCountsAll; 
  double Gall = 2.0/(avgfactorXall + avgfactorYall);
  
  TGraph *thresh_sample_all = new TGraph(nmodules);
  TGraph *thresh_strip_all = new TGraph(nmodules);
  TGraph *thresh_cluster_all = new TGraph(nmodules);

  thresh_sample_all->SetName("gThresh_sample_all");
  thresh_sample_all->SetTitle("Sample Threshold vs Module");
  thresh_sample_all->SetMarkerStyle(kFullDotLarge);
  thresh_sample_all->SetLineColor(0);
 
  thresh_strip_all->SetName("gThresh_strip_all");
  thresh_strip_all->SetTitle("Strip Threshold vs Module");
  thresh_strip_all->SetMarkerStyle(kFullDotLarge);
  thresh_strip_all->SetLineColor(0);

  thresh_cluster_all->SetName("gThresh_cluster_all");
  thresh_cluster_all->SetTitle("Cluster Threshold vs Module");
  thresh_cluster_all->SetMarkerStyle(kFullDotLarge);
  thresh_cluster_all->SetLineColor(0);

 


  for( int i=0; i<nmodules; i++ ){
    //double Gmod = RelativeGainByModule[i]*MPV_all/target_ADC;

    double Gmod = 1.0;

    TString hname;
    hname.Form("hADCavg_corrected_module%d",i);
    htemp = hADCavg_module_corrected->ProjectionY( hname.Data(), i+1, i+1 );

    if( htemp->GetEntries() >= 5000 ){
      TFitResultPtr ADCfit_module = htemp->Fit("landau","qS","",2000,25000.0);
      
      double MPV_mod = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("landau") ) )->GetParameter("MPV");
      
      // Gmod = MPV_mod / target_ADC;

      double avgfactorX = AppliedFactorsX[i] / HitCounts[i];
      double avgfactorY = AppliedFactorsY[i] / HitCounts[i]; 
      
      Gmod = 2.0 / (avgfactorX + avgfactorY );

      double thresh_sample_mod = thresh_sample / (Gall * MPV_all) * (Gmod * MPV_mod );
      double thresh_strip_mod = thresh_strip / (Gall * MPV_all) * (Gmod * MPV_mod );
      double thresh_cluster_mod = thresh_cluster / (Gall * MPV_all) * (Gmod * MPV_mod );
      
      thresh_sample_all->SetPoint(i,i+1,thresh_sample_mod);
      thresh_strip_all->SetPoint(i,i+1,thresh_strip_mod);
      thresh_cluster_all->SetPoint(i,i+1,thresh_cluster_mod);

      //cout << "# Module " << i << " average gain relative to target ADC of " << target_ADC << " = " << Gmod << endl;
      TString dbentry;
      dbentry.Form("%s.m%d.modulegain = %g",detname.Data(),i,Gmod);
      outfile_db << dbentry << endl;
      //cout << dbentry << endl;
      /*  Thresholds are no longer set by this script
      dbentry.Form( "%s.m%d.threshold_sample = %g", detname.Data(), i, thresh_sample_mod );

      // outfile_db << dbentry << endl;
      // cout << dbentry << endl;

      dbentry.Form( "%s.m%d.threshold_stripsum = %g", detname.Data(), i, thresh_strip_mod );
      // outfile_db << dbentry << endl;
      // cout << dbentry << endl;

      dbentry.Form( "%s.m%d.threshold_clustersum = %g", detname.Data(), i, thresh_cluster_mod );
      // // // // outfile_db << dbentry << endl;
      // // // // cout << dbentry << endl;

      // outfile_db << endl << endl;

      
      outfile_db << endl << endl;
      */
    } else {
      thresh_sample_all->SetPoint(i,i+1,thresh_sample);
      thresh_strip_all->SetPoint(i,i+1,thresh_strip);
      thresh_cluster_all->SetPoint(i,i+1,thresh_cluster);
      
      TString dbentry;
      dbentry.Form("%s.m%d.modulegain = %g",detname.Data(),i,1.0);
      outfile_db << dbentry << endl;
      //cout << dbentry << endl;
      /*  Thresholds are no longer set by this script
      dbentry.Form( "%s.m%d.threshold_sample = %g", detname.Data(), i, thresh_sample );

      // outfile_db << dbentry << endl;
      // cout << dbentry << endl;

      dbentry.Form( "%s.m%d.threshold_stripsum = %g", detname.Data(), i, thresh_strip );
      // outfile_db << dbentry << endl;
      // cout << dbentry << endl;

      dbentry.Form( "%s.m%d.threshold_clustersum = %g", detname.Data(), i, thresh_cluster );
      // outfile_db << dbentry << endl;
      // cout << dbentry << endl;

      // outfile_db << endl << endl; */
      
    }
  }
  thresh_sample_all ->Draw("AP");
  thresh_sample_all ->Write();
  thresh_strip_all ->Draw("AP");
  thresh_strip_all ->Write();
  thresh_cluster_all ->Draw("AP");
  thresh_cluster_all ->Write();
  outfile << endl;
  outfile_db << endl;
  fout->Write();

  TString outfilepdf(outfname);
  outfilepdf.ReplaceAll(".root",".pdf");

  // Plot mean and standard deviation for after gain matching
  TGraphErrors* g2 = new TGraphErrors();
  point = 0;
  for( int i=0; i<nmodules; i++ ){
    TString hnametemp;
    hnametemp.Form("ADCasym_module_corrected%d",i);
    htemp = hADCasym_module_corrected->ProjectionY(hnametemp.Data(), i+1, i+1 );

    if( htemp->GetEntries() >= 100 ){
      
      double max = htemp->GetMaximum();
      int binmax,binlo,binhi;
      
      binmax = htemp->GetMaximumBin();
      binlo = binmax;
      binhi = binmax;
      
      while( htemp->GetBinContent(binlo) >= 0.5*max && binlo > 1 ){binlo--;}
      while( htemp->GetBinContent(binhi) >= 0.5*max && binhi < 500 ){binhi++; }
      
      htemp->Fit("gaus","qS","",htemp->GetBinCenter(binlo),htemp->GetBinCenter(binhi) );
      
      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double x = hADCasym_module->GetXaxis()->GetBinCenter(i+1); // x-position
      double y = fitfunc->GetParameter(1);  // mean from Gaussian
      double ex = 0;           // no x error
      double ey = fitfunc->GetParameter(2); // sigma from fit

      g2->SetPoint(point, x, y);
      g2->SetPointError(point, ex, ey);
      point++;
    }
  }


  ////// Here we create all the plots for the output PDF ///////
  TCanvas *c2 = new TCanvas("c2","",1200,1000);
  c2->Divide(2,2);

  c2->cd(1);
  hADCasym_module->Draw("colz");
  g->SetMarkerStyle(21);
  g->SetMarkerColor(kRed);
  g->Draw("P SAME");

  c2->cd(2);
  hADCasym_module_corrected->Draw("colz");
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(kRed);
  g2->Draw("P SAME");

  c2->cd(3);
  hADCasym_APV->Draw("colz");

  c2->cd(4);
  hADCasym_APV_corrected->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  Gain_APV_all->Draw("AP");

  c2->Print(outfilepdf + "(");
  c3->Print(outfilepdf + ")");
}

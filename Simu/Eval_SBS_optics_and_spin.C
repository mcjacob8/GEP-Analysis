// ------------------------------------------------------------------------------------- //
// Drop-in replacement for Fit_SBS_optics_and_spin.C that does NOT fit. Instead it loads //
// coefficients from the text files written by a previous fit and then runs the same     //
// reconstruction loop, producing the same output ROOT file (Tout tree with identical    //
// branches + identical histograms) so existing comparison machinery works unchanged.    //
//                                                                                       //
// Usage: same as the fit script:                                                        //
//   Eval_SBS_optics_and_spin("config.cfg","outfile.root")                               //
//                                                         //                            //
// Config file: same layout as for the fit script (file list / endlist / cut / endcut /  //
// key-value block / endconfig). Extra keys understood in the key-value block:           //
//                                                                                       //
//   SBStheta_central   18.6127              # default is Kin. 3 value                   //
//   optics_coeff_file  /path/to/fitA_optics_coeff.txt          (required)               //
//   spin_coeff_file    /path/to/fitA_spin_coeff.txt            (optional)               //
//   forward_coeff_file /path/to/fitA_forward_optics_coeff.txt  (optional)               //
//                                                                                       //
// Fit-only keys (order_optics, order_spin, use_xtar_flag, trkrdist, fix_ytar, fix_yptar)//
// are accepted and ignored: the orders are inferred from the coefficient files, and the //
// fix_ytar/fix_yptar modifications are already baked into the written coefficients.     //
//                                                                                       //
// Terms are matched by their exponent columns (xfp ypfp... as printed), so a file of any//
// order can be loaded; any term not present in a file gets coefficient 0.               //
// ---------                                                                             //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 8-30-2026                                //
// ---------                                                                             //
// ** Do not tamper with this sticker! Log any updates to the script above.              //
// ------------------------------------------------------------------------------------- //

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "/work/halla/sbs/mcjacob/g4sbs/g4sbs/root_macros/gep_optics_tree.C"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTreeFormula.h"
#include "TCut.h"
#include "TObjArray.h"
#include "TObjString.h"

const double Mp = 0.938272046; //GeV/c^2
const double mu_p = 2.792847356; //mu_N
const double kappa_p = mu_p - 1.0; //anomalous magnetic moment
const double PI = TMath::Pi();
const double SBS_thetabend = 0.0*PI/180.0; //central bend angle (same as fit script)

// Number of terms for a given order (same nested loops as the fit script)
static int NTermsForOrder( int order ){
  int n = 0;
  for( int i=0; i<=order; i++ )
    for( int j=0; j<=order-i; j++ )
      for( int k=0; k<=order-i-j; k++ )
	for( int l=0; l<=order-i-j-k; l++ )
	  for( int m=0; m<=order-i-j-k-l; m++ ) n++;
  return n;
}

// Index of the term with printed exponents (m,l,k,j,i) within the fit script's loop ordering
static int TermIndex( int order, int em, int el, int ek, int ej, int ei ){
  int idx = 0;
  for( int i=0; i<=order; i++ )
    for( int j=0; j<=order-i; j++ )
      for( int k=0; k<=order-i-j; k++ )
	for( int l=0; l<=order-i-j-k; l++ )
	  for( int m=0; m<=order-i-j-k-l; m++ ){
	    if( i==ei && j==ej && k==ek && l==el && m==em ) return idx;
	    idx++;
	  }
  return -1;
}

// Read all data rows (ncoef coefficient columns + 5 exponent columns). Comment lines and
// rows with any other number of columns (e.g. "number of terms" lines) are skipped.
static bool ReadCoeffRows( const TString &fname, int ncoef, vector< vector<double> > &rows, int &order ){
  rows.clear();
  order = 0;
  ifstream in( fname.Data() );
  if( !in ){
    cout << "ERROR: cannot open coefficient file " << fname << endl;
    return false;
  }
  string line;
  while( getline(in,line) ){
    size_t first = line.find_first_not_of(" \t");
    if( first == string::npos || line[first] == '#' ) continue;
    istringstream ss(line);
    vector<double> vals;
    double x;
    while( ss >> x ) vals.push_back(x);
    if( (int)vals.size() != ncoef+5 ) continue;
    int tot = 0;
    for( int c=0; c<5; c++ ) tot += (int) TMath::Nint(vals[ncoef+c]);
    if( tot > order ) order = tot;
    rows.push_back(vals);
  }
  cout << "read " << rows.size() << " terms (max order " << order << ") from " << fname << endl;
  return !rows.empty();
}

// Fill arrays b[icol][icoeff] (icol<ncoef) from rows
static void FillCoeffs( const vector< vector<double> > &rows, int ncoef, int order, vector< vector<double> > &b ){
  b.assign( ncoef, vector<double>( NTermsForOrder(order), 0.0 ) );
  for( size_t r=0; r<rows.size(); r++ ){
    int em = (int) TMath::Nint(rows[r][ncoef+0]);
    int el = (int) TMath::Nint(rows[r][ncoef+1]);
    int ek = (int) TMath::Nint(rows[r][ncoef+2]);
    int ej = (int) TMath::Nint(rows[r][ncoef+3]);
    int ei = (int) TMath::Nint(rows[r][ncoef+4]);
    int idx = TermIndex( order, em, el, ek, ej, ei );
    if( idx < 0 ){
      cout << "WARNING: term with exponents " << em << el << ek << ej << ei << " not representable at order " << order << endl;
      continue;
    }
    for( int c=0; c<ncoef; c++ ) b[c][idx] = rocodews[r][c];
  }
}

void Eval_SBS_optics_and_spin( const char *configfilename, const char *outfilename="SBS_optics_fitresult.root" ){

  double SBStheta_central = 18.6127; //Kin. 3 value (same default as fit script)
  double trkrdist = 3.685; //only used in an unused first-order variable in the fit script

  TString opticsfile_in = "", spinfile_in = "", forwardfile_in = "";

  TChain *C = new TChain("Tout");

  ifstream infile(configfilename);
  if( !infile ) return;

  TString currentline;

  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" \t");
      if( tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();

	if( skey == "SBStheta_central" ) SBStheta_central = sval.Atof();
	if( skey == "trkrdist" ) trkrdist = sval.Atof();
	if( skey == "optics_coeff_file" ) opticsfile_in = sval;
	if( skey == "spin_coeff_file" ) spinfile_in = sval;
	if( skey == "forward_coeff_file" ) forwardfile_in = sval;
      }
      delete tokens;
    }
  }
  SBStheta_central *= PI/180.0;

  if( opticsfile_in == "" ){
    cout << "ERROR: config must give optics_coeff_file" << endl;
    return;
  }

  // ---------------- load coefficients (replaces all of the fitting) ----------------
  vector< vector<double> > rows_optics, rows_forward, rows_spin;
  int order_optics = 0, order_forward = 0, order_spin = 0;

  if( !ReadCoeffRows( opticsfile_in, 4, rows_optics, order_optics ) ) return;
  if( forwardfile_in != "" ){
    if( !ReadCoeffRows( forwardfile_in, 4, rows_forward, order_forward ) ) return;
  } else {
    cout << "WARNING: no forward_coeff_file given; forward xfp/yfp/xpfp/ypfp will be zero" << endl;
  }
  if( spinfile_in != "" ){
    if( !ReadCoeffRows( spinfile_in, 9, rows_spin, order_spin ) ) return;
  } else {
    cout << "WARNING: no spin_coeff_file given; spin matrices will be the pure dipole approximation" << endl;
  }

  // optics and forward arrays share one ordering (as in the fit script), so use a common order
  if( order_forward > order_optics ) order_optics = order_forward;
  int ncoeff_optics = NTermsForOrder( order_optics );
  int ncoeff_spin = NTermsForOrder( order_spin );
  cout << "order_optics = " << order_optics << " (" << ncoeff_optics << " terms), "
       << "order_spin = " << order_spin << " (" << ncoeff_spin << " terms)" << endl;

  vector< vector<double> > bo, bf, bs;
  FillCoeffs( rows_optics, 4, order_optics, bo );
  if( !rows_forward.empty() ) FillCoeffs( rows_forward, 4, order_optics, bf );
  else bf.assign( 4, vector<double>( ncoeff_optics, 0.0 ) );
  if( !rows_spin.empty() ) FillCoeffs( rows_spin, 9, order_spin, bs );
  else bs.assign( 9, vector<double>( ncoeff_spin, 0.0 ) );

  // same names as in the fit script so the reconstruction loop below is unchanged
  TVectorD b_xptar(ncoeff_optics), b_yptar(ncoeff_optics), b_ytar(ncoeff_optics), b_pinv(ncoeff_optics);
  TVectorD b_vz(ncoeff_optics);
  TVectorD b_xfp(ncoeff_optics), b_yfp(ncoeff_optics), b_xpfp(ncoeff_optics), b_ypfp(ncoeff_optics);
  for( int ic=0; ic<ncoeff_optics; ic++ ){
    b_xptar(ic) = bo[0][ic];
    b_yptar(ic) = bo[1][ic];
    b_ytar(ic)  = bo[2][ic];
    b_pinv(ic)  = bo[3][ic];
    b_vz(ic)    = 0.0; //vz coefficients are not written by the fit script
    b_xfp(ic)   = bf[0][ic];
    b_yfp(ic)   = bf[1][ic];
    b_xpfp(ic)  = bf[2][ic];
    b_ypfp(ic)  = bf[3][ic];
  }

  TVectorD b_Sxx(ncoeff_spin), b_Sxy(ncoeff_spin), b_Sxz(ncoeff_spin);
  TVectorD b_Syx(ncoeff_spin), b_Syy(ncoeff_spin), b_Syz(ncoeff_spin);
  TVectorD b_Szx(ncoeff_spin), b_Szy(ncoeff_spin), b_Szz(ncoeff_spin);
  for( int ic=0; ic<ncoeff_spin; ic++ ){
    b_Sxx(ic) = bs[0][ic]; b_Sxy(ic) = bs[1][ic]; b_Sxz(ic) = bs[2][ic];
    b_Syx(ic) = bs[3][ic]; b_Syy(ic) = bs[4][ic]; b_Syz(ic) = bs[5][ic];
    b_Szx(ic) = bs[6][ic]; b_Szy(ic) = bs[7][ic]; b_Szz(ic) = bs[8][ic];
  }

  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C );
  gep_optics_tree *T = new gep_optics_tree(C);

  // ---------------- output (identical to the fit script) ----------------
  TFile *fout = new TFile(outfilename,"RECREATE");
  TTree *Tout = new TTree("Tout","SBS optics fit result");

  double xfp,yfp,xpfp,ypfp,xfprecon,yfprecon,xpfprecon,ypfprecon;
  double xptar,yptar,ytar,xtar,p,thetabend,vz;
  double xptarrecon,yptarrecon,ytarrecon,xtarrecon,precon,thetabendrecon,vzrecon_ytar,vzrecon_direct;
  double xfpforward,yfpforward,xpfpforward,ypfpforward;

  Tout->Branch("xfp",&xfp);
  Tout->Branch("yfp",&yfp);
  Tout->Branch("xpfp",&xpfp);
  Tout->Branch("ypfp",&ypfp);
  Tout->Branch("xfprecon",&xfprecon);
  Tout->Branch("yfprecon",&yfprecon);
  Tout->Branch("xpfprecon",&xpfprecon);
  Tout->Branch("ypfprecon",&ypfprecon);
  Tout->Branch("xptar",&xptar);
  Tout->Branch("yptar",&yptar);
  Tout->Branch("ytar",&ytar);
  Tout->Branch("xtar",&xtar);
  Tout->Branch("p",&p);
  Tout->Branch("thetabend",&thetabend);
  Tout->Branch("vz", &vz);
  Tout->Branch("xptarrecon",&xptarrecon);
  Tout->Branch("yptarrecon",&yptarrecon);
  Tout->Branch("ytarrecon",&ytarrecon);
  Tout->Branch("xtarrecon",&xtarrecon);
  Tout->Branch("precon",&precon);
  Tout->Branch("thetabendrecon",&thetabendrecon);
  Tout->Branch("vzrecon_ytar", &vzrecon_ytar);
  Tout->Branch("vzrecon_direct", &vzrecon_direct);
  Tout->Branch("xfpforward",&xfpforward);
  Tout->Branch("yfpforward",&yfpforward);
  Tout->Branch("xpfpforward",&xpfpforward);
  Tout->Branch("ypfpforward",&ypfpforward);

  double Pxtg, Pytg, Pztg;
  double Pxfp, Pyfp, Pzfp;
  double Pxfpdipole,Pyfpdipole,Pzfpdipole;
  double Pxfprecon,Pyfprecon,Pzfprecon;
  Tout->Branch("Pxtg",&Pxtg);
  Tout->Branch("Pytg",&Pytg);
  Tout->Branch("Pztg",&Pztg);

  Tout->Branch("Pxfp",&Pxfp);
  Tout->Branch("Pyfp",&Pyfp);
  Tout->Branch("Pzfp",&Pzfp);

  Tout->Branch("Pxfprecon",&Pxfprecon);
  Tout->Branch("Pyfprecon",&Pyfprecon);
  Tout->Branch("Pzfprecon",&Pzfprecon);

  Tout->Branch("Pxfpdipole",&Pxfpdipole);
  Tout->Branch("Pyfpdipole",&Pyfpdipole);
  Tout->Branch("Pzfpdipole",&Pzfpdipole);
  double det;
  Tout->Branch("Determinant",&det);

  TH1D *hxptardiff = new TH1D("hxptardiff","",250,-0.01,0.01);
  TH1D *hyptardiff = new TH1D("hyptardiff","",250,-0.01,0.01);
  TH1D *hytardiff = new TH1D("hytardiff","",250,-0.02,0.02);
  TH1D *hpdiff = new TH1D("hpdiff","",250,-0.1,0.1);
  TH1D *hthetabend_diff = new TH1D("hthetabend_diff","",250,-.01,.01);

  TH1D *hxfpdiff = new TH1D("hxfpdiff","",250,-0.02,0.02);
  TH1D *hyfpdiff = new TH1D("hyfpdiff","",250,-0.02,0.02);
  TH1D *hxpfpdiff = new TH1D("hxpfpdiff","",250,-0.01,0.01);
  TH1D *hypfpdiff = new TH1D("hypfpdiff","",250,-0.01,0.01);

  TH2D *hxptar_recon_vs_true = new TH2D("hxptar_recon_vs_true","",250,-0.3,0.3,250,-0.3,0.3);
  TH2D *hyptar_recon_vs_true = new TH2D("hyptar_recon_vs_true","",250,-0.1,0.1,250,-0.1,0.1);
  TH2D *hytar_recon_vs_true = new TH2D("hytar_recon_vs_true","",250,-0.1,0.1,250,-0.1,0.1);
  TH2D *hp_recon_vs_true = new TH2D("hp_recon_vs_true","",250,0.0,12.5,250,0.0,12.5);
  TH2D *hthetabend_recon_vs_true = new TH2D("hthetabend_recon_vs_true","",250,0.0,0.5,250,0.0,0.5);
  TH2D *hpdiff_p = new TH2D("hpdiff_p","",250,0.0,12.5,250,-0.1,0.1);

  TH2D *hxfp_recon_vs_true = new TH2D("hxfp_recon_vs_true","",250,-0.8,0.8,250,-0.8,0.8);
  TH2D *hyfp_recon_vs_true = new TH2D("hyfp_recon_vs_true","",250,-0.25,0.25,250,-0.25,0.25);
  TH2D *hxpfp_recon_vs_true = new TH2D("hxpfp_recon_vs_true","",250,-0.4,0.4,250,-0.4,0.4);
  TH2D *hypfp_recon_vs_true = new TH2D("hypfp_recon_vs_true","",250,-0.1,0.1,250,-0.1,0.1);

  TH1D *hPxfpdiff = new TH1D("hPxfpdiff","",250,-0.05,0.05);
  TH1D *hPyfpdiff = new TH1D("hPyfpdiff","",250,-0.05,0.05);
  TH1D *hPzfpdiff = new TH1D("hPzfpdiff","",250,-0.05,0.05);

  TH2D *hPxfp_recon_vs_true = new TH2D("hPxfp_recon_vs_true","",250,-1.05,1.05,250,-1.05,1.05);
  TH2D *hPyfp_recon_vs_true = new TH2D("hPyfp_recon_vs_true","",250,-1.05,1.05,250,-1.05,1.05);
  TH2D *hPzfp_recon_vs_true = new TH2D("hPzfp_recon_vs_true","",250,-1.05,1.05,250,-1.05,1.05);

  TH1D *hSxx = new TH1D("hSxx","",250,-1.1,1.1);
  TH1D *hSxy = new TH1D("hSxy","",250,-1.1,1.1);
  TH1D *hSxz = new TH1D("hSxz","",250,-1.1,1.1);

  TH1D *hSyx = new TH1D("hSyx","",250,-1.1,1.1);
  TH1D *hSyy = new TH1D("hSyy","",250,-1.1,1.1);
  TH1D *hSyz = new TH1D("hSyz","",250,-1.1,1.1);

  TH1D *hSzx = new TH1D("hSzx","",250,-1.1,1.1);
  TH1D *hSzy = new TH1D("hSzy","",250,-1.1,1.1);
  TH1D *hSzz = new TH1D("hSzz","",250,-1.1,1.1);

  TH1D *hDeterminant = new TH1D("hDeterminant","",250,0.95,1.05);

  // ---------------- reconstruction loop (same as fit script second pass) ----------------
  // NOTE: like the fit script, this loop does not apply the global cut.
  long nevent=0;
  while( T->GetEntry(nevent++) ){
    if( nevent%100 == 0 ) cout << nevent << endl;
    xfp = T->xfp;
    yfp = T->yfp;
    xpfp = T->xpfp;
    ypfp = T->ypfp;
    xfprecon = T->xfprecon;
    yfprecon = T->yfprecon;
    xpfprecon = T->xpfprecon;
    ypfprecon = T->ypfprecon;

    xptar = T->xptar;
    yptar = T->yptar;
    ytar = T->ytar;
    xtar = T->xtar;
    p = T->p;

    vz = T->vz;

    Pxtg = T->Pxtg;
    Pytg = T->Pytg;
    Pztg = T->Pztg;
    Pxfp = T->Pxfp;
    Pyfp = T->Pyfp;
    Pzfp = T->Pzfp;

    TVector3 zaxis_FT(-sin(SBS_thetabend),0,cos(SBS_thetabend));
    TVector3 yaxis_FT(0,1,0);
    TVector3 xaxis_FT = yaxis_FT.Cross(zaxis_FT).Unit();

    TVector3 nhat_fp(xpfp,ypfp,1.0);
    nhat_fp = nhat_fp.Unit();
    TVector3 nhat_fp_global = nhat_fp.X() * xaxis_FT +
      nhat_fp.Y() * yaxis_FT +
      nhat_fp.Z() * zaxis_FT;

    TVector3 nhat_tgt(xptar,yptar,1.0);

    nhat_tgt = nhat_tgt.Unit();

    thetabend = acos(nhat_fp_global.Dot(nhat_tgt) );

    TVector3 bend_axis_track = nhat_tgt.Cross(nhat_fp_global).Unit();
    TRotation Rtrack;

    Rtrack.Rotate( thetabend, bend_axis_track );

    double pinv_sum = 0.0;
    double xptar_sum = 0.0;
    double yptar_sum = 0.0;
    double ytar_sum = 0.0;
    double vz_sum = 0.0;

    double xfp_sum = 0.0;
    double yfp_sum = 0.0;
    double xpfp_sum = 0.0;
    double ypfp_sum = 0.0;

    int icoeff=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      //use smeared fp quantities to account for effects of tracking resolution
	      double term = pow(xfprecon,m)*pow(yfprecon,l)*pow(xpfprecon,k)*pow(ypfprecon,j)*pow(T->xtar,i);
	      xptar_sum += term * b_xptar[icoeff];
	      yptar_sum += term * b_yptar[icoeff];
	      ytar_sum += term * b_ytar[icoeff];
	      pinv_sum += term * b_pinv[icoeff];

	      vz_sum += term * b_vz[icoeff];

	      double fterm = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->ytar,k)*pow(1.0/T->p,j)*pow(T->xtar,i);
	      xfp_sum += fterm * b_xfp[icoeff];
	      yfp_sum += fterm * b_yfp[icoeff];
	      xpfp_sum += fterm * b_xpfp[icoeff];
	      ypfp_sum += fterm * b_ypfp[icoeff];

	      icoeff++;
	    }
	  }
	}
      }
    }

    vzrecon_direct = vz_sum;

    xfpforward = xfp_sum;
    yfpforward = yfp_sum;
    xpfpforward = xpfp_sum;
    ypfpforward = ypfp_sum;

    hxfpdiff->Fill(xfpforward-xfp);
    hyfpdiff->Fill(yfpforward-yfp);
    hxpfpdiff->Fill(xpfpforward-xpfp);
    hypfpdiff->Fill(ypfpforward-ypfp);

    hxfp_recon_vs_true->Fill(xfp,xfpforward);
    hyfp_recon_vs_true->Fill(yfp,yfpforward);
    hxpfp_recon_vs_true->Fill(xpfp,xpfpforward);
    hypfp_recon_vs_true->Fill(ypfp,ypfpforward);

    xptarrecon = xptar_sum;
    yptarrecon = yptar_sum;

    ytarrecon = ytar_sum;
    yptarrecon = yptar_sum;

    vzrecon_ytar = ytarrecon/(sin(SBStheta_central)-yptarrecon*cos(SBStheta_central));

    //now pinv_sum = p*thetabend, so we need to compute the reconstructed and true thetabend!

    hxptardiff->Fill(xptarrecon-xptar);
    hyptardiff->Fill(yptarrecon-yptar);
    hytardiff->Fill(ytarrecon-ytar);

    hxptar_recon_vs_true->Fill(xptar,xptarrecon);
    hyptar_recon_vs_true->Fill(yptar,yptarrecon);
    hytar_recon_vs_true->Fill(ytar,ytarrecon);

    TVector3 nhat_fp_recon(xpfprecon,ypfprecon,1.);
    nhat_fp_recon = nhat_fp_recon.Unit();

    TVector3 nhat_fp_recon_global = nhat_fp_recon.X() * xaxis_FT +
      nhat_fp_recon.Y() * yaxis_FT +
      nhat_fp_recon.Z() * zaxis_FT;

    TVector3 nhat_tgt_recon(xptarrecon,yptarrecon,1.);
    nhat_tgt_recon = nhat_tgt_recon.Unit();

    thetabendrecon = acos( nhat_tgt_recon.Dot(nhat_fp_recon_global) );

    precon = pinv_sum/thetabendrecon;

    hpdiff->Fill(precon/T->p-1.0);
    hp_recon_vs_true->Fill(T->p,precon);
    hpdiff_p->Fill( T->p, precon/T->p-1.0);

    hthetabend_diff->Fill(thetabendrecon-thetabend);
    hthetabend_recon_vs_true->Fill(thetabend,thetabendrecon);

    TVector3 bend_axis_recon = nhat_tgt_recon.Cross( nhat_fp_recon_global ).Unit();

    double gammarecon = sqrt(1.0+pow(precon/Mp,2));
    double chirecon = gammarecon*kappa_p*thetabendrecon;

    double thetaspinrecon = chirecon + thetabendrecon;

    TRotation Rdipole;
    Rdipole.Rotate( thetaspinrecon, bend_axis_recon );
    Rdipole.RotateY( SBS_thetabend );

    double Sxx_sum = Rdipole.XX(), Sxy_sum = Rdipole.XY(), Sxz_sum = Rdipole.XZ();
    double Syx_sum = Rdipole.YX(), Syy_sum = Rdipole.YY(), Syz_sum = Rdipole.YZ();
    double Szx_sum = Rdipole.ZX(), Szy_sum = Rdipole.ZY(), Szz_sum = Rdipole.ZZ();

    icoeff=0;
    for( int i=0; i<=order_spin; i++ ){
      for( int j=0; j<=order_spin-i; j++ ){
	for( int k=0; k<=order_spin-i-j; k++ ){
	  for( int l=0; l<=order_spin-i-j-k; l++ ){
	    for( int m=0; m<=order_spin-i-j-k-l; m++ ){
	      double term = pow(xptarrecon,m)*pow(yptarrecon,l)*pow(T->xtar,i)*pow(T->ytar,k)*pow(1.0/precon,j);
	      Sxx_sum += b_Sxx[icoeff] * term;
	      Syx_sum += b_Syx[icoeff] * term;
	      Szx_sum += b_Szx[icoeff] * term;

	      Sxy_sum += b_Sxy[icoeff] * term;
	      Syy_sum += b_Syy[icoeff] * term;
	      Szy_sum += b_Szy[icoeff] * term;

	      Sxz_sum += b_Sxz[icoeff] * term;
	      Syz_sum += b_Syz[icoeff] * term;
	      Szz_sum += b_Szz[icoeff] * term;

	      icoeff++;
	    }
	  }
	}
      }
    }

    TMatrixD Rot(3,3);
    Rot(0,0) = Sxx_sum; Rot(0,1) = Sxy_sum; Rot(0,2) = Sxz_sum;
    Rot(1,0) = Syx_sum; Rot(1,1) = Syy_sum; Rot(1,2) = Syz_sum;
    Rot(2,0) = Szx_sum; Rot(2,1) = Szy_sum; Rot(2,2) = Szz_sum;

    TVectorD Ptg(3);
    Ptg(0) = Pxtg;
    Ptg(1) = Pytg;
    Ptg(2) = Pztg;

    TVectorD Pfp(3);
    Pfp = Rot * Ptg;
    Pxfprecon = Pfp(0);
    Pyfprecon = Pfp(1);
    Pzfprecon = Pfp(2);

    TVector3 Pfpdipole(Pxtg,Pytg,Pztg);
    Pfpdipole *= Rdipole;

    Pxfpdipole = Pfpdipole.X();
    Pyfpdipole = Pfpdipole.Y();
    Pzfpdipole = Pfpdipole.Z();

    hPxfpdiff->Fill(Pxfprecon-Pxfp);
    hPyfpdiff->Fill(Pyfprecon-Pyfp);
    hPzfpdiff->Fill(Pzfprecon-Pzfp);

    hPxfp_recon_vs_true->Fill(Pxfp,Pxfprecon);
    hPyfp_recon_vs_true->Fill(Pyfp,Pyfprecon);
    hPzfp_recon_vs_true->Fill(Pzfp,Pzfprecon);

    hSxx->Fill( Sxx_sum );
    hSxy->Fill( Sxy_sum );
    hSxz->Fill( Sxz_sum );

    hSyx->Fill( Syx_sum );
    hSyy->Fill( Syy_sum );
    hSyz->Fill( Syz_sum );

    hSzx->Fill( Szx_sum );
    hSzy->Fill( Szy_sum );
    hSzz->Fill( Szz_sum );

    det = Rot.Determinant();

    hDeterminant->Fill(det);

    Tout->Fill();
  }

  hxptardiff->GetXaxis()->SetTitle("#Deltax'_{tar}");
  hyptardiff->GetXaxis()->SetTitle("#Deltay'_{tar}");
  hytardiff->GetXaxis()->SetTitle("#Deltay_{tar} (m)");
  hpdiff->GetXaxis()->SetTitle("p_{recon}/p_{true}-1");

  hthetabend_diff->GetXaxis()->SetTitle("#Delta#theta_{bend} (rad)");

  hxfpdiff->GetXaxis()->SetTitle("#Deltax_{fp} (m)");
  hyfpdiff->GetXaxis()->SetTitle("#Deltay_{fp} (m)");
  hxpfpdiff->GetXaxis()->SetTitle("#Deltax'_{fp}" );
  hypfpdiff->GetXaxis()->SetTitle("#Deltay'_{fp}" );

  hxptar_recon_vs_true->GetXaxis()->SetTitle("x'_{tar} (true)");
  hxptar_recon_vs_true->GetYaxis()->SetTitle("x'_{tar} (reconstructed)");

  hyptar_recon_vs_true->GetXaxis()->SetTitle("y'_{tar} (true)");
  hyptar_recon_vs_true->GetYaxis()->SetTitle("y'_{tar} (reconstructed)");

  hytar_recon_vs_true->GetXaxis()->SetTitle("y_{tar}^{true} (m)");
  hytar_recon_vs_true->GetYaxis()->SetTitle("y_{tar}^{recon} (m)");

  hp_recon_vs_true->GetXaxis()->SetTitle("p_{true} (GeV)");
  hp_recon_vs_true->GetYaxis()->SetTitle("p_{recon} (GeV)");

  hthetabend_recon_vs_true->GetXaxis()->SetTitle("#theta_{bend}^{true} (rad)");
  hthetabend_recon_vs_true->GetYaxis()->SetTitle("#theta_{bend}^{recon} (rad)");

  hpdiff_p->GetXaxis()->SetTitle("p_{true} (GeV)");
  hpdiff_p->GetYaxis()->SetTitle("p_{recon}/p_{true}-1 (%)");

  hxfp_recon_vs_true->GetXaxis()->SetTitle("x_{fp}^{true} (m)");
  hxfp_recon_vs_true->GetYaxis()->SetTitle("x_{fp}^{recon} (m)");

  hyfp_recon_vs_true->GetXaxis()->SetTitle("y_{fp}^{true} (m)");
  hyfp_recon_vs_true->GetYaxis()->SetTitle("y_{fp}^{recon} (m)");

  hxpfp_recon_vs_true->GetXaxis()->SetTitle("(x'_{fp})_{true} (m)");
  hxpfp_recon_vs_true->GetYaxis()->SetTitle("(x'_{fp})_{recon} (m)");

  hypfp_recon_vs_true->GetXaxis()->SetTitle("(y'_{fp})_{true} (m)");
  hypfp_recon_vs_true->GetYaxis()->SetTitle("(y'_{fp})_{recon} (m)");

  hPxfpdiff->GetXaxis()->SetTitle("#DeltaP_{x}");
  hPyfpdiff->GetXaxis()->SetTitle("#DeltaP_{y}");
  hPzfpdiff->GetXaxis()->SetTitle("#DeltaP_{z}");

  hPxfp_recon_vs_true->GetXaxis()->SetTitle("P_{x}^{true}");
  hPyfp_recon_vs_true->GetXaxis()->SetTitle("P_{y}^{true}");
  hPzfp_recon_vs_true->GetXaxis()->SetTitle("P_{z}^{true}");

  hPxfp_recon_vs_true->GetYaxis()->SetTitle("P_{x}^{recon}");
  hPyfp_recon_vs_true->GetYaxis()->SetTitle("P_{y}^{recon}");
  hPzfp_recon_vs_true->GetYaxis()->SetTitle("P_{z}^{recon}");

  hSxx->GetXaxis()->SetTitle("S_{xx}");
  hSxy->GetXaxis()->SetTitle("S_{xy}");
  hSxz->GetXaxis()->SetTitle("S_{xz}");

  hSyx->GetXaxis()->SetTitle("S_{yx}");
  hSyy->GetXaxis()->SetTitle("S_{yy}");
  hSyz->GetXaxis()->SetTitle("S_{yz}");

  hSzx->GetXaxis()->SetTitle("S_{zx}");
  hSzy->GetXaxis()->SetTitle("S_{zy}");
  hSzz->GetXaxis()->SetTitle("S_{zz}");

  hDeterminant->GetXaxis()->SetTitle("Determinant");

  fout->Write();
}
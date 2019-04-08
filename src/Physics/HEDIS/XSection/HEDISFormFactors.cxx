#include "Physics/HEDIS/XSection/HEDISFormFactors.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"

#include <TSystem.h>

#include <iostream>
#include <fstream>

#ifdef __GENIE_APFEL_ENABLED__
#include "APFEL/APFEL.h"
#endif

#ifdef __GENIE_LHAPDF6_ENABLED__
#include "LHAPDF/LHAPDF.h"
LHAPDF::PDF* pdf;
#endif


using namespace genie;
using namespace genie::constants;

//constants for NC
double Sin2thw   = 1 - kMw2 / kMz2;
double gvu    =  1./2. - 4./3.*Sin2thw;
double gvd    = -1./2. + 2./3.*Sin2thw;
double gau    =  1./2.;
double gad    = -1./2.;

struct NCCouplings {
  double c2;
  double c3;
};

NCCouplings GetNCCouplings(double gv, double ga) {
  NCCouplings cp;
  cp.c2 = TMath::Power(gv, 2.) + TMath::Power(ga, 2.);
  cp.c3 = 2*gv*ga;
  return cp;
}


#ifdef __GENIE_APFEL_ENABLED__
double ReadPDFs(int* ii, double* xx, double* qq2, bool* fst) {

  int ipdf;
  if (*ii==0) ipdf = 21;
  else        ipdf = *ii;
  double x = *xx;
  double q2 = *qq2;

#ifdef __GENIE_LHAPDF5_ENABLED__
  return fmax( LHAPDF::xfx(x, TMath::Sqrt(q2), ipdf) , 0. );
#endif
#ifdef __GENIE_LHAPDF6_ENABLED__
  return fmax( pdf->xfxQ2( ipdf, x, q2 ) , 0. );
#endif

}
#endif

//_________________________________________________________________________
HEDISFormFactors * HEDISFormFactors::fgInstance = 0;
//_________________________________________________________________________
HEDISFormFactors::HEDISFormFactors(string LHAPDFmember, bool NLO, string Scheme, int QrkThrs, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax, const double CKM[9], double MZ, double MW, double Sin2ThW)
{

  //creating directory where form factors will be stored and check if it already exists
  string snlo = (NLO) ? "NLO" : "LO";
  fFormFactorsDir = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/formfactors/" + LHAPDFmember + "_" + snlo + "_" + Scheme + "_nx" + to_string(NX) + "_nq2" + to_string(NQ2) + "/";
  if ( gSystem->mkdir(fFormFactorsDir.c_str())==0 ) LOG("HEDISFormFactors", pINFO) << "Creating Form Factors directory: " << fFormFactorsDir;

  fQrkThrs = QrkThrs;

  //filling CKM matrix using values from configuration file
  fVud2 = TMath::Power(CKM[0],2); fVus2 = TMath::Power(CKM[1],2); fVub2 = TMath::Power(CKM[2],2);
  fVcd2 = TMath::Power(CKM[3],2); fVcs2 = TMath::Power(CKM[4],2); fVcb2 = TMath::Power(CKM[5],2);
  fVtd2 = TMath::Power(CKM[6],2); fVts2 = TMath::Power(CKM[7],2); fVtb2 = TMath::Power(CKM[8],2);

#ifdef __GENIE_LHAPDF6_ENABLED__
  LOG("HEDISFormFactors", pINFO) << "Initialising LHAPDF6...";
  pdf = LHAPDF::mkPDF(LHAPDFmember, 0);
  xPDFmin  = pdf->xMin();
  Q2PDFmin = pdf->q2Min();
  Q2PDFmax = pdf->q2Max();
  LOG("HEDISFormFactors", pINFO) << "PDF info:";
  LOG("HEDISFormFactors", pINFO) << "OrderQCD = " << pdf->orderQCD();
  LOG("HEDISFormFactors", pINFO) << "FlavorScheme = " << pdf->info().get_entry("FlavorScheme");
  LOG("HEDISFormFactors", pINFO) << "NumFlavors = " << pdf->info().get_entry("NumFlavors");
  LOG("HEDISFormFactors", pINFO) << "Xmin = " << xPDFmin << "  Xmax = " << pdf->xMax() << "  Q2min = " << Q2PDFmin << "  Q2max = " << Q2PDFmax;
  LOG("HEDISFormFactors", pINFO) << "MZ = " << pdf->info().get_entry("MZ");
  for (int i=1; i<7; i++) {
    mPDFQrk[i] = pdf->quarkMass(i);
    LOG("HEDISFormFactors", pINFO) << "M" << i << " = " << mPDFQrk[i];
  }
#endif
#ifdef __GENIE_LHAPDF5_ENABLED__
  LOG("HEDISFormFactors", pINFO) << "Initialising LHAPDF5...";
  LHAPDF::initPDFByName(LHAPDFmember, LHAPDF::LHGRID, 0);
  xPDFmin  = LHAPDF::getXmin(0);
  Q2PDFmin = LHAPDF::getQ2min(0);
  Q2PDFmax = LHAPDF::getQ2max(0);
  LOG("HEDISFormFactors", pINFO) << "PDF info:";
  LOG("HEDISFormFactors", pINFO) << "Xmin = " << xPDFmin << "  Xmax = " << pdf->xMax() << "  Q2min = " << Q2PDFmin << "  Q2max = " << Q2PDFmax;
  for (int i=1; i<7; i++) {
    mPDFQrk[i] = LHAPDF::getQMass(i);
    LOG("HEDISFormFactors", pINFO) << "M" << i << " = " << mPDFQrk[i];
  }
#endif

  if ( xGRIDmin < xPDFmin ) {
    LOG("HEDISFormFactors", pWARN) << "Lower boundary in X is smaller than input PDF";
    LOG("HEDISFormFactors", pWARN) << "xPDFmin = "  << xPDFmin;
    LOG("HEDISFormFactors", pWARN) << "xGRIDmin = " << xGRIDmin;
  }
  else if ( Q2GRIDmin < Q2PDFmin ) {
    LOG("HEDISFormFactors", pWARN) << "Lower boundary in Q2 is smaller than input PDF";
    LOG("HEDISFormFactors", pWARN) << "Q2PDFmin = "  << Q2PDFmin;
    LOG("HEDISFormFactors", pWARN) << "Q2GRIDmin = " << Q2GRIDmin;
  }
  else if ( Q2GRIDmax > Q2PDFmax ) {
    LOG("HEDISFormFactors", pWARN) << "Upper boundary in Q2 is bigger than input PDF";
    LOG("HEDISFormFactors", pWARN) << "Q2PDFmax = "  << Q2PDFmax;
    LOG("HEDISFormFactors", pWARN) << "Q2GRIDmax = " << Q2GRIDmax;
  }

  // define arrays to fill from data files
  double dlogq2 = TMath::Abs( TMath::Log10(Q2GRIDmin)-TMath::Log10(Q2GRIDmax) ) / NQ2;
  double dlogx  = TMath::Abs( TMath::Log10(xPDFmin)-TMath::Log10(1.) ) / NX;

  LOG("HEDISFormFactors", pINFO) << "Grid x,Q2 :" << NX << " , " << NQ2;

  for ( double logq2 = TMath::Log10(Q2GRIDmin); logq2<TMath::Log10(Q2GRIDmax); logq2+= dlogq2 ) {
    if ( TMath::Power( 10, logq2 + 0.5*dlogq2 )>Q2PDFmax ) continue;
    ff_logq2_array.push_back(TMath::Power( 10, logq2 + 0.5*dlogq2 ));
    LOG("HEDISFormFactors", pDEBUG) << "q2: " << ff_logq2_array.back();
  }
  for ( double logx = TMath::Log10(xGRIDmin); logx<TMath::Log10(1.); logx+= dlogx ) {
    if ( TMath::Power( 10, logx + 0.5*dlogx )>1. ) continue;
    ff_logx_array.push_back(TMath::Power( 10, logx + 0.5*dlogx ));
    LOG("HEDISFormFactors", pDEBUG) << "x: " << ff_logx_array.back();
  }

  //load form factors for each quark at LO
  for( int ch=1; ch<kHEDIS_numofchannels; ch++ ) LoadFormfactors( (HEDISChannel_t)ch );


#ifdef __GENIE_APFEL_ENABLED__
  //load form factors for each nucleon at LO
  if (NLO) {
    LOG("HEDISFormFactors", pINFO) << "Initialising APFEL..." ; 
    APFEL::SetPDFSet(LHAPDFmember);
    APFEL::SetReplica(0);
    if (Scheme=="FONLL") {
      APFEL::SetMassScheme("FONLL-B");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[6]);
    } 
    else if (Scheme=="ZM-VFNS") {
      APFEL::SetMassScheme("ZM-VFNS");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[5]+0.1);
    }
    APFEL::SetQLimits(TMath::Sqrt(Q2GRIDmin),TMath::Sqrt(Q2GRIDmax));
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1,90,3,xGRIDmin);
    APFEL::SetGridParameters(2,50,5,1e-1);
    APFEL::SetGridParameters(3,40,5,8e-1);
    APFEL::SetPerturbativeOrder(1);
    APFEL::SetAlphaQCDRef(pdf->alphasQ(MZ),MZ);    
    APFEL::SetProtonMass(kProtonMass);
    APFEL::SetWMass(MW);
    APFEL::SetZMass(MZ);
    APFEL::SetSin2ThetaW(Sin2thw);
    APFEL::SetCKM(CKM[0], CKM[1], CKM[2],
                  CKM[3], CKM[4], CKM[5],
                  CKM[6], CKM[7], CKM[8]);

    for( int ch=1; ch<kHEDISNucl_numofchannels; ch++ ) LoadNuclFormfactors( (HEDISNuclChannel_t)ch );
  }
#endif

  fgInstance = 0;

  delete pdf;

}
//_________________________________________________________________________
HEDISFormFactors::~HEDISFormFactors()
{

}
//_________________________________________________________________________
HEDISFormFactors * HEDISFormFactors::Instance(string LHAPDFmember, bool NLO, string Scheme, int QrkThrs, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax, const double CKM[9], double MZ, double MW, double Sin2ThW)
{
  if(fgInstance == 0) {
    LOG("HEDISFormFactors", pINFO) << "Late initialization";
    static HEDISFormFactors::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new HEDISFormFactors(LHAPDFmember, NLO, Scheme, QrkThrs, NX, xGRIDmin, NQ2, Q2GRIDmin, Q2GRIDmax, CKM, MZ, MW, Sin2ThW);
  }  
  return fgInstance;
}
//_________________________________________________________________________
void HEDISFormFactors::LoadFormfactors( HEDISChannel_t ch )
{

  string formfactorFile = fFormFactorsDir + "FFLO_" + HEDISChannel::AsString(ch) + ".dat";

  // make sure data files are available
  LOG("HEDISFormFactors", pINFO) << "Checking if file " << formfactorFile << " exists...";        
  if ( gSystem->AccessPathName( formfactorFile.c_str()) ) CreateFormFactorFile( ch, formfactorFile );


  for(int ff = 1; ff <= kFFnumber; ++ff) {
    fFormFactorsTables[ch].Table[(HEDISFormFactorType_t)ff] = ReadFormFactorFile( formfactorFile, (HEDISFormFactorType_t)ff );
  }

}
#ifdef __GENIE_APFEL_ENABLED__
//_________________________________________________________________________
void HEDISFormFactors::LoadNuclFormfactors( HEDISNuclChannel_t ch )
{

  string NLOformfactorFile = fFormFactorsDir + "NuclFFNLO_" + HEDISChannel::AsString(ch) + ".dat";
  LOG("HEDISFormFactors", pINFO) << "Checking if file " << NLOformfactorFile << " exists...";         
  if ( gSystem->AccessPathName( NLOformfactorFile.c_str()) ) CreateNLONuclFormFactorFile( ch, NLOformfactorFile );


  string LOformfactorFile = fFormFactorsDir + "NuclFFLO_" + HEDISChannel::AsString(ch) + ".dat";
  LOG("HEDISFormFactors", pINFO) << "Checking if file " << LOformfactorFile << " exists...";          
  if ( gSystem->AccessPathName( LOformfactorFile.c_str()) ) CreateLONuclFormFactorFile( ch, LOformfactorFile );


  for(int ff = 1; ff <= kFFnumber; ++ff) {
    fNLONuclFormFactorsTables[ch].Table[(HEDISFormFactorType_t)ff] = ReadFormFactorFile( NLOformfactorFile, (HEDISFormFactorType_t)ff );
    fLONuclFormFactorsTables[ch].Table[(HEDISFormFactorType_t)ff] = ReadFormFactorFile( LOformfactorFile, (HEDISFormFactorType_t)ff );
  }

}
#endif
//_________________________________________________________________________
BLI2DNonUnifGrid * HEDISFormFactors::ReadFormFactorFile( string filename, HEDISFormFactorType_t ffType )
{

  // open file
  std::ifstream ff_stream(filename.c_str(), std::ios::in);

  // check file exists
  if( !ff_stream.good() ){
    LOG("HEDISFormFactors", pERROR) << "Bad file name: " << filename;
    return NULL;
  }

  int nx = ff_logq2_array.size();
  int ny = ff_logx_array.size();

  double x[nx];
  double y[ny];
  double z[nx*ny];

  for (int i=0; i<nx; i++) x[i] = ff_logq2_array[i];
  for (int j=0; j<ny; j++) y[j] = ff_logx_array[j];

  for(int ff = 1; ff < kFFnumber; ++ff) {
    for (int ij=0; ij<nx*ny; ij++) {
      double aux;
      if ( (HEDISFormFactorType_t)ff == ffType ) ff_stream >> z[ij];
      else                                       ff_stream >> aux;
    }
  }
  
  return new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );

}
//____________________________________________________________________________
void HEDISFormFactors::CreateFormFactorFile( HEDISChannel_t ch, string filename ) {

  int pdg_nucl  = HEDISChannel::HitNuclPdg(ch);
  bool sea_iq   = HEDISChannel::HitQuarkSea(ch);
  int pdg_iq    = HEDISChannel::HitQuarkPdg(ch);
  int pdg_fq    = HEDISChannel::FnlQuarkPdg(ch);

  int qrkd = 0;
  int qrku = 0;
  if      ( pdg_nucl==2212 ) { qrkd = 1 ; qrku = 2; }
  else if ( pdg_nucl==2112 ) { qrkd = 2 ; qrku = 1; }

  //variables associated to the PDF and coupling of the quarks
  int qpdf1 = -999;
  int qpdf2 = -999;
  double Cp2 = -999;
  double Cp3 = -999;
  double sign3 = (HEDISChannel::IsNu(ch)) ? +1. : -1.;
  if ( HEDISChannel::InteractionType(ch) == kIntWeakCC ) {
    if ( HEDISChannel::IsNu(ch) ) {
      if      ( pdg_iq== 1 && !sea_iq && pdg_fq== 2 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*fVud2; Cp3 =  2*fVud2; }
      else if ( pdg_iq== 1 && !sea_iq && pdg_fq== 4 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*fVcd2; Cp3 =  2*fVcd2; }
      else if ( pdg_iq== 1 && !sea_iq && pdg_fq== 6 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*fVtd2; Cp3 =  2*fVtd2; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 2 ) { qpdf1 = -qrkd;                Cp2 = 2*fVud2; Cp3 =  2*fVud2; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 4 ) { qpdf1 = -qrkd;                Cp2 = 2*fVcd2; Cp3 =  2*fVcd2; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 6 ) { qpdf1 = -qrkd;                Cp2 = 2*fVtd2; Cp3 =  2*fVtd2; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 2 ) { qpdf1 =  3;                   Cp2 = 2*fVus2; Cp3 =  2*fVus2; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  3;                   Cp2 = 2*fVcs2; Cp3 =  2*fVcs2; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 6 ) { qpdf1 =  3;                   Cp2 = 2*fVts2; Cp3 =  2*fVts2; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 2 ) { qpdf1 =  5;                   Cp2 = 2*fVub2; Cp3 =  2*fVub2; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  5;                   Cp2 = 2*fVcb2; Cp3 =  2*fVcb2; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 6 ) { qpdf1 =  5;                   Cp2 = 2*fVtb2; Cp3 =  2*fVtb2; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -qrku;                Cp2 = 2*fVud2; Cp3 = -2*fVud2; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -qrku;                Cp2 = 2*fVus2; Cp3 = -2*fVus2; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -qrku;                Cp2 = 2*fVub2; Cp3 = -2*fVub2; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -4;                   Cp2 = 2*fVcd2; Cp3 = -2*fVcd2; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -4;                   Cp2 = 2*fVcs2; Cp3 = -2*fVcs2; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -4;                   Cp2 = 2*fVcb2; Cp3 = -2*fVcb2; }
    }
    else {
      if      ( pdg_iq== 2 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fVud2; Cp3 =  2*fVud2; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 3 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fVus2; Cp3 =  2*fVus2; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 5 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fVub2; Cp3 =  2*fVub2; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrku;                Cp2 = 2*fVud2; Cp3 =  2*fVud2; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 3 ) { qpdf1 = -qrku;                Cp2 = 2*fVus2; Cp3 =  2*fVus2; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 5 ) { qpdf1 = -qrku;                Cp2 = 2*fVub2; Cp3 =  2*fVub2; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 1 ) { qpdf1 =  4;                   Cp2 = 2*fVcd2; Cp3 =  2*fVcd2; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  4;                   Cp2 = 2*fVcs2; Cp3 =  2*fVcs2; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  4;                   Cp2 = 2*fVcb2; Cp3 =  2*fVcb2; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrkd;                Cp2 = 2*fVud2; Cp3 = -2*fVud2; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -qrkd;                Cp2 = 2*fVcd2; Cp3 = -2*fVcd2; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -qrkd;                Cp2 = 2*fVtd2; Cp3 = -2*fVtd2; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -3;                   Cp2 = 2*fVus2; Cp3 = -2*fVus2; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -3;                   Cp2 = 2*fVcs2; Cp3 = -2*fVcs2; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -3;                   Cp2 = 2*fVts2; Cp3 = -2*fVts2; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -5;                   Cp2 = 2*fVub2; Cp3 = -2*fVub2; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -5;                   Cp2 = 2*fVcb2; Cp3 = -2*fVcb2; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -5;                   Cp2 = 2*fVtb2; Cp3 = -2*fVtb2; }
    }
  }
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) {           
    NCCouplings Au, Ad;
    Ad = GetNCCouplings(gvd,gad);
    Au = GetNCCouplings(gvu,gau);
    if      ( pdg_iq== 1 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = Ad.c2; Cp3 =  Ad.c3; }
    else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 2 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = Au.c2; Cp3 =  Au.c3; }
    else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrkd;                Cp2 = Ad.c2; Cp3 =  Ad.c3; }
    else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 2 ) { qpdf1 = -qrku;                Cp2 = Au.c2; Cp3 =  Au.c3; }
    else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  3;                   Cp2 = Ad.c2; Cp3 =  Ad.c3; }
    else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  4;                   Cp2 = Au.c2; Cp3 =  Au.c3; }
    else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  5;                   Cp2 = Ad.c2; Cp3 =  Ad.c3; }
    else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -qrkd;                Cp2 = Ad.c2; Cp3 = -Ad.c3; }
    else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrku;                Cp2 = Au.c2; Cp3 = -Au.c3; }
    else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -3;                   Cp2 = Ad.c2; Cp3 = -Ad.c3; }
    else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -4;                   Cp2 = Au.c2; Cp3 = -Au.c3; }
    else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -5;                   Cp2 = Ad.c2; Cp3 = -Ad.c3; }
  }   

  //varaibles used for certain threshold options
  double mass_fq   = mPDFQrk[TMath::Abs(pdg_fq)];
  double mass_nucl = (pdg_nucl==2212) ? kProtonMass : kNeutronMass;

  // open file
  LOG("HEDISFormFactors", pINFO) << "Creating LO FormFactor file";
  std::ofstream ff_stream(filename.c_str());

  for(int ff = 1; ff < kFFnumber; ++ff) {
    for ( const auto& Q2 : ff_logq2_array ) {
      for ( const auto& x : ff_logx_array ) {

        double z = x; //in case you want to apply scaling

        if      (fQrkThrs==1) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { ff_stream << 0. << "  "; continue; }
        } 
        else if (fQrkThrs==2) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { ff_stream << 0. << "  "; continue; }
            z *= 1+mass_fq*mass_fq/Q2;
        }
        else if (fQrkThrs==3) {
            z *= 1+mass_fq*mass_fq/Q2;
        }

        double xPDF = TMath::Max( z, xPDFmin );
        double Q2PDF = TMath::Max( Q2, Q2PDFmin );
        Q2PDF = TMath::Min( Q2, Q2PDFmax  );

        double fPDF = fmax( pdf->xfxQ2(qpdf1, xPDF, Q2PDF)/z , 0.);
        if (qpdf2!= -999) fPDF -= fmax( pdf->xfxQ2(qpdf2, xPDF, Q2PDF)/z , 0.);
                    
        double tmp = -999;
        if      ( (HEDISFormFactorType_t)ff==kFFT1 ) tmp = fPDF*Cp2/2;
        else if ( (HEDISFormFactorType_t)ff==kFFT2 ) tmp = fPDF*Cp2*z;
        else if ( (HEDISFormFactorType_t)ff==kFFT3 ) tmp = fPDF*Cp3*sign3;

        LOG("HEDISFormFactors", pDEBUG) << "QrkLOF" << ff << "[x=" << x << "," << Q2 << "] = " << tmp;

        ff_stream << tmp << "  ";
        
      }
    }
  }

  ff_stream.close();

}
#ifdef __GENIE_APFEL_ENABLED__
//____________________________________________________________________________
void HEDISFormFactors::CreateNLONuclFormFactorFile( HEDISNuclChannel_t ch, string filename )
{

  if ( HEDISChannel::IsNu(ch) ) APFEL::SetProjectileDIS("neutrino");
  else                          APFEL::SetProjectileDIS("antineutrino");
  if      ( HEDISChannel::InteractionType(ch) == kIntWeakCC ) APFEL::SetProcessDIS("CC");
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) APFEL::SetProcessDIS("NC");
  if      ( HEDISChannel::HitNuclPdg(ch)==2212 ) APFEL::SetTargetDIS("proton");
  else if ( HEDISChannel::HitNuclPdg(ch)==2112 ) APFEL::SetTargetDIS("neutron");

  APFEL::InitializeAPFEL_DIS();

  int nx  = ff_logx_array.size();
  int nq2 = ff_logq2_array.size();
  double xlist[nx*nq2];
  double q2list[nx*nq2];
  double F2list[nx*nq2];
  double FLlist[nx*nq2];
  double xF3list[nx*nq2];

  int nlist = 0;
  for ( const auto& Q2 : ff_logq2_array ) {
    double Q = TMath::Sqrt(Q2);
    double norm = (HEDISChannel::InteractionType(ch)==kIntWeakCC) ? 1. : 2./TMath::Power( Q2/(Q2 + TMath::Power(APFEL::GetZMass(),2))/4/APFEL::GetSin2ThetaW()/(1-APFEL::GetSin2ThetaW()), 2 );
    APFEL::SetAlphaQCDRef(pdf->alphasQ(Q),Q);
    APFEL::ComputeStructureFunctionsAPFEL(Q,Q);
    for ( const auto& x : ff_logx_array ) {
      q2list[nlist]  = Q2;
      xlist[nlist]   = x;
      FLlist[nlist]  = norm*APFEL::FLtotal(x);
      F2list[nlist]  = norm*APFEL::F2total(x);
      xF3list[nlist] = norm*APFEL::F3total(x);
      nlist++;
    }
  }

  // open file
  LOG("HEDISFormFactors", pINFO) << "Creating NLO NuclFormFactor file";
  std::ofstream ff_stream(filename.c_str());

  double sign3 = (HEDISChannel::IsNu(ch)) ? +1. : -1.;
  for(int ff = 1; ff < kFFnumber; ++ff) {
    for (int i=0; i<nx*nq2; i++) {
      double tmp = 0;
      if      ( (HEDISFormFactorType_t)ff==kFFT1 ) tmp = (F2list[i]-FLlist[i])/2/xlist[i];
      else if ( (HEDISFormFactorType_t)ff==kFFT2 ) tmp = F2list[i];
      else if ( (HEDISFormFactorType_t)ff==kFFT3 ) tmp = sign3 * xF3list[i] / xlist[i];
      LOG("HEDISFormFactors", pDEBUG) << "NuclNLOF" << ff << "[x=" << xlist[i] << "," << q2list[i] << "] = " << tmp;
      ff_stream << tmp << "  ";
    }
  }
    
  ff_stream.close();

}
//____________________________________________________________________________
void HEDISFormFactors::CreateLONuclFormFactorFile( HEDISNuclChannel_t ch, string filename ) {

  int frstqrkch = HEDISChannel::GetFirstHEDISChannel(ch);
  int lastqrkch = HEDISChannel::GetLastHEDISChannel(ch);

  // open file
  LOG("HEDISFormFactors", pINFO) << "Creating LO NuclFormFactor file";
  std::ofstream ff_stream(filename.c_str());

  for(int ff = 1; ff < kFFnumber; ++ff) {
    for ( const auto& Q2 : ff_logq2_array ) {
      for ( const auto& x : ff_logx_array ) {
        double tmp = 0;
        for( int qch=frstqrkch; qch<=lastqrkch; qch++ ) tmp += fFormFactorsTables[(HEDISChannel_t)qch].Table[(HEDISFormFactorType_t)ff]->Evaluate(x,Q2);
        LOG("HEDISFormFactors", pDEBUG) << "NuclLOF" << ff << "[x=" << x << "," << Q2 << "] = " << tmp;
        ff_stream << tmp << "  ";
      }
    }
  }

  ff_stream.close();

}
#endif
//____________________________________________________________________________
FF_xQ2 HEDISFormFactors::EvalFFQrkLO( HEDISChannel_t ch, double x, double Q2 ) 
{
  FF_xQ2 ff;
  ff.F1 = fFormFactorsTables[ch].Table[kFFT1]->Evaluate(Q2,x);
  ff.F2 = fFormFactorsTables[ch].Table[kFFT2]->Evaluate(Q2,x);
  ff.F3 = fFormFactorsTables[ch].Table[kFFT3]->Evaluate(Q2,x);
  return ff;
}
#ifdef __GENIE_APFEL_ENABLED__
//____________________________________________________________________________
FF_xQ2 HEDISFormFactors::EvalNuclFFLO( HEDISNuclChannel_t ch, double x, double Q2 ) 
{
  FF_xQ2 ff;
  ff.F1 = fLONuclFormFactorsTables[ch].Table[kFFT1]->Evaluate(Q2,x);
  ff.F2 = fLONuclFormFactorsTables[ch].Table[kFFT2]->Evaluate(Q2,x);
  ff.F3 = fLONuclFormFactorsTables[ch].Table[kFFT3]->Evaluate(Q2,x);
  return ff;
}
//____________________________________________________________________________
FF_xQ2 HEDISFormFactors::EvalNuclFFNLO( HEDISNuclChannel_t ch, double x, double Q2 ) 
{
  FF_xQ2 ff;
  ff.F1 = fNLONuclFormFactorsTables[ch].Table[kFFT1]->Evaluate(Q2,x);
  ff.F2 = fNLONuclFormFactorsTables[ch].Table[kFFT2]->Evaluate(Q2,x);
  ff.F3 = fNLONuclFormFactorsTables[ch].Table[kFFT3]->Evaluate(Q2,x);
  return ff;
}
#endif

#include "Physics/HEDIS/XSection/HEDISFormFactors.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Conventions/Constants.h"

#include <TSystem.h>

#include <iostream>
#include <fstream>

#include "LHAPDF/LHAPDF.h"
#include "QCDNUM/QCDNUM.h"

using namespace genie;
using namespace genie::constants;

LHAPDF::PDF* pdf;

//constants for NC
double Sin2thw   = 1 - kMw2 / kMz2;
double gvu       =  1./2. - 4./3.*Sin2thw;
double gvd       = -1./2. + 2./3.*Sin2thw;
double gau_nu    =  1./2.;
double gau_nubar = -1./2.;
double gad_nu    = -1./2.;
double gad_nubar =  1./2.;

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


double ReadPDFs(int* ii, double* xx, double* qq2, bool* fst) {
  int ipdf;
  if (*ii==0) ipdf = 21;
  else        ipdf = *ii;
  double x = *xx;
  double q2 = *qq2;
  return fmax( pdf->xfxQ2( ipdf, x, q2 ) , 0. );
}

//_________________________________________________________________________
HEDISFormFactors * HEDISFormFactors::fgInstance = 0;
//_________________________________________________________________________
HEDISFormFactors::HEDISFormFactors(bool NLO, string LHAPDFmember, int NX, int NQ2, const double CKM[9])
{

  fVud2 = CKM[0]; fVus2 = CKM[1]; fVub2 = CKM[2];
  fVcd2 = CKM[3]; fVcs2 = CKM[4]; fVcb2 = CKM[5];
  fVtd2 = CKM[6]; fVts2 = CKM[7]; fVtb2 = CKM[8];

  LOG("HEDISFormFactors", pINFO) << "Initialising LHAPDF5...";

  pdf = LHAPDF::mkPDF(LHAPDFmember, 0);;

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

  // define arrays to fill from data files
  double dlogq2 = TMath::Abs( TMath::Log10(Q2PDFmin)-TMath::Log10(Q2PDFmax) ) / NQ2;
  double dlogx  = TMath::Abs( TMath::Log10(xPDFmin)-TMath::Log10(1.) ) / NX;

  LOG("HEDISFormFactors", pINFO) << "Grid x,Q2 :" << NX << " , " << NQ2;

  for ( double logq2 = TMath::Log10(Q2PDFmin); logq2<TMath::Log10(Q2PDFmax); logq2+= dlogq2 ) {
    if ( TMath::Power( 10, logq2 + 0.5*dlogq2 )>Q2PDFmax ) continue;
    ff_logq2_array.push_back(TMath::Power( 10, logq2 + 0.5*dlogq2 ));
    LOG("HEDISFormFactors", pDEBUG) << "q2: " << ff_logq2_array.back();
  }
  for ( double logx = TMath::Log10(xPDFmin); logx<TMath::Log10(1.); logx+= dlogx ) {
    if ( TMath::Power( 10, logx + 0.5*dlogx )>1. ) continue;
    ff_logx_array.push_back(TMath::Power( 10, logx + 0.5*dlogx ));
    LOG("HEDISFormFactors", pDEBUG) << "x: " << ff_logx_array.back();
  }


  for( int ch=1; ch<kHEDIS_numofchannels; ch++ ) LoadFormfactors( (HEDISChannel_t)ch );


  if (NLO) {
    for( int ch=1; ch<kHEDISNucl_numofchannels; ch++ ) LoadNuclFormfactors( (HEDISNuclChannel_t)ch );
  }

  fgInstance = 0;

  delete pdf;

}
//_________________________________________________________________________
HEDISFormFactors::~HEDISFormFactors()
{

}
//_________________________________________________________________________
HEDISFormFactors * HEDISFormFactors::Instance(bool NLO, string LHAPDFmenber, int NX, int NQ2, const double CKM[9])
{
  if(fgInstance == 0) {
    LOG("HEDISFormFactors", pINFO) << "Late initialization";
    static HEDISFormFactors::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new HEDISFormFactors(NLO,LHAPDFmenber,NX,NQ2,CKM);
  }  
  return fgInstance;
}
//_________________________________________________________________________
void HEDISFormFactors::InitQCDNUM(void) {

  LOG("HEDISFormFactors", pINFO) << "Initialising qcdnum v17..." ; 

  double massZ = atof(pdf->info().get_entry("MZ").c_str()); 

  QCDNUM::qcinit(6,"");
  QCDNUM::setord(2);
  QCDNUM::setalf(pdf->alphasQ(massZ),TMath::Power(massZ,2));
  QCDNUM::setval("elim",-999); //avoid error spline oscilation

  double xmin[] = {xPDFmin};
  int    iwt[]  = {1};
  int nxout;
  QCDNUM::gxmake(xmin,iwt,1,1000,nxout,3);  //2pol because 3pol returns FF<0 for x close to 1
  
  double qarr[] = {Q2PDFmin, Q2PDFmax};
  double warr[] = {1., 1.};
  int nqout;
  QCDNUM::gqmake(qarr,warr,2,1000,nqout);
  
  int iqc = QCDNUM::iqfrmq(mPDFQrk[4]*mPDFQrk[4]);
  int iqb = QCDNUM::iqfrmq(mPDFQrk[5]*mPDFQrk[5]);
  int iqt = QCDNUM::iqfrmq(mPDFQrk[6]*mPDFQrk[6]);    
  
  //QCDNUM::setcbt(0,iqc,iqb,iqt); 
  //QCDNUM::setcbt(0,iqc,iqb,0); //vfns
  QCDNUM::setcbt(6,0,0,0); //ffns

  int idmin = 0 ;
  int idmax = 0;
  int nwpdf = 0;
  QCDNUM::fillwt(1,idmin,idmax,nwpdf); 
  QCDNUM::dmpwgt(1,22,"unpol.wgt");    

  int nwords = 0;
  QCDNUM::zmfillw(nwords);
  QCDNUM::zmdumpw(22,"zm_weights.wgt");    

  double eps;
  QCDNUM::extpdf(ReadPDFs,5,0,0,eps);

  QCDNUM::zswitch(5);

  QCDNUMIsAlreadyInit = true;

}

//_________________________________________________________________________
void HEDISFormFactors::LoadFormfactors( HEDISChannel_t ch )
{

  string formfactorFile = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/formfactors/FFLO_" + HEDISChannel::AsString(ch) + ".dat";

  // make sure data files are available
  LOG("HEDISFormFactors", pINFO) << "Checking if file " << formfactorFile << " exists...";        
  if ( gSystem->AccessPathName( formfactorFile.c_str()) ) CreateFormFactorFile( ch, formfactorFile );


  for(int ff = 1; ff <= kFFnumber; ++ff) {
    fFormFactorsTables[ch].Table[(HEDISFormFactorType_t)ff] = ReadFormFactorFile( formfactorFile, (HEDISFormFactorType_t)ff );
  }

}
//_________________________________________________________________________
void HEDISFormFactors::LoadNuclFormfactors( HEDISNuclChannel_t ch )
{

  string NLOformfactorFile = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/formfactors/NuclFFNLO_" + HEDISChannel::AsString(ch) + ".dat";
  LOG("HEDISFormFactors", pINFO) << "Checking if file " << NLOformfactorFile << " exists...";         
  if ( gSystem->AccessPathName( NLOformfactorFile.c_str()) ) CreateNLONuclFormFactorFile( ch, NLOformfactorFile );


  string LOformfactorFile = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/formfactors/NuclFFLO_" + HEDISChannel::AsString(ch) + ".dat";
  LOG("HEDISFormFactors", pINFO) << "Checking if file " << LOformfactorFile << " exists...";          
  if ( gSystem->AccessPathName( LOformfactorFile.c_str()) ) CreateLONuclFormFactorFile( ch, LOformfactorFile );


  for(int ff = 1; ff <= kFFnumber; ++ff) {
    fNLONuclFormFactorsTables[ch].Table[(HEDISFormFactorType_t)ff] = ReadFormFactorFile( NLOformfactorFile, (HEDISFormFactorType_t)ff );
    fLONuclFormFactorsTables[ch].Table[(HEDISFormFactorType_t)ff] = ReadFormFactorFile( LOformfactorFile, (HEDISFormFactorType_t)ff );
  }

}
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

  int nx = ff_logx_array.size();
  int ny = ff_logq2_array.size();

  double x[nx];
  double y[ny];
  double z[nx*ny];

  for (int i=0; i<nx; i++) x[i] = ff_logx_array[i];
  for (int j=0; j<ny; j++) y[j] = ff_logq2_array[j];

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
      if      ( pdg_iq== 2 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fVud2; Cp3 = -2*fVud2; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 3 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fVus2; Cp3 = -2*fVus2; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 5 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fVub2; Cp3 = -2*fVub2; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrku;                Cp2 = 2*fVud2; Cp3 = -2*fVud2; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 3 ) { qpdf1 = -qrku;                Cp2 = 2*fVus2; Cp3 = -2*fVus2; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 5 ) { qpdf1 = -qrku;                Cp2 = 2*fVub2; Cp3 = -2*fVub2; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 1 ) { qpdf1 =  4;                   Cp2 = 2*fVcd2; Cp3 = -2*fVcd2; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  4;                   Cp2 = 2*fVcs2; Cp3 = -2*fVcs2; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  4;                   Cp2 = 2*fVcb2; Cp3 = -2*fVcb2; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrkd;                Cp2 = 2*fVud2; Cp3 =  2*fVud2; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -qrkd;                Cp2 = 2*fVcd2; Cp3 =  2*fVcd2; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -qrkd;                Cp2 = 2*fVtd2; Cp3 =  2*fVtd2; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -3;                   Cp2 = 2*fVus2; Cp3 =  2*fVus2; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -3;                   Cp2 = 2*fVcs2; Cp3 =  2*fVcs2; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -3;                   Cp2 = 2*fVts2; Cp3 =  2*fVts2; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -5;                   Cp2 = 2*fVub2; Cp3 =  2*fVub2; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -5;                   Cp2 = 2*fVcb2; Cp3 =  2*fVcb2; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -5;                   Cp2 = 2*fVtb2; Cp3 =  2*fVtb2; }
    }
  }
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) {           
    double gad = 0;
    double gau = 0;
    if ( HEDISChannel::IsNu(ch) ) {
      gad = gad_nu;
      gau = gau_nu;
    }
    else {
      gad = gad_nubar;
      gau = gau_nubar;
    }

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
  int thrs_type = 1;
  double mass_fq   = mPDFQrk[TMath::Abs(pdg_fq)];
  double mass_nucl = PDGLibrary::Instance()->Find( pdg_nucl )->Mass();

  // open file
  LOG("HEDISFormFactors", pINFO) << "Creating LO FormFactor file";
  std::ofstream ff_stream(filename.c_str());

  for(int ff = 1; ff < kFFnumber; ++ff) {
    for ( const auto& x : ff_logx_array ) {
      for ( const auto& Q2 : ff_logq2_array ) {
                      
        double z = x; //in case you want to apply scaling

        if      (thrs_type==1) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { ff_stream << 0. << "  "; continue; }
        } 
        else if (thrs_type==2) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { ff_stream << 0. << "  "; continue; }
            z *= 1+mass_fq*mass_fq/Q2;
        }
        else if (thrs_type==3) {
            z *= 1+mass_fq*mass_fq/Q2;
        }

        double fPDF = fmax( pdf->xfxQ2(qpdf1, z, Q2)/z , 0.);
        if (qpdf2!= -999) fPDF -= fmax( pdf->xfxQ2(qpdf2, z, Q2)/z , 0.);
                    
        double tmp = -999;
        if      ( (HEDISFormFactorType_t)ff==kFFT1 ) tmp = fPDF*Cp2/2;
        else if ( (HEDISFormFactorType_t)ff==kFFT2 ) tmp = fPDF*Cp2*z;
        else if ( (HEDISFormFactorType_t)ff==kFFT3 ) tmp = fPDF*Cp3;

        LOG("HEDISFormFactors", pDEBUG) << "QrkLOF" << ff << "[x=" << x << "," << Q2 << "] = " << tmp;

        ff_stream << tmp << "  ";
        
      }
    }
  }

  ff_stream.close();

}
//____________________________________________________________________________
void HEDISFormFactors::CreateNLONuclFormFactorFile( HEDISNuclChannel_t ch, string filename )
{

  if (!QCDNUMIsAlreadyInit) InitQCDNUM();

    // open file
    LOG("HEDISFormFactors", pINFO) << "Creating NLO NuclFormFactor file";
    std::ofstream ff_stream(filename.c_str());

  double wf2[13],wf3[13];

  if      (ch==kHEDISNucl_v_cc_p) {
    //                     tb   bb   cb   sb   ub   db    g    d    u    s    c    b    t
    double wf2aux[13] = {  2.,  0.,  2.,  0.,  2.,  0.,  0.,  2.,  0.,  2.,  0.,  2.,  0. };
    double wf3aux[13] = { -2.,  0., -2.,  0., -2.,  0.,  0.,  2.,  0.,  2.,  0.,  2.,  0. };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  else if (ch==kHEDISNucl_v_cc_n) {
    double wf2aux[13] = {  2.,  0.,  2.,  0.,  0.,  2.,  0.,  0.,  2.,  2.,  0.,  2.,  0. };
    double wf3aux[13] = { -2.,  0., -2.,  0.,  0., -2.,  0.,  0.,  2.,  2.,  0.,  2.,  0. };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  if      (ch==kHEDISNucl_vbar_cc_p) {
    double wf2aux[13] = {  0.,  2.,  0.,  2.,  0.,  2.,  0.,  0.,  2.,  0.,  2.,  0.,  2. };
    double wf3aux[13] = {  0.,  2.,  0.,  2.,  0.,  2.,  0.,  0., -2.,  0., -2.,  0., -2. };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  else if (ch==kHEDISNucl_vbar_cc_n) {
    double wf2aux[13] = {  0.,  2.,  0.,  2.,  2.,  0.,  0.,  2.,  0.,  0.,  2.,  0.,  2. };
    double wf3aux[13] = {  0.,  2.,  0.,  2.,  2.,  0.,  0., -2.,  0.,  0., -2.,  0., -2. };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  else if (ch==kHEDISNucl_v_nc_p) {
    NCCouplings Au, Ad;
    Ad = GetNCCouplings(gvd,gad_nu);
    Au = GetNCCouplings(gvu,gau_nu);
    double wf2aux[13] = {  Au.c2,  Ad.c2,  Au.c2,  Ad.c2,  Au.c2,  Ad.c2, 0., Ad.c2, Au.c2, Ad.c2, Au.c2, Ad.c2, Au.c2 };
    double wf3aux[13] = { -Au.c3, -Ad.c3, -Au.c3, -Ad.c3, -Au.c3, -Ad.c3, 0., Ad.c3, Au.c3, Ad.c3, Au.c3, Ad.c3, Au.c3 };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  else if (ch==kHEDISNucl_v_nc_n) {
    NCCouplings Au, Ad;
    Ad = GetNCCouplings(gvd,gad_nu);
    Au = GetNCCouplings(gvu,gau_nu);
    double wf2aux[13] = {  Au.c2,  Ad.c2,  Au.c2,  Ad.c2,  Ad.c2,  Au.c2, 0., Au.c2, Ad.c2, Ad.c2, Au.c2, Ad.c2, Au.c2 };
    double wf3aux[13] = { -Au.c3, -Ad.c3, -Au.c3, -Ad.c3, -Ad.c3, -Au.c3, 0., Au.c3, Ad.c3, Ad.c3, Au.c3, Ad.c3, Au.c3 };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  else if (ch==kHEDISNucl_vbar_nc_p) {
    NCCouplings Au, Ad;
    Ad = GetNCCouplings(gvd,gad_nubar);
    Au = GetNCCouplings(gvu,gau_nubar);
    double wf2aux[13] = {  Au.c2,  Ad.c2,  Au.c2,  Ad.c2,  Au.c2,  Ad.c2, 0., Ad.c2, Au.c2, Ad.c2, Au.c2, Ad.c2, Au.c2 };
    double wf3aux[13] = { -Au.c3, -Ad.c3, -Au.c3, -Ad.c3, -Au.c3, -Ad.c3, 0., Ad.c3, Au.c3, Ad.c3, Au.c3, Ad.c3, Au.c3 };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }
  else if (ch==kHEDISNucl_vbar_nc_n) {
    NCCouplings Au, Ad;
    Ad = GetNCCouplings(gvd,gad_nubar);
    Au = GetNCCouplings(gvu,gau_nubar);
    double wf2aux[13] = {  Au.c2,  Ad.c2,  Au.c2,  Ad.c2,  Ad.c2,  Au.c2, 0., Au.c2, Ad.c2, Ad.c2, Au.c2, Ad.c2, Au.c2 };
    double wf3aux[13] = { -Au.c3, -Ad.c3, -Au.c3, -Ad.c3, -Ad.c3, -Au.c3, 0., Au.c3, Ad.c3, Ad.c3, Au.c3, Ad.c3, Au.c3 };
    std::copy( std::begin(wf2aux), std::end(wf2aux), std::begin(wf2) );
    std::copy( std::begin(wf3aux), std::end(wf3aux), std::begin(wf3) );
  }


  int nx  = ff_logx_array.size();
  int nq2 = ff_logq2_array.size();
  double xlist[nx*nq2];
  double q2list[nx*nq2];

  int nlist = 0;
  for ( const auto& x : ff_logx_array ) {
    for ( const auto& Q2 : ff_logq2_array ) {
      xlist[nlist]  = x;
      q2list[nlist] = Q2;
      nlist++;
    }
  }

  double F2list[nx*nq2];
  double FLlist[nx*nq2];
  double xF3list[nx*nq2];
  QCDNUM::zmstfun(1, wf2, xlist, q2list,  FLlist, nx*nq2, 0);
  QCDNUM::zmstfun(2, wf2, xlist, q2list,  F2list, nx*nq2, 0);
  QCDNUM::zmstfun(3, wf3, xlist, q2list, xF3list, nx*nq2, 0);
  for(int ff = 1; ff < kFFnumber; ++ff) {
    for (int i=0; i<nx*nq2; i++) {
      double tmp = 0;
      if      ( (HEDISFormFactorType_t)ff==kFFT1 ) tmp = (F2list[i]-FLlist[i])/2/xlist[i];
      else if ( (HEDISFormFactorType_t)ff==kFFT2 ) tmp = F2list[i];
      else if ( (HEDISFormFactorType_t)ff==kFFT3 ) tmp = xF3list[i] / xlist[i];
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
    for ( const auto& x : ff_logx_array ) {
      for ( const auto& Q2 : ff_logq2_array ) {
        double tmp = 0;
        for( int qch=frstqrkch; qch<=lastqrkch; qch++ ) tmp += fFormFactorsTables[(HEDISChannel_t)qch].Table[(HEDISFormFactorType_t)ff]->Evaluate(x,Q2);
        LOG("HEDISFormFactors", pDEBUG) << "NuclLOF" << ff << "[x=" << x << "," << Q2 << "] = " << tmp;
        ff_stream << tmp << "  ";
      }
    }
  }

  ff_stream.close();

}
//____________________________________________________________________________
FF_xQ2 HEDISFormFactors::EvalFFQrkLO( HEDISChannel_t ch, double x, double Q2 ) 
{
  FF_xQ2 ff;
  ff.F1 = fFormFactorsTables[ch].Table[kFFT1]->Evaluate(x,Q2);
  ff.F2 = fFormFactorsTables[ch].Table[kFFT2]->Evaluate(x,Q2);
  ff.F3 = fFormFactorsTables[ch].Table[kFFT3]->Evaluate(x,Q2);
  return ff;
}
//____________________________________________________________________________
FF_xQ2 HEDISFormFactors::EvalNuclFFLO( HEDISNuclChannel_t ch, double x, double Q2 ) 
{
  FF_xQ2 ff;
  ff.F1 = fLONuclFormFactorsTables[ch].Table[kFFT1]->Evaluate(x,Q2);
  ff.F2 = fLONuclFormFactorsTables[ch].Table[kFFT2]->Evaluate(x,Q2);
  ff.F3 = fLONuclFormFactorsTables[ch].Table[kFFT3]->Evaluate(x,Q2);
  return ff;
}
//____________________________________________________________________________
FF_xQ2 HEDISFormFactors::EvalNuclFFNLO( HEDISNuclChannel_t ch, double x, double Q2 ) 
{
  FF_xQ2 ff;
  ff.F1 = fNLONuclFormFactorsTables[ch].Table[kFFT1]->Evaluate(x,Q2);
  ff.F2 = fNLONuclFormFactorsTables[ch].Table[kFFT2]->Evaluate(x,Q2);
  ff.F3 = fNLONuclFormFactorsTables[ch].Table[kFFT3]->Evaluate(x,Q2);
  return ff;
}


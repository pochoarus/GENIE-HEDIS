
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

#include <TSystem.h>
#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include <string>
#include <iostream>
#include <fstream>

#ifdef __GENIE_APFEL_ENABLED__
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
LHAPDF::PDF* pdf;
#endif

using namespace std;

using namespace genie;
using namespace genie::constants;

int fNucPdg = 0;

#ifdef __GENIE_APFEL_ENABLED__
class PhotonConv: public ROOT::Math::IBaseFunctionOneDim
{
  public:
    PhotonConv(double x, double q) : ROOT::Math::IBaseFunctionOneDim(),xmin(x),Qin(q){};
   ~PhotonConv(){};
    ROOT::Math::IBaseFunctionOneDim * Clone  (void)     const { return new PhotonConv(xmin,Qin); }
    unsigned int                      NDim   (void)     const { return 1; }
    double                            DoEval (double y) const { return 2. * ( TMath::Power(xmin/y,2)+TMath::Power( 1.-xmin/y, 2) ) * pdf->xfxQ(22, y, Qin); }
    void                              SetPar (double x, double q) { xmin=x; Qin=q; }

  private:
    double xmin,Qin;
};

// External functions accessed by apfel to parameterise xf(x0,Q0)
extern "C" void externalsetapfellept_(double* x, double* q, int* irep, double* xf, double* xl) {
  if (*x >= 1 || *x < 0) {
    for ( int i=0; i<13; i++ ) xf[i] = 0.;
    for ( int i=0; i<7; i++ )  xl[i] = 0.;
    return;
  }
  else{
    for ( int i=0; i<13; i++ ) {
      xf[i] = pdf->xfxQ(i-6, *x, *q);
      if ( pdg::IsNeutron(fNucPdg) ) {
        if     ( i==4 ) xf[i] = pdf->xfxQ(-1, *x, *q); 
        else if( i==5 ) xf[i] = pdf->xfxQ(-2, *x, *q);
        else if( i==7 ) xf[i] = pdf->xfxQ( 2, *x, *q);
        else if( i==8 ) xf[i] = pdf->xfxQ( 1, *x, *q);
      }
    }
    for ( int i=0; i<7; i++ ) {
      if ( pdg::IsProton (fNucPdg) ) {
        if      (i==0 || i==6) xl[i] = 0.;
        else if (i==3 )        xl[i] = pdf->xfxQ(22, *x, *q);
        else{       
          double mlep;
          if      (i==1 || i==5) mlep = kMuonMass;
          else if (i==2 || i==4) mlep = kElectronMass;
          ROOT::Math::IBaseFunctionOneDim * func = new PhotonConv(*x,*q);
          ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;
          ROOT::Math::Integrator ig(*func,ig_type,1,0.01,100000);
          double res = ig.Integral(*x,1.);
          xl[i] = APFEL::AlphaQED(*q)/2./kPi*TMath::Log( *q/mlep ) * res;
        }
      }
      else if ( pdg::IsNeutron(fNucPdg) ) xl[i] = 0.0;
    }
  }
  return;
}
#endif

//____________________________________________________________________________
int main(int argc, char ** argv)
{

  const int nx = 1000.;

  int nucs[2] = { kPdgProton, kPdgNeutron };
  int pdgs[6] = { kPdgNuE, kPdgAntiNuE, kPdgNuMu, kPdgAntiNuMu, kPdgNuTau, kPdgAntiNuTau };
  
#ifdef __GENIE_APFEL_ENABLED__
  for (int k=0; k<2; k++) {
    
    fNucPdg = nucs[k];

    // initialising APFEL framework
    LOG("gmkglrersf", pINFO) << "Initialising APFEL..." ; 
    string pdfset;
    if      (pdg::IsProton (fNucPdg) ) pdfset = "NNPDF31_nnlo_as_0118_luxqed";
    else if (pdg::IsNeutron(fNucPdg) ) pdfset = "NNPDF31_nnlo_as_0118";

    pdf = LHAPDF::mkPDF(pdfset,0);

    double xPDFmin,QPDFmin,QPDFmax,mc,mb,mt;
    xPDFmin = pdf->xMin();
    QPDFmin = pdf->qMin();
    QPDFmax = pdf->qMax();
    mc      = pdf->quarkMass(4);
    mb      = pdf->quarkMass(5);
    mt      = pdf->quarkMass(6);
    LOG("gmkglrersf", pINFO) << "xPDFmin = " << xPDFmin;
    LOG("gmkglrersf", pINFO) << "QPDFmin = " << QPDFmin;
    LOG("gmkglrersf", pINFO) << "QPDFmax = " << QPDFmax;
    LOG("gmkglrersf", pINFO) << "mc = " << mc;
    LOG("gmkglrersf", pINFO) << "mb = " << mb;
    LOG("gmkglrersf", pINFO) << "mt = " << mt;

    APFEL::CleanUp();

    APFEL::SetPDFSet(pdfset);
    APFEL::SetReplica(0);
    APFEL::SetPerturbativeOrder(2);
    APFEL::SetQLimits(QPDFmin,QPDFmax);
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1,100,5,1e-9);
    APFEL::SetGridParameters(2,40,3,1e-1);
    APFEL::SetGridParameters(3,60,3,8e-1);
    APFEL::SetGFermi(kGF);
    APFEL::SetPoleMasses(mc,mb,mt);
    APFEL::SetTheory("QUniD");
    APFEL::EnableLeptonEvolution(true);
    APFEL::SetFastEvolution(true);
    APFEL::SetPDFSet("leptexternal");
    APFEL::InitializeAPFEL();

    LOG("gmkglrersf", pWARN) << "Init EvolveAPFEL";        
    APFEL::EvolveAPFEL(QPDFmin,kMw);
    LOG("gmkglrersf", pWARN) << "End EvolveAPFEL";        

    // open file in which SF will be stored

    double x[nx];
    for ( int i=0; i<nx; i++ ) x[i] = TMath::Power( 10, TMath::Log10(xPDFmin) + i*(TMath::Log10(1.)-TMath::Log10(xPDFmin))/(1000.-1) );

    for(int j=0; j<6; j++) {
      string SFname = string(gSystem->Getenv("GENIE")) + "/data/evgen/glres/PhotonSF_hitnuc"+to_string(fNucPdg)+"_hitlep"+to_string(pdgs[j])+".dat";
      std::ofstream sf_stream(SFname);
      for ( int i=0; i<nx; i++ ) {
        double tmp = 0;
        if      ( pdg::IsNuE      (pdgs[j]) ) tmp = APFEL::xLepton( 1,x[i]);
        else if ( pdg::IsAntiNuE  (pdgs[j]) ) tmp = APFEL::xLepton(-1,x[i]);
        else if ( pdg::IsNuMu     (pdgs[j]) ) tmp = APFEL::xLepton( 2,x[i]);
        else if ( pdg::IsAntiNuMu (pdgs[j]) ) tmp = APFEL::xLepton(-2,x[i]);
        else if ( pdg::IsNuTau    (pdgs[j]) ) tmp = APFEL::xLepton( 3,x[i]);
        else if ( pdg::IsAntiNuTau(pdgs[j]) ) tmp = APFEL::xLepton(-3,x[i]);
        LOG("gmkglrersf", pWARN) << "SF " << pdgs[j] << " [x=" << x[i] << "] = " << tmp;
        sf_stream << x[i] << " " << tmp << endl;
      }
      // Close file in which SF are stored
      sf_stream.close();
    }        
    

  }
#endif

}



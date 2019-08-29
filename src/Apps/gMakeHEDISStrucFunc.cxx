//____________________________________________________________________________
/*!
\program gmkspl
\brief   GENIE utility program building Structure Functions needed for HEDIS
         package.
         Syntax :
           gmksf [-h]
                 --outdir output_sf_dir_name
                 --lhapdf lhapdf_name
                 --mz mass_boson_z
                 --mw mass_boson_w
                 --rho weinberg_angle_em_correction
                 [--sin2thw sin_weinberg_angle]
                 --ckm ckm_matrix
                 --qrkthrs quark_prod_threshold
                 [--nlo]
                 [--scheme nlo_mass_scheme]
                 --x x_log10_grid
                 --q2 q2_log10_grid
         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.
         Options :
           --outdir
               Name of output directory where structure functions will be 
               stored.
           --lhapdf
               Name of the pdf set used to compute the SF using LHAPDF 
           --mz
               Mass of the Z boson in GeV
           --mw
               Mass of the W boson in GeV
           --rho
               EM correction for the weinberg angle [Sin2ThW = 1.-MW2/MZ2/(1+rho)]
               If it is set to 0, then sin2thw must be provided by user
           --sin2thw
               Weinberg angle. Only input when rho=0
           --ckm
               CKM matrix (9 elements). Input them with comma separated values.
           --qrkthrs
               Quark production threshold. Four options are available:
                 0: No threshold           & z=x
                 1: W>M_nucleon+M_finalqrk & z=x
                 2: W>M_nucleon+M_finalqrk & z=x(1+M_finalqrk2/Q2) 
                 3: No threshold           & z=x(1+M_finalqrk2/Q2) 
           --nlo
               The SF will be computed using NLO coefficients. This option can be
               used if APFEL is available.
           --scheme
               Mass scheme used to account for heavy quark effects. Two options:
                 BRG: Mass scheme used in arXiv:1808.02034
                 CSMS: Mass scherme used in arXiv:1106.3723
           --x
               Number of bins in log10x and min value of x for the SF grid.
               Two values must be comma separated.
           --q2
               Number of bins in log10q2 and min/max value of Q2 for the SF grid.
               Three values must be comma separated.
        ***  See the User Manual for more details and examples. ***
\author  Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)
\created July 8, 2019
\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Interaction/HEDISChannel.h"

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

void   GetCommandLineArgs (int argc, char ** argv);
void   PrintSyntax        (void);

void   CreateQrkStrucFunc (HEDISQrkChannel_t ch, vector<double> sf_x_array, vector<double> sf_q2_array);

#ifdef __GENIE_APFEL_ENABLED__
void   CreateNucStrucFunc (HEDISNucChannel_t ch, vector<double> sf_x_array, vector<double> sf_q2_array);
#endif

// User-specified options:
string gSFname = "";             // Name output directory where SF are stored
string gLHAPDFmember;            // Name PDF set in LHAPDF
double gMZ;                      // Mass Z boson
double gMW;                      // Mass W boson
double gRho;                     // EM correction Weinberg angle
double gSin2ThW;                 // Weinberg angle
double gVud,gVus,gVub;           // CKM matrix elements
double gVcd,gVcs,gVcb;
double gVtd,gVts,gVtb;
double gQrkThrs;                 // Quark production threshold
bool gNLO;                       // Use NLO coefficients
string gMassScheme;              // Mass scheme in NLO computation
int gNX;                         // Number of bins in x for SF grid
double gxGRIDmin;                // Minimum values of x for SF grid
int gNQ2;                        // Number of bins in Q2 for SF grid
double gQ2GRIDmin;               // Minimum values of Q2 for SF grid
double gQ2GRIDmax;               // Maximum values of Q2 for SF grid

// values from LHAPDF set
double xPDFmin;                  // Minimum values of x in grid from LHPADF set
double Q2PDFmin;                 // Minimum values of Q2 in grid from LHPADF set
double Q2PDFmax;                 // Maximum values of Q2 in grid from LHPADF set
std::map<int, double> mPDFQrk;   // Mass of the quark from LHAPDF set

//____________________________________________________________________________
int main(int argc, char ** argv)
{

  GetCommandLineArgs (argc,argv);   // Get command line arguments

  // creating directory where structure functions will be stored and check if it already exists
  if ( gSystem->mkdir(gSFname.c_str())==0 ) LOG("gmkhedissf", pINFO) << "Creating structure functions directory: " << gSFname;
  else {
    LOG("gmkhedissf", pERROR) << "Directory already exists.";
    return 0;
  }

  // initialising lhapdf
#ifdef __GENIE_LHAPDF6_ENABLED__
  LOG("gmkhedissf", pINFO) << "Initialising LHAPDF6...";
  pdf = LHAPDF::mkPDF(gLHAPDFmember, 0);
  xPDFmin  = pdf->xMin();
  Q2PDFmin = pdf->q2Min();
  Q2PDFmax = pdf->q2Max();
  LOG("gmkhedissf", pINFO) << "PDF info:";
  LOG("gmkhedissf", pINFO) << "OrderQCD = " << pdf->orderQCD();
  LOG("gmkhedissf", pINFO) << "FlavorScheme = " << pdf->info().get_entry("FlavorScheme");
  LOG("gmkhedissf", pINFO) << "NumFlavors = " << pdf->info().get_entry("NumFlavors");
  LOG("gmkhedissf", pINFO) << "Xmin = " << xPDFmin << "  Xmax = " << pdf->xMax() << "  Q2min = " << Q2PDFmin << "  Q2max = " << Q2PDFmax;
  LOG("gmkhedissf", pINFO) << "MZ = " << pdf->info().get_entry("MZ");
  for (int i=1; i<7; i++) {
    mPDFQrk[i] = pdf->quarkMass(i);
    LOG("gmkhedissf", pINFO) << "M" << i << " = " << mPDFQrk[i];
  }
#endif
#ifdef __GENIE_LHAPDF5_ENABLED__
  LOG("gmkhedissf", pINFO) << "Initialising LHAPDF5...";
  LHAPDF::initPDFByName(fLHAPDFmember, LHAPDF::LHGRID, 0);
  xPDFmin  = LHAPDF::getXmin(0);
  Q2PDFmin = LHAPDF::getQ2min(0);
  Q2PDFmax = LHAPDF::getQ2max(0);
  LOG("gmkhedissf", pINFO) << "PDF info:";
  LOG("gmkhedissf", pINFO) << "Xmin = " << xPDFmin << "  Xmax = " << pdf->xMax() << "  Q2min = " << Q2PDFmin << "  Q2max = " << Q2PDFmax;
  for (int i=1; i<7; i++) {
    mPDFQrk[i] = LHAPDF::getQMass(i);
    LOG("gmkhedissf", pINFO) << "M" << i << " = " << mPDFQrk[i];
  }
#endif

  // checking if LHAPDF boundaries matches boundaries defined by user
  if ( gxGRIDmin < xPDFmin ) {
    LOG("gmkhedissf", pWARN) << "Lower boundary in X is smaller than input PDF";
    LOG("gmkhedissf", pWARN) << "xPDFmin = "  << xPDFmin;
    LOG("gmkhedissf", pWARN) << "xGRIDmin = " << gxGRIDmin;
  }
  else if ( gQ2GRIDmin < Q2PDFmin ) {
    LOG("gmkhedissf", pWARN) << "Lower boundary in Q2 is smaller than input PDF";
    LOG("gmkhedissf", pWARN) << "Q2PDFmin = "  << Q2PDFmin;
    LOG("gmkhedissf", pWARN) << "Q2GRIDmin = " << gQ2GRIDmin;
  }
  else if ( gQ2GRIDmax > Q2PDFmax ) {
    LOG("gmkhedissf", pWARN) << "Upper boundary in Q2 is bigger than input PDF";
    LOG("gmkhedissf", pWARN) << "Q2PDFmax = "  << Q2PDFmax;
    LOG("gmkhedissf", pWARN) << "Q2GRIDmax = " << gQ2GRIDmax;
  }

  // create binning for the SF grid in log10 space for x and Q2
  LOG("gmkhedissf", pINFO) << "Grid x,Q2 :" << gNX << " , " << gNQ2;
  vector<double> sf_x_array;
  vector<double> sf_q2_array;
  double dlogq2 = TMath::Abs( TMath::Log10(gQ2GRIDmin)-TMath::Log10(gQ2GRIDmax) ) / gNQ2;
  double dlogx  = TMath::Abs( TMath::Log10(gxGRIDmin)-TMath::Log10(1.) ) / gNX;
  for ( double logq2 = TMath::Log10(gQ2GRIDmin); logq2<TMath::Log10(gQ2GRIDmax); logq2+= dlogq2 ) {
    double q2 = TMath::Power( 10, logq2 + 0.5*dlogq2 );
    if (q2>gQ2GRIDmax) continue;
    sf_q2_array.push_back(q2);
    LOG("gmkhedissf", pDEBUG) << "q2: " << sf_q2_array.back();
  }
  for ( double logx = TMath::Log10(gxGRIDmin); logx<TMath::Log10(1.); logx+= dlogx ) {
    double x = TMath::Power( 10, logx + 0.5*dlogx );
    if ( x>1. ) continue;
    sf_x_array.push_back(x);
    LOG("gmkhedissf", pDEBUG) << "x: " << sf_x_array.back();
  }


  // create file in which all the information used to compute the SF is saved
  LOG("gmkhedissf", pINFO) << "Creating Metafile with input information";
  string filename = gSFname + "/Inputs.txt";
  std::ofstream meta_stream(filename.c_str());
  meta_stream << "# NX" << std::endl;
  meta_stream << gNX << std::endl;
  meta_stream << "# Xmin" << std::endl;
  meta_stream << gxGRIDmin << std::endl;
  meta_stream << "# NQ2" << std::endl;
  meta_stream << gNQ2 << std::endl;
  meta_stream << "# Q2min" << std::endl;
  meta_stream << gQ2GRIDmin << std::endl;
  meta_stream << "# Q2max" << std::endl;
  meta_stream << gQ2GRIDmax << std::endl;
  meta_stream << "# LHAPDF member" << std::endl;
  meta_stream << gLHAPDFmember << std::endl;
  meta_stream << "# Mass Z" << std::endl;
  meta_stream << gMZ << std::endl;
  meta_stream << "# Mass W" << std::endl;
  meta_stream << gMW << std::endl;
  meta_stream << "# Rho" << std::endl;
  meta_stream << gRho << std::endl;
  meta_stream << "# Sin2ThW" << std::endl;
  meta_stream << gSin2ThW << std::endl;
  meta_stream << "# Quark threshold" << std::endl;
  meta_stream << gQrkThrs << std::endl;
  meta_stream << "# CKM" << std::endl;
  meta_stream << gVud << std::endl;
  meta_stream << gVus << std::endl;
  meta_stream << gVub << std::endl;
  meta_stream << gVcd << std::endl;
  meta_stream << gVcs << std::endl;
  meta_stream << gVcb << std::endl;
  meta_stream << gVtd << std::endl;
  meta_stream << gVts << std::endl;
  meta_stream << gVtb << std::endl;
  meta_stream << "# Mass quarks" << std::endl;
  meta_stream << mPDFQrk[0] << std::endl;
  meta_stream << mPDFQrk[1] << std::endl;
  meta_stream << mPDFQrk[2] << std::endl;
  meta_stream << mPDFQrk[3] << std::endl;
  meta_stream << mPDFQrk[4] << std::endl;
  meta_stream << mPDFQrk[5] << std::endl;
  meta_stream << "# NLO" << std::endl;
  meta_stream << gNLO << std::endl;
  if (gNLO) {
    meta_stream << "# Mass Scheme" << std::endl;
    meta_stream << gMassScheme << std::endl;
  }
  meta_stream.close();

  // load structure functions for each quark channel at LO
  for( int ch=1; ch<kHEDISQrk_numofchannels; ch++ ) CreateQrkStrucFunc( (HEDISQrkChannel_t)ch, sf_x_array, sf_q2_array );


#ifdef __GENIE_APFEL_ENABLED__
  // load structure functions for each nucleon at NLO
  if (gNLO) {
    // initialising APFEL framework
    LOG("gmkhedissf", pINFO) << "Initialising APFEL..." ; 
    APFEL::SetPDFSet(gLHAPDFmember);
    APFEL::SetReplica(0);
    if (gMassScheme=="BRG") {
      APFEL::SetMassScheme("FONLL-B");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[6]);
    } 
    else if (gMassScheme=="CSMS") {
      APFEL::SetMassScheme("ZM-VFNS");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[5]+0.1);
    }
    APFEL::SetQLimits(TMath::Sqrt(gQ2GRIDmin),TMath::Sqrt(gQ2GRIDmax));
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1,90,3,gxGRIDmin);
    APFEL::SetGridParameters(2,50,5,1e-1);
    APFEL::SetGridParameters(3,40,5,8e-1);
    APFEL::SetPerturbativeOrder(1);
    APFEL::SetAlphaQCDRef(pdf->alphasQ(gMZ),gMZ);    
    APFEL::SetProtonMass(kProtonMass);
    APFEL::SetWMass(gMW);
    APFEL::SetZMass(gMZ);
    APFEL::SetSin2ThetaW(gSin2ThW);
    APFEL::SetCKM(gVud, gVus, gVub,
                  gVcd, gVcs, gVcb,
                  gVtd, gVts, gVtb);

    // load structure functions for each nucleon channel at NLO
    for( int ch=1; ch<kHEDISNuc_numofchannels; ch++ ) CreateNucStrucFunc( (HEDISNucChannel_t)ch, sf_x_array, sf_q2_array );
  }
#endif

  delete pdf;

}
//____________________________________________________________________________
void CreateQrkStrucFunc( HEDISQrkChannel_t ch, vector<double> sf_x_array, vector<double> sf_q2_array ) 
{

  // variables used to tag the SF for particular channel
  int pdg_nucl  = HEDISChannel::HitNuclPdg(ch);   // PDG of target nucleon
  bool sea_iq   = HEDISChannel::HitQuarkSea(ch);  // Valence/sea hit quark
  int pdg_iq    = HEDISChannel::HitQuarkPdg(ch);  // PDG hit quark
  int pdg_fq    = HEDISChannel::FnlQuarkPdg(ch);  // PDG produced quark

  // variables used for certain quark production threshold options
  double mass_fq   = mPDFQrk[TMath::Abs(pdg_fq)];                     // mass produced quark
  double mass_nucl = (pdg_nucl==2212) ? kProtonMass : kNeutronMass;   // mass target nucleon

  // up and down quark swicth depending on proton or neutron interaction  
  int qrkd = 0;
  int qrku = 0;
  if      ( pdg_nucl==2212 ) { qrkd = 1 ; qrku = 2; }
  else if ( pdg_nucl==2112 ) { qrkd = 2 ; qrku = 1; }

  // variables associated to the PDF and coupling of the quarks
  int qpdf1 = -999;                                     // number used to compute PDF
  int qpdf2 = -999;                                     // auxiliary number used to compute PDF for valence quarks
  double Cp2 = -999;                                    // couping for F1,F2
  double Cp3 = -999;                                    // couping for F3
  double sign3 = (HEDISChannel::IsNu(ch)) ? +1. : -1.;  // sign change for nu/nubar in F3
  if ( HEDISChannel::InteractionType(ch) == kIntWeakCC ) {
    if ( HEDISChannel::IsNu(ch) ) {
      if      ( pdg_iq== 1 && !sea_iq && pdg_fq== 2 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*gVud*gVud; Cp3 =  2*gVud*gVud; }
      else if ( pdg_iq== 1 && !sea_iq && pdg_fq== 4 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*gVcd*gVcd; Cp3 =  2*gVcd*gVcd; }
      else if ( pdg_iq== 1 && !sea_iq && pdg_fq== 6 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*gVtd*gVtd; Cp3 =  2*gVtd*gVtd; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 2 ) { qpdf1 = -qrkd;                Cp2 = 2*gVud*gVud; Cp3 =  2*gVud*gVud; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 4 ) { qpdf1 = -qrkd;                Cp2 = 2*gVcd*gVcd; Cp3 =  2*gVcd*gVcd; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 6 ) { qpdf1 = -qrkd;                Cp2 = 2*gVtd*gVtd; Cp3 =  2*gVtd*gVtd; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 2 ) { qpdf1 =  3;                   Cp2 = 2*gVus*gVus; Cp3 =  2*gVus*gVus; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  3;                   Cp2 = 2*gVcs*gVcs; Cp3 =  2*gVcs*gVcs; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 6 ) { qpdf1 =  3;                   Cp2 = 2*gVts*gVts; Cp3 =  2*gVts*gVts; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 2 ) { qpdf1 =  5;                   Cp2 = 2*gVub*gVub; Cp3 =  2*gVub*gVub; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  5;                   Cp2 = 2*gVcb*gVcb; Cp3 =  2*gVcb*gVcb; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 6 ) { qpdf1 =  5;                   Cp2 = 2*gVtb*gVtb; Cp3 =  2*gVtb*gVtb; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -qrku;                Cp2 = 2*gVud*gVud; Cp3 = -2*gVud*gVud; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -qrku;                Cp2 = 2*gVus*gVus; Cp3 = -2*gVus*gVus; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -qrku;                Cp2 = 2*gVub*gVub; Cp3 = -2*gVub*gVub; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -4;                   Cp2 = 2*gVcd*gVcd; Cp3 = -2*gVcd*gVcd; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -4;                   Cp2 = 2*gVcs*gVcs; Cp3 = -2*gVcs*gVcs; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -4;                   Cp2 = 2*gVcb*gVcb; Cp3 = -2*gVcb*gVcb; }
    }
    else {
      if      ( pdg_iq== 2 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*gVud*gVud; Cp3 =  2*gVud*gVud; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 3 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*gVus*gVus; Cp3 =  2*gVus*gVus; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 5 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*gVub*gVub; Cp3 =  2*gVub*gVub; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrku;                Cp2 = 2*gVud*gVud; Cp3 =  2*gVud*gVud; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 3 ) { qpdf1 = -qrku;                Cp2 = 2*gVus*gVus; Cp3 =  2*gVus*gVus; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 5 ) { qpdf1 = -qrku;                Cp2 = 2*gVub*gVub; Cp3 =  2*gVub*gVub; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 1 ) { qpdf1 =  4;                   Cp2 = 2*gVcd*gVcd; Cp3 =  2*gVcd*gVcd; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  4;                   Cp2 = 2*gVcs*gVcs; Cp3 =  2*gVcs*gVcs; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  4;                   Cp2 = 2*gVcb*gVcb; Cp3 =  2*gVcb*gVcb; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrkd;                Cp2 = 2*gVud*gVud; Cp3 = -2*gVud*gVud; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -qrkd;                Cp2 = 2*gVcd*gVcd; Cp3 = -2*gVcd*gVcd; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -qrkd;                Cp2 = 2*gVtd*gVtd; Cp3 = -2*gVtd*gVtd; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -3;                   Cp2 = 2*gVus*gVus; Cp3 = -2*gVus*gVus; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -3;                   Cp2 = 2*gVcs*gVcs; Cp3 = -2*gVcs*gVcs; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -3;                   Cp2 = 2*gVts*gVts; Cp3 = -2*gVts*gVts; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -5;                   Cp2 = 2*gVub*gVub; Cp3 = -2*gVub*gVub; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -5;                   Cp2 = 2*gVcb*gVcb; Cp3 = -2*gVcb*gVcb; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -5;                   Cp2 = 2*gVtb*gVtb; Cp3 = -2*gVtb*gVtb; }
    }
  }
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) {           
    double c2u = TMath::Power( 1./2. - 4./3.*gSin2ThW,2) + 1./4.;
    double c2d = TMath::Power(-1./2. + 2./3.*gSin2ThW,2) + 1./4.;
    double c3u = 1./2. - 4./3.*gSin2ThW;
    double c3d = 1./2. - 2./3.*gSin2ThW;
    if      ( pdg_iq== 1 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 2 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = c2u; Cp3 =  c3u; }
    else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrkd;                Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 2 ) { qpdf1 = -qrku;                Cp2 = c2u; Cp3 =  c3u; }
    else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  3;                   Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  4;                   Cp2 = c2u; Cp3 =  c3u; }
    else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  5;                   Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -qrkd;                Cp2 = c2d; Cp3 = -c3d; }
    else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrku;                Cp2 = c2u; Cp3 = -c3u; }
    else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -3;                   Cp2 = c2d; Cp3 = -c3d; }
    else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -4;                   Cp2 = c2u; Cp3 = -c3u; }
    else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -5;                   Cp2 = c2d; Cp3 = -c3d; }
  }   

  // open file in which SF will be stored
  LOG("gmkhedissf", pINFO) << "Creating file QrkSF at LO -> " << HEDISChannel::AsString(ch);
  string filename = gSFname + "/QrkSF_LO_" + HEDISChannel::AsString(ch) + ".dat";
  std::ofstream sf_stream(filename.c_str());

  // loop over 3 different SF: F1,F2,F3
  for(int sf = 1; sf < 4; sf++) {
    for ( unsigned int i=0; i<sf_q2_array.size(); i++ ) {
      double Q2 = sf_q2_array[i];
      for ( unsigned int j=0; j<sf_x_array.size(); j++ ) {
        double x = sf_x_array[j];

        double z = x; // this variable is introduce in case you want to apply scaling

        // W threshold
        if      (gQrkThrs==1) { 
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { sf_stream << 0. << "  "; continue; }
        } 
        // W threshold and slow rescaling
        else if (gQrkThrs==2) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { sf_stream << 0. << "  "; continue; }
            z *= 1+mass_fq*mass_fq/Q2;
        }
        // Slow rescaling
        else if (gQrkThrs==3) {
            z *= 1+mass_fq*mass_fq/Q2;
        }

        // Fill x,Q2 used to extract PDF. If values outside boundaries then freeze them.
        double xPDF = TMath::Max( z, xPDFmin );
        double Q2PDF = TMath::Max( Q2, Q2PDFmin );
        Q2PDF = TMath::Min( Q2PDF, Q2PDFmax  );

        // Extract PDF requiring then to be higher than zero
        double fPDF = fmax( pdf->xfxQ2(qpdf1, xPDF, Q2PDF)/z , 0.);
        if (qpdf2!= -999) fPDF -= fmax( pdf->xfxQ2(qpdf2, xPDF, Q2PDF)/z , 0.);
                    
        // Compute SF
        double tmp = -999;
        if      ( sf==1 ) tmp = fPDF*Cp2/2;
        else if ( sf==2 ) tmp = fPDF*Cp2*z;
        else if ( sf==3 ) tmp = fPDF*Cp3*sign3;

        // Save SF for particular x and Q2 in file
        LOG("gmkhedissf", pDEBUG) << "QrkSFLO" << sf << "[x=" << x << "," << Q2 << "] = " << tmp;
        sf_stream << tmp << "  ";
        
      }
    }
  }

  // Close file in which SF are stored
  sf_stream.close();

}
#ifdef __GENIE_APFEL_ENABLED__
//____________________________________________________________________________
void CreateNucStrucFunc( HEDISNucChannel_t ch, vector<double> sf_x_array, vector<double> sf_q2_array )
{

  // Define the channel that is used in APFEL
  if ( HEDISChannel::IsNu(ch) ) APFEL::SetProjectileDIS("neutrino");
  else                          APFEL::SetProjectileDIS("antineutrino");
  if      ( HEDISChannel::InteractionType(ch) == kIntWeakCC ) APFEL::SetProcessDIS("CC");
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) APFEL::SetProcessDIS("NC");
  if      ( HEDISChannel::HitNuclPdg(ch)==2212 ) APFEL::SetTargetDIS("proton");
  else if ( HEDISChannel::HitNuclPdg(ch)==2112 ) APFEL::SetTargetDIS("neutron");

  APFEL::InitializeAPFEL_DIS();

  // Using APFEL format to store the SF grid
  int nx  = sf_x_array.size();
  int nq2 = sf_q2_array.size();
  double xlist[nx*nq2];
  double q2list[nx*nq2];
  double F2list[nx*nq2];
  double FLlist[nx*nq2];
  double xF3list[nx*nq2];

  int nlist = 0;
  for ( unsigned int i=0; i<sf_q2_array.size(); i++ ) {
    double Q2 = sf_q2_array[i];
    double Q  = TMath::Sqrt(Q2);
    // SF from APFEL are multiplied by a prefactor in NC. We dont want that prefactor
    double norm = (HEDISChannel::InteractionType(ch)==kIntWeakCC) ? 1. : 2./TMath::Power( Q2/(Q2 + TMath::Power(APFEL::GetZMass(),2))/4/APFEL::GetSin2ThetaW()/(1-APFEL::GetSin2ThetaW()), 2 );
    APFEL::SetAlphaQCDRef(pdf->alphasQ(Q),Q);
    APFEL::ComputeStructureFunctionsAPFEL(Q,Q);
    for ( unsigned int j=0; j<sf_x_array.size(); j++ ) {
      double x = sf_x_array[j];
      q2list[nlist]  = Q2;
      xlist[nlist]   = x;
      FLlist[nlist]  = norm*APFEL::FLtotal(x);
      F2list[nlist]  = norm*APFEL::F2total(x);
      xF3list[nlist] = norm*APFEL::F3total(x);
      nlist++;
    }
  }

  // open file in which SF will be stored
  LOG("gmkhedissf", pINFO) << "Creating file Nucleon Structure Function at NLO -> " << HEDISChannel::AsString(ch);
  string filename = gSFname + "/NucSF_NLO_" + HEDISChannel::AsString(ch) + ".dat";
  std::ofstream sf_stream(filename.c_str());

  double sign3 = (HEDISChannel::IsNu(ch)) ? +1. : -1.;  // sign change for nu/nubar in F3
  // loop over 3 different SF: F1,F2,F3
  for(int sf = 1; sf < 4; sf++) {
    for (int i=0; i<nx*nq2; i++) {
      double tmp = 0;
      if      ( sf==1 ) tmp = (F2list[i]-FLlist[i])/2/xlist[i];
      else if ( sf==2 ) tmp = F2list[i];
      else if ( sf==3 ) tmp = sign3 * xF3list[i] / xlist[i];
      // Save SF for particular x and Q2 in file
      LOG("gmkhedissf", pDEBUG) << "NucSFNLO" << sf << "[x=" << xlist[i] << "," << q2list[i] << "] = " << tmp;
      sf_stream << tmp << "  ";
    }
  }
    
  // Close file in which SF are stored
  sf_stream.close();

}
#endif
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  CmdLnArgParser parser(argc,argv);

  //help
  if(parser.OptionExists('h')) {
      PrintSyntax();
      exit(0);
  }

  //output directory for SF
  if     (parser.OptionExists("outdir")){
    gSFname = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/sf/" + parser.ArgAsString("outdir");
    LOG("gmkhedissf", pDEBUG) << "gSFname: " << gSFname;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified SF name - Exiting";
    PrintSyntax();
    exit(1);
  }

  //lhapdf set name
  if(parser.OptionExists("lhapdf")){
    gLHAPDFmember = parser.ArgAsString("lhapdf");
    LOG("gmkhedissf", pDEBUG) << "gLHAPDFmember: " << gLHAPDFmember;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified LHAPDF name - Exiting";
    PrintSyntax();
    exit(1);
  }

  //mass Z boson
  if(parser.OptionExists("mz")){
    gMZ = parser.ArgAsDouble("mz");
    LOG("gmkhedissf", pDEBUG) << "gMZ: " << gMZ;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified mass Z boson - Exiting";
    PrintSyntax();
    exit(1);
  }

  //mass W boson
  if(parser.OptionExists("mw")){
    gMW = parser.ArgAsDouble("mw");
    LOG("gmkhedissf", pDEBUG) << "gMW: " << gMW;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified mass W boson - Exiting";
    PrintSyntax();
    exit(1);
  }

  //EM correction Weinberg angle
  if(parser.OptionExists("rho")){
    gRho = parser.ArgAsDouble("rho");
    LOG("gmkhedissf", pDEBUG) << "gRho: " << gRho;
    if (gRho==0) {
      //Weinberg angle when rho=0
      if(parser.OptionExists("sin2thw")){
        gSin2ThW = parser.ArgAsDouble("sin2thw");
      }
      else {
        LOG("gmkhedissf", pFATAL) << "Unspecified Weignber angle - Exiting";
        PrintSyntax();
        exit(1);
      }
    } 
    else gSin2ThW = 1.-gMW*gMW/gMZ/gMZ/(1+gRho);
    LOG("gmkhedissf", pDEBUG) << "gSin2ThW: " << gSin2ThW;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified Raditive correction Rho - Exiting";
    PrintSyntax();
    exit(1);
  }


  //CKM matrix elements
  if(parser.OptionExists("ckm")){
    string sckm = parser.ArgAsString("ckm");
    // split the comma separated list
    vector<string> vsckm = utils::str::Split(sckm, ",");
    assert(vsckm.size() == 9);
    gVud = atof(vsckm[0].c_str());
    gVus = atof(vsckm[1].c_str());
    gVub = atof(vsckm[2].c_str());
    gVcd = atof(vsckm[3].c_str());
    gVcs = atof(vsckm[4].c_str());
    gVcb = atof(vsckm[5].c_str());
    gVtd = atof(vsckm[6].c_str());
    gVts = atof(vsckm[7].c_str());
    gVtb = atof(vsckm[8].c_str());
    LOG("gmkhedissf", pDEBUG) << "CKM: ";
    LOG("gmkhedissf", pDEBUG) << gVud << "  " << gVus << "  " << gVub;
    LOG("gmkhedissf", pDEBUG) << gVcd << "  " << gVcs << "  " << gVcb;
    LOG("gmkhedissf", pDEBUG) << gVtd << "  " << gVts << "  " << gVtb;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified CKM matrix - Exiting";
    PrintSyntax();
    exit(1);
  }

  //Quark production threshold
  if(parser.OptionExists("qrkthrs")){
    gQrkThrs = parser.ArgAsInt("qrkthrs");
    LOG("gmkhedissf", pDEBUG) << "gQrkThrs: " << gQrkThrs;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified Quark threshold - Exiting";
    PrintSyntax();
    exit(1);
  }

  //Compute SF with NLO coefficients.
  gNLO = parser.OptionExists("nlo");
  LOG("gmkhedissf", pDEBUG) << "gNLO: " << gNLO;
  if (gNLO) {
    if(parser.OptionExists("scheme")){
      gMassScheme = parser.ArgAsString("scheme");
      LOG("gmkhedissf", pDEBUG) << "gMassScheme: " << gMassScheme;
    }
    else {
      LOG("gmkhedissf", pFATAL) << "Unspecified Mass Scheme - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //Number of bins and boundaries in log10x axis
  if(parser.OptionExists("x")){
    string sx = parser.ArgAsString("x");
    // split the comma separated list
    vector<string> vsx = utils::str::Split(sx, ",");
    assert(vsx.size() == 2);
    gNX       = atoi(vsx[0].c_str());
    gxGRIDmin = atof(vsx[1].c_str());
    LOG("gmkhedissf", pDEBUG) << "gNX: " << gNX;
    LOG("gmkhedissf", pDEBUG) << "gxGRIDmin: " << gxGRIDmin;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified x GRID - Exiting";
    PrintSyntax();
    exit(1);
  }

  //Number of bins and boundaries in log10Q2 axis
  if(parser.OptionExists("q2")){
    string sq2 = parser.ArgAsString("q2");
    // split the comma separated list
    vector<string> vsq2 = utils::str::Split(sq2, ",");
    assert(vsq2.size() == 3);
    gNQ2       = atoi(vsq2[0].c_str());
    gQ2GRIDmin = atof(vsq2[1].c_str());
    gQ2GRIDmax = atof(vsq2[2].c_str());
    LOG("gmkhedissf", pDEBUG) << "gNQ2: " << gNQ2;
    LOG("gmkhedissf", pDEBUG) << "gQ2GRIDmin: " << gQ2GRIDmin;
    LOG("gmkhedissf", pDEBUG) << "gQ2GRIDmax: " << gQ2GRIDmax;
  }
  else {
    LOG("gmkhedissf", pFATAL) << "Unspecified Q2 GRID - Exiting";
    PrintSyntax();
    exit(1);
  }

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmkhedissf", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gmkhedissf [-h]"
    << "\n                  --outdir output_sf_dir_name"
    << "\n                  --lhapdf lhapdf_name"
    << "\n                  --mz mass_boson_z"
    << "\n                  --mw mass_boson_w"
    << "\n                  --rho weinberg_angle_em_correction"
    << "\n                  [--sin2thw sin_weinberg_angle]"
    << "\n                  --ckm ckm_matrix"
    << "\n                  --qrkthrs quark_prod_threshold"
    << "\n                  [--nlo]"
    << "\n                  [--scheme nlo_mass_scheme]"
    << "\n                  --x x_log10_grid"
    << "\n                  --q2 q2_log10_grid"
    << "\n";
}

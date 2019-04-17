#include "Physics/HEDIS/XSection/HEDISPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

#include <TMath.h>
#include <TSystem.h>


#include <fstream>

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
HEDISPXSec::HEDISPXSec() :
XSecAlgorithmI("genie::HEDISPXSec")
{

}
//____________________________________________________________________________
HEDISPXSec::HEDISPXSec(string config) :
XSecAlgorithmI("genie::HEDISPXSec", config)
{

}
//____________________________________________________________________________
HEDISPXSec::~HEDISPXSec()
{

}
//____________________________________________________________________________
double HEDISPXSec::XSec(
     const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidKinematics (interaction) ) return 0.;

  HEDISFormFactors * formfactors = HEDISFormFactors::Instance(fLHAPDFmember,fNLO,fScheme,fQrkThres,fNX,fXmin,fNQ2,fQ2min,fQ2max,fCKM,fMassZ,fMassW,fSin2ThW);

  const Kinematics   & kinematics = interaction -> Kine();
  double x     = kinematics.x();
  double Q2    = kinematics.Q2();

  if ( x<fXmin   ) return 0.;  
  if ( Q2<fQ2min ) return 0.;
  if ( Q2>fQ2max ) return 0.;

  const InitialState & init_state = interaction -> InitState();

  HEDISChannel_t ch = interaction->ExclTag().HEDISChannel();
  double E     = init_state.ProbeE(kRfLab);
  double Mnuc  = init_state.Tgt().HitNucMass();
  double Mlep2 = TMath::Power(interaction->FSPrimLepton()->Mass(),2);
  double y     = kinematics.y();
  FF_xQ2 ff = formfactors->EvalFFQrkLO( ch, x, Q2 );

  double xsec = ds_dxdy( ff, E, Mnuc, Mlep2, x, y );

#ifdef __GENIE_APFEL_ENABLED__
  if (fNLO && xsec>0.) {
    HEDISNuclChannel_t nuclch  = HEDISChannel::HEDISNuclChannel(ch);
    FF_xQ2 fflo  = formfactors->EvalNuclFFLO(nuclch,x,Q2);
    FF_xQ2 ffnlo = formfactors->EvalNuclFFNLO(nuclch,x,Q2);
    double lo  = ds_dxdy( fflo, E, Mnuc, Mlep2, x, y );
    if (lo>0.) {
      double nlo = ds_dxdy( ffnlo, E, Mnuc, Mlep2, x, y );
      xsec *= nlo / lo;
    }
  }
#endif

  // Compute the front factor
  double propagator = 0;
  if (interaction -> ProcInfo().IsWeakCC()) propagator = TMath::Power( fMassW*fMassW/(Q2+fMassW*fMassW), 2);
  else                                      propagator = TMath::Power( fMassZ*fMassZ/(Q2+fMassZ*fMassZ)/(1.-fRho), 2);

  xsec *= kGF2 * Mnuc* E / kPi * propagator;

  LOG("HEDISPXSec", pINFO) << "d2xsec/dxdy[FreeN] (E= " << E << ", x= " << x  << ", y= " << y << ", Q2= " << Q2 << ") = " << xsec;

  // The algorithm computes d^2xsec/dxdy. Check whether variable tranformation is needed
  if( kps!=kPSxyfE ) xsec *= utils::kinematics::Jacobian(interaction,kPSxyfE,kps);

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must have been included in the structure functions)
  int NNucl = (pdg::IsProton(init_state.Tgt().HitNucPdg())) ? init_state.Tgt().Z() : init_state.Tgt().N(); 
  xsec *= NNucl; 

  return xsec;

}
//____________________________________________________________________________
double HEDISPXSec::ds_dxdy(FF_xQ2 ff, double e, double mt, double ml2, double x, double y ) const
{

    //double term1 = y * ( x*y + ml2/2/e/mt );
    //double term2 = ( 1 - y - mt*x*y/2/e - ml2/4/e/e );
    //double term3 = sgn_t3 * (x*y*(1-y/2) - y*ml2/4/mt/e);
    //double term4 = x*y*ml2/2/mt/e + ml2*ml2/4/mt/mt/e/e;
    //double term5 = -1.*ml2/2/mt/e;

    double term1 = y * ( x*y );
    double term2 = ( 1 - y );
    double term3 = ( x*y*(1-y/2) );

    LOG("HEDISPXSec", pDEBUG) << ff.F1 << "  " << ff.F2 << "  " << ff.F3;
    LOG("HEDISPXSec", pDEBUG) << term1*ff.F1 + term2*ff.F2 + term3*ff.F3;

    return fmax( term1*ff.F1 + term2*ff.F2 + term3*ff.F3 , 0.);

}
//____________________________________________________________________________
double HEDISPXSec::Integral(const Interaction * interaction) const
{

  return fXSecIntegrator->Integrate(this,interaction);

}
//____________________________________________________________________________
bool HEDISPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsHEDIS()) return false;

  const InitialState & init_state = interaction -> InitState();
  int probe_pdg = init_state.ProbePdg();
  if(!pdg::IsLepton(probe_pdg)) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;

  return true;
}
//____________________________________________________________________________
void HEDISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  GetParam("CKM-Vud", fCKM[0]);
  GetParam("CKM-Vus", fCKM[1]);
  GetParam("CKM-Vub", fCKM[2]);
  GetParam("CKM-Vcd", fCKM[3]);
  GetParam("CKM-Vcs", fCKM[4]);
  GetParam("CKM-Vcb", fCKM[5]);
  GetParam("CKM-Vtd", fCKM[6]);
  GetParam("CKM-Vts", fCKM[7]);
  GetParam("CKM-Vtb", fCKM[8]);

  for (int i=0; i<9; i++) LOG("HEDISPXSec", pERROR) << "CKM" << i << " = " << fCKM[i];

  GetParamDef("LHAPDF-member", fLHAPDFmember, string("") ) ;
  GetParamDef("Is-NLO", fNLO, false ) ;
  GetParamDef("Scheme", fScheme, string("") ) ;
  GetParamDef("Quark-Threshold", fQrkThres, 0 ) ;
  GetParamDef("MassW", fMassW, kMw ) ;
  GetParamDef("MassZ", fMassZ, kMz ) ;
  GetParam("WeinbergAngle", fSin2ThW) ;
  GetParamDef("Rho", fRho, 0. ) ;

  if (fSin2ThW==0.) fSin2ThW = 1.-fMassW*fMassW/fMassZ/fMassZ/(1+fRho) ;
  else              fSin2ThW = TMath::Power(TMath::Sin(fSin2ThW),2);

  LOG("HEDISPXSec", pERROR) << "fSin2ThW = " << fSin2ThW;


  GetParamDef("NGridX",  fNX,  100 );
  GetParamDef("NGridQ2", fNQ2, 100 );
  GetParamDef("XGrid-Min",    fXmin, 1e-10 );
  GetParamDef("Q2Grid-Min",  fQ2min,   0.1 );
  GetParamDef("Q2Grid-Max",  fQ2max,  1e10 );

}

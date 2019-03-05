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

  HEDISFormFactors * formfactors = HEDISFormFactors::Instance(fNLO,fLHAPDFmember,fNX,fNQ2,fCKM);

  const Kinematics   & kinematics = interaction -> Kine();
  double x     = kinematics.x();
  double Q2    = kinematics.Q2();

  if( !formfactors->ValidKinematicsPDF(x,Q2) ) return 0.;  

  const InitialState & init_state = interaction -> InitState();

  HEDISChannel_t ch = interaction->ExclTag().HEDISChannel();
  double E     = init_state.ProbeE(kRfLab);
  double Mnuc  = init_state.Tgt().HitNucMass();
  double Mlep2 = TMath::Power(interaction->FSPrimLepton()->Mass(),2);
  double y     = kinematics.y();
  FF_xQ2 ff = formfactors->EvalFFQrkLO( ch, x, Q2 );

  double xsec = ds_dxdy( ff, E, Mnuc, Mlep2, x, y );

  if (fNLO && xsec>0) {
    HEDISNuclChannel_t nuclch  = HEDISChannel::HEDISNuclChannel(ch);
    FF_xQ2 fflo  = formfactors->EvalNuclFFLO(nuclch,x,Q2);
    FF_xQ2 ffnlo = formfactors->EvalNuclFFNLO(nuclch,x,Q2);
    double lo  = ds_dxdy( fflo, E, Mnuc, Mlep2, x, y );
    double nlo = ds_dxdy( ffnlo, E, Mnuc, Mlep2, x, y );
    xsec *= nlo / lo;
  }

  // Compute the front factor
  double front_factor = 0;
  if (interaction -> ProcInfo().IsWeakCC()) front_factor = kGF2 * kMw2 * kMw2 / TMath::Power((Q2 + kMw2), 2) * Mnuc* E / kPi;
  else                                      front_factor = kGF2 * kMz2 * kMz2 / TMath::Power((Q2 + kMz2), 2) * Mnuc* E / kPi;

  xsec *= front_factor;

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

    double term1 = y * ( x*y );                //genhen
    double term2 = ( 1 - y );                  //genhen
    double term3 = ( x*y*(1-y/2) );   //genhen

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

  GetParamDef("Is-NLO", fNLO, false ) ;

  GetParamDef("CKM-Vud", fCKM[0], 1. );
  GetParamDef("CKM-Vus", fCKM[1], 1. );
  GetParamDef("CKM-Vub", fCKM[2], 1. );
  GetParamDef("CKM-Vcd", fCKM[3], 1. );
  GetParamDef("CKM-Vcs", fCKM[4], 1. );
  GetParamDef("CKM-Vcb", fCKM[5], 1. );
  GetParamDef("CKM-Vtd", fCKM[6], 1. );
  GetParamDef("CKM-Vts", fCKM[7], 1. );
  GetParamDef("CKM-Vtb", fCKM[8], 1. );
  for (int i=0; i<9; i++) fCKM[i] = TMath::Power(fCKM[i],2);

  GetParamDef("LHAPDF-member", fLHAPDFmember, string("") ) ;

  GetParamDef("NX",  fNX,  100 );
  GetParamDef("NQ2", fNQ2, 100 );

}

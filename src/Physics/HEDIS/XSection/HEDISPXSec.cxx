//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/XSection/HEDISPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

#include <TMath.h>

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

  // Load SF tables 
  HEDISStrucFunc * sf_tbl = HEDISStrucFunc::Instance(fSFinfo);

  // W limits are computed using kinematics assumption.
  // The lower limit is tuneable because hadronization might have issues at low W (as in PYTHIA6).
  // To be consistent the cross section must be computed in the W range where the events are generated.
  const Kinematics  & kinematics = interaction -> Kine();
  const KPhaseSpace & ps         = interaction -> PhaseSpace();
  double W = kinematics.W();
  Range1D_t Wl  = ps.WLim();
  Wl.min = TMath::Max(Wl.min,fWmin);
  if      ( W<Wl.min ) return 0.;
  else if ( W>Wl.max ) return 0.;

  const InitialState & init_state = interaction -> InitState();

  HEDISQrkChannel_t ch = interaction->ExclTag().HEDISQrkChannel();
  double y     = kinematics.y();
  double Q2    = kinematics.Q2();
  double x     = kinematics.x();

  // Get F1,F2,F3 for particular quark channel and compute differential xsec
  SF_xQ2 sf = sf_tbl->EvalQrkSFLO( ch, x, Q2 );
  double xsec = ds_dxdy( sf, x, y );

  // If NLO is enable we compute sigma_NLO/sigma_LO. Then the quark xsec 
  // is multiplied by this ratio.
  // This is done because at NLO we can only compute the nucleon xsec. But
  // for the hadronization we need the different quark contributions.
  // This could be avoid if a NLO parton showering is introduced.
  if (fSFinfo.IsNLO && xsec>0.) {
    HEDISNucChannel_t nucch  = HEDISChannel::HEDISNucChannel(ch);
    SF_xQ2 sflo  = sf_tbl->EvalNucSFLO(nucch,x,Q2);
    SF_xQ2 sfnlo = sf_tbl->EvalNucSFNLO(nucch,x,Q2);
    double lo  = ds_dxdy( sflo, x, y );
    if (lo>0.) {
      double nlo = ds_dxdy( sfnlo, x, y );
      xsec *= nlo / lo;
    }
  }

  // Compute the front factor
  double propagator = 0;
  if (interaction -> ProcInfo().IsWeakCC()) propagator = TMath::Power( fSFinfo.MassW*fSFinfo.MassW/(Q2+fSFinfo.MassW*fSFinfo.MassW), 2);
  else                                      propagator = TMath::Power( fSFinfo.MassZ*fSFinfo.MassZ/(Q2+fSFinfo.MassZ*fSFinfo.MassZ)/(1.-fSFinfo.Rho), 2);

  xsec *= kGF2/(2*kPi*x) * propagator;

  LOG("HEDISPXSec", pINFO) << "d2xsec/dxdy[FreeN] (x= " << x  << ", y= " << y << ", Q2= " << Q2 << ") = " << xsec;

  // The algorithm computes d^2xsec/dxdy. Check whether variable tranformation is needed
  if( kps!=kPSxQ2fE ) xsec *= utils::kinematics::Jacobian(interaction,kPSxQ2fE,kps);

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must have been included in the structure functions)
  int NNucl = (pdg::IsProton(init_state.Tgt().HitNucPdg())) ? init_state.Tgt().Z() : init_state.Tgt().N(); 
  xsec *= NNucl; 

  return xsec;

}
//____________________________________________________________________________
double HEDISPXSec::ds_dxdy(SF_xQ2 sf, double x, double y ) const
{

    // We neglect F4 and F5 and higher order terms.
    double term1 = y * ( x*y );
    double term2 = ( 1 - y );
    double term3 = ( x*y*(1-y/2) );

    LOG("HEDISPXSec", pDEBUG) << sf.F1 << "  " << sf.F2 << "  " << sf.F3;
    LOG("HEDISPXSec", pDEBUG) << term1*sf.F1 + term2*sf.F2 + term3*sf.F3;

    return fmax( term1*sf.F1 + term2*sf.F2 + term3*sf.F3 , 0.);

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

  // Minimum value of W (typically driven by hadronization limitation)
  GetParam("Wmin",     fWmin);
  // Information about Structure Functions
  GetParam("HEDIS-LHAPDF-set",      fSFinfo.LHAPDFset    );
  GetParam("HEDIS-LHAPDF-member",   fSFinfo.LHAPDFmember );
  GetParam("HEDIS-Is-NLO",          fSFinfo.IsNLO        );
  GetParam("HEDIS-Scheme",          fSFinfo.Scheme       );
  GetParam("HEDIS-Quark-Threshold", fSFinfo.QrkThrs      );
  GetParam("HEDIS-NGridX",          fSFinfo.NGridX       );
  GetParam("HEDIS-NGridQ2",         fSFinfo.NGridQ2      );
  GetParam("HEDIS-XGrid-Min",       fSFinfo.XGridMin     );
  GetParam("HEDIS-Q2Grid-Min",      fSFinfo.Q2GridMin    );
  GetParam("HEDIS-Q2Grid-Max",      fSFinfo.Q2GridMax    );
  GetParam("HEDIS-MassW",           fSFinfo.MassW        );
  GetParam("HEDIS-MassZ",           fSFinfo.MassZ        );
  GetParam("HEDIS-Rho",             fSFinfo.Rho          );
  GetParam("HEDIS-Sin2ThW",         fSFinfo.Sin2ThW      );
  GetParam("HEDIS-CKM-Vud",         fSFinfo.Vud          );
  GetParam("HEDIS-CKM-Vus",         fSFinfo.Vus          );
  GetParam("HEDIS-CKM-Vub",         fSFinfo.Vub          );
  GetParam("HEDIS-CKM-Vcd",         fSFinfo.Vcd          );
  GetParam("HEDIS-CKM-Vcs",         fSFinfo.Vcs          );
  GetParam("HEDIS-CKM-Vcb",         fSFinfo.Vcb          );
  GetParam("HEDIS-CKM-Vtd",         fSFinfo.Vtd          );
  GetParam("HEDIS-CKM-Vts",         fSFinfo.Vts          );
  GetParam("HEDIS-CKM-Vtb",         fSFinfo.Vtb          );

}
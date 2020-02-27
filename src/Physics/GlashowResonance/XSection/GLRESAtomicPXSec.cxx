//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ November 08, 2019 - Alfonso Garcia
   Modified to generate the kinematics of outgoing lepton properly.
   Phys. Rev. D 22, 2122 – Published 1 November 1980

*/
//____________________________________________________________________________

#include "Physics/GlashowResonance/XSection/GLRESAtomicPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

#include <TMath.h>

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GLRESAtomicPXSec::GLRESAtomicPXSec() :
XSecAlgorithmI("genie::GLRESAtomicPXSec")
{
  born = new GLRESBornPXSec();
}
//____________________________________________________________________________
GLRESAtomicPXSec::GLRESAtomicPXSec(string config) :
XSecAlgorithmI("genie::GLRESAtomicPXSec", config)
{

}
//____________________________________________________________________________
GLRESAtomicPXSec::~GLRESAtomicPXSec()
{

}
//____________________________________________________________________________
double GLRESAtomicPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const XclsTag &      xclstag    = interaction -> ExclTag();

  int probepdg = init_state.ProbePdg();
  int lout     = xclstag.FinalLeptonPdg();
  double mlout = interaction->FSPrimLepton()->Mass();
  double mlin  = 0.;
  if      (pdg::IsNuE  (TMath::Abs(probepdg))) mlin = kElectronMass;
  else if (pdg::IsNuMu (TMath::Abs(probepdg))) mlin = kMuonMass;
  else if (pdg::IsNuTau(TMath::Abs(probepdg))) mlin = kTauMass;

  double E = init_state.ProbeE(kRfLab);
  double s = 2 * mlin * E + mlin*mlin;

  double x1 = kinematics.GetKV(kKVGLx1);
  double x2 = kinematics.GetKV(kKVGLx2);
  double t  = born->GetT(0.,mlin,mlout,0.,s,x1);

  double s_r = s;
  double t_r = t;
  double pdf_soft = 1.;

  //nlo correction
  if (fIsNLO) {
    double zeta = born->GetReAlpha()/kPi*(2.*TMath::Log(TMath::Sqrt(-t)/kElectronMass)-1.);
    double omx  = TMath::Power(x2, 1./zeta );
    pdf_soft = TMath::Exp(zeta*(3./4.-TMath::EulerGamma()))/TMath::Gamma(1.+zeta) + omx*(omx-2.)/2./x2;
    s_r *= (1. - omx);
    t_r *= (1. - omx);
  }

  double xsec = kPi/4./(s-mlin*mlin) * pdf_soft ;
  
  if ( pdg::IsPion(lout) ) {
    if ( TMath::Sqrt(s_r)<fWmin ) return 0.; // The W limit is because hadronization might have issues at low W (as in PYTHIA6).
    xsec *= 64.41/10.63;    
  }

  xsec *= born->PXSecLepton(s_r,t_r,probepdg,lout,mlin*mlin,mlout*mlout);
   
  //----- If requested return the free electron xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeElectron) ) return xsec;
   
  //----- Scale for the number of scattering centers at the target
  int Ne = init_state.Tgt().Z(); // num of scattering centers
  xsec *= Ne;

  if(kps!=kPSGLx1x2fE) {
      LOG("GLRESAtomicPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSGLx1x2fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  LOG("GLRESAtomicPXSec", pINFO) << "dxsec/dy (E= " << E << ", x1= " << x1 << ", x2=" << x2 << ") = " << xsec;

  return xsec;

}
//____________________________________________________________________________
double GLRESAtomicPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool GLRESAtomicPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsGlashowResonanceAtomic()) return false;
  if(!proc_info.IsWeakCC()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(!pdg::IsAntiNuE(init_state.ProbePdg())) return false;

  if(init_state.Tgt().HitNucIsSet()) return false;
 
  return true;
}
//____________________________________________________________________________
void GLRESAtomicPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESAtomicPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESAtomicPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  GetParam( "Wmin", fWmin ) ;
  GetParam( "GLRESAtomic-Is-NLO", fIsNLO );

}

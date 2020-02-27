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

#include "Physics/GlashowResonance/XSection/GLRESCoherentPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

const double a  = 0.523 * (units::fermi);
const double r0 = 1.126 * (units::fermi);

//____________________________________________________________________________
GLRESCoherentPXSec::GLRESCoherentPXSec() :
XSecAlgorithmI("genie::GLRESCoherentPXSec")
{
  born = new GLRESBornPXSec();
}
//____________________________________________________________________________
GLRESCoherentPXSec::GLRESCoherentPXSec(string config) :
XSecAlgorithmI("genie::GLRESCoherentPXSec", config)
{

}
//____________________________________________________________________________
GLRESCoherentPXSec::~GLRESCoherentPXSec()
{

}
//____________________________________________________________________________
double GLRESCoherentPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();

  int probepdg = init_state.ProbePdg();
  double E     = init_state.ProbeE(kRfLab);

  double mlout = 0;
  if      (pdg::IsNuE  (TMath::Abs(probepdg))) mlout = kElectronMass;
  else if (pdg::IsNuMu (TMath::Abs(probepdg))) mlout = kMuonMass;
  else if (pdg::IsNuTau(TMath::Abs(probepdg))) mlout = kTauMass;

  int A        = init_state.Tgt().A(); 
  int Z        = init_state.Tgt().Z(); 
  double Mtgt = Z*kProtonMass + (A-Z)*kNeutronMass;

  double x1 = kinematics.GetKV(kKVGLx1);
  double x2 = kinematics.GetKV(kKVGLx2);
  double x3 = kinematics.GetKV(kKVGLx3);

  double ME,dps2,dP;
  if (fIsEPA) {

    double s = 2.*Mtgt*E + Mtgt*Mtgt;
    double y0 = TMath::Power(kMw+mlout,2)/s;
    double y = TMath::Exp( TMath::Log(y0)+(TMath::Log(1)-TMath::Log(y0))*x2 );

    double Q2min = TMath::Power(y*Mtgt,2)/(1.-y);
    double Q2max = s-TMath::Power(Mtgt,2);
    double Q2    = TMath::Exp( TMath::Log(Q2min)+(TMath::Log(Q2max)-TMath::Log(Q2min))*x3 );
   
    double s_r = s * y;
    double t_r = born->GetT(0.,0.,kMw,mlout,s_r,x1); 
   
    ME = born->PXSecPhoton(s_r,t_r,mlout*mlout);
    dps2 = 1./64./kPi2/s_r * TMath::Sqrt( born->Lambda(1.,kMw2/s_r,mlout*mlout/s_r) ) * (TMath::Log(Q2max)-TMath::Log(Q2min)) * (TMath::Log(1)-TMath::Log(y0));
    dP = born->GetReAlpha()*TMath::Power(Z,2) * (2.-2.*y+y*y) * F2_Q( TMath::Sqrt(Q2), r0*TMath::Power(A,1./3.) ) * (Q2-Q2min)/Q2;

  }
  else {

    double mL = mlout+kMw;
    double Delta = TMath::Sqrt( TMath::Power(2.*E*Mtgt-mL*mL,2)-4.*TMath::Power(Mtgt*mL,2) );
    
    double s12_min = E/(2.*E+Mtgt)*(mL*mL+2.*E*Mtgt-Delta);
    double s12_max = E/(2.*E+Mtgt)*(mL*mL+2.*E*Mtgt+Delta);
    double s12 = TMath::Exp( TMath::Log(s12_min)+(TMath::Log(s12_max)-TMath::Log(s12_min))*x2);

    double Q2_min = TMath::Power(s12,2)*Mtgt/2./E/(2.*E*Mtgt-s12);
    double Q2_max = s12 - mL*mL;
    double Q2 = TMath::Exp( TMath::Log(Q2_min) + (TMath::Log(Q2_max)-TMath::Log(Q2_min))*x3 );

    double s = s12 - Q2;
    double s13 = s12/2.*((1.+kMw2/s-mlout*mlout/s)-TMath::Sqrt(born->Lambda(1.,kMw2/s,mlout*mlout/s))*x1);

    double ME_T = born->PXSecPhoton_T(s12,s13,Q2,mlout*mlout) * (1.-s12/2./E/Mtgt-TMath::Power(s12/2./E,2)/Q2);
    double ME_L = born->PXSecPhoton_L(s12,s13,Q2,mlout*mlout) * TMath::Power(1.-s12/4./E/Mtgt,2);
    
    ME   = ME_T+ME_L;
    dps2 = 1/32./kPi2/s12 * TMath::Sqrt( born->Lambda(1.,kMw2/s,mlout*mlout/s) ) * (TMath::Log(Q2_max)-TMath::Log(Q2_min)) * (TMath::Log(s12_max)-TMath::Log(s12_min));    
    dP   = born->GetReAlpha()*TMath::Power(Z,2)*F2_Q( TMath::Sqrt(Q2), r0*TMath::Power(A,1./3.) );

  }

  double xsec = ME * dps2 * dP;
   
  if(kps!=kPSGLx1x2x3fE) {
      LOG("GLRESCoherentPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSGLx1x2x3fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  LOG("GLRESCoherentPXSec", pINFO) << "dxsec/dy (E= " << E << ", x1= " << x1 << ", x2=" << x2 << ", x3=" << x3 << ") = " << xsec;

  return xsec;
}
//____________________________________________________________________________
double GLRESCoherentPXSec::F2_Q(double Q, double r0) const
{
  // Analytic Woods-Saxon, A.3 of https://arxiv.org/pdf/1807.10973.pdf
  double FF = 0.0;
  double coth = 1./TMath::TanH(kPi*Q*a);
  FF = 3.*kPi*a/(TMath::Power(r0,2)+TMath::Power(kPi*a,2)) / (Q*r0* TMath::SinH(kPi*Q*a));
  FF *= (kPi*a*coth*TMath::Sin(Q*r0) - r0 * TMath::Cos(Q*r0));        
  return TMath::Power(FF,2);
}
//____________________________________________________________________________
double GLRESCoherentPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool GLRESCoherentPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsGlashowResonanceCoh()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(!pdg::IsLepton(init_state.ProbePdg())) return false;

  if(init_state.Tgt().HitNucIsSet()) return false;
  
  return true;
}
//____________________________________________________________________________


//____________________________________________________________________________
void GLRESCoherentPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESCoherentPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESCoherentPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  GetParam( "GLRESCoh-Is-EPA", fIsEPA );

}

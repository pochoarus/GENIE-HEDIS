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

#include "Physics/GlashowResonance/XSection/GLRESInelasticStrucFunc.h"
#include "Physics/GlashowResonance/XSection/GLRESInelasticPXSec.h"
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

//____________________________________________________________________________
GLRESInelasticPXSec::GLRESInelasticPXSec() :
XSecAlgorithmI("genie::GLRESInelasticPXSec")
{
  born = new GLRESBornPXSec();
}
//____________________________________________________________________________
GLRESInelasticPXSec::GLRESInelasticPXSec(string config) :
XSecAlgorithmI("genie::GLRESInelasticPXSec", config)
{

}
//____________________________________________________________________________
GLRESInelasticPXSec::~GLRESInelasticPXSec()
{

}
//____________________________________________________________________________
double GLRESInelasticPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  // Load SF tables 
  GLRESStrucFunc * sf_tbl = GLRESStrucFunc::Instance();

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const XclsTag &      xclstag    = interaction -> ExclTag();

  int probepdg = init_state.ProbePdg();
  int lout     = xclstag.FinalLeptonPdg();
  double mlout = interaction->FSPrimLepton()->Mass();

  int tgtpdg        = init_state.Tgt().HitNucPdg();
  double Mnuc = init_state.Tgt().HitNucMass();

  double E = init_state.ProbeE(kRfLab);
  double s = (2 * Mnuc * E + Mnuc*Mnuc);

  double x1 = kinematics.GetKV(kKVGLx1);
  double x2 = kinematics.GetKV(kKVGLx2);

  double xmin = fQ2PDFmin/2./E/Mnuc;
  double x = TMath::Exp( TMath::Log(xmin) + (TMath::Log(1.0)-TMath::Log(xmin))*x2 );

  if (x<fxPDFmin) return 0.;

  double s_r = x*s;
  double t_r = born->GetT(0.,0.,mlout,0.,s_r,x1);

  double xsec = kPi/4./(s_r-Mnuc*Mnuc) * sf_tbl->EvalSF(tgtpdg,probepdg,x) * (TMath::Log(1.0)-TMath::Log(xmin)) ;
  
  if ( pdg::IsPion(lout) ) {
    if ( TMath::Sqrt(s_r)<fWmin ) return 0.; // The W limit is because hadronization might have issues at low W (as in PYTHIA6).
    xsec *= 64.41/10.63;    
  }

  xsec *= born->PXSecLepton(s_r,t_r,probepdg,lout,0.,mlout*mlout);
   
  if(kps!=kPSGLx1x2fE) {
      LOG("GLRESInelasticPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSGLx1x2fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must have been included in the structure functions)
  int NNucl = (pdg::IsProton(tgtpdg)) ? init_state.Tgt().Z() : init_state.Tgt().N(); 
  xsec *= NNucl; 

  LOG("GLRESInelasticPXSec", pINFO) << "dxsec/dy (E= " << E << ", x1= " << x1 << ", x2=" << x2 << ") = " << xsec;

  return xsec;

}
//____________________________________________________________________________
double GLRESInelasticPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool GLRESInelasticPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsGlashowResonanceInel()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(!pdg::IsLepton(init_state.ProbePdg())) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;
  
  return true;
}
//____________________________________________________________________________


//____________________________________________________________________________
void GLRESInelasticPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESInelasticPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESInelasticPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  GetParam( "Wmin", fWmin ) ;

  GetParam("HEDIS-Q2Grid-Min", fQ2PDFmin );
  GetParam("HEDIS-XGrid-Min", fxPDFmin );

}

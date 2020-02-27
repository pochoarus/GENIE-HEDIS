//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/Integrator.h>

#include "Physics/GlashowResonance/XSection/GLRESXSec.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
GLRESXSec::GLRESXSec() :
XSecIntegratorI("genie::GLRESXSec")
{

}
//____________________________________________________________________________
GLRESXSec::GLRESXSec(string config) :
XSecIntegratorI("genie::GLRESXSec", config)
{

}
//____________________________________________________________________________
GLRESXSec::~GLRESXSec()
{

}
//____________________________________________________________________________
double GLRESXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("GLRESXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }


  const ProcessInfo &  proc_info  = in->ProcInfo();
  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  
  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if( !xsl->IsEmpty() ) {

    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();

    int NNucl = 0;
    if ( proc_info.IsGlashowResonanceAtomic() ) {
      NNucl = init_state.Tgt().Z();
      target->SetId(kPdgTgtFreeP);
    }
    else if ( proc_info.IsGlashowResonanceInel() ) {
      int nucpdgc = init_state.Tgt().HitNucPdg();
      if (pdg::IsProton(nucpdgc)) {
        NNucl = init_state.Tgt().Z();
        target->SetId(kPdgTgtFreeP);
      }
      else {
        NNucl = init_state.Tgt().N();
        target->SetId(kPdgTgtFreeN);
      }
    }

    if(xsl->SplineExists(model,interaction)) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("GLRESXSec", pINFO)  << "From XSecSplineList: XSec["<<init_state.ProbePdg()<<","<<proc_info.ScatteringTypeAsString()<<","<<interaction->ExclTag().FinalLeptonPdg()<< ",free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if( !interaction->TestBit(kIAssumeFreeNucleon) ) { 
        xsec *= NNucl; 
        LOG("GLRESXSec", pINFO)  << "XSec["<<init_state.ProbePdg()<<","<<proc_info.ScatteringTypeAsString()<<","<<interaction->ExclTag().FinalLeptonPdg()<< "] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }
   
  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);

  double xsec = 0.;
  double kine_min[2] = { -1.,  0. }; 
  double kine_max[2] = {  1.,  1. }; 
  ROOT::Math::IBaseFunctionMultiDim * func = new utils::gsl::d2Xsec_GLRES(model, interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  ROOT::Math::IntegratorMultiDim ig(*func,ig_type,1,fGSLRelTol,fGSLMaxEval);
  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
  delete func;

  LOG("GLRESXSec", pDEBUG) << "*** XSec["<<init_state.ProbePdg()<<","<<proc_info.ScatteringTypeAsString()<<","<<interaction->ExclTag().FinalLeptonPdg()<< "] (E=" << interaction->InitState().ProbeE(kRfLab) << ") = " << xsec * (1E+38/units::cm2) << " x 1E-38 cm^2";

  delete interaction;
  return xsec;
}
//____________________________________________________________________________
void GLRESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef("gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef("gsl-relative-tolerance", fGSLRelTol, 1E-2 ) ;
  int max_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  fGSLMaxEval  = (unsigned int) max_eval ;

}
//____________________________________________________________________________

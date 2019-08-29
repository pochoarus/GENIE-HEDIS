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

#include "Physics/HEDIS/XSection/HEDISXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Numerical/GSLUtils.h"

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"
#include <TSystem.h>

#include <fstream>

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
HEDISXSec::HEDISXSec() :
XSecIntegratorI("genie::HEDISXSec")
{

}
//____________________________________________________________________________
HEDISXSec::HEDISXSec(string config) :
XSecIntegratorI("genie::HEDISXSec", config)
{

}
//____________________________________________________________________________
HEDISXSec::~HEDISXSec()
{

}
//____________________________________________________________________________
double HEDISXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{

  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("DISXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
 
  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty() ) {

    int nucpdgc = init_state.Tgt().HitNucPdg();
    int NNucl   = (pdg::IsProton(nucpdgc)) ? init_state.Tgt().Z() : init_state.Tgt().N();

    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();
    if(pdg::IsProton(nucpdgc)) { target->SetId(kPdgTgtFreeP); }
    else                       { target->SetId(kPdgTgtFreeN); }
    if(xsl->SplineExists(model,interaction)) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("HEDISXSec", pINFO) << "From XSecSplineList: XSec[HEDIS,free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if( !interaction->TestBit(kIAssumeFreeNucleon) ) { 
          xsec *= NNucl; 
          LOG("HEDISXSec", pINFO)  << "XSec[HEDIS] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }

  Range1D_t Wl  = kps.WLim();
  LOG("HEDISXSec", pDEBUG) << "W integration range = [" << Wl.min << ", " << Wl.max << "]";
  Range1D_t Q2l = kps.Q2Lim();
  LOG("HEDISXSec", pDEBUG) << "Q2 integration range = [" << Q2l.min << ", " << Q2l.max << "]";

  bool phsp_ok = 
      (Q2l.min >= 0. && Q2l.max >= 0. && Q2l.max >= Q2l.min &&
        Wl.min >= 0. &&  Wl.max >= 0. &&  Wl.max >=  Wl.min);

  if (!phsp_ok) return 0.;

  // Just go ahead and integrate the input differential cross section for the 
  // specified interaction.
  //
  double xsec = 0.;
  double d2xsec_max = 0.;

  Interaction * interaction = new Interaction(*in);
  
  // First we compute the integral using the grid.
  // This is always needed because we need to get the Max Xsec.
  // TODO: This might be change in the future if Max Xsec are computed with minimization.
  for ( const auto& y : fVy ) {
    interaction->KinePtr()->Sety( y );
    double dxsecdy = 0.;
    for ( const auto& x : fVx ) {
      interaction->KinePtr()->Setx(x);
      utils::kinematics::UpdateWQ2FromXY(interaction);
      double d2xsecdxdy = model->XSec(interaction, kPSxyfE);     
      dxsecdy += d2xsecdxdy * x;
      double d2xsec = d2xsecdxdy * x * y;
      xsec += d2xsec;
      if ( d2xsec > d2xsec_max ) d2xsec_max = d2xsec;
    }
  }
  d2xsec_max *= TMath::Log(10.) * TMath::Log(10.);
  xsec       *= TMath::Log(10.) * TMath::Log(10.) * fdlogx * fdlogy;
  
  // If a GSL option has been chosen, then the total xsec is recomptued
  if (fGSLIntgType!="") {
       ROOT::Math::IBaseFunctionMultiDim * func = new utils::gsl::d2XSec_dWdQ2_E(model, interaction);
       ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
       double abstol = 1; //We mostly care about relative tolerance.
       ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
       double kine_min[2] = { Wl.min, Q2l.min };
       double kine_max[2] = { Wl.max, Q2l.max };
       xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
       delete func;
  }

  delete interaction;

  LOG("HEDISXSec", pINFO)  << "XSec[HEDIS] (E = " << Ev << " GeV) = " << xsec * (1E+38/units::cm2) << " x 1E-38 cm^2  (max = " << d2xsec_max << ")";

  // File for the Max Xsec is created here and Max Xsec for each neutrino energy.
  // The program checks if the file already exists. If so then it crashes becuase
  // we dont want to overwrite a previously produce max xsec.
  string sabsnupdg = std::to_string( TMath::Abs(init_state.ProbePdg()) );
  HEDISQrkChannel_t chID = in->ExclTag().HEDISQrkChannel();
  string maxxsecfilename = fMaxXsecDirName + "/Flvr" + sabsnupdg + "_" + HEDISChannel::AsString(chID) + ".dat";

  if (!fMaxXsecFileExists) {
    if ( gSystem->AccessPathName( maxxsecfilename.c_str()) ) {
      LOG("HEDISXSec", pINFO) << "Creating Max Xsec file: " << maxxsecfilename;        
      fMaxXsecFileExists = true;
    }
    else {
      LOG("HEDISXSec", pERROR) << "Max Xsec file already exists: " << maxxsecfilename;        
      LOG("HEDISXSec", pERROR) << "Stop to avoid overwriting!!!";        
      assert(0);
    }
  }

  std::ofstream maxxsec_stream(maxxsecfilename.c_str(),std::fstream::app);
  maxxsec_stream << Ev << "  " << d2xsec_max << std::endl;
  maxxsec_stream.close();

  return xsec;

}
//____________________________________________________________________________
void HEDISXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISXSec::LoadConfig(void)
{

  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1E-2 ) ;

  int max_eval, min_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  GetParamDef( "gsl-min-eval", min_eval, 10000 ) ;
  fGSLMaxEval  = (unsigned int) max_eval ;
  fGSLMinEval  = (unsigned int) min_eval ;

  // Binning size (in log10) used to intergrate over x and y
  GetParamDef("DlogY", fdlogy, 0.01 ) ;
  GetParamDef("DlogX", fdlogx, 0.01 ) ;
  std::ostringstream sdlogy,sdlogx; 
  sdlogy << fdlogy;
  sdlogx << fdlogx;

  // The minimum value for x and y is 1e-10. This value is safe for energies below 1e15 GeV.
  // Binning is generated for x and y
  double flogymin = -10;
  double flogxmin = -10;
  for ( int iy=0; iy<TMath::Abs(flogymin)/fdlogy; iy++) fVy.push_back( TMath::Power( 10., flogymin + (iy+0.5)*fdlogy ) );
  for ( int ix=0; ix<TMath::Abs(flogxmin)/fdlogx; ix++) fVx.push_back( TMath::Power( 10., flogxmin + (ix+0.5)*fdlogx ) );

  // Name of the directory where SF tables are saved. We need it here to define a name for the Max Xsec directory.
  // TODO: This might be removed in the future if Max Xsec are included in Splines.
  string SFname;
  GetParamDef("SF-name", SFname, string("") ) ;

  // Name Max Xsec directory. It will be created if it doesnt exist.
  fMaxXsecDirName = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/maxxsec/" + SFname + "_dx" + sdlogx.str() + "dy" + sdlogy.str();
  if ( gSystem->mkdir(fMaxXsecDirName.c_str())==0 ) {
    LOG("HEDISXSec", pINFO) << "Creating Max Xsec directory: " << fMaxXsecDirName;
  }
  else {
    LOG("HEDISXSec", pINFO) << "Max Xsec directory already exist: " << fMaxXsecDirName;
  }


}

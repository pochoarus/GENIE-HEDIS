#include "Physics/HEDIS/XSection/HEDISXSec.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Numerical/Spline.h"

#include <TMath.h>
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
      LOG("DISXSec", pINFO) << "From XSecSplineList: XSec[HEDIS,free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if( !interaction->TestBit(kIAssumeFreeNucleon) ) { 
          xsec *= NNucl; 
          LOG("DISXSec", pINFO)  << "XSec[HEDIS] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }

  Range1D_t Wl  = kps.WLim();
  Range1D_t Q2l = kps.Q2Lim();
  LOG("HEDISXSec", pDEBUG) << "W integration range = [" << Wl.min << ", " << Wl.max << "]";
  LOG("HEDISXSec", pDEBUG) << "Q2 integration range = [" << Q2l.min << ", " << Q2l.max << "]";

  bool phsp_ok = 
      (Q2l.min >= 0. && Q2l.max >= 0. && Q2l.max >= Q2l.min &&
        Wl.min >= 0. &&  Wl.max >= 0. &&  Wl.max >=  Wl.min);

  if (!phsp_ok) return 0.;

  string sabsnupdg = std::to_string( TMath::Abs(init_state.ProbePdg()) );
  HEDISChannel_t chID = in->ExclTag().HEDISChannel();

  string filename = "Flvr" + sabsnupdg + "_" + HEDISChannel::AsString(chID) + ".dat";

  // Just go ahead and integrate the input differential cross section for the 
  // specified interaction.
  //
  double xsec = 0.;
  double d2xsec_max = 0.;

  Interaction * interaction = new Interaction(*in);
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
  delete interaction;

  d2xsec_max *= TMath::Log(10.) * TMath::Log(10.);
  xsec       *= TMath::Log(10.) * TMath::Log(10.) * fdlogx * fdlogy;

  LOG("HEDISXSec", pINFO)  << "XSec[HEDIS] (E = " << Ev << " GeV) = " << xsec * (1E+38/units::cm2) << " x 1E-38 cm^2  (max = " << d2xsec_max << ")";

  string maxxsecfilename = fMaxXsecDirName + "/" + filename;
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

  GetParamDef("MaxXSec-DirName", fMaxXsecDirName, string("") ) ;
  fMaxXsecDirName = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/maxxsec/" + fMaxXsecDirName;
  if ( gSystem->mkdir(fMaxXsecDirName.c_str())==0 ) LOG("HEDISFormFactors", pINFO) << "Creating Max Xsec directory: " << fMaxXsecDirName;

  GetParamDef("DlogY", fdlogy, 0.01 ) ;
  GetParamDef("DlogX", fdlogx, 0.01 ) ;

  double flogymin = -10;
  double flogxmin = -10;
  for ( int iy=0; iy<TMath::Abs(flogymin)/fdlogy; iy++) fVy.push_back( TMath::Power( 10., flogymin + (iy+0.5)*fdlogy ) );
  for ( int ix=0; ix<TMath::Abs(flogxmin)/fdlogx; ix++) fVx.push_back( TMath::Power( 10., flogxmin + (ix+0.5)*fdlogx ) );

}

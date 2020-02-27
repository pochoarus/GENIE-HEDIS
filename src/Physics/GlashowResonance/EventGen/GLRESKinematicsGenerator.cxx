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

#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "Physics/GlashowResonance/EventGen/GLRESKinematicsGenerator.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
GLRESKinematicsGenerator::GLRESKinematicsGenerator() :
KineGeneratorWithCache("genie::GLRESKinematicsGenerator")
{

}
//___________________________________________________________________________
GLRESKinematicsGenerator::GLRESKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::GLRESKinematicsGenerator", config)
{

}
//___________________________________________________________________________
GLRESKinematicsGenerator::~GLRESKinematicsGenerator()
{

}
//___________________________________________________________________________
void GLRESKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  Interaction * interaction = evrec->Summary();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = this->MaxXSec(evrec);

  double x1min = -1.;
  double x1max =  1.;
  double x2min =  0.;
  double x2max =  1.;
  double dx1   = x1max-x1min;
  double dx2   = x2max-x2min;

  //-- Try to select a valid inelastisity y
  double xsec = -1;
  unsigned int iter = 0;
  bool accept = false;
  
  while(1) {
    iter++;
    if(iter > 1000000) {
      LOG("GLRESKinematics", pWARN)
            << "*** Could not select a valid y after "
                                            << iter << " iterations";
      evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
      genie::exceptions::EVGThreadException exception;
      exception.SetReason("Couldn't select kinematics");
      exception.SwitchOnFastForward();
      throw exception;
    }

    double x1 = x1min + dx1 * rnd->RndKine().Rndm();
    double x2 = x2min + dx2 * rnd->RndKine().Rndm();
    interaction->KinePtr()->SetKV(kKVGLx1,x1);
    interaction->KinePtr()->SetKV(kKVGLx2,x2);

    LOG("GLRESKinematics", pINFO) << "Trying: x1 = " << x1 << ", x2 = " << x2;

    //-- computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction, kPSGLx1x2fE);

    this->AssertXSecLimits(interaction, xsec, xsec_max);

    double t = xsec_max * rnd->RndKine().Rndm();
    LOG("GLRESKinematics", pDEBUG) << "xsec= "<< xsec<< ", J= 1, Rnd= "<< t;

    accept = (t<xsec);

    //-- If the generated kinematics are accepted, finish-up module's job
    if(accept) {
      LOG("GLRESKinematics", pINFO) << "Selected: x1 = " << x1 << ", x2 = " << x2;

      // set the cross section for the selected kinematics
      evrec->SetDiffXSec(xsec,kPSGLx1x2fE);

      // lock selected kinematics & clear running values
      interaction->KinePtr()->ClearRunningValues();

      LOG("GLRESKinematics", pINFO) << "Bye";

      return;
    }
  }// iterations

}
//___________________________________________________________________________
double GLRESKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  double max_xsec = -1.0;

  const int Nscan = 100;
  const int x1min = -1.;
  const int x1max =  1.;
  const int x2min =  0.;
  const int x2max =  1.;
  const double dx1 = (x1max-x1min)/(double)Nscan;
  const double dx2 = (x2max-x2min)/(double)
  Nscan;

  double scan_x1 = 0.;
  double scan_x2 = 0.;
  for (int i=0; i<Nscan; i++) {
    double x1 = x1min + dx1*i;
    for (int j=0; j<Nscan; j++) {
      double x2 = x2min + dx2*j;
      interaction->KinePtr()->SetKV(kKVGLx1,x1);
      interaction->KinePtr()->SetKV(kKVGLx2,x2);
      double dxsec = fXSecModel->XSec(interaction, kPSGLx1x2fE);
      if ( dxsec > max_xsec ) {
        scan_x1 = x1;
        scan_x2 = x2;
        max_xsec = dxsec;
      }    
    }    
  }

  utils::gsl::d2Xsec_GLRES f(fXSecModel,interaction,-1);
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit");
  min->SetFunction( f );
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(0.05);
  min->SetLimitedVariable ( 0, "x1", scan_x1, 0.001, TMath::Max(-1.,scan_x1-0.1), TMath::Min(1.,scan_x1+0.1));
  min->SetLimitedVariable ( 1, "x2", scan_x2,   0.1, TMath::Max(-0.,scan_x2-0.1), TMath::Min(1.,scan_x2+0.1));
  min->Minimize();
  interaction->KinePtr()->SetKV(kKVGLx1,min->X()[0]);
  interaction->KinePtr()->SetKV(kKVGLx2,min->X()[1]);
  max_xsec = fXSecModel->XSec(interaction, kPSGLx1x2fE);
  SLOG("GLRESKinematics", pDEBUG) << "Minimum found -> x1: " << min->X()[0] << ", x2: " << min->X()[1];

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("GLRESKinematics", pDEBUG) << interaction->AsString();
  SLOG("GLRESKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("GLRESKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double GLRESKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void GLRESKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESKinematicsGenerator::LoadConfig(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Safety factor for the maximum differential cross section
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor,  2. ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

}
//____________________________________________________________________________

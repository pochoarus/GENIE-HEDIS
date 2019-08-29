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

#include "Physics/HEDIS/EventGen/HEDISKinematicsGenerator.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
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
#include "Framework/Utils/RunOpt.h"
#include "Framework/ParticleData/PDGUtils.h"

#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

#include <string>

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
HEDISKinematicsGenerator::HEDISKinematicsGenerator() :
KineGeneratorWithCache("genie::HEDISKinematicsGenerator")
{

}
//___________________________________________________________________________
HEDISKinematicsGenerator::HEDISKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::HEDISKinematicsGenerator", config)
{

}
//___________________________________________________________________________
HEDISKinematicsGenerator::~HEDISKinematicsGenerator()
{

}
//___________________________________________________________________________
void HEDISKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction 
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  //-- Get neutrino energy and hit 'nucleon mass' 
  const InitialState & init_state = interaction->InitState();
  double Ev  = init_state.ProbeE(kRfLab);
  double M   = init_state.Tgt().HitNucP4().M(); // can be off m-shell

  //-- Get the physical W range 
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t W  = kps.Limits(kKVW);
  if(W.max <=0 || W.min>=W.max) {
     LOG("HEDISKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }


  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  //-- x and y lower limit restrict by limits in SF tables 
  double xymin = fQ2min/2/M/Ev;   // x*y = Q2/(2ME)
  if ( fXmin > xl.min ) xl.min = fXmin;
  if ( xymin > xl.min ) xl.min = xymin;  //assuming y=1 -> xmin = Q2min/(2ME)
  if ( xymin > yl.min ) yl.min = xymin;  //assuming x=1 -> xmin = Q2min/(2ME)

  double log10xmin = TMath::Log10(xl.min);  
  double log10xmax = TMath::Log10(xl.max); 
  double log10ymin = TMath::Log10(yl.min);  
  double log10ymax = TMath::Log10(yl.max); 

  LOG("HEDISKinematics", pNOTICE) << "log10x: [" << log10xmin << ", " << log10xmax << "]"; 
  LOG("HEDISKinematics", pNOTICE) << "log10y: [" << log10ymin << ", " << log10ymax << "]"; 

  //-- For the subsequent kinematic selection with the rejection method:
  double xsec_max = this->ComputeMaxXSec(interaction);

  //-- Try to select a valid (x,y) pair using the rejection method
  double dlog10x = log10xmax - log10xmin; 
  double dlog10y = log10ymax - log10ymin; 
  
  double gx=-1, gy=-1, gW=-1, gQ2=-1, xsec=-1;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
       LOG("HEDISKinematics", pWARN)
         << " Couldn't select kinematics after " << iter << " iterations";
       evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }
    
     gx = TMath::Power( 10., log10xmin + dlog10x * rnd->RndKine().Rndm() ); 
     gy = TMath::Power( 10., log10ymin + dlog10y * rnd->RndKine().Rndm() ); 

     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->Sety(gy);
     kinematics::UpdateWQ2FromXY(interaction);

     LOG("HEDISKinematics", pDEBUG) 
        << "Trying: x = " << gx << ", y = " << gy 
        << " (W  = " << interaction->KinePtr()->W()  << ","
        << "  Q2 = " << interaction->KinePtr()->Q2() << ")";

     //-- compute the cross section for current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- jacobian account for the log10 sampling
     double J = TMath::Log(10.)*gx*TMath::Log(10.)*gy; 

     //-- decide whether to accept the current kinematics
     this->AssertXSecLimits(interaction, J*xsec, xsec_max);

     double t = xsec_max * rnd->RndKine().Rndm();

     LOG("HEDISKinematics", pDEBUG) << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;

     accept = (t < J*xsec);

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
         LOG("HEDISKinematics", pNOTICE) 
            << "Selected:  x = " << gx << ", y = " << gy
            << " (W  = " << interaction->KinePtr()->W()  << ","
            << " (Q2 = " << interaction->KinePtr()->Q2() << ")";

         // reset trust bits
         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec,kPSxyfE);

         // compute W,Q2 for selected x,y
         kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

         LOG("HEDISKinematics", pNOTICE) 
                        << "Selected x,y => W = " << gW << ", Q2 = " << gQ2;

         // lock selected kinematics & clear running values
         interaction->KinePtr()->SetW (gW,  true);
         interaction->KinePtr()->SetQ2(gQ2, true);
         interaction->KinePtr()->Setx (gx,  true);
         interaction->KinePtr()->Sety (gy,  true);
         interaction->KinePtr()->ClearRunningValues();
         return;
     }
  } // iterations
}
//___________________________________________________________________________
double HEDISKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction ) const
{

  // Max Xsec must be load the first time we call this method and not anymore.
  // TODO: This might be removed in the future if Max Xsec are included in Splines.
  if (!fMaxXsecIsLoad) LoadMaxXsecFromAscii();

  LOG("HEDISKinematics", pINFO)<< "Computing max xsec in allowed phase space";
  double max_xsec = 0.0;

  const InitialState & init_state = interaction->InitState();
  double Ev  = init_state.ProbeE(kRfLab);

  int absnupdg = TMath::Abs(interaction->InitState().ProbePdg());
  HEDISQrkChannel_t chID = interaction->ExclTag().HEDISQrkChannel();
  max_xsec = fspl_max.at(chID).Spline.at(absnupdg)->Evaluate(Ev);

  int  nuc_pdgc = interaction->InitState().Tgt().HitNucPdg();
  if      ( pdg::IsProton  (nuc_pdgc) ) { max_xsec *= interaction->InitState().Tgt().Z(); }
  else if ( pdg::IsNeutron (nuc_pdgc) ) { max_xsec *= interaction->InitState().Tgt().N(); }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HEDISKinematics", pINFO) << interaction->AsString();
  LOG("HEDISKinematics", pINFO) << "Max xsec in phase space ( E= " << Ev << " GeV) = " << max_xsec;
  LOG("HEDISKinematics", pINFO) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//____________________________________________________________________________
void HEDISKinematicsGenerator::LoadMaxXsecFromAscii() const
{

  // Check that Max Xsec directory exists
  if ( gSystem->AccessPathName( fMaxXsecDirName.c_str()) ) {
    LOG("HEDISKinematics", pERROR) << "Max Xsec directory does not exist...";
    LOG("HEDISKinematics", pERROR) << fMaxXsecDirName;
    assert(0);
  }

  // Load Max Xsec to Splines
  int absnupdg[3] = { 12, 14, 16 };
  for( int n=0; n<3; n++ ) {
    for( int ch=1; ch<kHEDISQrk_numofchannels; ch++ ) {      
      string filename = fMaxXsecDirName + "/Flvr" + std::to_string(absnupdg[n]) + "_" + HEDISChannel::AsString((HEDISQrkChannel_t)ch) + ".dat";      
      if ( !gSystem->AccessPathName(filename.c_str()) ) {
        LOG("HEDISKinematics", pINFO)<< "Loading splines from: " << filename;
        fspl_max[(HEDISQrkChannel_t)ch].Spline[absnupdg[n]] = new Spline(filename,"","",false);
        fspl_max[(HEDISQrkChannel_t)ch].Spline[absnupdg[n]]->Multiply(fSafetyFactor);
      }
    
    }
  }

  // Change to true to make sure this method is called only once.
  fMaxXsecIsLoad = true;

}
//___________________________________________________________________________
void HEDISKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISKinematicsGenerator::LoadConfig(void)
{

  //-- Safety factor for the maximum differential cross section
  GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 2. ) ;
  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef("MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
  assert(fMaxXSecDiffTolerance>=0);

  // The name of the directory in which SF tables are stored is need for two reasons.
  // 1) To extract some information from the metafile about the limits of the tables
  // 2) To get the name of the Max Xsec directory
  string SFname;
  GetParamDef("SF-name", SFname, string("") ) ;

  string metaFile = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/sf/" + SFname + "/Inputs.txt";
  // make sure meta files are available
  LOG("HEDISKinematics", pDEBUG) << "Checking if file " << metaFile << " exists...";        
  if ( gSystem->AccessPathName( metaFile.c_str()) ) {
    LOG("HEDISKinematics", pERROR) << "File doesnt exist";        
    LOG("HEDISKinematics", pERROR) << "HEDIS package requires precomputation of SF using gMakeStrucFunc";        
    assert(0);
  }

  // From the metafile we need the minimum value of x and Q2 in the SF tables.
  // They are used to restrict the sampling phase space for the kinematics.
  std::ifstream meta_stream(metaFile.c_str(), std::ios::in);
  string saux;
  std::getline (meta_stream,saux); //# NX
  std::getline (meta_stream,saux);
  std::getline (meta_stream,saux); //# Xmin
  std::getline (meta_stream,saux); fXmin = atof(saux.c_str());
  std::getline (meta_stream,saux); //# NQ2
  std::getline (meta_stream,saux);
  std::getline (meta_stream,saux); //# Q2min
  std::getline (meta_stream,saux); fQ2min = atof(saux.c_str());
  meta_stream.close();

  // The following information is needed to get the name of the Max Xsec directory.
  // TODO: This might be removed in the future if Max Xsec are included in Splines.
  double dlogy;
  double dlogx;
  GetParamDef("DlogY", dlogy, 0.01 ) ;
  GetParamDef("DlogX", dlogx, 0.01 ) ;
  std::ostringstream sdlogy,sdlogx; 
  sdlogy << dlogy;
  sdlogx << dlogx;

  // Name Max Xsec directory. This shoud match the one generated by HEDISXSec.
  fMaxXsecDirName = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis/maxxsec/" + SFname + "_dx" + sdlogx.str() + "dy" + sdlogy.str();



}

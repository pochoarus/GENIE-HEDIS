//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Was first added in v2.5.1 
 @ Apr 24, 2010 - CA
   Add code to decay the off-the-mass-shell W- using PYTHIA6. 
   First complete version of the GLRES event thread.
 @ November 08, 2019 - Alfonso Garcia
   Modified to generate the kinematics of outgoing lepton properly.
   Phys. Rev. D 22, 2122 – Published 1 November 1980
*/
//____________________________________________________________________________

#include <cstring>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TClonesArray.h>
#include <TMath.h>

#include "Physics/GlashowResonance/EventGen/GLRESCohGenerator.h"
#include "Physics/HEDIS/EventGen/HEDISGenerator.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;
using namespace genie::constants;

const double kOffShellDm = 0.002; // 2 MeV

//___________________________________________________________________________
GLRESCohGenerator::GLRESCohGenerator() :
EventRecordVisitorI("genie::GLRESCohGenerator")
{
  this->Initialize();
}
//___________________________________________________________________________
GLRESCohGenerator::GLRESCohGenerator(string config) :
EventRecordVisitorI("genie::GLRESCohGenerator", config)
{
  this->Initialize();
}
//___________________________________________________________________________
GLRESCohGenerator::~GLRESCohGenerator()
{

}
//____________________________________________________________________________                                                                                        
void GLRESCohGenerator::Initialize(void) const
{
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number                                                                                                                                   
  RandomGen::Instance();
}
//___________________________________________________________________________
void GLRESCohGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu = evrec->Probe();

  GHepParticle * target = evrec -> TargetNucleus();
  if(target) evrec->AddParticle(target->Pdg(), kIStStableFinalState, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  LongLorentzVector p4_nu (nu->P4());

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  long double Ev = init_state.ProbeE(kRfLab); 

  long double Mtgt = init_state.Tgt().Z()*kProtonMass + init_state.Tgt().N()*kNeutronMass;

  long double x1 = interaction->Kine().GetKV(kKVGLx1);
  long double x2 = interaction->Kine().GetKV(kKVGLx2);
  long double x3 = interaction->Kine().GetKV(kKVGLx3);

  long double costh = x1;
  long double sinth = sqrtl(1-costh*costh);
  
  long double mlout = 0;
  if      (pdg::IsNuE  (TMath::Abs(nu->Pdg()))) mlout = kElectronMass;
  else if (pdg::IsNuMu (TMath::Abs(nu->Pdg()))) mlout = kMuonMass;
  else if (pdg::IsNuTau(TMath::Abs(nu->Pdg()))) mlout = kTauMass;

  long double s_r = 0.;
  long double Q2 = 0.;
  if (fIsEPA) {
    long double s  = 2 * Mtgt * Ev + Mtgt*Mtgt;
    long double y0 = powl(kMw+mlout,2)/s;
    long double y  = expl( logl(y0)+(logl(1)-logl(y0))*x2 );
    s_r = s * y;
  }
  else {
    long double mL      = mlout+kMw;
    long double Delta   = sqrtl( powl(2.*Ev*Mtgt-mL*mL,2)-4.*powl(Mtgt*mL,2) );
    long double s12_min = Ev/(2.*Ev+Mtgt)*(mL*mL+2.*Ev*Mtgt-Delta);
    long double s12_max = Ev/(2.*Ev+Mtgt)*(mL*mL+2.*Ev*Mtgt+Delta);
    long double s12     = expl( logl(s12_min)+(logl(s12_max)-logl(s12_min))*x2);
    long double Q2_min  = powl(s12,2)*Mtgt/(2.*Ev*(2.*Ev*Mtgt-s12));
    long double Q2_max  = s12 - mL*mL;
    Q2  = expl( logl(Q2_min) + (logl(Q2_max)-logl(Q2_min))*x3 );
    s_r = s12 - Q2;
  }
  
  long double EW = (s_r+kMw2-mlout*mlout)/sqrtl(s_r)/2.;
  long double El = (s_r-kMw2+mlout*mlout)/sqrtl(s_r)/2.;
  long double p  = sqrtl( EW*EW - kMw2 );
  LongLorentzVector p4_lout( 0., -p*sinth, -p*costh, El );

  TLorentzVector p4l_o;
  long double bz = 0.;
  long double by = 0.;
  if (fIsEPA) {
    long double EvCM = sqrtl(s_r)/2.;
    bz = (powl(Ev,2)-powl(EvCM,2))/(powl(Ev,2)+powl(EvCM,2));
    p4_lout.BoostZ(bz);
    p4l_o.SetE ((double)p4_lout.E ());
    p4l_o.SetPx((double)p4_lout.Px());
    p4l_o.SetPy((double)p4_lout.Py());
    p4l_o.SetPz((double)p4_lout.Pz());
  }
  else {
    bz = 4.*Ev*Mtgt*Q2/(Q2+s_r)/(2.*Ev*Mtgt-Q2) - (2.*Ev*Mtgt+Q2)/(2.*Ev*Mtgt-Q2);
    by = sqrtl(Mtgt*powl(Q2+s_r,2)/(2.*Ev*Q2*(s_r+Q2-2.*Ev*Mtgt))+1.);
    p4_lout.BoostZ(-bz);
    p4_lout.BoostY(-by);
    p4l_o.SetE ((double)p4_lout.E ());
    p4l_o.SetPx((double)p4_lout.Px());
    p4l_o.SetPy((double)p4_lout.Py());
    p4l_o.SetPz((double)p4_lout.Pz());
    p4l_o.RotateX((double)acosl(by)-kPi/2.);
  }

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2 * kPi * rnd->RndLep().Rndm();
  p4l_o.RotateZ(phi);

  //rotate from LAB=[0,0,Ev,Ev]->[px,py,pz,E]
  p4l_o.RotateUz(unit_nu);

  int pdglout = 0;
  if      (pdg::IsAntiNuE  (nu->Pdg())) pdglout = kPdgPositron;
  else if (pdg::IsNuE      (nu->Pdg())) pdglout = kPdgElectron;
  else if (pdg::IsAntiNuMu (nu->Pdg())) pdglout = kPdgAntiMuon;
  else if (pdg::IsNuMu     (nu->Pdg())) pdglout = kPdgMuon;
  else if (pdg::IsAntiNuTau(nu->Pdg())) pdglout = kPdgAntiTau;
  else if (pdg::IsNuTau    (nu->Pdg())) pdglout = kPdgTau;

  // Create a GHepParticle and add it to the event record
  evrec->AddParticle( pdglout,  kIStStableFinalState, 0, -1, -1, -1, p4l_o, *(nu->X4()) );

  int pdgboson = pdg::IsNeutrino(init_state.ProbePdg()) ? kPdgWP : kPdgWM;

  double width = fPythia->GetPMAS(24,2);
  fPythia->SetPMAS(24,2,0.); //set to 0 the width of the W boson to avoid problems with energy conservation
  fPythia->Py1ent( -1, pdgboson, EW, acosl(costh), kPi/2. ); //k(1,2) = 2
  fPythia->SetPMAS(24,2,width);

  fPythia->Pyexec();
  //fPythia->Pylist(3);

  fPythia->GetPrimaries();
  TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
  int np = pythia_particles->GetEntries();
  assert(np>0);

  TMCParticle * particle = 0;
  TIter piter(pythia_particles);
  while( (particle = (TMCParticle *) piter.Next()) ) {
    
    int pdgc = particle->GetKF();
    int ks = particle->GetKS();
   
    LongLorentzVector p4longo(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetEnergy());
    
    TLorentzVector p4o;
    if (fIsEPA) {
      p4longo.BoostZ(bz);
      p4o.SetE ((double)p4longo.E ());
      p4o.SetPx((double)p4longo.Px());
      p4o.SetPy((double)p4longo.Py());
      p4o.SetPz((double)p4longo.Pz());
    }
    else {
      p4longo.BoostZ(-bz);
      p4longo.BoostY(-by);
      p4o.SetE ((double)p4longo.E ());
      p4o.SetPx((double)p4longo.Px());
      p4o.SetPy((double)p4longo.Py());
      p4o.SetPz((double)p4longo.Pz());
      p4o.RotateX((double)acosl(by)-kPi/2.);
    }

    //p4o.RotateZ(phi);
    p4o.RotateUz(unit_nu);

    TParticlePDG * part = PDGLibrary::Instance()->Find(pdgc);
    if ( (ks==1 || ks==4) && p4o.E() < part->Mass() ) {
      LOG("GLRESCohGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("GLRESCohGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
      LOG("GLRESCohGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P()); 
      LOG("GLRESCohGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
      LOG("GLRESCohGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << part->Mass();
      p4o.SetXYZT(0,0,0,part->Mass());
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

    // fix numbering scheme used for mother/daughter assignments
    int firstmother = -1;
    int lastmother  = -1;
    int firstchild  = -1;
    int lastchild   = -1;
    
    if (particle->GetParent()==0) {
      firstmother = 0;
    }
    else {
      firstmother = particle->GetParent() + 3;       
      if (particle->GetFirstChild()!=0) firstchild  = particle->GetFirstChild() + 3;
      if (particle->GetLastChild() !=0) lastchild   = particle->GetLastChild() + 3;

    }

    double lightspeed = 299792458e3; //c in mm/s. Used for time in PYTHIA t[s]=t_pythia[mm]/c[mm/s]
    double vx = nu->X4()->X() + particle->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
    double vy = nu->X4()->Y() + particle->GetVy()*1e12;
    double vz = nu->X4()->Z() + particle->GetVz()*1e12;
    double vt = nu->X4()->T() + particle->GetTime()/lightspeed;
    TLorentzVector pos( vx, vy, vz, vt );

    evrec->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );

  }

  delete particle;
  pythia_particles->Clear("C");

}
//___________________________________________________________________________
void GLRESCohGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESCohGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESCohGenerator::LoadConfig(void)
{

  // PYTHIA parameters only valid for HEDIS
  double wmin;        GetParam( "Wmin",          wmin ) ;
  int warnings;       GetParam( "Warnings",      warnings ) ;
  int errors;         GetParam( "Errors",        errors ) ;
  int qrk_mass;       GetParam( "QuarkMass",     qrk_mass ) ;
  int decaycut;       GetParam( "DecayCutOff",   decaycut ) ;
  double decaylength; GetParam( "DecayLength",   decaylength ) ;
  fPythia->SetPARP(2,  wmin);         //(D = 10. GeV) lowest c.m. energy for the event as a whole that the program will accept to simulate. (bellow 2GeV pythia crashes)
  fPythia->SetMSTU(26, warnings);     // (Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, errors);       // (Default=10) maximum number of errors that are printed
  fPythia->SetMSTJ(93, qrk_mass);     // light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.
  fPythia->SetMSTJ(22, decaycut);     // (Default=1) cut-off on decay length for a particle that is allowed to decay according to MSTJ(21) and the MDCY value
  fPythia->SetPARJ(71, decaylength);  // (Default=10. mm) maximum average proper lifetime cτ for particles allowed to decay

  fPythia->SetMDME(192,1,0);          //swicht off W decay to top
  fPythia->SetMDME(196,1,0); 
  fPythia->SetMDME(200,1,0); 

  fPythia->SetPMAS(24,1,kMw); //mass of the W boson (pythia=80.450 // genie=80.385)

  GetParam( "GLRESCoh-Is-EPA",    fIsEPA );

}
//____________________________________________________________________________


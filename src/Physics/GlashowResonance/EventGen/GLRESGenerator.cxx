#include <cstring>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TClonesArray.h>
#include <TMath.h>

#include "Physics/GlashowResonance/EventGen/GLRESGenerator.h"
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
GLRESGenerator::GLRESGenerator() :
EventRecordVisitorI("genie::GLRESGenerator")
{
  this->Initialize();
}
//___________________________________________________________________________
GLRESGenerator::GLRESGenerator(string config) :
EventRecordVisitorI("genie::GLRESGenerator", config)
{
  this->Initialize();
}
//___________________________________________________________________________
GLRESGenerator::~GLRESGenerator()
{

}
//____________________________________________________________________________                                                                                        
void GLRESGenerator::Initialize(void) const
{
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number                                                                                                                                   
  RandomGen::Instance();
}
//___________________________________________________________________________
void GLRESGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu     = evrec->Probe();
  GHepParticle * el     = evrec->HitElectron();
  GHepParticle * target = evrec -> TargetNucleus();
  if(target) evrec->AddParticle(target->Pdg(), kIStStableFinalState, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  TLorentzVector p4_nu (* nu->P4());
  TLorentzVector p4_el (* el->P4());

  TLorentzVector p4_W = p4_nu + p4_el;

  double Wmass = p4_W.M();
  LOG("GLRESGenerator", pINFO) << "Wmass : " << Wmass;

  if( Wmass <= 2. ) {
    LOG("GLRESGenerator", pWARN) << "Low invariant mass, W = " << Wmass << " GeV!!";
    
    evrec->EventFlags()->SetBitNumber(kHadroSysGenErr, true);

    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    
    return;
  }

  // Final state primary lepton PDG code
  int pdgl = interaction->FSPrimLeptonPdg();
  assert(pdgl!=0);

  LOG("GLRESGenerator", pINFO) << "Channel : " << interaction->FSPrimLeptonPdg();

  if ( pdg::IsElectron(pdgl) || pdg::IsMuon(pdgl) || pdg::IsTau(pdgl) ) {

    // Get selected kinematics
    double y = interaction->Kine().y(true);
    assert(y>0 && y<1);
  
    // Compute the neutrino and muon energy
    double Ev  = (long double) init_state.ProbeE(kRfLab); 
    double El  = y*Ev;

    LOG("GLRESGenerator", pINFO) << "Ev = " << Ev << ", y = " << y << ", -> El = " << El;
    
    // Compute the momentum transfer and scattering angle
    double El2   = TMath::Power(El,2);
    double me    = (long double) kElectronMass;
    double ml    = (long double) interaction->FSPrimLepton()->Mass();
    double ml2   = TMath::Power(ml,2);
    double pl    = TMath::Sqrt(El2-ml2);   
    
    double Q2    = 2*(Ev-El)*me + me*me;
    double costh = (El-0.5*(Q2+ml2)/Ev)/pl;
    double sinth = TMath::Sqrt( 1-TMath::Power(costh,2.) );
    // Randomize transverse components
    RandomGen * rnd = RandomGen::Instance();
    double phi  = 2* M_PIl * rnd->RndLep().Rndm();
    
    LOG("GLRESGenerator", pNOTICE) << "Q2 = " << Q2 << ", cos(theta) = " << costh << ", phi = " << phi;
   
    double plx = pl*sinth*TMath::Cos(phi);
    double ply = pl*sinth*TMath::Sin(phi);
    double plz = pl*costh;

    //rotate from LAB (where neutrino is [0,0,E,E]) to LAB' (where neutrino is [px,py,pz,E])
    TVector3 unit_nudir = p4_nu.Vect().Unit(); 
    TVector3 p3l(plx,ply,plz);
    p3l.RotateUz(unit_nudir);

    TLorentzVector p4_lpout( p3l, El );
    TLorentzVector p4_nuout = p4_nu + p4_el - p4_lpout;

    int pdgvout = 0;
    if      ( pdg::IsElectron(pdgl) ) pdgvout = kPdgAntiNuE;
    else if ( pdg::IsMuon(pdgl)     ) pdgvout = kPdgAntiNuMu;
    else if ( pdg::IsTau(pdgl)      ) pdgvout = kPdgAntiNuTau;

    // Create a GHepParticle and add it to the event record
    evrec->AddParticle( kPdgWM,  kIStDecayedState,     0, -1,  5,  6, p4_W,     *(nu->X4()) ); //W [mothers: nuebar_in,e_in][daugthers: nulbar_out,lp_out]
    evrec->AddParticle( pdgvout, kIStStableFinalState, 4, -1, -1, -1, p4_nuout, *(nu->X4()) );
    evrec->AddParticle( pdgl,    kIStStableFinalState, 4, -1, -1, -1, p4_lpout, *(nu->X4()) );
    evrec->Summary()->KinePtr()->SetFSLeptonP4(p4_lpout);

  }
  else {

    char p6frame[10], p6nu[10], p6tgt[10];
    strcpy(p6frame, "CMS"    );
    strcpy(p6nu,    "nu_ebar");
    strcpy(p6tgt,   "e-"     );
    fPythia->Pyinit(p6frame, p6nu, p6tgt, (double)Wmass);
    
    fPythia->SetMDME(206,1,0); //swicht off W decay leptonic modes
    fPythia->SetMDME(207,1,0); 
    fPythia->SetMDME(208,1,0); 

    fPythia->Pyevnt();
    //fPythia->Pylist(2);

    fPythia->SetMDME(206,1,1); //swicht them back on
    fPythia->SetMDME(207,1,1); 
    fPythia->SetMDME(208,1,1); 

    // get LUJETS record
    fPythia->GetPrimaries();
    TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
    int np = pythia_particles->GetEntries();
    assert(np>0);

    // Vector defining rotation from LAB to LAB' (z:= \vec{resonance momentum})
    TVector3 unit_Wdir = p4_W.Vect().Unit();

    // Boost velocity LAB' -> Resonance rest frame
    TVector3 beta(0,0,p4_W.P()/p4_W.Energy());

    TMCParticle * particle = 0;
    TIter piter(pythia_particles);
    while( (particle = (TMCParticle *) piter.Next()) ) {
      
      if ( particle->GetKS()==21 ) { continue; } //we dont want to save first particles from pythia (init states)

      TLorentzVector p4o(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetEnergy());
      p4o.Boost(beta); 
      TVector3 p3 = p4o.Vect();
      p3.RotateUz(unit_Wdir); 
      TLorentzVector p4(p3,p4o.Energy());

      TParticlePDG * part = PDGLibrary::Instance()->Find(particle->GetKF());
      if ( particle->GetKS()==1 && p4.E() < part->Mass() ) {
        LOG("GLRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("GLRESGenerator", pWARN) << "PDG = " << particle->GetKF() << " // State = " << particle->GetKS();
        LOG("GLRESGenerator", pWARN) << "E = " << p4.E() << " // |p| = " << TMath::Sqrt(p4.P()); 
        LOG("GLRESGenerator", pWARN) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
        LOG("GLRESGenerator", pWARN) << "m    = " << p4.M() << " // mpdg = " << part->Mass();
        p4.SetXYZT(0,0,0,part->Mass());
      }

      // copy final state particles to the event record
      GHepStatus_t ist = (particle->GetKS()==1) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

      // fix numbering scheme used for mother/daughter assignments
      int firstmother = -1;
      int lastmother  = -1;
      int firstchild  = -1;
      int lastchild   = -1;

      if ( particle->GetParent() < 10 ) {
        if      ( TMath::Abs(particle->GetKF())<7 ) {   //outgoing quarks: mother will be the boson (saved in position 4)
          firstmother = 4;       
          firstchild  = particle->GetFirstChild() - 6;
          lastchild   = particle->GetLastChild()  - 6;
        }
        else if ( TMath::Abs(particle->GetKF())==24 ) { //produced W boson: mother will be the incoming neutrino
          firstmother = 0;
          firstchild  = particle->GetFirstChild() - 6;
          lastchild   = particle->GetLastChild()  - 6;
        }
        else if ( particle->GetKF()==22 ) {             //radiative photons: mother will be the incoming electron
          firstmother = 2; 
        }
      }
      else { //rest
        firstmother = particle->GetParent()     - 6; //shift to match boson position
        firstchild  = (particle->GetFirstChild()==0) ? particle->GetFirstChild() - 1 : particle->GetFirstChild() - 6;
        lastchild   = (particle->GetLastChild()==0)  ? particle->GetLastChild()  - 1 : particle->GetLastChild()  - 6;
      }

      double lightspeed = 299792458e3; //c in mm/s. Used for time in PYTHIA t[s]=t_pythia[mm]/c[mm/s]
      double vx = nu->X4()->X() + particle->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
      double vy = nu->X4()->Y() + particle->GetVy()*1e12;
      double vz = nu->X4()->Z() + particle->GetVz()*1e12;
      double vt = nu->X4()->T() + particle->GetTime()/lightspeed;
      TLorentzVector pos( vx, vy, vz, vt );

      LOG("GLRESGenerator", pDEBUG) << particle->GetKF() << "  " << ist << "  " << firstmother << "  " << firstchild;

      evrec->AddParticle(particle->GetKF(), ist, firstmother, lastmother, firstchild, lastchild, p4, pos );

    }
  
    delete particle;
    pythia_particles->Clear("C");


  }
}
//___________________________________________________________________________
void GLRESGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::LoadConfig(void)
{

  fPythia->SetMSTJ(93, 1);   //light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.
  fPythia->SetPARJ(22, 2);   //(Default=1) cut-off on decay length for a particle that is allowed to decay according to MSTJ(21) and the MDCY value
  fPythia->SetPARJ(32, 0.5); //(Default=1. GeV) is, with quark masses added, used to define the minimum allowable energy of a colour-singlet parton system.
  fPythia->SetPARJ(71, 2.);  //(Default=10. mm) maximum average proper lifetime cτ for particles allowed to decay
  fPythia->SetMSTU(26, 1000000);  //(Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, 1000000);  //(Default=10) maximum number of errors that are printed
  fPythia->SetPARP(2, 2.);  //(D = 10. GeV) lowest c.m. energy for the event as a whole that the program will accept to simulate. (bellow 2GeV pythia crashes)

}
//____________________________________________________________________________


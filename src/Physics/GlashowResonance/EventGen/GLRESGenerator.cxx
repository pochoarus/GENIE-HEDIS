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

#include "math.h"

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

  // Neutrino 4p
  long double p4v[4];
  p4v[0] = (long double)nu->Px();
  p4v[1] = (long double)nu->Py();
  p4v[2] = (long double)nu->Pz();
  p4v[3] = (long double)nu->E();

  // electron 4p
  long double p4e[4];
  p4e[0] = (long double)el->Px();
  p4e[1] = (long double)el->Py();
  p4e[2] = (long double)el->Pz();
  p4e[3] = (long double)el->E();

  long double p4Wlong[4];
  for (int i=0; i<4; i++) p4Wlong[i] = p4v[i] + p4e[i];  

  TLorentzVector p4_W( (double)p4Wlong[0], (double)p4Wlong[1], (double)p4Wlong[2], (double)p4Wlong[3] );

  long double Wmass = sqrtl( p4Wlong[3]*p4Wlong[3] - p4Wlong[0]*p4Wlong[0] - p4Wlong[1]*p4Wlong[1] - p4Wlong[2]*p4Wlong[2] );
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
    long double y = interaction->Kine().y(true);
    assert(y>0 && y<1);

    
    // Compute the neutrino and muon energy
    long double Ev  = (long double) init_state.ProbeE(kRfLab); 
    long double El  = y*Ev;

    LOG("GLRESGenerator", pINFO) << "Ev = " << Ev << ", y = " << y << ", -> El = " << El;
    
    // Compute the momentum transfer and scattering angle
    long double El2   = powl(El,2);
    long double me    = (long double) kElectronMass;
    long double ml    = (long double) interaction->FSPrimLepton()->Mass();
    long double ml2   = powl(ml,2);
    long double pl    = sqrtl(El2-ml2);   
    assert(El2>=ml2);
    
    long double Q2    = 2*(Ev-El)*me + me*me;
    long double costh = (El-0.5*(Q2+ml2)/Ev)/pl;
    long double sinth = sqrtl( fmaxl(0., 1-powl(costh,2.)) );
    
    LOG("GLRESGenerator", pNOTICE) << "Q2 = " << Q2 << ", cos(theta) = " << costh;
    
    //warn about overflow in costheta and ignore it if it is small or abort
    if( abs(costh)>1 ) LOG("GLRESGenerator", pWARN) << "El = " << El << ", Ev = " << Ev << ", cos(theta) = " << costh;
    assert(abs(costh)<=1);
   
    // Randomize transverse components
    RandomGen * rnd = RandomGen::Instance();
    long double phi  = 2* M_PIl * rnd->RndLep().Rndm();
    long double p4loutlong[4];
    p4loutlong[0] = pl * sinth * cosl(phi);
    p4loutlong[1] = pl * sinth * sinl(phi);
    p4loutlong[2] = pl * costh; // p(//)
    p4loutlong[3] = El;

    // Rotate lepton momentum vector from the reference frame (x'y'z') where 
    // {z':(neutrino direction), z'x':(theta plane)} to the LAB
    long double tot2 = p4v[0]*p4v[0] + p4v[1]*p4v[1] + p4v[2]*p4v[2];
    long double tot = (tot2 > 0) ?  1.0/sqrtl(tot2) : 1.0;
    long double unit_nudir[3] = { p4v[0]*tot, p4v[1]*tot, p4v[2]*tot };

    long double up = unit_nudir[0]*unit_nudir[0] + unit_nudir[1]*unit_nudir[1];
    if (up) {
      up = sqrtl(up);
      long double px = p4loutlong[0],  py = p4loutlong[1],  pz = p4loutlong[2];
      p4loutlong[0] = (unit_nudir[0]*unit_nudir[2]*px - unit_nudir[1]*py + unit_nudir[0]*up*pz)/up;
      p4loutlong[1] = (unit_nudir[1]*unit_nudir[2]*px + unit_nudir[0]*py + unit_nudir[1]*up*pz)/up;
      p4loutlong[2] = (unit_nudir[2]*unit_nudir[2]*px -               px + unit_nudir[2]*up*pz)/up;
    } 
    else if (unit_nudir[2] < 0.) { // phi=0  teta=pi
      p4loutlong[0] = -p4loutlong[0]; 
      p4loutlong[2] = -p4loutlong[2]; 
    }
  
    // Figure out the Final State Lepton PDG Code
    TParticlePDG * p = PDGLibrary::Instance()->Find(pdgl);
    double Mpdg = p->Mass();
    if( abs(sqrtl(p4loutlong[3]*p4loutlong[3] - p4loutlong[0]*p4loutlong[0] - p4loutlong[1]*p4loutlong[1] - p4loutlong[2]*p4loutlong[2])-Mpdg) > kOffShellDm ) {
      LOG("GLRESGenerator", pERROR)
                 << "*** Selected kinematics lead to off mass shell lepton!";
       evrec->EventFlags()->SetBitNumber(kLeptoGenErr, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("E<m for final state lepton");
       exception.SwitchOnFastForward();
       throw exception;
    }

    TLorentzVector p4lout( (double)p4loutlong[0], (double)p4loutlong[1], (double)p4loutlong[2], (double)p4loutlong[3] );
    p4lout.SetE( TMath::Sqrt( p4lout.Vect().Mag2() + Mpdg*Mpdg ) );

    long double p4voutlong[4];
    for (int i=0; i<4; i++) p4voutlong[i] = p4v[i] + p4e[i] - p4loutlong[i];  

    TLorentzVector p4vout( (double)p4voutlong[0], (double)p4voutlong[1], (double)p4voutlong[2], (double)p4voutlong[3] );
    p4vout.SetE( TMath::Sqrt( p4vout.Vect().Mag2() ) );

    int pdgvout = 0;
    if      ( pdg::IsElectron(pdgl) ) pdgvout = kPdgAntiNuE;
    else if ( pdg::IsMuon(pdgl)     ) pdgvout = kPdgAntiNuMu;
    else if ( pdg::IsTau(pdgl)      ) pdgvout = kPdgAntiNuTau;
    assert(pdgvout!=0);

    // Create a GHepParticle and add it to the event record
    evrec->AddParticle( kPdgWM,  kIStDecayedState,     0, -1,  5,  6, p4_W,   *(nu->X4()) ); //W [mothers: nuebar_in,e_in][daugthers: nulbar_out,l_out]
    evrec->AddParticle( pdgvout, kIStStableFinalState, 4, -1, -1, -1, p4vout, *(nu->X4()) );

    // Create a GHepParticle and add it to the event record
    if ( pdg::IsTau(abs(pdgl)) ) {

      evrec->AddParticle( pdgl, kIStDecayedState, 4, -1, -1, -1, p4lout, *(nu->X4()) );
      evrec->Summary()->KinePtr()->SetFSLeptonP4(p4lout);
      int mom = evrec->GetEntries() - 1;

      fPythia->Py1ent( -1, pdgl, p4lout.E(), acos(p4lout.Pz()/p4lout.E()), atan(p4lout.Py()/p4lout.Px()) ); //k(1,2) = 2
      fPythia->Pyexec();
      // get LUJETS record
      fPythia->GetPrimaries();
      TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
      int np = pythia_particles->GetEntries();
      assert(np>0);

      TMCParticle * particle = 0;
      TIter piter(pythia_particles);
      while( (particle = (TMCParticle *) piter.Next()) ) {
        
        if ( particle->GetParent()==0 ) continue; //we dont want to save tau, already saved

        int pdgc = particle->GetKF();
        int ks   = particle->GetKS();

        TLorentzVector p4(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetEnergy());

        // copy final state particles to the event record
        GHepStatus_t ist = (ks==1) ? kIStStableFinalState : kIStDecayedState;

        int im  = mom - 1 + particle->GetParent();
        int ifc = (particle->GetFirstChild() <= 0) ? -1 : mom - 1 + particle->GetFirstChild();
        int ilc = (particle->GetLastChild()  <= 0) ? -1 : mom - 1 + particle->GetLastChild();

        TParticlePDG * part = PDGLibrary::Instance()->Find(pdgc);
        if ( ks==1 && p4.E() < part->Mass() ) {
          LOG("GLRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
          LOG("GLRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
          LOG("GLRESGenerator", pWARN) << "E = " << p4.E() << " // |p| = " << TMath::Sqrt(p4.P()); 
          LOG("GLRESGenerator", pWARN) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
          LOG("GLRESGenerator", pWARN) << "m    = " << p4.M() << " // mpdg = " << part->Mass();
          p4.SetXYZT(0,0,0,part->Mass());
        }

        double lightspeed = 299792458e3; //c in mm/s. Used for time in PYTHIA t[s]=t_pythia[mm]/c[mm/s]
        double vx = evrec->Probe()->X4()->X() + particle->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
        double vy = evrec->Probe()->X4()->Y() + particle->GetVy()*1e12;
        double vz = evrec->Probe()->X4()->Z() + particle->GetVz()*1e12;
        double vt = evrec->Probe()->X4()->T() + particle->GetTime()/lightspeed;
        TLorentzVector pos( vx, vy, vz, vt );

        LOG("GLRESGenerator", pDEBUG) << pdgc << "  " << ist << "  " << im << "  " << ifc;
        evrec->AddParticle( pdgc, ist, im, -1, ifc, ilc, p4, pos );

      }
      delete particle;
      pythia_particles->Clear("C");

    }     
    else {
      evrec->AddParticle( pdgl,    kIStStableFinalState, 4, -1, -1, -1, p4lout, *(nu->X4()) );
      evrec->Summary()->KinePtr()->SetFSLeptonP4(p4lout);
    }

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
    long double tot2 = p4Wlong[0]*p4Wlong[0] + p4Wlong[1]*p4Wlong[1] + p4Wlong[2]*p4Wlong[2];
    long double tot = (tot2 > 0) ?  1.0/sqrtl(tot2) : 1.0;
    long double unit_Wdir[3] = { p4Wlong[0]*tot, p4Wlong[1]*tot, p4Wlong[2]*tot };

    // Boost velocity LAB' -> Resonance rest frame
    long double beta[3] = { p4Wlong[0]/p4Wlong[3], p4Wlong[1]/p4Wlong[3], p4Wlong[2]/p4Wlong[3] };

    TMCParticle * particle = 0;
    TIter piter(pythia_particles);
    while( (particle = (TMCParticle *) piter.Next()) ) {
      
      if ( particle->GetKS()==21 ) { continue; } //we dont want to save first particles from pythia (init states)

      long double p4o[4];
      p4o[0] = (long double)particle->GetPx();
      p4o[1] = (long double)particle->GetPy();
      p4o[2] = (long double)particle->GetPz();
      p4o[3] = (long double)particle->GetEnergy();

      LongBoost( beta[0], beta[1], beta[2], p4o[0], p4o[1], p4o[2], p4o[3]);

      long double up = unit_Wdir[0]*unit_Wdir[0] + unit_Wdir[1]*unit_Wdir[1];
      if (up) {
        up = sqrtl(up);
        long double px = p4o[0],  py = p4o[1],  pz = p4o[2];
        p4o[0] = (unit_Wdir[0]*unit_Wdir[2]*px - unit_Wdir[1]*py + unit_Wdir[0]*up*pz)/up;
        p4o[1] = (unit_Wdir[1]*unit_Wdir[2]*px + unit_Wdir[0]*py + unit_Wdir[1]*up*pz)/up;
        p4o[2] = (unit_Wdir[2]*unit_Wdir[2]*px -              px + unit_Wdir[2]*up*pz)/up;
      } 
      else if (unit_Wdir[2] < 0.) { // phi=0  teta=pi
        p4o[0] = -p4o[0]; 
        p4o[2] = -p4o[2]; 
      }

      TLorentzVector p4((double)p4o[0],(double)p4o[1],(double)p4o[2],(double)p4o[3]);

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
void GLRESGenerator::LongBoost(long double bx, long double by, long double bz, long double &px, long double &py, long double &pz, long double &E) const
{

  long double b2 = bx*bx + by*by + bz*bz;
  long double gamma = 1.0 / sqrtl(1.0 - b2);
  long double bp = bx*px + by*py + bz*pz;
  long double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
  px = px + gamma2*bp*bx + gamma*bx*E;
  py = py + gamma2*bp*by + gamma*by*E;
  pz = pz + gamma2*bp*bz + gamma*bz*E;
  E = gamma*(E + bp);

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


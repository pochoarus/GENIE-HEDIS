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

#include "Physics/GlashowResonance/EventGen/GLRESGenerator.h"
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
GLRESGenerator::GLRESGenerator() :
EventRecordVisitorI("genie::GLRESGenerator")
{
  this->Initialize();
  born = new GLRESBornPXSec();
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
  const ProcessInfo & proc_info   = interaction->ProcInfo();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu     = evrec->Probe();
  GHepParticle * el;
  if      ( proc_info.IsGlashowResonanceAtomic() ) el = evrec->HitElectron();
  else if ( proc_info.IsGlashowResonanceInel()   ) el = evrec->HitNucleon();

  GHepParticle * target = evrec -> TargetNucleus();
  if(target) evrec->AddParticle(target->Pdg(), kIStStableFinalState, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  LongLorentzVector p4_nu (nu->P4());
  LongLorentzVector p4_el (el->P4());

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  long double Ev = init_state.ProbeE(kRfLab); 

  long double Mtarget = 0.;
  long double mlin    = 0.;
  long double mlout   = interaction->FSPrimLepton()->Mass();
  if      ( proc_info.IsGlashowResonanceAtomic() ) {
    if      (pdg::IsNuE  (TMath::Abs(nu->Pdg()))) mlin = kElectronMass;
    else if (pdg::IsNuMu (TMath::Abs(nu->Pdg()))) mlin = kMuonMass;
    else if (pdg::IsNuTau(TMath::Abs(nu->Pdg()))) mlin = kTauMass;
    Mtarget = mlin;
  }
  else if ( proc_info.IsGlashowResonanceInel()) {
    Mtarget = init_state.Tgt().HitNucMass();
  }
  
  long double s = 2 * Mtarget * Ev + Mtarget*Mtarget;

  long double x1 = interaction->Kine().GetKV(kKVGLx1);
  long double x2 = interaction->Kine().GetKV(kKVGLx2);

  long double costh = x1;
  long double sinth = sqrtl(1-costh*costh);
  
  long double s_r = s;
  if (proc_info.IsGlashowResonanceAtomic()) {
    if (fIsNLO) {
      long double t = born->GetT(0.,mlin,mlout,0.,s,x1);
      long double zeta  = born->GetReAlpha()/kPi*(2.0*logl(sqrtl(-t)/kElectronMass)-1.0);
      long double omx   = powl(x2, 1.0/zeta );
      s_r *= ( 1.-omx );
    }
  }
  else if (proc_info.IsGlashowResonanceInel()) {
    long double xmin = fQ2PDFmin/2./Ev/Mtarget;
    long double x = expl( logl(xmin) + (logl(1.0)-logl(xmin))*x2 );
    s_r *= x; 
  }

  // Boost velocity CM -> LAB
  long double EvCM = (s_r-mlin*mlin)/sqrtl(s_r)/2.;
  long double beta = (powl(Ev,2)-powl(EvCM,2))/(powl(Ev,2)+powl(EvCM,2));

  // Final state primary lepton PDG code
  int pdgl = interaction->FSPrimLeptonPdg();
  assert(pdgl!=0);

  if ( pdg::IsElectron(TMath::Abs(pdgl)) || pdg::IsMuon(TMath::Abs(pdgl)) || pdg::IsTau(TMath::Abs(pdgl)) ) {

    long double Elpout = (s_r+mlout*mlout)/sqrtl(s_r)/2.;
    long double Enuout = (s_r-mlout*mlout)/sqrtl(s_r)/2.;
    LongLorentzVector p4_lpout( 0.,  Enuout*sinth,  Enuout*costh, Elpout );
    LongLorentzVector p4_nuout( 0., -Enuout*sinth, -Enuout*costh, Enuout );

    p4_lpout.BoostZ(beta);
    p4_nuout.BoostZ(beta);

    TLorentzVector p4lp_o( (double)p4_lpout.Px(), (double)p4_lpout.Py(), (double)p4_lpout.Pz(), (double)p4_lpout.E() );
    TLorentzVector p4nu_o( (double)p4_nuout.Px(), (double)p4_nuout.Py(), (double)p4_nuout.Pz(), (double)p4_nuout.E() );

    // Randomize transverse components
    RandomGen * rnd = RandomGen::Instance();
    double phi  = 2* kPi * rnd->RndLep().Rndm();
    p4lp_o.RotateZ(phi);
    p4nu_o.RotateZ(phi);

    //rotate from LAB=[0,0,Ev,Ev]->[px,py,pz,E]
    p4lp_o.RotateUz(unit_nu);
    p4nu_o.RotateUz(unit_nu);

    int pdgvout = 0;
    if      ( pdg::IsElectron(pdgl) ) pdgvout = kPdgAntiNuE;
    else if ( pdg::IsPositron(pdgl) ) pdgvout = kPdgNuE;
    else if ( pdg::IsMuon(pdgl)     ) pdgvout = kPdgAntiNuMu;
    else if ( pdg::IsAntiMuon(pdgl) ) pdgvout = kPdgNuMu;
    else if ( pdg::IsTau(pdgl)      ) pdgvout = kPdgAntiNuTau;
    else if ( pdg::IsAntiTau(pdgl)  ) pdgvout = kPdgNuTau;

    int pdgboson = pdg::IsNeutrino(init_state.ProbePdg()) ? kPdgWP : kPdgWM;

    // Create a GHepParticle and add it to the event record
    evrec->AddParticle( pdgboson, kIStDecayedState,     0, -1,  5,  6, *nu->P4()+*el->P4(), *(nu->X4()) ); //W [mothers: nuebar_in,e_in][daugthers: nulbar_out,lp_out]
    evrec->AddParticle( pdgl,     kIStStableFinalState, 4, -1, -1, -1, p4lp_o,              *(nu->X4()) );
    evrec->AddParticle( pdgvout,  kIStStableFinalState, 4, -1, -1, -1, p4nu_o,              *(nu->X4()) );
    evrec->Summary()->KinePtr()->SetFSLeptonP4(p4lp_o);

  }
  else {

    char p6frame[10];
    strcpy(p6frame, "CMS"    );

    char p6nu[10], p6tgt[10];
    if      (pdg::IsAntiNeutrino(nu->Pdg())) { strcpy(p6nu, "nu_ebar");   strcpy(p6tgt, "e-");   }
    else if (pdg::IsNeutrino    (nu->Pdg())) { strcpy(p6nu, "nu_e");      strcpy(p6tgt, "e+");   }

    fPythia->Pyinit(p6frame, p6nu, p6tgt, sqrtl(s_r));
    
    fPythia->Pyevnt();
    //fPythia->Pylist(3);

    // get LUJETS record
    fPythia->GetPrimaries();
    TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
    int np = pythia_particles->GetEntries();
    assert(np>0);

    TMCParticle * particle = 0;
    TIter piter(pythia_particles);
    while( (particle = (TMCParticle *) piter.Next()) ) {
      
      int pdgc = particle->GetKF();
      int ks = particle->GetKS();

      if ( ks==21 ) { continue; } //we dont want to save first particles from pythia (init states)

      LongLorentzVector p4longo(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetEnergy());
      p4longo.BoostZ(beta); 

      TLorentzVector p4o( (double)p4longo.Px(), (double)p4longo.Py(), (double)p4longo.Pz(), (double)p4longo.E() );
      p4o.RotateUz(unit_nu); 
     
      TParticlePDG * part = PDGLibrary::Instance()->Find(pdgc);
      if ( (ks==1 || ks==4) && p4o.E() < part->Mass() ) {
        LOG("GLRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("GLRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
        LOG("GLRESGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P()); 
        LOG("GLRESGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
        LOG("GLRESGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << part->Mass();
        p4o.SetXYZT(0,0,0,part->Mass());
      }

      // copy final state particles to the event record
      GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

      // fix numbering scheme used for mother/daughter assignments
      int firstmother = -1;
      int lastmother  = -1;
      int firstchild  = -1;
      int lastchild   = -1;

      if ( particle->GetParent() < 10 ) {
        if      ( TMath::Abs(pdgc)<7 ) {   //outgoing quarks: mother will be the boson (saved in position 4)
          firstmother = 4;       
          firstchild  = particle->GetFirstChild() - 6;
          lastchild   = particle->GetLastChild()  - 6;
        }
        else if ( TMath::Abs(pdgc)==24 ) { //produced W boson: mother will be the incoming neutrino
          firstmother = 0;
          firstchild  = particle->GetFirstChild() - 6;
          lastchild   = particle->GetLastChild()  - 6;
        }
        else if ( pdgc==22 ) {             //radiative photons: mother will be the incoming electron
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

      evrec->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );

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

  fPythia->SetMSTP(61,0);             // (Default=2) master switch for initial-state QCD and QED radiation.
  fPythia->SetMSTP(71,0);             // (Default=2) master switch for initial-state QCD and QED radiation.
  fPythia->SetMDME(206,1,0);          //swicht off W decay leptonic modes
  fPythia->SetMDME(207,1,0); 
  fPythia->SetMDME(208,1,0); 
  fPythia->SetMDME(192,1,0);          //swicht off W decay to top
  fPythia->SetMDME(196,1,0); 
  fPythia->SetMDME(200,1,0); 

  fPythia->SetPMAS(24,1,kMw); //mass of the W boson (pythia=80.450 // genie=80.385)

  GetParam( "GLRESAtomic-Is-NLO", fIsNLO );
  GetParam( "HEDIS-Q2Grid-Min", fQ2PDFmin );

}
//____________________________________________________________________________


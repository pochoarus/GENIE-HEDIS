#include "Physics/HEDIS/EventGen/HEDISGenerator.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <TPythia6.h>

using namespace genie;
using namespace genie::constants;

// the actual PYTHIA call
extern "C" {
  double pyangl_( double *,  double * );
  void   pykfdi_( int *,  int *, int *, int * );
  void   pyzdis_( int *,  int *, double *, double * );
  void   pyrobo_( int *,  int *, double *, double *, double *, double *, double * );
  void   pydecy_( int * );
  void   py2ent_( int *,  int *, int *, double * );
} 

//___________________________________________________________________________
HEDISGenerator::HEDISGenerator() :
HadronicSystemGenerator("genie::HEDISGenerator")
{
  this->Initialize();
}
//___________________________________________________________________________
HEDISGenerator::HEDISGenerator(string config) :
HadronicSystemGenerator("genie::HEDISGenerator", config)
{
  this->Initialize();
}
//___________________________________________________________________________
HEDISGenerator::~HEDISGenerator()
{

}
//____________________________________________________________________________                                                                                        
void HEDISGenerator::Initialize(void) const
{
  fPythia = TPythia6::Instance();
  fPythia->SetMSTJ(93, 1);   //light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.
  fPythia->SetPARJ(22, 2);   //(Default=1) cut-off on decay length for a particle that is allowed to decay according to MSTJ(21) and the MDCY value
  fPythia->SetPARJ(32, 0.5); //(Default=1. GeV) is, with quark masses added, used to define the minimum allowable energy of a colour-singlet parton system.
  fPythia->SetPARJ(71, 2.);  //(Default=10. mm) maximum average proper lifetime cτ for particles allowed to decay
  fPythia->SetMSTU(26, 1000000);  //(Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, 1000000);  //(Default=10) maximum number of errors that are printed

  // sync GENIE/PYTHIA6 seed number                                                                                                                                   
  RandomGen::Instance();
}
//___________________________________________________________________________
void HEDISGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- Add the target remnant
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the primary lepton
  this->AddPrimaryLepton(evrec);

  //-- Add an entry for the DIS Pre-Fragm. Hadronic State
  this->AddFinalHadronicSyst(evrec);

  //-- Add the fragmentation products
  this->AddFragmentationProducts(evrec);

}
//___________________________________________________________________________
void HEDISGenerator::AddPrimaryLepton(GHepRecord * evrec) const
{

  Interaction * interaction = evrec->Summary();

  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB

  // Look-up selected kinematics & other needed kinematical params
  double Q2  = interaction->Kine().Q2(true);
  double y   = interaction->Kine().y(true);
  double Ev  = p4v->E(); 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(fmax(0.,El*El-plp*plp-ml2)); // p(-|)
  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2 * M_PIl * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  //rotate from LAB (where neutrino is [0,0,E,E]) to LAB' (where neutrino is [px,py,pz,E])
  TVector3 unit_nudir = p4v->Vect().Unit(); 
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  TLorentzVector p4l(p3l,El);
  LOG("HEDISGenerator", pINFO) << "LEPTON @ LAB' =>   E = " << p4l.E() << " // p = " << p4l.Px() << " , "  << p4l.Py() << " , "  << p4l.Pz() << " // p = " << p4l.P();
  LOG("HEDISGenerator", pINFO) << "                         m = " << p4l.M();

  // Decay tau
  int pdgl = interaction->FSPrimLepton()->PdgCode();
  evrec->AddParticle(pdgl, kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, p4l, *(evrec->Probe()->X4()));
  evrec->Summary()->KinePtr()->SetFSLeptonP4(p4l);

}
//___________________________________________________________________________
void HEDISGenerator::AddFragmentationProducts(GHepRecord * evrec) const
{

  //-- Compute the hadronic system invariant mass
  TLorentzVector p4Had = this->Hadronic4pLAB(evrec);
  double W = p4Had.M();

  if( W < 2 ) {
    LOG("HEDISGenerator", pWARN) << "Low invariant mass, W = " << W << " GeV!!";
    evrec->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    return;
  }

  //-- Run the hadronization model and get the fragmentation products:
  //   A collection of ROOT TMCParticles (equiv. to a LUJETS record)

  TClonesArray * plist;
  plist = this->Hadronize(evrec,W);
  if(!plist) {
     LOG("HEDISGenerator", pWARN) << "Got an empty particle list. Hadronizer failed!";
     LOG("HEDISGenerator", pWARN) << "Quitting the current event generation thread";
     evrec->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Could not simulate the hadronic system");
     exception.SwitchOnFastForward();
     throw exception;
     return;
  }

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.

  int mom = evrec->FinalStateHadronicSystemPosition();
  assert(mom!=-1);
 
  TMCParticle * p = 0;
  TIter particle_iter(plist);


  const TLorentzVector & vtx = *(evrec->Probe()->X4());

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());


  while( (p = (TMCParticle *) particle_iter.Next()) ) {

    int pdgc = p->GetKF();
    int ks   = p->GetKS();

    // The fragmentation products are generated in the hadronic CM frame
    // where the z>0 axis is the \vec{phad} direction. For each particle 
    // returned by the hadronizer:
    // - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
    // - rotate its 3-momentum from LAB' to LAB

    TLorentzVector p4o(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());
    p4o.Boost(beta); 
    TVector3 p3 = p4o.Vect();
    p3.RotateUz(unitvq); 
    TLorentzVector p4(p3,p4o.Energy());

    double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
    if ( ks==1 && p4.E() < massPDG ) {
      LOG("HEDISGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("HEDISGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
      LOG("HEDISGenerator", pWARN) << "E = " << p4.E() << " // |p| = " << TMath::Sqrt(p4.P()); 
      LOG("HEDISGenerator", pWARN) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
      LOG("HEDISGenerator", pWARN) << "m    = " << p4.M() << " // mpdg = " << massPDG;
      p4.SetXYZT(0,0,0,massPDG);
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks==1) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

    int im  = mom + 1 + p->GetParent();
    int ifc = (p->GetFirstChild() <= -1) ? -1 : mom + 1 + p->GetFirstChild();
    int ilc = (p->GetLastChild()  <= -1) ? -1 : mom + 1 + p->GetLastChild();

    double lightspeed = 299792458e3; //c in mm/s. Used for time in PYTHIA t[s]=t_pythia[mm]/c[mm/s]
    double vx = vtx.X() + p->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
    double vy = vtx.Y() + p->GetVy()*1e12;
    double vz = vtx.Z() + p->GetVz()*1e12;
    double vt = vtx.T() + p->GetTime()/lightspeed;
    TLorentzVector pos( vx, vy, vz, vt );

    LOG("HEDISGenerator", pDEBUG) << pdgc << "  " << ist << "  " << im << "  " << ifc;

    evrec->AddParticle( pdgc, ist, im,-1, ifc, ilc, p4, pos );

  } // fragmentation-products-iterator

  plist->Delete();
  delete plist;

}

//____________________________________________________________________________
TClonesArray * HEDISGenerator::Hadronize(GHepRecord * evrec, double W) const
{
  LOG("HEDISGenerator", pDEBUG) << "Running HEPYTHIA hadronizer";

  // get kinematics / init-state / process-info

  Interaction * interaction = evrec->Summary();
  interaction->KinePtr()->SetW(W);

  const XclsTag &      xclstag    = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();
  const Target &       target     = init_state.Tgt();

  assert(target.HitQrkIsSet()); 

  bool isp        = pdg::IsProton(target.HitNucPdg());
  int  hit_quark  = target.HitQrkPdg();
  int  frag_quark = xclstag.FinalQuarkPdg();

  LOG("HEDISGenerator", pDEBUG) << "Hit nucleon pdgc = " << target.HitNucPdg() << ", W = " << W;
  LOG("HEDISGenerator", pDEBUG) << "Selected hit quark pdgc = " << hit_quark << " // Fragmentation quark = " << frag_quark;

  // check hit-nucleon assignment, input neutrino & interaction type
  RandomGen * rnd = RandomGen::Instance();

  if ( pdg::IsDQuark(hit_quark) ) {
    int diquark = 0;
    if (isp) diquark = kPdgUUDiquarkS1;
    else     diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
    double m_frag    = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_diquark = PDGLibrary::Instance()->Find(diquark)->Mass();
    if( W <= m_frag + m_diquark + fPythia->GetPARJ(32) ) {
      LOG("HEDISGenerator", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("HEDISGenerator", pWARN) << "frag_quark = " << frag_quark << "    -> m = " << m_frag;
      LOG("HEDISGenerator", pWARN) << "diquark    = " << diquark    << " -> m = "    << m_diquark;
      return 0;
    }
    fPythia->Py1ent( -1, frag_quark, (W*W - m_diquark*m_diquark + m_frag*m_frag)/2./W, 0., 0. ); //k(1,2) = 2
    if ( pdg::IsTQuark(frag_quark) ) { //decay the top quark
      int ip = 1;
      pydecy_(&ip);
    }
    fPythia->Py1ent( fPythia->GetN()+1,    diquark,  (W*W + m_diquark*m_diquark - m_frag*m_frag)/2./W, fPythia->GetPARU(1), 0. ); //k(2,2) = 1
  }
  else if ( pdg::IsUQuark(hit_quark) ) {
    int diquark = 0;
    if (isp) diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
    else     diquark = kPdgDDDiquarkS1;
    double m_frag    = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_diquark = PDGLibrary::Instance()->Find(diquark)->Mass();
    if( W <= m_frag + m_diquark + fPythia->GetPARJ(32) ) {
      LOG("HEDISGenerator", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("HEDISGenerator", pWARN) << "frag_quark = " << frag_quark << "    -> m = " << m_frag;
      LOG("HEDISGenerator", pWARN) << "diquark    = " << diquark    << " -> m = "    << m_diquark;
      return 0;
    }
    fPythia->Py1ent( -1, frag_quark, (W*W - m_diquark*m_diquark + m_frag*m_frag)/2./W, 0., 0. ); //k(1,2) = 2
    fPythia->Py1ent( fPythia->GetN()+1,    diquark,  (W*W + m_diquark*m_diquark - m_frag*m_frag)/2./W, fPythia->GetPARU(1), 0. ); //k(2,2) = 1
  }
  else {

    int rema_hit_quark = -hit_quark;

    double m_frag     = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_rema_hit = PDGLibrary::Instance()->Find(rema_hit_quark)->Mass();
    if (W <= m_frag + m_rema_hit + 0.9 + fPythia->GetPARJ(32) ) {
      LOG("HEDISGenerator", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("HEDISGenerator", pWARN) << " frag_quark     = " << frag_quark     << " -> m = " << m_frag;
      LOG("HEDISGenerator", pWARN) << " rema_hit_quark = " << rema_hit_quark << " -> m = " << m_rema_hit;
      return 0;
    }

    int hadron     = 0;
    int rema_quark = 0;
    int ntwoq = isp ? 2 : 1; //proton two ups & neutron two downs

    int counter = 0;
    while( counter<100 ) {

      while(hadron==0) {
        int valquark = int(1.+ntwoq/3.+rnd->RndHadro().Rndm());
        int diquark  = 0;
        if ( valquark==ntwoq ) diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
        else                   diquark = 1000*ntwoq+100*ntwoq+3;

        int idum;
        if ( rema_hit_quark>0 ) { //create a baryon (qqq)
          //to generate a new quark or diquark flavour and to combine it with an existing flavour to give a hadron.
          pykfdi_(&diquark,&rema_hit_quark,&idum,&hadron);
          rema_quark = valquark;
        }
        else {                    //create a meson (qqbar)
          pykfdi_(&valquark,&rema_hit_quark,&idum,&hadron);
          rema_quark = diquark;
        }
      }

      double m_hadron = PDGLibrary::Instance()->Find(hadron)->Mass();
      double m_rema   = PDGLibrary::Instance()->Find(rema_quark)->Mass();

      double pT  = 0.35 * TMath::Sqrt( -1*TMath::Log( rnd->RndHadro().Rndm() ) );
      double pT2 = TMath::Power(pT,2);
      double pr  = TMath::Power(m_hadron,2)+pT2;
      int kfl1 = 1;
      int kfl3 = 0;
      double z;
      //to generate the longitudinal scaling variable z in jet fragmentation, either according to the Lund symmetric fragmentation function, or according to a choice of other shapes.
      pyzdis_(&kfl1,&kfl3,&pr,&z);  

      int krema_quark = fPythia->Pycomp(rema_quark);
      if      ( krema_quark==90 )                m_rema -= 0.2;
      else if ( krema_quark>0 && krema_quark<7 ) m_rema -= 0.1;
      int khadron = fPythia->Pycomp(hadron);
      if      ( khadron==90 )            m_hadron -= 0.2;
      else if ( khadron>0 && khadron<7 ) m_hadron -= 0.1;

      double tm_hadron = ( TMath::Power(m_hadron,2) + pT2 ) / z / W;
      double E_hadron   = 0.5 * ( z*W + tm_hadron );
      double E_pz       = W - tm_hadron;
      double WT         = (1-z) * W * E_pz - pT2;
      
      if ( WT > TMath::Power(m_frag+m_rema+fPythia->GetPARJ(32),2) ) {

        WT = TMath::Sqrt( WT + pT2 );
        double tm_rema   = TMath::Power(m_rema,2) + pT2;
        double E_frag    = 0.5 * ( WT + ( TMath::Power(m_frag,2) - tm_rema)/WT );
        double E_rema    = 0.5 * ( WT + (-TMath::Power(m_frag,2) + tm_rema)/WT );
        double x_rema    = -1 * TMath::Sqrt( TMath::Power(E_rema,2) - tm_rema );
        double theta_rema = pyangl_(&x_rema,&pT);

        double phi = 2*kPi*rnd->RndHadro().Rndm();

        fPythia->Py1ent( -1, frag_quark, E_frag, 0.,         0. ); //k(1,2) = 2
        if (TMath::Abs(frag_quark) > 5 ) {
          int ip = 1;
          pydecy_(&ip);
        }
        fPythia->Py1ent( fPythia->GetN()+1, rema_quark, E_rema, theta_rema, phi ); //k(2,2) = 1

        int imin     = 0;
        int imax     = 0;
        double the  = 0.; double ph   = 0.;
        double dbex = 0.; double dbey = 0.; double dbez = (E_pz-(1-z)*W)/(E_pz+(1-z)*W);
        pyrobo_( &imin , &imax, &the, &ph, &dbex, &dbey , &dbez );
      
        double pz_hadron  = -0.5 * ( z*W - tm_hadron );
        double theta_hadron = pyangl_(&pz_hadron,&pT);
        fPythia->SetMSTU( 10, 1 ); //keep the mass value stored in P(I,5), whatever it is.
        fPythia->SetP( fPythia->GetN()+1, 5, m_hadron );
        fPythia->Py1ent( fPythia->GetN()+1, hadron, E_hadron, theta_hadron, phi + kPi );
        fPythia->SetMSTU( 10, 2 ); //find masses according to mass tables as usual.

        if ( fPythia->GetP(fPythia->GetN()-1,3)<0 && fPythia->GetP(fPythia->GetN(),3)<0 ) break;
        
        LOG("HEDISGenerator", pWARN) << "Not backward hadron or rema_quark";
        LOG("HEDISGenerator", pWARN) << "hadron     = " << hadron     << " -> Pz = " << fPythia->GetP(fPythia->GetN(),3) ;
        LOG("HEDISGenerator", pWARN) << "rema_quark = " << rema_quark << " -> Pz = " << fPythia->GetP(fPythia->GetN()-1,3) ;

      }
      else {
        LOG("HEDISGenerator", pWARN) << "Low WT value ... ";
        LOG("HEDISGenerator", pWARN) << "WT = " << TMath::Sqrt(WT) << " // m_frag = " << m_frag << " // m_rema = " << m_rema;
      }

      LOG("HEDISGenerator", pWARN) << "Hadronization paricles not suitable. Trying again... " << counter;
      counter++;
      if (counter==100) {
        LOG("HEDISGenerator", pWARN) << "Hadronization particles failed after " << counter << " iterations";
        return 0;
      }

    }

  }

  int imin     = 0;
  int imax     = 0;
  double dbex = 0.; double dbey = 0.; double dbez = 0;
  double phi   = -2*kPi*rnd->RndHadro().Rndm();
  double theta = 0.;
  pyrobo_( &imin , &imax, &theta, &phi, &dbex, &dbey , &dbez );

  double pT  = fpT * TMath::Sqrt( -1*TMath::Log( rnd->RndHadro().Rndm() ) );
  phi   = -1 * phi;
  theta = TMath::ATan(2.*pT/W);
  pyrobo_( &imin , &imax, &theta, &phi, &dbex, &dbey , &dbez );

  fPythia->Pyexec();
  
  if (fPromptPythiaList) fPythia->Pylist(3);

  // get LUJETS record
  fPythia->GetPrimaries();

  TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("TMCParticle", np);
  particle_list->SetOwner(true);

  bool isTop = false;
  register unsigned int i = 0;
  TMCParticle * particle = 0;
  TIter particle_iter(pythia_particles);
  while( (particle = (TMCParticle *) particle_iter.Next()) ) {
    
    //LOG("HEDISGenerator", pDEBUG) << "Adding final state particle pdgc = " << particle->GetKF() << " with status = " << particle->GetKS();

    if(particle->GetKS() == 1) {
      if( pdg::IsQuark(particle->GetKF()) || pdg::IsDiQuark(particle->GetKF()) ) {
        LOG("HEDISGenerator", pERROR) << "Hadronization failed! Bare quark/di-quarks appear in final state!";
        particle_list->Delete();
        delete particle_list;
        return 0;            
      }  
    }

    if ( i==0 && pdg::IsTQuark( TMath::Abs(particle->GetKF()) ) ) { isTop=true; continue; } //we dont safe the top but its daughters (to avoid problem with parent/daughter)

    // fix numbering scheme used for mother/daughter assignments
    if ( isTop ) {
      (particle->GetParent()==0) ? particle->SetParent(particle->GetParent() - 1) : particle->SetParent(particle->GetParent() - 2);
      particle->SetFirstChild (particle->GetFirstChild() - 2);
      particle->SetLastChild  (particle->GetLastChild()  - 2);
    }
    else  {
      particle->SetParent(particle->GetParent() - 1);
      particle->SetFirstChild (particle->GetFirstChild() - 1);
      particle->SetLastChild  (particle->GetLastChild()  - 1);
    }                                  
    // insert the particle in the list
    new ( (*particle_list)[i++] ) TMCParticle(*particle);
  }

  delete particle;
  pythia_particles->Clear("C");

  return particle_list;

}
//___________________________________________________________________________
void HEDISGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISGenerator::LoadConfig(void)
{

  GetParamDef("PromptPythiaList", fPromptPythiaList, false ) ;
  GetParamDef("pT", fpT, 0.44 ) ;

}

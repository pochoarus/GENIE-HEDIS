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

#include "Physics/HEDIS/EventGen/HEDISGenerator.h"
#include "Physics/Hadronization/LeptoHadronization.h"
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

using namespace genie;
using namespace genie::constants;

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

}
//___________________________________________________________________________
void HEDISGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- Add the target remnant
  this->AddTargetNucleusRemnant(evrec);
  GHepParticle * target = evrec -> TargetNucleus();
  if(target) evrec->Particle(evrec->RemnantNucleusPosition())->SetStatus(kIStFinalStateNuclearRemnant);

  //-- Add the primary lepton
  this->AddPrimaryLepton(evrec);

  //-- Add the fragmentation products
  this->AddFragmentationProducts(evrec);

}
//___________________________________________________________________________
void HEDISGenerator::AddPrimaryLepton(GHepRecord * evrec) const
{

  Interaction * interaction = evrec->Summary();

  // Neutrino 4p
  LongLorentzVector p4v( evrec->Probe()->P4() );
  LOG("HEDISGenerator", pINFO) << "NEUTRINO @ LAB' =>  E = " << p4v.E() << " //  m = " << p4v.M() << " // p = " << p4v.P();
  LOG("HEDISGenerator", pINFO) << "                  dir = " << p4v.Dx() << " , "  << p4v.Dy() << " , "  << p4v.Dz();

  // Look-up selected kinematics & other needed kinematical params
  long double Q2  = interaction->Kine().Q2(true);
  long double y   = interaction->Kine().y(true);
  long double Ev  = p4v.E(); 
  long double ml  = interaction->FSPrimLepton()->Mass();
  long double ml2 = powl(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  long double El  = (1-y)*Ev;
  long double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  long double plt = sqrtl(fmaxl(0.,El*El-plp*plp-ml2)); // p(-|)
  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  long double phi  = 2 * M_PIl * rnd->RndLep().Rndm();
  long double pltx = plt * cosl(phi);
  long double plty = plt * sinl(phi);

  // Lepton 4-momentum in the LAB frame
  LongLorentzVector p4llong( pltx, plty, plp, El );
  p4llong.Rotate(p4v);
  LOG("HEDISGenerator", pINFO) << "LEPTON     @ LAB' =>  E = " << p4llong.E() << " //  m = " << p4llong.M() << " // p = " << p4llong.P();
  LOG("HEDISGenerator", pINFO) << "                    dir = " << p4llong.Dx() << " , "  << p4llong.Dy() << " , "  << p4llong.Dz();
 
  // Translate from long double to double
  TLorentzVector p4l( (double)p4llong.Px(), (double)p4llong.Py(), (double)p4llong.Pz(), (double)p4llong.E() );

  // Add lepton to EventRecord
  int pdgl = interaction->FSPrimLepton()->PdgCode();
  evrec->AddParticle(pdgl, kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, p4l, *(evrec->Probe()->X4()));
  evrec->Summary()->KinePtr()->SetFSLeptonP4(p4l);

}
//___________________________________________________________________________
void HEDISGenerator::AddFragmentationProducts(GHepRecord * evrec) const
{

  // Compute kinematics of hadronic system with energy/momentum conservation
  LongLorentzVector p4v( evrec->Probe()->P4()                   );
  LongLorentzVector p4N( evrec->HitNucleon()->P4()              );
  LongLorentzVector p4l( evrec->FinalStatePrimaryLepton()->P4() );
  LongLorentzVector p4Hadlong( p4v.Px()+p4N.Px()-p4l.Px(), p4v.Py()+p4N.Py()-p4l.Py(), p4v.Pz()+p4N.Pz()-p4l.Pz(), p4v.E()+p4N.E()-p4l.E() );

  LOG("HEDISGenerator", pDEBUG) << "v [LAB']: " << p4v.E() << " // " << p4v.M2() << " // [ " << p4v.Dx() << " , " << p4v.Dy() << " , " << p4v.Dz() << " ]";
  LOG("HEDISGenerator", pDEBUG) << "N [LAB']: " << p4N.E() << " // " << p4N.M2() << " // [ " << p4N.Dx() << " , " << p4N.Dy() << " , " << p4N.Dz() << " ]";
  LOG("HEDISGenerator", pDEBUG) << "l [LAB']: " << p4l.E() << " // " << p4l.M2() << " // [ " << p4l.Dx() << " , " << p4l.Dy() << " , " << p4l.Dz() << " ]";
  LOG("HEDISGenerator", pDEBUG) << "H [LAB']: " << p4Hadlong.E() << " // " << p4Hadlong.M2() << " // [ " << p4Hadlong.Dx() << " , " << p4Hadlong.Dy() << " , " << p4Hadlong.Dz() << " ]";

  // Transferred energy computed with long double for precission issue with double
  double W = p4Hadlong.M();

  // Translate from long double to double
  const TLorentzVector & vtx = *( evrec->Probe()->X4());
  TLorentzVector p4Had( (double)p4Hadlong.Px(), (double)p4Hadlong.Py(), (double)p4Hadlong.Pz(), (double)p4Hadlong.E() );
  evrec->AddParticle(kPdgHadronicSyst, kIStDISPreFragmHadronicState, evrec->HitNucleonPosition(),-1,-1,-1, p4Had, vtx);
  

  Interaction * interaction = evrec->Summary();
  interaction->KinePtr()->SetHadSystP4(p4Had);
  interaction->KinePtr()->SetW(W);

  //-- Run the hadronization model and get the fragmentation products:
  //   A collection of ROOT TMCParticles (equiv. to a LUJETS record)
  TClonesArray * plist;
  plist = fHadronizationModel->Hadronize(interaction);
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

  //-- Take the hadronic system weight to handle cases that the hadronizer
  //   was asked to produce weighted events
  double wght = fHadronizationModel->Weight();

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.
  int mom = evrec->FinalStateHadronicSystemPosition();
  assert(mom!=-1);
 
  TMCParticle * p = 0;
  TIter particle_iter(plist);

  // Boost velocity HCM -> LAB
  long double beta = p4Hadlong.P()/p4Hadlong.E();

  while( (p = (TMCParticle *) particle_iter.Next()) ) {

    int pdgc = p->GetKF();
    int ks   = p->GetKS();

    LongLorentzVector p4long( p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy()  );
    p4long.BoostZ(beta);
    p4long.Rotate(p4Hadlong);

    // Translate from long double to double
    TLorentzVector p4( (double)p4long.Px(), (double)p4long.Py(), (double)p4long.Pz(), (double)p4long.E() );

    // Somtimes PYTHIA output particles with E smaller than its mass. This is wrong,
    // so we assume that the are at rest.
    double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
    if ( (ks==1 || ks==4) && p4.E() < massPDG ) {
      LOG("HEDISGenerator", pINFO) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("HEDISGenerator", pINFO) << "PDG = " << pdgc << " // State = " << ks;
      LOG("HEDISGenerator", pINFO) << "E = " << p4.E() << " // |p| = " << p4.P(); 
      LOG("HEDISGenerator", pINFO) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
      LOG("HEDISGenerator", pINFO) << "m    = " << p4.M() << " // mpdg = " << massPDG;
      p4.SetXYZT(0,0,0,massPDG);
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

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

  //-- Handle the case that the hadronizer produced weighted events and
  //   take into account that the current event might be already weighted
  evrec->SetWeight (wght * evrec->Weight());

  plist->Delete();
  delete plist;

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
  fHadronizationModel = 0;

  //-- Get the requested hadronization model
  fHadronizationModel = 
     dynamic_cast<const HadronizationModelI *> (this->SubAlg("Hadronizer"));
  assert(fHadronizationModel);

}

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
 @ November 08, 2019 - Alfonso Garcia
   Modified to generate the kinematics of outgoing lepton properly.
   Phys. Rev. D 22, 2122 – Published 1 November 1980

*/
//____________________________________________________________________________

#include "Physics/GlashowResonance/EventGen/GLRESAtomicInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
GLRESAtomicInteractionListGenerator::GLRESAtomicInteractionListGenerator() :
InteractionListGeneratorI("genie::GLRESAtomicInteractionListGenerator")
{

}
//___________________________________________________________________________
GLRESAtomicInteractionListGenerator::GLRESAtomicInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::GLRESAtomicInteractionListGenerator", config)
{

}
//___________________________________________________________________________
GLRESAtomicInteractionListGenerator::~GLRESAtomicInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
   GLRESAtomicInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nuebar + e- -> W- -> nuebar + e-
// nuebar + e- -> W- -> nuebar + mu-
// nuebar + e- -> W- -> nuebar + tau-
// nuebar + e- -> W- -> hadrons

  if(init_state.ProbePdg() != kPdgAntiNuE) {
     LOG("IntLst", pDEBUG) 
          << "Return *null* interaction list";
     return 0;
  }

  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);  

  ProcessInfo   proc_info(kScGlashowResonanceAtomic, kIntWeakCC);

  InteractionList * intlist = new InteractionList;

  Interaction * interaction = new Interaction(init_state, proc_info);
  XclsTag exclusive_tag;
  if      (fIsMu)  exclusive_tag.SetFinalLepton(kPdgMuon);
  else if (fIsTau) exclusive_tag.SetFinalLepton(kPdgTau);
  else if (fIsEle) exclusive_tag.SetFinalLepton(kPdgElectron);
  else if (fIsHad) exclusive_tag.SetFinalLepton(kPdgPiP);
  interaction->SetExclTag(exclusive_tag);
  intlist->push_back(interaction);  

  if(intlist->size() == 0) {
     LOG("IntLst", pERROR)
         << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }

  return intlist;

}
//___________________________________________________________________________
void GLRESAtomicInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESAtomicInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESAtomicInteractionListGenerator::LoadConfigData(void)
{

  GetParamDef("is-Mu",  fIsMu,  false ) ;
  GetParamDef("is-Tau", fIsTau, false ) ;
  GetParamDef("is-Ele", fIsEle, false ) ;
  GetParamDef("is-Had", fIsHad, false ) ;

}

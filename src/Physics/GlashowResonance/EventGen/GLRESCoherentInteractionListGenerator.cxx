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

#include "Physics/GlashowResonance/EventGen/GLRESCoherentInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
GLRESCoherentInteractionListGenerator::GLRESCoherentInteractionListGenerator() :
InteractionListGeneratorI("genie::GLRESCoherentInteractionListGenerator")
{

}
//___________________________________________________________________________
GLRESCoherentInteractionListGenerator::GLRESCoherentInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::GLRESCoherentInteractionListGenerator", config)
{

}
//___________________________________________________________________________
GLRESCoherentInteractionListGenerator::~GLRESCoherentInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
   GLRESCoherentInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nuebar   + N[e-]   -> W- -> nuebar   + e-/mu-/tau-/hadrons
// numubar  + N[mu-]  -> W- -> numubar  + e-/mu-/tau-/hadrons
// nutaubar + N[tau-] -> W- -> nutaubar + e-/mu-/tau-/hadrons
// nue      + N[e+]   -> W+ -> nue      + e+/mu+/tau+/hadrons
// numu     + N[mu+]  -> W+ -> numu     + e+/mu+/tau+/hadrons
// nutau    + N[tau+] -> W+ -> nutau    + e+/mu+/tau+/hadrons


  ProcessInfo   proc_info(kScGlashowResonanceCoh, kIntWeakCC);

  InteractionList * intlist = new InteractionList;

  Interaction * interaction = new Interaction(init_state, proc_info);
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
void GLRESCoherentInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESCoherentInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESCoherentInteractionListGenerator::LoadConfigData(void)
{

}

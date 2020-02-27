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

#include "Physics/GlashowResonance/EventGen/GLRESInelasticInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
GLRESInelasticInteractionListGenerator::GLRESInelasticInteractionListGenerator() :
InteractionListGeneratorI("genie::GLRESInelasticInteractionListGenerator")
{

}
//___________________________________________________________________________
GLRESInelasticInteractionListGenerator::GLRESInelasticInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::GLRESInelasticInteractionListGenerator", config)
{

}
//___________________________________________________________________________
GLRESInelasticInteractionListGenerator::~GLRESInelasticInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * 
   GLRESInelasticInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nuebar   + p/n[e-]   -> W- -> nuebar   + e-/mu-/tau-/hadrons
// numubar  + p/n[mu-]  -> W- -> numubar  + e-/mu-/tau-/hadrons
// nutaubar + p/n[tau-] -> W- -> nutaubar + e-/mu-/tau-/hadrons
// nue      + p/n[e+]   -> W+ -> nue      + e+/mu+/tau+/hadrons
// numu     + p/n[mu+]  -> W+ -> numu     + e+/mu+/tau+/hadrons
// nutau    + p/n[tau+] -> W+ -> nutau    + e+/mu+/tau+/hadrons

  int probepdg = init_state.ProbePdg();

  ProcessInfo   proc_info(kScGlashowResonanceInel, kIntWeakCC);

  InteractionList * intlist = new InteractionList;

  bool hasP = (init_state.Tgt().Z() > 0);
  bool hasN = (init_state.Tgt().N() > 0);

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };

  for(int inucl=0; inucl<2; inucl++) {
    int struck_nucleon = nuclpdg[inucl];
    if( (struck_nucleon == kPdgProton  && hasP) || (struck_nucleon == kPdgNeutron && hasN) ) {
      Interaction * interaction = new Interaction(init_state, proc_info);
      Target * target = interaction->InitStatePtr()->TgtPtr();
      target->SetHitNucPdg(struck_nucleon);
      if (fIsMu) {
        XclsTag exclusive_tag;
        exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgAntiMuon : kPdgMuon );
        interaction->SetExclTag(exclusive_tag);
      }
      else if (fIsTau) {
        XclsTag exclusive_tag;
        exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgAntiTau : kPdgTau );
        interaction->SetExclTag(exclusive_tag);
      }
      else if (fIsEle) {
        XclsTag exclusive_tag;
        exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgPositron : kPdgElectron );
        interaction->SetExclTag(exclusive_tag);
      }
      else if (fIsHad) {
        XclsTag exclusive_tag;
        exclusive_tag.SetFinalLepton( (probepdg>0) ? kPdgPiP : kPdgPiM );
        interaction->SetExclTag(exclusive_tag);
      }
      intlist->push_back(interaction);
    }
  }

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
void GLRESInelasticInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESInelasticInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void GLRESInelasticInteractionListGenerator::LoadConfigData(void)
{

  GetParamDef("is-Mu",  fIsMu,  false ) ;
  GetParamDef("is-Tau", fIsTau, false ) ;
  GetParamDef("is-Ele", fIsEle, false ) ;
  GetParamDef("is-Had", fIsHad, false ) ;

}

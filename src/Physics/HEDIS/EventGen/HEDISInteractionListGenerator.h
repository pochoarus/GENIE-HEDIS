#ifndef _HEDIS_INTERACTION_LIST_GENERATOR_H_
#define _HEDIS_INTERACTION_LIST_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/Interaction/HEDISChannel.h"

namespace genie {

class Interaction;

class HEDISInteractionListGenerator : public InteractionListGeneratorI {

public :
  HEDISInteractionListGenerator();
  HEDISInteractionListGenerator(string config);
 ~HEDISInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void AddFinalStateInfo(Interaction * interaction, HEDISChannel_t hedischan) const;
  void LoadConfigData(void);
  
  bool fIsCC;
  bool fIsNC;

};

}      // genie namespace

#endif // _HEDIS_INTERACTION_LIST_GENERATOR_H_

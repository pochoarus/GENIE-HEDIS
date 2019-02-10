#ifndef _GLASHOW_RESONANCE_INTERACTION_GENERATOR_H_
#define _GLASHOW_RESONANCE_INTERACTION_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class GLRESInteractionListGenerator : public InteractionListGeneratorI {

public :
  GLRESInteractionListGenerator();
  GLRESInteractionListGenerator(string config);
 ~GLRESInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace

#endif // _GLASHOW_RESONANCE_INTERACTION_GENERATOR_H_

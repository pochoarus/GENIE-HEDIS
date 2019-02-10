#ifndef _GLASHOW_RESONANCE_GENERATOR_H_
#define _GLASHOW_RESONANCE_GENERATOR_H_

#include <TPythia6.h>

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class GLRESGenerator : public EventRecordVisitorI {

public :
  GLRESGenerator();
  GLRESGenerator(string config);
 ~GLRESGenerator();

  // implement the EventRecordVisitorI interface
  void Initialize         (void)               const;
  void ProcessEventRecord (GHepRecord * evrec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);

  void Configure(string config);

private:

  void LongBoost                (long double bx, long double by, long double bz, long double &px, long double &py, long double &pz, long double &E) const;
  void LoadConfig(void);

  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class
};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_GENERATOR_H_

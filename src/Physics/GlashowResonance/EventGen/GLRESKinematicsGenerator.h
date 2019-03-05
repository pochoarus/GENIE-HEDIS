#ifndef _GLASHOW_RESONANCE_KINEMATICS_GENERATOR_H_
#define _GLASHOW_RESONANCE_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"

namespace genie {

class GLRESKinematicsGenerator : public KineGeneratorWithCache {

public :
  GLRESKinematicsGenerator();
  GLRESKinematicsGenerator(string config);
 ~GLRESKinematicsGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

public:

  //-- methods to load sub-algorithms and config data from the Registry
  void   LoadConfig (void);

  //-- overload KineGeneratorWithCache methods
  double ComputeMaxXSec (const Interaction * in) const;
  double Energy         (const Interaction * in) const;
};

}      // genie namespace
#endif // _GLASHOW_RESONANCE_KINEMATICS_GENERATOR_H_

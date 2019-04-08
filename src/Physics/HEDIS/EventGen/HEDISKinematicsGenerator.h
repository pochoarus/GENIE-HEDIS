#ifndef _HEDIS_KINEMATICS_GENERATOR_H_
#define _HEDIS_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Numerical/Spline.h"

#include <map>

namespace genie {

class HEDISKinematicsGenerator : public KineGeneratorWithCache {

public :
  HEDISKinematicsGenerator();
  HEDISKinematicsGenerator(string config);
  ~HEDISKinematicsGenerator();

  class HEDISMaxXsecSpline 
  {
  public:
    HEDISMaxXsecSpline() { }
    ~HEDISMaxXsecSpline() { }
    map< int, genie::Spline * > Spline;
  };

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

	string fMaxXsecDirName;
	
  mutable bool fMaxXsecIsAlreadyLoaded = false;
  mutable map<HEDISChannel_t, HEDISMaxXsecSpline> fspl_max;

  double fXmin;
  double fQ2min;

  double ComputeMaxXSec       (const Interaction * interaction) const;
  void   LoadMaxXsecFromAscii (void) const;

  void   LoadConfig           (void);

};

}      // genie namespace

#endif // _HEDIS_KINEMATICS_GENERATOR_H_

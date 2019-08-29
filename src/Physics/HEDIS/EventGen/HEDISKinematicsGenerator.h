//____________________________________________________________________________
/*!

\class    genie::HEDISKinematicsGenerator

\brief    Generates values for the kinematic variables describing HEDIS v
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

          Max Xsec are precomputed as the total xsec and they are stored in 
          ascii files. This saves a lot of computational time.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

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

  double ComputeMaxXSec       (const Interaction * interaction) const;
  void   LoadMaxXsecFromAscii (void) const;

  void   LoadConfig           (void);

  string fMaxXsecDirName;                ///< name  Max Xsec direcotry
  mutable bool fMaxXsecIsLoad = false;   ///< check Max Xsec spl are loaded
  mutable map<HEDISQrkChannel_t, HEDISMaxXsecSpline> 
                fspl_max;  ///< splines to store max xsec for each channel


  double fXmin;   ///< minimum value of x for which SF tables are computed
  double fQ2min;  ///< minimum value of Q2 for which SF tables are computed 

};

}      // genie namespace

#endif // _HEDIS_KINEMATICS_GENERATOR_H_

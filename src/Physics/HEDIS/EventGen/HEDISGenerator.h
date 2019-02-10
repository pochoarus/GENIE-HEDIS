#ifndef _HEDIS_GENERATOR_H_
#define _HEDIS_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

class TPythia6;

namespace genie {

class HadronizationModelI;

class HEDISGenerator : public HadronicSystemGenerator {

public :
  HEDISGenerator();
  HEDISGenerator(string config);
 ~HEDISGenerator();

  // implement the EventRecordVisitorI interface
  void Initialize        (void)               const;
  void ProcessEventRecord(GHepRecord * evrec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void AddPrimaryLepton         (GHepRecord * evrec) const;
  void AddFragmentationProducts (GHepRecord * evrec) const;
  TClonesArray * Hadronize      (GHepRecord * evrec, double W) const;

  void LoadConfig (void);

  mutable TPythia6 *             fPythia;     

  bool   fPromptPythiaList;
  double fpT;

};

}      // genie namespace

#endif // _HEDIS_HADRONIC_SYSTEM_GENERATOR_H_

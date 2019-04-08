#ifndef _HEDIS_XSEC_H_
#define _HEDIS_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

#include <vector>

namespace genie {

class HEDISXSec : public XSecIntegratorI {

public:
  HEDISXSec();
  HEDISXSec(string config);
  virtual ~HEDISXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig (void);

  string fMaxXsecDirName;

  double fdlogy;
  double fdlogx;
  std::vector<double> fVy;
  std::vector<double> fVx;

};

}       // genie namespace
#endif  // _DIS_XSEC_H_

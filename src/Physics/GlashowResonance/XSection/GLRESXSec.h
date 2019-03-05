#ifndef _GLASHOW_RESONANCE_XSEC_H_
#define _GLASHOW_RESONANCE_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class GLRESXSec : public XSecIntegratorI {

public:
  GLRESXSec();
  GLRESXSec(string config);
  virtual ~GLRESXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _GLASHOW_RESONANCE_XSEC_H_

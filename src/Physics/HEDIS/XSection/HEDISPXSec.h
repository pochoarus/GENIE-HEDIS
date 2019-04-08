#ifndef _HEDIS_PXSEC_H_
#define _HEDIS_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HEDIS/XSection/HEDISFormFactors.h"

namespace genie {

  class XSecIntegratorI;

  class HEDISPXSec : public XSecAlgorithmI {

    public:
      HEDISPXSec();
      HEDISPXSec(string config);
      virtual ~HEDISPXSec();

      // XSecAlgorithmI interface implementation
      double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
      double Integral        (const Interaction * i) const;
      bool   ValidProcess    (const Interaction * i) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:
      void   LoadConfig (void);
      double ds_dxdy  (FF_xQ2 ff, double e, double mt, double ml2, double x, double y) const;

      const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator
    
      double fCKM[9];
      string fLHAPDFmember;
      bool fNLO;
      string fScheme;
      int fQrkThres;
      double fMassW;
      double fMassZ;
      double fRho;
      double fSin2ThW;
      int fNX;
      int fNQ2;
      double fXmin;
      double fQ2min;
      double fQ2max;
  };

}       // genie namespace

#endif  // _HEDIS_PXSEC_H_

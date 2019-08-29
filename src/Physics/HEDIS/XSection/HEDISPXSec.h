//____________________________________________________________________________
/*!

\class    genie::HEDISPXSec

\brief    Computes the double differential Cross Section for HEDIS. \n
          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_PXSEC_H_
#define _HEDIS_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HEDIS/XSection/HEDISStrucFunc.h"

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
      double ds_dxdy  (SF_xQ2 sf, double x, double y) const;

      const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator
    
      double fWmin;            ///< Minimum value of W

      string fSFname;          ///< name of the SF directory
      bool fSFNLO;             ///< SF is computed using LO or NLO
      int fSFNX;               ///< Number of bins in x for SF grid
      double fSFXmin;          ///< Minimum value of x for SF grid
      int fSFNQ2;              ///< Number of bins in Q2 for SF grid
      double fSFQ2min;         ///< Minimum value of Q2 for SF grid
      double fSFQ2max;         ///< Maximum value of Q2 for SF grid
      double fMassW;           ///< Mass of W boson used to compute SF
      double fMassZ;           ///< Mass of Z boson used to compute SF
      double fRho;             ///< EM correction for Weinberg angle  used to compute SF

  };

}       // genie namespace

#endif  // _HEDIS_PXSEC_H_

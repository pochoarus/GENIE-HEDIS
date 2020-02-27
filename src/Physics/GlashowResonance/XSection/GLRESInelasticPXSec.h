//____________________________________________________________________________
/*!

\class    genie::GLRESInelasticPXSec

\brief    Inelastic cross section at the Glashow resonance (nu + p/n -> W-).
          Is a concrete implementation of the XSecAlgorithmI interface. 

\ref      R.Gauld, Phys. Rev. D 100, 091301(R) (2019)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_INELASTIC_PXSEC_H_
#define _GLASHOW_RESONANCE_INELASTIC_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/GlashowResonance/XSection/GLRESBornPXSec.h"

namespace genie {

class XSecIntegratorI;

class GLRESInelasticPXSec : public XSecAlgorithmI {

public:
  GLRESInelasticPXSec ();
  GLRESInelasticPXSec (string config);
  virtual ~GLRESInelasticPXSec ();

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

  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  double fWmin;            ///< Minimum value of W

  double fQ2PDFmin;
  double fxPDFmin;

  GLRESBornPXSec * born;


};

}       // genie namespace

#endif  // _GLASHOW_RESONANCE_INELASTIC_PXSEC_H_

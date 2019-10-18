//____________________________________________________________________________
/*!

\class    genie::LeptoHadronization

\brief    Provides access to the LEPTO hadronization models. \n
          Is a concrete implementation of the HadronizationModelI interface.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF (Amsterdam)

\created  October 18, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LEPTO_HADRONIZATION_H_
#define _LEPTO_HADRONIZATION_H_

#include <TPythia6.h>

#include "Physics/Hadronization/HadronizationModelBase.h"

namespace genie {

class DecayModelI;
class LeptoHadronization : public HadronizationModelBase {

public:
  LeptoHadronization();
  LeptoHadronization(string config);
  virtual ~LeptoHadronization();

  //-- implement the HadronizationModelI interface
  void           Initialize       (void)                                  const;
  TClonesArray * Hadronize        (const Interaction*)                    const;
  double         Weight           (void)                                  const;
  PDGCodeList *  SelectParticles  (const Interaction*)                    const;
  TH1D *         MultiplicityProb (const Interaction*, Option_t* opt="")  const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig     (void);

  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class

  //-- configuration parameters
  int    fMaxIterHad;         // Maxmium number of iterations to look for a combination of hadrons
  double fPrimordialKT;       // Width of Gaussian distribution for the primordial transverse momentum kT of partons in the nucleon.
  double fRemnantPT;          // Width of Gaussian distribution in transverse momentum when a non-trivial target remnant is split into two particles
  double fMinESinglet;        // It is, with quark masses added, used to define the minimum allowable energy of a colour-singlet parton system.
  bool   fPromptPythiaList;   // Print the list of particles from PYTHIA

};

}         // genie namespace

#endif    // _LEPTO_HADRONIZATION__H_


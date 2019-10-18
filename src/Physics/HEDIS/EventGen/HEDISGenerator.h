//____________________________________________________________________________
/*!

\class    genie::HEDISGenerator

\brief    Generates the final state leptonic and hadronic system in v HEDIS 
          interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_GENERATOR_H_
#define _HEDIS_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

#include "math.h"

namespace genie {

class HadronizationModelI;

// This class has been created to perform several operations with long 
// doubles. It is needed in HEDIS because the kinematics of the outgoing
// particles can be so large that the on-shell feature is not fulfilled 
// many times due to the precission of double. 
class LongLorentzVector {

  public :
    LongLorentzVector(long double px, long double py, long double pz, long double e) : fPx(px),fPy(py),fPz(pz),fE(e) {}
    LongLorentzVector(TLorentzVector * p4) { 
      fPx = (long double)p4->Px();  
      fPy = (long double)p4->Py();  
      fPz = (long double)p4->Pz();  
      fE  = (long double)p4->E();  
    }
   ~LongLorentzVector() {}

    long double Px (void) { return fPx; }
    long double Py (void) { return fPy; }
    long double Pz (void) { return fPz; }
    long double E  (void) { return fE;  }
    long double P  (void) { return sqrtl(fPx*fPx+fPy*fPy+fPz*fPz);     }
    long double M  (void) { return sqrtl(fE*fE-fPx*fPx-fPy*fPy-fPz*fPz); }
    long double M2 (void) { return fE*fE-fPx*fPx-fPy*fPy-fPz*fPz; }
    long double Dx (void) { return fPx/sqrtl(fPx*fPx+fPy*fPy+fPz*fPz); }
    long double Dy (void) { return fPy/sqrtl(fPx*fPx+fPy*fPy+fPz*fPz); }
    long double Dz (void) { return fPz/sqrtl(fPx*fPx+fPy*fPy+fPz*fPz); }

    void Rotate    (LongLorentzVector axis) {
      long double up = axis.Dx()*axis.Dx() + axis.Dy()*axis.Dy();
      if (up) {
        up = sqrtl(up);
        long double pxaux = fPx,  pyaux = fPy,  pzaux = fPz;
        fPx = (axis.Dx()*axis.Dz()*pxaux - axis.Dy()*pyaux + axis.Dx()*up*pzaux)/up;
        fPy = (axis.Dy()*axis.Dz()*pxaux + axis.Dx()*pyaux + axis.Dy()*up*pzaux)/up;
        fPz = (axis.Dz()*axis.Dz()*pxaux -           pxaux + axis.Dz()*up*pzaux)/up;
      } 
      else if (axis.Dz() < 0.) { // phi=0  teta=pi
        fPx = -fPx; 
        fPz = -fPz; 
      }
    }    

    void Boost    (long double bz) {
      long double b2 = bz*bz;
      long double gamma = 1.0 / sqrtl(1.0 - b2);
      long double bp = bz*fPz;
      long double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
      fPz = fPz + gamma2*bp*bz + gamma*bz*fE;
      fE  = gamma*(fE + bp);    
    }    


  private :

    long double fPx;
    long double fPy;
    long double fPz;
    long double fE;
};

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

  void LoadConfig (void);

  const HadronizationModelI * fHadronizationModel;

};

}      // genie namespace

#endif // _HEDIS_HADRONIC_SYSTEM_GENERATOR_H_

//____________________________________________________________________________
/*!

\class    genie::GeomAnalyzerI

\brief    Defines the GENIE Geometry Analyzer Interface

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  July 13, 2005

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GEOMETRY_ANALYZER_I_H_
#define _GEOMETRY_ANALYZER_I_H_

#include <utility>
#include <vector>

class TLorentzVector;
class TVector3;
class TGeoMaterial;
class TGeoMixture;

namespace genie {

class PDGCodeList;
class PathLengthList;
class MatLengthList;

class GeomAnalyzerI {

public :

  virtual ~GeomAnalyzerI();

  // define the GeomAnalyzerI interface

  virtual const PDGCodeList &
            ListOfTargetNuclei (void) = 0;

  virtual const PathLengthList & 
            ComputeMaxPathLengths (void) = 0;
  virtual const PathLengthList &
            ComputePathLengths (
              const TLorentzVector & x, const TLorentzVector & p) = 0;
  virtual const TVector3 &
            GenerateVertex (
              const TLorentzVector & x, const TLorentzVector & p, int tgtpdg) = 0;

  virtual std::vector< std::pair<double, const TGeoMaterial*> >
            ComputeMatLengths (
              const TLorentzVector & x, const TLorentzVector & p) = 0;

  virtual int    GetTargetPdgCode        (const TGeoMaterial * const m) const = 0;
  virtual int    GetTargetPdgCode        (const TGeoMixture * const m, int ielement) const = 0;



protected:

  GeomAnalyzerI();
};

}      // genie namespace

#endif // _GEOMETRY_ANALYZER_I_H_
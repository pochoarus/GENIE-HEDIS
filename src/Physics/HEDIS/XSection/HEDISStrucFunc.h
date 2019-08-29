//____________________________________________________________________________
/*!

\class    genie::HEDISStrcuFunc

\brief    Singleton class to load Structure Functions used in HEDIS.

\author   Alfonso Garcia <alfonsog \at nikhef.nl>
          NIKHEF

\created  August 28, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HEDIS_STRUC_FUNC_H_
#define _HEDIS_STRUC_FUNC_H_

#include "Framework/Numerical/BLI2D.h"
#include "Framework/Interaction/HEDISChannel.h"

#include <map>
#include <vector>
#include <string>

using std::map;
using std::vector;
using std::string;
using std::to_string;

namespace genie {

  struct SF_xQ2 {

    double F1;
    double F2;
    double F3;

  };

  class HEDISStrucFunc
  {
    public:

      // ................................................................
      // HEDIS structure functions type
      //

      typedef enum HEDISStrucFuncType {  
        kMHTUndefined = 0,
        kSFT1, 
        kSFT2, 
        kSFT3, 
        kSFnumber, 
      } 
      HEDISStrucFuncType_t;

      // ................................................................
      // HEDIS form factor type
      //

      class HEDISStrucFuncTable 
      {
        public:
          HEDISStrucFuncTable() { }
          ~HEDISStrucFuncTable() { /* note: should delete the grids! */ }
          map< HEDISStrucFunc::HEDISStrucFuncType_t, genie::BLI2DNonUnifGrid * > Table;
      };

      // ................................................................

      static HEDISStrucFunc * Instance(string SFdirectory, bool NLO, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax);

      // method to return values of the SF for a particular channel in x and Q2
      SF_xQ2 EvalQrkSFLO  ( HEDISQrkChannel_t ch, double x, double Q2 );
      SF_xQ2 EvalNucSFLO  ( HEDISNucChannel_t ch, double x, double Q2 ); 
      SF_xQ2 EvalNucSFNLO ( HEDISNucChannel_t ch, double x, double Q2 );


    private:

      // Ctors & dtor
      HEDISStrucFunc(string SFdirectory, bool NLO, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax);
      HEDISStrucFunc(const HEDISStrucFunc &);
     ~HEDISStrucFunc();

      // Self
      static HEDISStrucFunc * fgInstance;

      // These map holds all SF tables (interaction channel is the key)
      map<HEDISQrkChannel_t, HEDISStrucFuncTable> fQrkSFLOTables;
      map<HEDISNucChannel_t, HEDISStrucFuncTable> fNucSFLOTables;
      map<HEDISNucChannel_t, HEDISStrucFuncTable> fNucSFNLOTables;

      // singleton cleaner
      struct Cleaner {
        void DummyMethodAndSilentCompiler(){}
          ~Cleaner(){
          if (HEDISStrucFunc::fgInstance !=0){
            delete HEDISStrucFunc::fgInstance;
            HEDISStrucFunc::fgInstance = 0;
          }
        }
      };
      friend struct Cleaner;
  };

} // genie namespace

#endif // _HEDIS_STRUC_FUNC_H_

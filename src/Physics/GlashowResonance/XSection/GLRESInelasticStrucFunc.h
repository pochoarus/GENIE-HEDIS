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

#ifndef _GLRES_STRUC_FUNC_H_
#define _GLRES_STRUC_FUNC_H_

#include "Framework/Numerical/Spline.h"

#include <map>

using std::map;

namespace genie {

  struct SF_x {
    double Fp1,Fm1;
    double Fp2,Fm2;
    double Fp3,Fm3;
  };

  class GLRESStrucFunc
  {
    public:

      // ................................................................
      // GLRES form factor type
      //

      class GLRESStrucFuncTable 
      {
        public:
          GLRESStrucFuncTable() { }
          ~GLRESStrucFuncTable() { /* note: should delete the grids! */ }
          map< int, genie::Spline * > Table;
      };

      // ................................................................

      static GLRESStrucFunc * Instance(void);

      double EvalSF  ( int hitnuc, int hitlep, double x ) { return fSFTables[hitnuc].Table[hitlep]->Evaluate(x); }

    private:

      // Ctors & dtor
      GLRESStrucFunc();
      GLRESStrucFunc(const GLRESStrucFunc &);
     ~GLRESStrucFunc();

      // Self
      static GLRESStrucFunc * fgInstance;

      // These map holds all SF tables (interaction channel is the key)
      map<int, GLRESStrucFuncTable> fSFTables;


      // singleton cleaner
      struct Cleaner {
        void DummyMethodAndSilentCompiler(){}
          ~Cleaner(){
          if (GLRESStrucFunc::fgInstance !=0){
            delete GLRESStrucFunc::fgInstance;
            GLRESStrucFunc::fgInstance = 0;
          }
        }
      };
      friend struct Cleaner;
  };

} // genie namespace

#endif // _GLRES_STRUC_FUNC_H_

#ifndef _HEDIS_FORM_FACTORS_H_
#define _HEDIS_FORM_FACTORS_H_

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

  struct FF_xQ2 {

    double F1;
    double F2;
    double F3;

  };

  class HEDISFormFactors
  {
    public:

      // ................................................................
      // HEDIS form factor type
      //

      typedef enum HEDISFormFactorType {  
        kMHTUndefined = 0,
        kFFT1, 
        kFFT2, 
        kFFT3, 
        kFFnumber, 
      } 
      HEDISFormFactorType_t;

      // ................................................................
      // HEDIS form factor type
      //

      class HEDISFormFactorTable 
      {
        public:
          HEDISFormFactorTable() { }
          ~HEDISFormFactorTable() { /* note: should delete the grids! */ }
          map< HEDISFormFactors::HEDISFormFactorType_t, genie::BLI2DNonUnifGrid * > Table;
      };

      // ................................................................

      static HEDISFormFactors * Instance(string LHAPDFmember, bool NLO, string Scheme, int QrkThrs, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax, const double CKM[9], double MZ, double MW, double Sin2ThW);

      FF_xQ2 EvalFFQrkLO   ( HEDISChannel_t ch, double x, double Q2 );
      FF_xQ2 EvalNuclFFLO  ( HEDISNuclChannel_t ch, double x, double Q2 ); 
      FF_xQ2 EvalNuclFFNLO ( HEDISNuclChannel_t ch, double x, double Q2 );


    private:

      // Ctors & dtor
      HEDISFormFactors(string LHAPDFmember, bool NLO, string Scheme, int QrkThrs, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax, const double CKM[9], double MZ, double MW, double Sin2ThW);
      HEDISFormFactors(const HEDISFormFactors &);
     ~HEDISFormFactors();

      // Self
      static HEDISFormFactors * fgInstance;

      // Load available form factor tables.
      void LoadFormfactors             ( HEDISChannel_t     ch );
      void LoadNuclFormfactors         ( HEDISNuclChannel_t ch );
      void CreateFormFactorFile        ( HEDISChannel_t     ch, string filename );
      void CreateLONuclFormFactorFile  ( HEDISNuclChannel_t ch, string filename );
      void CreateNLONuclFormFactorFile ( HEDISNuclChannel_t ch, string filename );
      BLI2DNonUnifGrid * ReadFormFactorFile ( string filename, HEDISFormFactorType_t fftype );

      double fQrkThrs;
      double fVud2; double fVus2; double fVub2;
      double fVcd2; double fVcs2; double fVcb2;
      double fVtd2; double fVts2; double fVtb2;
      double fSin2ThW;

      string fLHAPDFmember;
      map<int, double> mPDFQrk;

      double xPDFmin; double Q2PDFmin; double Q2PDFmax;

      vector<double> ff_logx_array;
      vector<double> ff_logq2_array;

      string fFormFactorsDir;

      bool APFELIsAlreadyInit = false;

      // This map holds all known tensor tables (target PDG code is the key)
      map<HEDISChannel_t, HEDISFormFactorTable> fFormFactorsTables;
      map<HEDISNuclChannel_t, HEDISFormFactorTable> fNLONuclFormFactorsTables;
      map<HEDISNuclChannel_t, HEDISFormFactorTable> fLONuclFormFactorsTables;

      // singleton cleaner
      struct Cleaner {
        void DummyMethodAndSilentCompiler(){}
          ~Cleaner(){
          if (HEDISFormFactors::fgInstance !=0){
            delete HEDISFormFactors::fgInstance;
            HEDISFormFactors::fgInstance = 0;
          }
        }
      };
      friend struct Cleaner;
  };

} // genie namespace

#endif // _HEDIS_FORM_FACTORS_H_

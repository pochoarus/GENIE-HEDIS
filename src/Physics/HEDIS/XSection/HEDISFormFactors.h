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

      static HEDISFormFactors * Instance(bool NLO, string LHAPDFmember, int NX, int NQ2, const double CKM[9]);

      // method to return whether the kinematics are valid for the PDF
      bool ValidKinematicsPDF (double x, double Q2) { return x>xPDFmin && Q2>Q2PDFmin && Q2<Q2PDFmax; }
      double GetXPDFmin       (void)                { return xPDFmin; }
      double GetQ2PDFmin      (void)                { return Q2PDFmin; }

      FF_xQ2 EvalFFQrkLO   ( HEDISChannel_t ch, double x, double Q2 );
      FF_xQ2 EvalNuclFFLO  ( HEDISNuclChannel_t ch, double x, double Q2 ); 
      FF_xQ2 EvalNuclFFNLO ( HEDISNuclChannel_t ch, double x, double Q2 );


    private:

      // Ctors & dtor
      HEDISFormFactors(bool NLO, string LHAPDFmember, int NX, int NQ2, const double CKM[9]);
      HEDISFormFactors(const HEDISFormFactors &);
     ~HEDISFormFactors();

      // Self
      static HEDISFormFactors * fgInstance;

      // Load available form factor tables.
      void InitQCDNUM                  ( void );
      void LoadFormfactors             ( HEDISChannel_t     ch );
      void LoadNuclFormfactors         ( HEDISNuclChannel_t ch );
      void CreateFormFactorFile        ( HEDISChannel_t     ch, string filename );
      void CreateLONuclFormFactorFile  ( HEDISNuclChannel_t ch, string filename );
      void CreateNLONuclFormFactorFile ( HEDISNuclChannel_t ch, string filename );
      BLI2DNonUnifGrid * ReadFormFactorFile ( string filename, HEDISFormFactorType_t fftype );

      vector<double> ff_logx_array;
      vector<double> ff_logq2_array;

      bool QCDNUMIsAlreadyInit = false;

      double xPDFmin; double Q2PDFmin; double Q2PDFmax;
      
      map<int, double> mPDFQrk;

      double fVud2; double fVus2; double fVub2;
      double fVcd2; double fVcs2; double fVcb2;
      double fVtd2; double fVts2; double fVtb2;

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

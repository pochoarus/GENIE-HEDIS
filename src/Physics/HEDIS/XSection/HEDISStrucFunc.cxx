//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/XSection/HEDISStrucFunc.h"
#include "Framework/Messenger/Messenger.h"

#include <TSystem.h>
#include <TMath.h>

#include <iostream>
#include <fstream>

using namespace genie;

//_________________________________________________________________________
HEDISStrucFunc * HEDISStrucFunc::fgInstance = 0;
//_________________________________________________________________________
HEDISStrucFunc::HEDISStrucFunc(string SFname, bool NLO, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax)
{

  // Check that the directory where SF tables are stored exists
  LOG("HEDISStrucFunc", pINFO) << "SF stored in following directory: " << SFname;
  if ( gSystem->mkdir(SFname.c_str())==0 ) {
    LOG("HEDISStrucFunc", pERROR) << "Directory doesnt exists.";
    LOG("HEDISStrucFunc", pERROR) << "HEDIS package requires precomputation of SF using gMakeStrucFunc";        
    assert(0);
  }
 
  // define arrays to fill from data files
  vector<double> sf_x_array;
  vector<double> sf_q2_array;
  double dlogq2 = TMath::Abs( TMath::Log10(Q2GRIDmin)-TMath::Log10(Q2GRIDmax) ) / NQ2;
  double dlogx  = TMath::Abs( TMath::Log10(xGRIDmin)-TMath::Log10(1.) ) / NX;
  LOG("HEDISStrucFunc", pINFO) << "Grid x,Q2 :" << NX << " , " << NQ2;
  for ( double logq2 = TMath::Log10(Q2GRIDmin); logq2<TMath::Log10(Q2GRIDmax); logq2+= dlogq2 ) {
    double q2 = TMath::Power( 10, logq2 + 0.5*dlogq2 );
    if (q2>Q2GRIDmax) continue;
    sf_q2_array.push_back(q2);
    LOG("HEDISStrucFunc", pDEBUG) << "q2: " << sf_q2_array.back();
  }
  for ( double logx = TMath::Log10(xGRIDmin); logx<TMath::Log10(1.); logx+= dlogx ) {
    double x = TMath::Power( 10, logx + 0.5*dlogx );
    if ( x>1. ) continue;
    sf_x_array.push_back(x);
    LOG("HEDISStrucFunc", pDEBUG) << "x: " << sf_x_array.back();
  }

  // Change to variables that are suitable for BLI2DNonUnifGrid
  int nx = sf_q2_array.size();
  int ny = sf_x_array.size();
  double x[nx];
  double y[ny];
  double z[nx*ny];
  for (int i=0; i<nx; i++) x[i] = sf_q2_array[i];
  for (int j=0; j<ny; j++) y[j] = sf_x_array[j];

  // Load structure functions for each quark at LO
  for( int ch=1; ch<kHEDISQrk_numofchannels; ch++ ) {
    string sfFile = SFname + "/QrkSF_LO_" + HEDISChannel::AsString((HEDISQrkChannel_t)ch) + ".dat";
    // Make sure data files are available
    LOG("HEDISStrucFunc", pINFO) << "Checking if file " << sfFile << " exists...";        
    if ( gSystem->AccessPathName( sfFile.c_str()) ) {
      LOG("HEDISStrucFunc", pERROR) << "File doesnt exist";        
      LOG("HEDISStrucFunc", pERROR) << "HEDIS package requires precomputation of SF using gMakeStrucFunc";        
      assert(0);
    }
    std::ifstream sf_stream(sfFile.c_str(), std::ios::in);
    // Loop over F1,F2,F3
    for(int sf = 1; sf <= kSFnumber; ++sf) {
      // Loop over x/Q2 bins
      for ( int ij=0; ij<nx*ny; ij++ ) sf_stream >> z[ij];
      // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
      fQrkSFLOTables[(HEDISQrkChannel_t)ch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
    }
  }
  if (NLO) {
    for( int ch=1; ch<kHEDISNuc_numofchannels; ch++ ) {
      // Load structure functions for each nucleon at NLO
      string sfFile = SFname + "/NucSF_NLO_" + HEDISChannel::AsString((HEDISNucChannel_t)ch) + ".dat";
      // Make sure data files are available
      LOG("HEDISStrucFunc", pINFO) << "Checking if file " << sfFile << " exists...";        
      if ( gSystem->AccessPathName( sfFile.c_str()) ) {
        LOG("HEDISStrucFunc", pERROR) << "File doesnt exist";        
		    LOG("HEDISStrucFunc", pERROR) << "HEDIS package requires precomputation of SF using gMakeStrucFunc";        
        assert(0);
      }
      std::ifstream sf_stream(sfFile.c_str(), std::ios::in);
      // Loop over F1,F2,F3
      for(int sf = 1; sf <= kSFnumber; ++sf) {
        // Loop over x/Q2 bins
        for ( int ij=0; ij<nx*ny; ij++ ) sf_stream >> z[ij];
        // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
        fNucSFNLOTables[(HEDISNucChannel_t)ch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
      }        

      //compute structure functions for each nucleon at LO using quark grids
      LOG("HEDISStrucFunc", pDEBUG) << "Creating LO " << sfFile;        
      int frstqrkch = HEDISChannel::GetFirstHEDISQrkChannel((HEDISNucChannel_t)ch);
      int lastqrkch = HEDISChannel::GetLastHEDISQrkChannel((HEDISNucChannel_t)ch);
      // Loop over F1,F2,F3
      for(int sf = 1; sf <= kSFnumber; ++sf) {
        int ij = 0;
        // Loop over Q2 bins
        for (int i=0; i<nx; i++) {
          // Loop over x bins
          for (int j=0; j<ny; j++) {
            double sum = 0.;
            // NucSF = sum_qrks QrkSF
            for( int qch=frstqrkch; qch<=lastqrkch; qch++ ) sum += fQrkSFLOTables[(HEDISQrkChannel_t)qch].Table[(HEDISStrucFuncType_t)sf]->Evaluate(x[i],y[j]);
            z[ij] = sum;
            ij++;
          }
        }
        // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
        fNucSFLOTables[(HEDISNucChannel_t)ch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
      }
    }
  }

  fgInstance = 0;

}
//_________________________________________________________________________
HEDISStrucFunc::~HEDISStrucFunc()
{

}
//_________________________________________________________________________
HEDISStrucFunc * HEDISStrucFunc::Instance(string SFname, bool NLO, int NX, double xGRIDmin, int NQ2, double Q2GRIDmin, double Q2GRIDmax)
{
  if(fgInstance == 0) {
    LOG("HEDISStrucFunc", pINFO) << "Late initialization";
    static HEDISStrucFunc::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new HEDISStrucFunc(SFname,NLO,NX,xGRIDmin,NQ2,Q2GRIDmin,Q2GRIDmax);
  }  
  return fgInstance;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalQrkSFLO( HEDISQrkChannel_t ch, double x, double Q2 ) 
{
  SF_xQ2 sf;
  sf.F1 = fQrkSFLOTables[ch].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fQrkSFLOTables[ch].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fQrkSFLOTables[ch].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalNucSFLO( HEDISNucChannel_t ch, double x, double Q2 ) 
{
  SF_xQ2 sf;
  sf.F1 = fNucSFLOTables[ch].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fNucSFLOTables[ch].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fNucSFLOTables[ch].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalNucSFNLO( HEDISNucChannel_t ch, double x, double Q2 ) 
{
  SF_xQ2 sf;
  sf.F1 = fNucSFNLOTables[ch].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fNucSFNLOTables[ch].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fNucSFNLOTables[ch].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
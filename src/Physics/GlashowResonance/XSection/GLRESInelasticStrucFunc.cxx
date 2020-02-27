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

#include "Physics/GlashowResonance/XSection/GLRESInelasticStrucFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "TSystem.h"

using namespace genie;

//_________________________________________________________________________
GLRESStrucFunc * GLRESStrucFunc::fgInstance = 0;
//_________________________________________________________________________
GLRESStrucFunc::GLRESStrucFunc()
{

  int nucs[2] = { kPdgProton, kPdgNeutron };
  int pdgs[6] = { kPdgNuE, kPdgAntiNuE, kPdgNuMu, kPdgAntiNuMu, kPdgNuTau, kPdgAntiNuTau };
  
  for (int k=0; k<2; k++) {
    for(int j=0; j<6; j++) {
      string SFname = string(gSystem->Getenv("GENIE")) + "/data/evgen/glres/PhotonSF_hitnuc"+std::to_string(nucs[k])+"_hitlep"+std::to_string(pdgs[j])+".dat";
      if ( gSystem->AccessPathName( SFname.c_str()) ) {
        LOG("GLRESStrucFunc", pWARN) << "File doesnt exist. SF table must be compute with gmkglressf.";        
        assert(0);
      }
      fSFTables[nucs[k]].Table[pdgs[j]] = new genie::Spline();
      fSFTables[nucs[k]].Table[pdgs[j]]->LoadFromAsciiFile(SFname);
    }        
  }

  fgInstance = 0;

}
//_________________________________________________________________________
GLRESStrucFunc::~GLRESStrucFunc()
{

}
//_________________________________________________________________________
GLRESStrucFunc * GLRESStrucFunc::Instance()
{
  if(fgInstance == 0) {
    LOG("GLRESStrucFunc", pINFO) << "Late initialization";
    static GLRESStrucFunc::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new GLRESStrucFunc();
  }  
  return fgInstance;
}
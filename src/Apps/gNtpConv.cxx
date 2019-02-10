#include <cassert>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Conventions/Units.h"

using namespace genie;

//input options (from command line arguments):
std::string     gOptInpFileName;         ///< input file name
std::string     gOptOutFileName;         ///< output file name

//consts
const int kNPmax = 500;

void GetCommandLineArgs(int argc, char ** argv)  {

  // Parse run options for this app
  CmdLnArgParser parser(argc,argv);

  // get input ROOT file (containing a GENIE GHEP event tree)
  if( parser.OptionExists('i') ) {
    LOG("gntpc", pINFO) << "Reading input filename";
    gOptInpFileName = parser.ArgAsString('i');
  } 
  else {
    LOG("gntpc", pFATAL) << "Unspecified input filename - Exiting";
    gAbortingInErr = true;
    exit(1);
  }
  // check input GENIE ROOT file
  bool inpok = !(gSystem->AccessPathName(gOptInpFileName.c_str()));
  if (!inpok) {
    LOG("gntpc", pFATAL) << "The input ROOT file [" << gOptInpFileName << "] is not accessible";
    gAbortingInErr = true;
    exit(2);
  }


  // get output file name
  if( parser.OptionExists('o') ) {
    LOG("gntpc", pINFO) << "Reading output filename";
    gOptOutFileName = parser.ArgAsString('o');
  } else {
    LOG("gntpc", pFATAL) << "Unspecified output filename - Exiting";
    gAbortingInErr = true;
    exit(1);
  }

  LOG("gntpc", pNOTICE) << "Input filename  = " << gOptInpFileName;
  LOG("gntpc", pNOTICE) << "Output filename = " << gOptOutFileName;

} 



//____________________________________________________________________________________
int main(int argc, char ** argv) {

  GetCommandLineArgs(argc, argv);

  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());

  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  //-- define the output rootracker tree branches

  // event info

  int         brEvtNum;                   // Event num.
  int         brEvtReac[2];         
  double      brEvtXSec;                  // Cross section for selected event (1E-38 cm2)
  double      brEvtDXSec;                 // Cross section for selected event kinematics (1E-38 cm2 /{K^n})
  double      brEvtWLim[2];  // W = [ kNeutronMass + kPhotontest + kASmallNum , TMath::Sqrt(s) - ml - kASmallNum ]
  double      brEvtQ2Lim[2]; // Q = -ml2 + 0.5*(s-M2)/s * { s + ml2 - W2 +- sqrt[ (s + ml2 - W2)^2 - 4 * s * ml2 ] } or Qmin = kMinQ2Limit
  double      brEvtXLim[2];  // x = [ ml2/(2*Ev_NRF*M_offshell) + kASmallNum , 1. - kASmallNum ]
  double      brEvtYLim[2];  // y = { 1 - ml2/(2*M*Ev*x) - ml2/(2*Ev2) +- sqrt[ ( 1 - ml2/(2*M*Ev*x) )^2 - b ] } / (2 + x*M/Ev) +- kASmallNum
  double      brEvtW;
  double      brEvtQ2;
  double      brEvtX;
  double      brEvtY;

  int         brProbeIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brProbeM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brProbePdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brProbeP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brTgtIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brTgtM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brTgtPdg[2];    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brTgtP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brTgtNucIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brTgtNucM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brTgtNucPdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brTgtNucP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brTgtQrkPdg[2];    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brOutQrkPdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brRecTgtIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brRecTgtM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brRecTgtPdg[2];    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brRecTgtP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brRecTgtNucPdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brLeptonIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brLeptonM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brLeptonPdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brLeptonP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brHadSystIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brHadSystM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brHadSystPdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brHadSystP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brHadBlbIdx;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brHadBlbM;    // Pdg codes (& generator specific codes for pseudoparticles)
  int         brHadBlbPdg;    // Pdg codes (& generator specific codes for pseudoparticles)
  double      brHadBlbP4[4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)

  int         brHadCluN;                  // Number of particles in particle array 
  double      brHadCluE;                  // Number of particles in particle array 
  int         brHadCluIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brHadCluM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brHadCluPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brHadCluP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brPreFragHadN;                  // Number of particles in particle array 
  double      brPreFragHadE;                  // Number of particles in particle array 
  int         brPreFragHadIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brPreFragHadM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brPreFragHadPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brPreFragHadP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFragHadN;                  // Number of particles in particle array 
  double      brFragHadE;                  // Number of particles in particle array 
  int         brFragHadIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFragHadM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFragHadPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFragHadP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFSIHadN;                  // Number of particles in particle array 
  double      brFSIHadE;                  // Number of particles in particle array 
  int         brFSIHadIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFSIHadM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFSIHadPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFSIHadRescat[kNPmax];                  // Number of particles in particle array 
  double      brFSIHadP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalHadPN;                  // Number of particles in particle array 
  double      brFinalHadPE;                  // Number of particles in particle array 
  int         brFinalHadPIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalHadPM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalHadPPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFinalHadPP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalHadNN;                  // Number of particles in particle array 
  double      brFinalHadNE;                  // Number of particles in particle array 
  int         brFinalHadNIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalHadNM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalHadNPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFinalHadNP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalHad0N;                  // Number of particles in particle array 
  double      brFinalHad0E;                  // Number of particles in particle array 
  int         brFinalHad0Idx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalHad0M     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalHad0Pdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFinalHad0P4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalMuN;                  // Number of particles in particle array 
  double      brFinalMuE;                  // Number of particles in particle array 
  int         brFinalMuIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalMuM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalMuPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFinalMuP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalEleN;                  // Number of particles in particle array 
  double      brFinalEleE;                  // Number of particles in particle array 
  int         brFinalEleIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalEleM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalElePdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFinalEleP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalOtherN;                  // Number of particles in particle array 
  double      brFinalOtherE;                  // Number of particles in particle array 
  int         brFinalOtherIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalOtherM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brFinalOtherPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  double      brFinalOtherP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  int         brFinalDeExtN;                  // Number of DeExticles in DeExticle array 
  double      brFinalDeExtE;                  // Number of DeExticles in DeExticle array 
  int         brFinalDeExtIdx   [kNPmax];     // Pdg codes (& generator specific codes for pseudoDeExticles)
  int         brFinalDeExtM     [kNPmax];     // Pdg codes (& generator specific codes for pseudoDeExticles)
  int         brFinalDeExtPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoDeExticles)
  double      brFinalDeExtP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)


  //-- open the output ROOT file
  TFile fout(gOptOutFileName.c_str(), "RECREATE");

  //-- create the output ROOT tree
  TTree * rootracker_tree = new TTree("gRooTracker","GENIE event tree rootracker format");
  rootracker_tree->Branch("EvtNum",          &brEvtNum,          "EvtNum/I");             
  rootracker_tree->Branch("EvtReac",          brEvtReac,         "EvtReac[2]/I");             
  rootracker_tree->Branch("EvtXSec",         &brEvtXSec,         "EvtXSec/D");            
  rootracker_tree->Branch("EvtDXSec",        &brEvtDXSec,        "EvtDXSec/D");           
  rootracker_tree->Branch("EvtWLim",          brEvtWLim,         "EvtWLim[2]/D");             
  rootracker_tree->Branch("EvtQ2Lim",         brEvtQ2Lim,        "EvtQ2Lim[2]/D");             
  rootracker_tree->Branch("EvtXLim",          brEvtXLim,         "EvtXLim[2]/D");             
  rootracker_tree->Branch("EvtYLim",          brEvtYLim,         "EvtYLim[2]/D");             
  rootracker_tree->Branch("EvtW",            &brEvtW,            "EvtW/D");             
  rootracker_tree->Branch("EvtQ2",           &brEvtQ2,           "EvtQ2/D");             
  rootracker_tree->Branch("EvtX",            &brEvtX,            "EvtX/D");             
  rootracker_tree->Branch("EvtY",            &brEvtY,            "EvtY/D");             
  rootracker_tree->Branch("ProbeIdx",        &brProbeIdx,        "ProbeIdx/I");              
  rootracker_tree->Branch("ProbeM",          &brProbeM,          "ProbeM/I");              
  rootracker_tree->Branch("ProbePdg",        &brProbePdg,        "ProbePdg/I");              
  rootracker_tree->Branch("ProbeP4",          brProbeP4,         "ProbeP4[4]/D");              
  rootracker_tree->Branch("TgtIdx",          &brTgtIdx,          "TgtIdx/I");              
  rootracker_tree->Branch("TgtM",            &brTgtM,            "TgtM/I");              
  rootracker_tree->Branch("TgtPdg",          &brTgtPdg,          "TgtPdg[2]/I");              
  rootracker_tree->Branch("TgtP4",            brTgtP4,           "TgtP4[4]/D");              
  rootracker_tree->Branch("TgtNucIdx",       &brTgtNucIdx,       "TgtNucIdx/I");              
  rootracker_tree->Branch("TgtNucM",         &brTgtNucM,         "TgtNucM/I");              
  rootracker_tree->Branch("TgtNucPdg",       &brTgtNucPdg,       "TgtNucPdg/I");              
  rootracker_tree->Branch("TgtNucP4",         brTgtNucP4,        "TgtNucP4[4]/D");              
  rootracker_tree->Branch("TgtQrkPdg",        brTgtQrkPdg,       "TgtQrkPdg[2]/I");              
  rootracker_tree->Branch("OutQrkPdg",        brOutQrkPdg,       "OutQrkPdg/I");              
  rootracker_tree->Branch("RecTgtIdx",       &brRecTgtIdx,       "RecTgtIdx/I");              
  rootracker_tree->Branch("RecTgtM",         &brRecTgtM,         "RecTgtM/I");              
  rootracker_tree->Branch("RecTgtPdg",        brRecTgtPdg,       "RecTgtPdg[2]/I");              
  rootracker_tree->Branch("RecTgtP4",         brRecTgtP4,        "RecTgtP4[4]/D");              
  rootracker_tree->Branch("RecTgtNucPdg",    &brRecTgtNucPdg,    "RecTgtNucPdg/I");              
  rootracker_tree->Branch("LeptonIdx",       &brLeptonIdx,       "LeptonIdx/I");              
  rootracker_tree->Branch("LeptonM",         &brLeptonM,         "LeptonM/I");              
  rootracker_tree->Branch("LeptonPdg",       &brLeptonPdg,       "LeptonPdg/I");              
  rootracker_tree->Branch("LeptonP4",         brLeptonP4,        "LeptonP4[4]/D");              
  rootracker_tree->Branch("HadSystIdx",      &brHadSystIdx,      "HadSystIdx/I");              
  rootracker_tree->Branch("HadSystM",        &brHadSystM,        "HadSystM/I");              
  rootracker_tree->Branch("HadSystPdg",      &brHadSystPdg,      "HadSystPdg/I");              
  rootracker_tree->Branch("HadSystP4",        brHadSystP4,       "HadSystP4[4]/D");              
  rootracker_tree->Branch("HadBlbIdx",       &brHadBlbIdx,       "HadBlbIdx/I");              
  rootracker_tree->Branch("HadBlbM",         &brHadBlbM,         "HadBlbM/I");              
  rootracker_tree->Branch("HadBlbPdg",       &brHadBlbPdg,       "HadBlbPdg/I");              
  rootracker_tree->Branch("HadBlbP4",         brHadBlbP4,        "HadBlbP4[4]/D");              
  rootracker_tree->Branch("HadCluN",         &brHadCluN,         "HadCluN/I");  
  rootracker_tree->Branch("HadCluE",         &brHadCluE,         "HadCluE/D");  
  rootracker_tree->Branch("HadCluIdx",        brHadCluIdx,       "HadCluIdx[HadCluN]/I");  
  rootracker_tree->Branch("HadCluM"  ,        brHadCluM,         "HadCluM[HadCluN]/I");  
  rootracker_tree->Branch("HadCluPdg",        brHadCluPdg,       "HadCluPdg[HadCluN]/I");  
  rootracker_tree->Branch("HadCluP4",         brHadCluP4,        "HadCluP4[HadCluN][4]/D"); 
  rootracker_tree->Branch("PreFragHadN",     &brPreFragHadN,     "PreFragHadN/I");  
  rootracker_tree->Branch("PreFragHadE",     &brPreFragHadE,     "PreFragHadE/D");  
  rootracker_tree->Branch("PreFragHadIdx",    brPreFragHadIdx,   "PreFragHadIdx[PreFragHadN]/I");  
  rootracker_tree->Branch("PreFragHadM"  ,    brPreFragHadM,     "PreFragHadM[PreFragHadN]/I");  
  rootracker_tree->Branch("PreFragHadPdg",    brPreFragHadPdg,   "PreFragHadPdg[PreFragHadN]/I");  
  rootracker_tree->Branch("PreFragHadP4",     brPreFragHadP4,    "PreFragHadP4[PreFragHadN][4]/D"); 
  rootracker_tree->Branch("FragHadN",        &brFragHadN,        "FragHadN/I");  
  rootracker_tree->Branch("FragHadE",        &brFragHadE,        "FragHadE/D");  
  rootracker_tree->Branch("FragHadIdx",       brFragHadIdx,      "FragHadIdx[FragHadN]/I");  
  rootracker_tree->Branch("FragHadM"  ,       brFragHadM,        "FragHadM[FragHadN]/I");  
  rootracker_tree->Branch("FragHadPdg",       brFragHadPdg,      "FragHadPdg[FragHadN]/I");  
  rootracker_tree->Branch("FragHadP4",        brFragHadP4,       "FragHadP4[FragHadN][4]/D"); 
  rootracker_tree->Branch("FSIHadN",         &brFSIHadN,         "FSIHadN/I");  
  rootracker_tree->Branch("FSIHadE",         &brFSIHadE,         "FSIHadE/D");  
  rootracker_tree->Branch("FSIHadIdx",        brFSIHadIdx,       "FSIHadIdx[FSIHadN]/I");  
  rootracker_tree->Branch("FSIHadM"  ,        brFSIHadM,         "FSIHadM[FSIHadN]/I");  
  rootracker_tree->Branch("FSIHadPdg",        brFSIHadPdg,       "FSIHadPdg[FSIHadN]/I");  
  rootracker_tree->Branch("FSIHadRescat",     brFSIHadRescat,    "FSIHadRescat[FSIHadN]/I");  
  rootracker_tree->Branch("FSIHadP4",         brFSIHadP4,        "FSIHadP4[FSIHadN][4]/D"); 
  rootracker_tree->Branch("FinalHadPN",      &brFinalHadPN,      "FinalHadPN/I");  
  rootracker_tree->Branch("FinalHadPE",      &brFinalHadPE,      "FinalHadPE/D");  
  rootracker_tree->Branch("FinalHadPIdx",     brFinalHadPIdx,    "FinalHadPIdx[FinalHadPN]/I");  
  rootracker_tree->Branch("FinalHadPM"  ,     brFinalHadPM,      "FinalHadPM[FinalHadPN]/I");  
  rootracker_tree->Branch("FinalHadPPdg",     brFinalHadPPdg,    "FinalHadPPdg[FinalHadPN]/I");  
  rootracker_tree->Branch("FinalHadPP4",      brFinalHadPP4,     "FinalHadPP4[FinalHadPN][4]/D"); 
  rootracker_tree->Branch("FinalHadNN",      &brFinalHadNN,      "FinalHadNN/I");  
  rootracker_tree->Branch("FinalHadNE",      &brFinalHadNE,      "FinalHadNE/D");  
  rootracker_tree->Branch("FinalHadNIdx",     brFinalHadNIdx,    "FinalHadNIdx[FinalHadNN]/I");  
  rootracker_tree->Branch("FinalHadNM"  ,     brFinalHadNM,      "FinalHadNM[FinalHadNN]/I");  
  rootracker_tree->Branch("FinalHadNPdg",     brFinalHadNPdg,    "FinalHadNPdg[FinalHadNN]/I");  
  rootracker_tree->Branch("FinalHadNP4",      brFinalHadNP4,     "FinalHadNP4[FinalHadNN][4]/D"); 
  rootracker_tree->Branch("FinalHad0N",      &brFinalHad0N,      "FinalHad0N/I");  
  rootracker_tree->Branch("FinalHad0E",      &brFinalHad0E,      "FinalHad0E/D");  
  rootracker_tree->Branch("FinalHad0Idx",     brFinalHad0Idx,    "FinalHad0Idx[FinalHad0N]/I");  
  rootracker_tree->Branch("FinalHad0M"  ,     brFinalHad0M,      "FinalHad0M[FinalHad0N]/I");  
  rootracker_tree->Branch("FinalHad0Pdg",     brFinalHad0Pdg,    "FinalHad0Pdg[FinalHad0N]/I");  
  rootracker_tree->Branch("FinalHad0P4",      brFinalHad0P4,     "FinalHad0P4[FinalHad0N][4]/D"); 
  rootracker_tree->Branch("FinalMuN",        &brFinalMuN,        "FinalMuN/I");  
  rootracker_tree->Branch("FinalMuE",        &brFinalMuE,        "FinalMuE/D");  
  rootracker_tree->Branch("FinalMuIdx",       brFinalMuIdx,      "FinalMuIdx[FinalMuN]/I");  
  rootracker_tree->Branch("FinalMuM"  ,       brFinalMuM,        "FinalMuM[FinalMuN]/I");  
  rootracker_tree->Branch("FinalMuPdg",       brFinalMuPdg,      "FinalMuPdg[FinalMuN]/I");  
  rootracker_tree->Branch("FinalMuP4",        brFinalMuP4,       "FinalMuP4[FinalMuN][4]/D"); 
  rootracker_tree->Branch("FinalEleN",        &brFinalEleN,      "FinalEleN/I");  
  rootracker_tree->Branch("FinalEleE",        &brFinalEleE,      "FinalEleE/D");  
  rootracker_tree->Branch("FinalEleIdx",       brFinalEleIdx,    "FinalEleIdx[FinalEleN]/I");  
  rootracker_tree->Branch("FinalEleM"  ,       brFinalEleM,      "FinalEleM[FinalEleN]/I");  
  rootracker_tree->Branch("FinalElePdg",       brFinalElePdg,    "FinalElePdg[FinalEleN]/I");  
  rootracker_tree->Branch("FinalEleP4",        brFinalEleP4,     "FinalEleP4[FinalEleN][4]/D"); 
  rootracker_tree->Branch("FinalOtherN",      &brFinalOtherN,    "FinalOtherN/I");  
  rootracker_tree->Branch("FinalOtherE",      &brFinalOtherE,    "FinalOtherE/D");  
  rootracker_tree->Branch("FinalOtherIdx",     brFinalOtherIdx,  "FinalOtherIdx[FinalOtherN]/I");  
  rootracker_tree->Branch("FinalOtherM"  ,     brFinalOtherM,    "FinalOtherM[FinalOtherN]/I");  
  rootracker_tree->Branch("FinalOtherPdg",     brFinalOtherPdg,  "FinalOtherPdg[FinalOtherN]/I");  
  rootracker_tree->Branch("FinalOtherP4",      brFinalOtherP4,   "FinalOtherP4[FinalOtherN][4]/D"); 
  rootracker_tree->Branch("FinalDeExtN",      &brFinalDeExtN,    "FinalDeExtN/I");  
  rootracker_tree->Branch("FinalDeExtE",      &brFinalDeExtE,    "FinalDeExtE/D");  
  rootracker_tree->Branch("FinalDeExtIdx",     brFinalDeExtIdx,  "FinalDeExtIdx[FinalDeExtN]/I");  
  rootracker_tree->Branch("FinalDeExtM"  ,     brFinalDeExtM,    "FinalDeExtM[FinalDeExtN]/I");  
  rootracker_tree->Branch("FinalDeExtPdg",     brFinalDeExtPdg,  "FinalDeExtPdg[FinalDeExtN]/I");  
  rootracker_tree->Branch("FinalDeExtP4",      brFinalDeExtP4,   "FinalDeExtP4[FinalDeExtN][4]/D"); 

  //-- open the input GENIE ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           gtree = 0;
  NtpMCTreeHeader * thdr  = 0;
  gtree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr  = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  gtree->SetBranchAddress("gmcrec", &mcrec);

  //-- figure out how many events to analyze
  Long64_t nmax = gtree->GetEntries();
  if (nmax<=0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return 0;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  //-- event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    gtree->GetEntry(iev);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);
    Interaction * interaction = event.Summary();

    //LOG("gntpc", pINFO) << rec_header;
    //LOG("gntpc", pINFO) << event;
    //LOG("gntpc", pINFO) << *interaction;

    //
    // clear output tree branches
    //
    brEvtNum      = -999;    
    brEvtReac[0]  = -999;
    brEvtReac[1]  = -999;
    brEvtXSec     = -999;
    brEvtDXSec    = -999;
    brEvtWLim[0]  = -999;
    brEvtWLim[1]  = -999;
    brEvtQ2Lim[0] = -999;
    brEvtQ2Lim[1] = -999;
    brEvtXLim[0]  = -999;
    brEvtXLim[1]  = -999;
    brEvtYLim[0]  = -999;
    brEvtYLim[1]  = -999;
    brEvtW  = -999;
    brEvtQ2 = -999;
    brEvtX  = -999;
    brEvtY  = -999;
    brProbeIdx     = -999;
    brProbeM       = -999;
    brProbePdg     = 0;
    brTgtIdx       = -999;
    brTgtM        = -999;
    brTgtPdg[0]    = -999; //Z
    brTgtPdg[1]    = -999; //A
    brRecTgtIdx    = -999;
    brRecTgtM     = -999;
    brRecTgtPdg[0] = -999; //Z
    brRecTgtPdg[1] = -999; //A
    brRecTgtNucPdg = 0;
    brTgtNucIdx    = -999;
    brTgtNucM      = -999;
    brTgtNucPdg    = 0;
    brTgtQrkPdg[0] = 0;
    brTgtQrkPdg[1] = -999; //isquarksea
    brOutQrkPdg    = -999;
    brLeptonIdx    = -999;
    brLeptonM      = -999;
    brLeptonPdg    = 0;
    brHadSystIdx   = -999;
    brHadSystM     = -999;
    brHadSystPdg   = 0;
    brHadBlbIdx    = -999;
    brHadBlbM      = -999;
    brHadBlbPdg    = 0;
    for (int k=0; k<4; k++) {
      brProbeP4[k]   = -9e-10;
      brTgtP4[k]     = -9e-10;
      brRecTgtP4[k]  = -9e-10;
      brTgtNucP4[k]  = -9e-10;
      brLeptonP4[k]  = -9e-10;
      brHadSystP4[k] = -9e-10;
      brHadBlbP4[k] = -9e-10;
    }
    brHadCluN = -999; 
    brHadCluE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brHadCluIdx   [i] = -999;  
       brHadCluM     [i] = -999;  
       brHadCluPdg   [i] = 0;  
       for(int k=0; k<4; k++) brHadCluP4 [i][k] = -9e-10;  
    }
    brPreFragHadN = -999; 
    brPreFragHadE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brPreFragHadIdx   [i] = -999;  
       brPreFragHadM     [i] = -999;  
       brPreFragHadPdg   [i] = 0;  
       for(int k=0; k<4; k++) brPreFragHadP4 [i][k] = -9e-10;  
    }
    brFragHadN = -999; 
    brFragHadE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFragHadIdx   [i] = -999;  
       brFragHadM     [i] = -999;  
       brFragHadPdg   [i] = 0;  
       for(int k=0; k<4; k++) brFragHadP4 [i][k] = -9e-10;  
    }
    brFSIHadN = -999; 
    brFSIHadE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFSIHadIdx   [i] = -999;  
       brFSIHadM     [i] = -999;  
       brFSIHadPdg   [i] = 0;  
       brFSIHadRescat[i] = 0;  
       for(int k=0; k<4; k++) brFSIHadP4 [i][k] = -9e-10;  
    }
    brFinalHadPN = -999; 
    brFinalHadPE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalHadPIdx   [i] = -999;  
       brFinalHadPM     [i] = -999;  
       brFinalHadPPdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalHadPP4 [i][k] = -9e-10;  
    }
    brFinalHadNN = -999; 
    brFinalHadNE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalHadNIdx   [i] = -999;  
       brFinalHadNM     [i] = -999;  
       brFinalHadNPdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalHadNP4 [i][k] = -9e-10;  
    }
    brFinalHad0N = -999; 
    brFinalHad0E = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalHad0Idx   [i] = -999;  
       brFinalHad0M     [i] = -999;  
       brFinalHad0Pdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalHad0P4 [i][k] = -9e-10;  
    }
    brFinalMuN = -999; 
    brFinalMuE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalMuIdx   [i] = -999;  
       brFinalMuM     [i] = -999;  
       brFinalMuPdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalMuP4 [i][k] = -9e-10;  
    }
    brFinalEleN = -999; 
    brFinalEleE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalEleIdx   [i] = -999;  
       brFinalEleM     [i] = -999;  
       brFinalElePdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalEleP4 [i][k] = -9e-10;  
    }
    brFinalOtherN = -999; 
    brFinalOtherE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalOtherIdx   [i] = -999;  
       brFinalOtherM     [i] = -999;  
       brFinalOtherPdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalOtherP4 [i][k] = -9e-10;  
    }
    brFinalDeExtN = -999; 
    brFinalDeExtE = -999; 
    for(int i=0; i<kNPmax; i++) {
       brFinalDeExtIdx   [i] = -999;  
       brFinalDeExtM     [i] = -999;  
       brFinalDeExtPdg   [i] = 0;  
       for(int k=0; k<4; k++) brFinalDeExtP4 [i][k] = -9e-10;  
    }

    //
    // copy current event info to output tree
    //
    brEvtNum    = (int) iev;    

    brEvtReac[0]  = interaction->ProcInfoPtr()->InteractionTypeId();
    brEvtReac[1]  = interaction->ProcInfoPtr()->ScatteringTypeId();

    brEvtXSec   = (1E+38/units::cm2) * event.XSec();    
    brEvtDXSec  = (1E+38/units::cm2) * event.DiffXSec();    
    
    brEvtWLim[0]  = interaction->PhaseSpacePtr()->WLim().min;
    brEvtWLim[1]  = interaction->PhaseSpacePtr()->WLim().max;
    brEvtQ2Lim[0] = interaction->PhaseSpacePtr()->Q2Lim().min;
    brEvtQ2Lim[1] = interaction->PhaseSpacePtr()->Q2Lim().max;
    brEvtXLim[0]  = interaction->PhaseSpacePtr()->XLim().min;
    brEvtXLim[1]  = interaction->PhaseSpacePtr()->XLim().max;
    brEvtYLim[0]  = interaction->PhaseSpacePtr()->YLim().min;
    brEvtYLim[1]  = interaction->PhaseSpacePtr()->YLim().max;

    brEvtW  = interaction->KinePtr()->W(true);
    brEvtQ2 = interaction->KinePtr()->Q2(true);
    brEvtX  = interaction->KinePtr()->x(true);
    brEvtY  = interaction->KinePtr()->y(true);

    brRecTgtNucPdg = interaction->RecoilNucleonPdg();
    brTgtQrkPdg[0] = interaction->InitStatePtr()->TgtPtr()->HitQrkPdg();
    brTgtQrkPdg[1] = interaction->InitStatePtr()->TgtPtr()->HitSeaQrk();
    brOutQrkPdg    = interaction->ExclTag().FinalQuarkPdg();

    int iHadClu=0;
    int iPreFragHad=0;
    int iFragHad=0;
    int iFSIHad=0;
    int iFinalHadP=0;
    int iFinalHadN=0;
    int iFinalHad0=0;
    int iFinalMu=0;
    int iFinalEle=0;
    int iFinalOther=0;
    int iFinalDeExt=0;
    double eneHadClu=0;
    double enePreFragHad=0;
    double eneFragHad=0;
    double eneFSIHad=0;
    double eneFinalHadP=0;
    double eneFinalHadN=0;
    double eneFinalHad0=0;
    double eneFinalMu=0;
    double eneFinalEle=0;
    double eneFinalOther=0;
    double eneFinalDeExt=0;
    GHepParticle * prt = 0;
    TIter event_iter(&event);
    while ( (prt = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
        assert(prt);
        
        if      ( event.ParticlePosition(prt)==event.ProbePosition() ) {
          brProbeIdx   = event.ParticlePosition(prt);
          brProbeM     = prt->FirstMother();
          brProbePdg   = prt->Pdg();
          brProbeP4[0] = prt->P4()->Px();
          brProbeP4[1] = prt->P4()->Py();
          brProbeP4[2] = prt->P4()->Pz();
          brProbeP4[3] = prt->P4()->E();
        }
        else if ( event.ParticlePosition(prt)==event.TargetNucleusPosition() ) {
          brTgtIdx    = event.ParticlePosition(prt);
          brTgtM      = prt->FirstMother();
          brTgtPdg[0] = prt->Z();
          brTgtPdg[1] = prt->A();
          brTgtP4[0]  = prt->P4()->Px();
          brTgtP4[1]  = prt->P4()->Py();
          brTgtP4[2]  = prt->P4()->Pz();
          brTgtP4[3]  = prt->P4()->E();
        }
        else if ( event.ParticlePosition(prt)==event.HitNucleonPosition() ) {
          brTgtNucIdx   = event.ParticlePosition(prt);
          brTgtNucM     = prt->FirstMother();
          brTgtNucPdg   = prt->Pdg();
          brTgtNucP4[0] = prt->P4()->Px();
          brTgtNucP4[1] = prt->P4()->Py();
          brTgtNucP4[2] = prt->P4()->Pz();
          brTgtNucP4[3] = prt->P4()->E();
        }
        else if ( event.ParticlePosition(prt)==event.FinalStatePrimaryLeptonPosition() ) {
          brLeptonIdx   = event.ParticlePosition(prt);
          brLeptonM     = prt->FirstMother();
          brLeptonPdg   = prt->Pdg();
          brLeptonP4[0] = prt->P4()->Px();
          brLeptonP4[1] = prt->P4()->Py();
          brLeptonP4[2] = prt->P4()->Pz();
          brLeptonP4[3] = prt->P4()->E();
        }
        else if ( event.ParticlePosition(prt)==event.FinalStateHadronicSystemPosition() ) {
          brHadSystIdx   = event.ParticlePosition(prt);
          brHadSystM     = prt->FirstMother();
          brHadSystPdg   = prt->Pdg();
          brHadSystP4[0] = prt->P4()->Px();
          brHadSystP4[1] = prt->P4()->Py();
          brHadSystP4[2] = prt->P4()->Pz();
          brHadSystP4[3] = prt->P4()->E();
        }
        else if ( prt->Status()==2 ) {
          brRecTgtIdx    = event.ParticlePosition(prt);
          brRecTgtM   = prt->FirstMother();
          brRecTgtPdg[0] = prt->Z();
          brRecTgtPdg[1] = prt->A();
          brRecTgtP4[0]  = prt->P4()->Px();
          brRecTgtP4[1]  = prt->P4()->Py();
          brRecTgtP4[2]  = prt->P4()->Pz();
          brRecTgtP4[3]  = prt->P4()->E();
        }
        else if ( prt->Status()==15 ) {
          brHadBlbIdx    = event.ParticlePosition(prt);
          brHadBlbM      = prt->FirstMother();
          brHadBlbPdg    = prt->Pdg();
          brHadBlbP4[0]  = prt->P4()->Px();
          brHadBlbP4[1]  = prt->P4()->Py();
          brHadBlbP4[2]  = prt->P4()->Pz();
          brHadBlbP4[3]  = prt->P4()->E();
        }
        else if ( prt->Status()==16 ) {
          brHadCluIdx   [iHadClu]    = event.ParticlePosition(prt);
          brHadCluM     [iHadClu]    = prt->FirstMother();
          brHadCluPdg   [iHadClu]    = prt->Pdg();
          brHadCluP4    [iHadClu][0] = prt->P4()->Px();
          brHadCluP4    [iHadClu][1] = prt->P4()->Py();
          brHadCluP4    [iHadClu][2] = prt->P4()->Pz();
          brHadCluP4    [iHadClu][3] = prt->P4()->E();
          eneHadClu += prt->P4()->E();
          iHadClu++;
        }
        else if ( prt->Status()==12 || prt->Status()==14 ) {
          
          if ( prt->Pdg()==92 || pdg::IsDiQuark(prt->Pdg()) || pdg::IsQuark(prt->Pdg()) ) {
            brPreFragHadIdx   [iPreFragHad]    = event.ParticlePosition(prt);
            brPreFragHadM     [iPreFragHad]    = prt->FirstMother();
            brPreFragHadPdg   [iPreFragHad]    = prt->Pdg(); 
            brPreFragHadP4    [iPreFragHad][0] = prt->P4()->Px(); 
            brPreFragHadP4    [iPreFragHad][1] = prt->P4()->Py(); 
            brPreFragHadP4    [iPreFragHad][2] = prt->P4()->Pz(); 
            brPreFragHadP4    [iPreFragHad][3] = prt->P4()->E(); 
            enePreFragHad += prt->P4()->E();
            iPreFragHad++;
          }
          if ( event.Particle(prt->FirstMother())->Pdg()==92 ) {
            brFragHadIdx   [iFragHad]    = event.ParticlePosition(prt);
            brFragHadM     [iFragHad]    = prt->FirstMother();
            brFragHadPdg   [iFragHad]    = prt->Pdg(); 
            brFragHadP4    [iFragHad][0] = prt->P4()->Px(); 
            brFragHadP4    [iFragHad][1] = prt->P4()->Py(); 
            brFragHadP4    [iFragHad][2] = prt->P4()->Pz(); 
            brFragHadP4    [iFragHad][3] = prt->P4()->E(); 
            eneFragHad += prt->P4()->E();
            iFragHad++;
          }
          if ( prt->Status()==14 ) {
            brFSIHadIdx   [iFSIHad]    = event.ParticlePosition(prt);
            brFSIHadM     [iFSIHad]    = prt->FirstMother();
            brFSIHadPdg   [iFSIHad]    = prt->Pdg(); 
            brFSIHadRescat[iFSIHad]    = prt->RescatterCode(); 
            brFSIHadP4    [iFSIHad][0] = prt->P4()->Px(); 
            brFSIHadP4    [iFSIHad][1] = prt->P4()->Py(); 
            brFSIHadP4    [iFSIHad][2] = prt->P4()->Pz(); 
            brFSIHadP4    [iFSIHad][3] = prt->P4()->E(); 
            eneFSIHad += prt->P4()->E();
            iFSIHad++;
          }
        }
        else if ( prt->Status()==1 ) {
          if ( pdg::IsHadron(prt->Pdg()) ) {
            if (prt->Charge()>0) {
              brFinalHadPIdx   [iFinalHadP]    = event.ParticlePosition(prt);
              brFinalHadPM     [iFinalHadP]    = prt->FirstMother();
              brFinalHadPPdg   [iFinalHadP]    = prt->Pdg(); 
              brFinalHadPP4    [iFinalHadP][0] = prt->P4()->Px(); 
              brFinalHadPP4    [iFinalHadP][1] = prt->P4()->Py(); 
              brFinalHadPP4    [iFinalHadP][2] = prt->P4()->Pz(); 
              brFinalHadPP4    [iFinalHadP][3] = prt->P4()->E(); 
              eneFinalHadP += prt->P4()->E();
              iFinalHadP++;
            }
            else if (prt->Charge()<0) {
              brFinalHadNIdx   [iFinalHadN]    = event.ParticlePosition(prt);
              brFinalHadNM     [iFinalHadN]    = prt->FirstMother();
              brFinalHadNPdg   [iFinalHadN]    = prt->Pdg(); 
              brFinalHadNP4    [iFinalHadN][0] = prt->P4()->Px(); 
              brFinalHadNP4    [iFinalHadN][1] = prt->P4()->Py(); 
              brFinalHadNP4    [iFinalHadN][2] = prt->P4()->Pz(); 
              brFinalHadNP4    [iFinalHadN][3] = prt->P4()->E(); 
              eneFinalHadN += prt->P4()->E();
              iFinalHadN++;
            }
            else {
              brFinalHad0Idx   [iFinalHad0]    = event.ParticlePosition(prt);
              brFinalHad0M     [iFinalHad0]    = prt->FirstMother();
              brFinalHad0Pdg   [iFinalHad0]    = prt->Pdg(); 
              brFinalHad0P4    [iFinalHad0][0] = prt->P4()->Px(); 
              brFinalHad0P4    [iFinalHad0][1] = prt->P4()->Py(); 
              brFinalHad0P4    [iFinalHad0][2] = prt->P4()->Pz(); 
              brFinalHad0P4    [iFinalHad0][3] = prt->P4()->E(); 
              eneFinalHad0 += prt->P4()->E();
              iFinalHad0++;
            }

          }
          else if (prt->FirstMother()==1) {
            brFinalDeExtIdx   [iFinalDeExt]   = event.ParticlePosition(prt);
            brFinalDeExtM     [iFinalDeExt]    = prt->FirstMother();
            brFinalDeExtPdg   [iFinalDeExt]    = prt->Pdg(); 
            brFinalDeExtP4    [iFinalDeExt][0] = prt->P4()->Px(); 
            brFinalDeExtP4    [iFinalDeExt][1] = prt->P4()->Py(); 
            brFinalDeExtP4    [iFinalDeExt][2] = prt->P4()->Pz(); 
            brFinalDeExtP4    [iFinalDeExt][3] = prt->P4()->E(); 
            eneFinalDeExt += prt->P4()->E();
            iFinalDeExt++;
          }
          else {
            if (TMath::Abs(prt->Pdg()==13)) {
              brFinalMuIdx   [iFinalMu]    = event.ParticlePosition(prt);
              brFinalMuM     [iFinalMu]    = prt->FirstMother();
              brFinalMuPdg   [iFinalMu]    = prt->Pdg(); 
              brFinalMuP4    [iFinalMu][0] = prt->P4()->Px(); 
              brFinalMuP4    [iFinalMu][1] = prt->P4()->Py(); 
              brFinalMuP4    [iFinalMu][2] = prt->P4()->Pz(); 
              brFinalMuP4    [iFinalMu][3] = prt->P4()->E(); 
              eneFinalMu += prt->P4()->E();
              iFinalMu++;
            }
            else if (TMath::Abs(prt->Pdg()==11)) {
              brFinalEleIdx   [iFinalEle]    = event.ParticlePosition(prt);
              brFinalEleM     [iFinalEle]    = prt->FirstMother();
              brFinalElePdg   [iFinalEle]    = prt->Pdg(); 
              brFinalEleP4    [iFinalEle][0] = prt->P4()->Px(); 
              brFinalEleP4    [iFinalEle][1] = prt->P4()->Py(); 
              brFinalEleP4    [iFinalEle][2] = prt->P4()->Pz(); 
              brFinalEleP4    [iFinalEle][3] = prt->P4()->E(); 
              eneFinalEle += prt->P4()->E();
              iFinalEle++;
            }
            else {
              brFinalOtherIdx   [iFinalOther]    = event.ParticlePosition(prt);
              brFinalOtherM     [iFinalOther]    = prt->FirstMother();
              brFinalOtherPdg   [iFinalOther]    = prt->Pdg(); 
              brFinalOtherP4    [iFinalOther][0] = prt->P4()->Px(); 
              brFinalOtherP4    [iFinalOther][1] = prt->P4()->Py(); 
              brFinalOtherP4    [iFinalOther][2] = prt->P4()->Pz(); 
              brFinalOtherP4    [iFinalOther][3] = prt->P4()->E(); 
              eneFinalOther += prt->P4()->E();
              iFinalOther++;
            }
          }
        }
        else std::cout << "Missing particle w/ state " << prt->Status() << std::endl;

    }
    brHadCluN     = iHadClu; 
    brPreFragHadN = iPreFragHad; 
    brFragHadN    = iFragHad; 
    brFSIHadN     = iFSIHad; 
    brFinalHadPN  = iFinalHadP; 
    brFinalHadNN  = iFinalHadN; 
    brFinalHad0N  = iFinalHad0; 
    brFinalMuN    = iFinalMu; 
    brFinalEleN   = iFinalEle; 
    brFinalOtherN = iFinalOther; 
    brFinalDeExtN = iFinalDeExt; 
    brHadCluE     = eneHadClu; 
    brPreFragHadE = enePreFragHad; 
    brFragHadE    = eneFragHad; 
    brFSIHadE     = eneFSIHad; 
    brFinalHadPE  = eneFinalHadP; 
    brFinalHadNE  = eneFinalHadN; 
    brFinalHad0E  = eneFinalHad0; 
    brFinalMuE    = eneFinalMu; 
    brFinalEleE   = eneFinalEle; 
    brFinalOtherE = eneFinalOther; 
    brFinalDeExtE = eneFinalDeExt; 

    // fill tree
    rootracker_tree->Fill();
    mcrec->Clear();

  } // event loop

  fin.Close();

  fout.Write();
  fout.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";

  return 0;

}
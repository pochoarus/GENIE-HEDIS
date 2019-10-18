//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TClonesArray.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Physics/Hadronization/LeptoHadronization.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

// the actual PYTHIA call
extern "C" {
  double pyangl_( double *,  double * );
  void   pykfdi_( int *,  int *, int *, int * );
  void   pyzdis_( int *,  int *, double *, double * );
  void   pyrobo_( int *,  int *, double *, double *, double *, double *, double * );
  void   pydecy_( int * );
  void   py2ent_( int *,  int *, int *, double * );
} 

//____________________________________________________________________________
LeptoHadronization::LeptoHadronization() :
HadronizationModelBase("genie::LeptoHadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
LeptoHadronization::LeptoHadronization(string config) :
HadronizationModelBase("genie::LeptoHadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
LeptoHadronization::~LeptoHadronization()
{

}
//____________________________________________________________________________
void LeptoHadronization::Initialize(void) const
{
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();
}
//____________________________________________________________________________
TClonesArray * 
  LeptoHadronization::Hadronize(
         const Interaction * interaction) const
{

  LOG("LeptoHad", pDEBUG) << "Running HEDIS PYTHIA hadronizer";

  // get kinematics / init-state / process-info
  const Kinematics &   kinematics = interaction->Kine();
  const XclsTag &      xclstag    = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();
  const Target &       target     = init_state.Tgt();

  assert(target.HitQrkIsSet()); 

  double W = kinematics.W();

  bool isp        = pdg::IsProton(target.HitNucPdg());
  int  hit_quark  = target.HitQrkPdg();
  int  frag_quark = xclstag.FinalQuarkPdg();

  LOG("LeptoHad", pDEBUG) << "Hit nucleon pdgc = " << target.HitNucPdg() << ", W = " << W;
  LOG("LeptoHad", pDEBUG) << "Selected hit quark pdgc = " << hit_quark << " // Fragmentation quark = " << frag_quark;

  RandomGen * rnd = RandomGen::Instance();

  //
  // Generate the hadron combination to input PYTHIA
  //

  //If the hit quark is a d we have these options:
  /* uud(->q)     => uu + q */
  /* uud d(->q)db => uu + q (d valence and db sea annihilates)*/
  /* udd(->q)     => ud + q */
  /* udd d(->q)db => ud + q (d valence and db sea annihilates)*/
  if ( pdg::IsDQuark(hit_quark) ) {
    // choose diquark system depending on proton or neutron
    int diquark = 0;
    if (isp) diquark = kPdgUUDiquarkS1;
    else     diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
    // Check that the trasnferred energy is higher than the mass of the produced quarks
    double m_frag    = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_diquark = PDGLibrary::Instance()->Find(diquark)->Mass();
    if( W <= m_frag + m_diquark + fMinESinglet ) {
      LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("LeptoHad", pWARN) << "frag_quark = " << frag_quark << "    -> m = " << m_frag;
      LOG("LeptoHad", pWARN) << "diquark    = " << diquark    << " -> m = "    << m_diquark;
      return 0;
    }
    // Input the two particles to PYTHIA back to back in the CM frame
    // If a top quark is produced we decay it because it does not hadronize
    fPythia->Py1ent( -1, frag_quark, (W*W - m_diquark*m_diquark + m_frag*m_frag)/2./W, 0., 0. ); //k(1,2) = 2
    if ( pdg::IsTQuark(frag_quark) ) {
      int ip = 1;
      pydecy_(&ip);
    }
    fPythia->Py1ent( fPythia->GetN()+1,    diquark,  (W*W + m_diquark*m_diquark - m_frag*m_frag)/2./W, fPythia->GetPARU(1), 0. ); //k(2,2) = 1
  }

  //If the hit quark is a u we have these options:
  /* u(->q)ud     => ud + q */
  /* uud u(->q)ub => ud + q (u valence and ub sea annihilates)*/
  /* u(->q)dd     => dd + q */
  /* udd u(->q)ub => dd + q (u valence and ub sea annihilates)*/
  else if ( pdg::IsUQuark(hit_quark) ) {
    // choose diquark system depending on proton or neutron
    int diquark = 0;
    if (isp) diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
    else     diquark = kPdgDDDiquarkS1;
    // Check that the trasnferred energy is higher than the mass of the produced quarks.
    double m_frag    = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_diquark = PDGLibrary::Instance()->Find(diquark)->Mass();
    if( W <= m_frag + m_diquark + fMinESinglet ) {
      LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("LeptoHad", pWARN) << "frag_quark = " << frag_quark << "    -> m = " << m_frag;
      LOG("LeptoHad", pWARN) << "diquark    = " << diquark    << " -> m = "    << m_diquark;
      return 0;
    }
    // Input the two particles to PYTHIA back to back in the CM frame
    fPythia->Py1ent( -1, frag_quark, (W*W - m_diquark*m_diquark + m_frag*m_frag)/2./W, 0., 0. ); //k(1,2) = 2
    fPythia->Py1ent( fPythia->GetN()+1,    diquark,  (W*W + m_diquark*m_diquark - m_frag*m_frag)/2./W, fPythia->GetPARU(1), 0. ); //k(2,2) = 1
  }
  else {

    // If the hit quark is not u or d then is more complicated.
    // We are using the same procedure use in LEPTO (see lqev.F)
    // Our initial systemt will look like this          ->  qqq + hit_q(->frag_q) + rema_q
    // And we have to input PYTHIA something like this  ->  frag_q + rema  + hadron
    // These are the posible combinations               ->  frag_q[q] + meson [qqb]  + diquark [qq]
    //                                                  ->  frag_q[qb] + baryon [qqq] + quark [q]

    // Remnant of the hit quark (which is from the sea) will be of opposite charge
    int rema_hit_quark = -hit_quark;

    // Check that the trasnfered energy is higher than the mass of the produce quarks plus remnant quark and nucleon
    double m_frag     = PDGLibrary::Instance()->Find(frag_quark)->Mass();
    double m_rema_hit = PDGLibrary::Instance()->Find(rema_hit_quark)->Mass();
    if (W <= m_frag + m_rema_hit + 0.9 + fMinESinglet ) {
      LOG("LeptoHad", pWARN) << "Low invariant mass, W = " << W << " GeV! Returning a null list";
      LOG("LeptoHad", pWARN) << " frag_quark     = " << frag_quark     << " -> m = " << m_frag;
      LOG("LeptoHad", pWARN) << " rema_hit_quark = " << rema_hit_quark << " -> m = " << m_rema_hit;
      return 0;
    }

    //PDG of the two hadronic particles for the final state
    int hadron = 0;      
    int rema   = 0;

    int ntwoq = isp ? 2 : 1; //proton two ups & neutron one up
    int counter = 0;

    // Here we select the id and kinematics of the hadron and rema particles
    // Some combinations can be kinematically forbiden so we repeat this process
    // up to 100 times before the event is discarded.
    while( counter<fMaxIterHad ) {

      // Loop to create a combination of hadron + rema. Two options are possible:
      // 1) diquark [qq] + meson [qqb]
      // 2) quark [q] + baryon [qqq]
      while(hadron==0) {
        //choose a valence quark and the remaining will be a diquark system
        int valquark = int(1.+ntwoq/3.+rnd->RndHadro().Rndm());
        int diquark  = 0;
        if ( valquark==ntwoq ) diquark = rnd->RndHadro().Rndm()>0.75 ? kPdgUDDiquarkS1 : kPdgUDDiquarkS0;
        else                   diquark = 1000*ntwoq+100*ntwoq+3;

        // Choose flavours using PYTHIA tool
        int idum;
        if ( rema_hit_quark>0 ) { //create a baryon (qqq)
          pykfdi_(&diquark,&rema_hit_quark,&idum,&hadron);
          rema = valquark;
        }
        else {                    //create a meson (qqbar)
          pykfdi_(&valquark,&rema_hit_quark,&idum,&hadron);
          rema = diquark;
        }
      }

      double m_hadron = PDGLibrary::Instance()->Find(hadron)->Mass();
      double m_rema   = PDGLibrary::Instance()->Find(rema)->Mass();

      // Give balancing pT to hadron and rema particles
      double pT  = fRemnantPT * TMath::Sqrt( -1*TMath::Log( rnd->RndHadro().Rndm() ) );
      double pT2 = TMath::Power(pT,2);
      double pr  = TMath::Power(m_hadron,2)+pT2;
      int kfl1 = 1;
      int kfl3 = 0;
      double z;
      //to generate the longitudinal scaling variable z in jet fragmentation using PYTHIA function
      // Split energy-momentum of remnant using PYTHIA function
      // z=E-pz fraction for rema forming jet-system with frag_q
      // 1-z=E-pz fraction for hadron
      pyzdis_(&kfl1,&kfl3,&pr,&z);  

      // Energy of trasnfered to the hadron
      double tm_hadron = pr / z / W;
      double E_hadron   = 0.5 * ( z*W + tm_hadron );  //E_hadron - pz = zW
      double E_pz       = W - tm_hadron;
      double WT         = (1-z) * W * E_pz - pT2;
      
      // Check if energy in jet system is enough for fragmentation.
      if ( WT > TMath::Power(m_frag+m_rema+fMinESinglet,2) ) {

        // Energy of transfered to the fragmented quark and rema system
        // Applying energy conservation
        WT = TMath::Sqrt( WT + pT2 );
        double tm_rema   = TMath::Power(m_rema,2) + pT2;
        double E_frag    = 0.5 * ( WT + ( TMath::Power(m_frag,2) - tm_rema)/WT ); //E_frag + E_rema = WT
        double E_rema    = 0.5 * ( WT + (-TMath::Power(m_frag,2) + tm_rema)/WT );
        double x_rema    = -1 * TMath::Sqrt( TMath::Power(E_rema,2) - tm_rema );
        double theta_rema = pyangl_(&x_rema,&pT);

        // Select a phi angle between between particles randomly
        double phi = 2*kPi*rnd->RndHadro().Rndm();

        // Input the three particles to PYTHIA in the CM frame
        // If a top quark is produced we decay it because it does not hadronize
        fPythia->Py1ent( -1, frag_quark, E_frag, 0.,         0. );           //k(1,2) = 2
        if (TMath::Abs(frag_quark) > 5 ) {
          int ip = 1;
          pydecy_(&ip);
        }
        fPythia->Py1ent( fPythia->GetN()+1, rema, E_rema, theta_rema, phi ); //k(2,2) = 1

        int imin     = 0;
        int imax     = 0;
        double the  = 0.; double ph   = 0.;
        double dbex = 0.; double dbey = 0.; double dbez = (E_pz-(1-z)*W)/(E_pz+(1-z)*W);
        pyrobo_( &imin , &imax, &the, &ph, &dbex, &dbey , &dbez );
      
        double pz_hadron  = -0.5 * ( z*W - tm_hadron );
        double theta_hadron = pyangl_(&pz_hadron,&pT);
        fPythia->SetMSTU( 10, 1 ); //keep the mass value stored in P(I,5), whatever it is.
        fPythia->SetP( fPythia->GetN()+1, 5, m_hadron );
        fPythia->Py1ent( fPythia->GetN()+1, hadron, E_hadron, theta_hadron, phi + kPi );
        fPythia->SetMSTU( 10, 2 ); //find masses according to mass tables as usual.

        // Target remnants required to go backwards in hadronic cms
        if ( fPythia->GetP(fPythia->GetN()-1,3)<0 && fPythia->GetP(fPythia->GetN(),3)<0 ) break; //quit the while from line 368
        
        LOG("LeptoHad", pINFO) << "Not backward hadron or rema";
        LOG("LeptoHad", pINFO) << "hadron     = " << hadron     << " -> Pz = " << fPythia->GetP(fPythia->GetN(),3) ;
        LOG("LeptoHad", pINFO) << "rema = " << rema << " -> Pz = " << fPythia->GetP(fPythia->GetN()-1,3) ;

      }
      else {
        LOG("LeptoHad", pINFO) << "Low WT value ... ";
        LOG("LeptoHad", pINFO) << "WT = " << TMath::Sqrt(WT) << " // m_frag = " << m_frag << " // m_rema = " << m_rema;
      }

      LOG("LeptoHad", pINFO) << "Hadronization paricles not suitable. Trying again... " << counter;
      counter++;
      if (counter==100) {
        LOG("LeptoHad", pWARN) << "Hadronization particles failed after " << counter << " iterations! Returning a null list";
        return 0;
      }

    }

  }

  // Introduce a primordial kT system
  double pT  = fPrimordialKT * TMath::Sqrt( -1*TMath::Log( rnd->RndHadro().Rndm() ) );
  double phi   = -2*kPi*rnd->RndHadro().Rndm();
  double theta = 0.;
  int imin     = 0;
  int imax     = 0;
  double dbex = 0.; double dbey = 0.; double dbez = 0;
  pyrobo_( &imin , &imax, &theta, &phi, &dbex, &dbey , &dbez );
  phi   = -1 * phi;
  theta = TMath::ATan(2.*pT/W);
  pyrobo_( &imin , &imax, &theta, &phi, &dbex, &dbey , &dbez );

  // Run PYTHIA with the input particles
  fPythia->Pyexec();
  
  // Use for debugging purposes
  if (fPromptPythiaList) fPythia->Pylist(3);

  // get LUJETS record
  fPythia->GetPrimaries();

  TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("TMCParticle", np);
  particle_list->SetOwner(true);

  bool isTop = false;
  register unsigned int i = 0;
  TMCParticle * particle = 0;
  TIter particle_iter(pythia_particles);
  while( (particle = (TMCParticle *) particle_iter.Next()) ) {
    
    // Final state particles can not be quarks or diquarks but colorless 
    if(particle->GetKS() == 1) {
      if( pdg::IsQuark(particle->GetKF()) || pdg::IsDiQuark(particle->GetKF()) ) {
        LOG("LeptoHad", pERROR) << "Hadronization failed! Bare quark/di-quarks appear in final state!";
        particle_list->Delete();
        delete particle_list;
        return 0;            
      }  
    }

    // When top quark is produced, it is immidiately decay before hadronization. Then the decayed
    // products are hadronized with the hadron remnants. Therefore, we remove the top quark from 
    // the list of particles so that the mother/daugher assigments is at the same level for decayed
    // products and hadron remnants.
    if ( i==0 && pdg::IsTQuark( TMath::Abs(particle->GetKF()) ) ) { isTop=true; continue; }

    // fix numbering scheme used for mother/daughter assignments
    if ( isTop ) {
      (particle->GetParent()==0) ? particle->SetParent(particle->GetParent() - 1) : particle->SetParent(particle->GetParent() - 2);
      particle->SetFirstChild (particle->GetFirstChild() - 2);
      particle->SetLastChild  (particle->GetLastChild()  - 2);
    }
    else  {
      particle->SetParent(particle->GetParent() - 1);
      particle->SetFirstChild (particle->GetFirstChild() - 1);
      particle->SetLastChild  (particle->GetLastChild()  - 1);
    }                                  
    // insert the particle in the list
    new ( (*particle_list)[i++] ) TMCParticle(*particle);
  }

  delete particle;
  pythia_particles->Clear("C");

  return particle_list;

}
//____________________________________________________________________________
PDGCodeList * LeptoHadronization::SelectParticles(
                                        const Interaction * interaction) const
{
  return 0;
}
//____________________________________________________________________________
TH1D * LeptoHadronization::MultiplicityProb(
           const Interaction * interaction, Option_t * opt) const
{
  return 0;
}
//____________________________________________________________________________
double LeptoHadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
void LeptoHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LeptoHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LeptoHadronization::LoadConfig(void)
{


  GetParamDef("MaxIter-Had", fMaxIterHad, 100 ) ;

  // Width of Gaussian distribution for transverse momentums
  // Define in LEPTO with PARL(3) and PARL(14)
  GetParamDef("Primordial-kT", fPrimordialKT, 0.44 ) ;
  GetParamDef("Remnant-pT",    fRemnantPT,    0.35 ) ;

  // It is, with quark masses added, used to define the minimum allowable energy of a colour-singlet parton system.
  GetParam( "Energy-Singlet", fMinESinglet ) ; 

  // PYTHIA parameters only valid for HEDIS
  int warnings;       GetParam( "PYTHIA-Warnings",      warnings    , 100000 ) ;
  int errors;         GetParam( "PYTHIA-Errors",        errors      , 100000 ) ;
  int qrk_mass;       GetParam( "PYTHIA-QuarkMass",     qrk_mass    ,      1 ) ;
  int decaycut;       GetParam( "PYTHIA-DecayCutOff",   decaycut    ,      2 ) ;
  double decaylength; GetParam( "PYTHIA-DecayLength",   decaylength , 0.0003 ) ;
  fPythia->SetMSTU(26, warnings);     // (Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, errors);       // (Default=10) maximum number of errors that are printed
  fPythia->SetMSTJ(93, qrk_mass);     // light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.
  fPythia->SetMSTJ(22, decaycut);     // (Default=1) cut-off on decay length for a particle that is allowed to decay according to MSTJ(21) and the MDCY value
  fPythia->SetPARJ(71, decaylength);  // (Default=10. mm) maximum average proper lifetime cτ for particles allowed to decay

  GetParamDef("PromptPythiaList", fPromptPythiaList, false ) ;

}
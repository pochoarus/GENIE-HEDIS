#include "Physics/HEDIS/EventGen/HEDISInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
HEDISInteractionListGenerator::HEDISInteractionListGenerator() :
InteractionListGeneratorI("genie::HEDISInteractionListGenerator")
{

}
//___________________________________________________________________________
HEDISInteractionListGenerator::HEDISInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::HEDISInteractionListGenerator", config)
{

}
//___________________________________________________________________________
HEDISInteractionListGenerator::~HEDISInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * HEDISInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  int nupdg = init_state.ProbePdg();
  if( !pdg::IsLepton(nupdg) ) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }


  const int n_nucc_channels = 36;
  const int n_nunc_channels = 24;
  HEDISChannel_t nucc_channels[n_nucc_channels] = {kHEDISNull};
  HEDISChannel_t nunc_channels[n_nunc_channels] = {kHEDISNull};

  if( pdg::IsNeutrino(nupdg) ) {    
    nucc_channels[0]  = kHEDIS_v_cc_p_dval_u,
    nucc_channels[1]  = kHEDIS_v_cc_p_dval_c,
    nucc_channels[2]  = kHEDIS_v_cc_p_dval_t;
    nucc_channels[3]  = kHEDIS_v_cc_p_dsea_u;
    nucc_channels[4]  = kHEDIS_v_cc_p_dsea_c;
    nucc_channels[5]  = kHEDIS_v_cc_p_dsea_t;
    nucc_channels[6]  = kHEDIS_v_cc_p_ssea_u;
    nucc_channels[7]  = kHEDIS_v_cc_p_ssea_c;
    nucc_channels[8]  = kHEDIS_v_cc_p_ssea_t;
    nucc_channels[9]  = kHEDIS_v_cc_p_bsea_u;
    nucc_channels[10] = kHEDIS_v_cc_p_bsea_c;
    nucc_channels[11] = kHEDIS_v_cc_p_bsea_t;
    nucc_channels[12] = kHEDIS_v_cc_p_ubarsea_dbar;
    nucc_channels[13] = kHEDIS_v_cc_p_ubarsea_sbar;
    nucc_channels[14] = kHEDIS_v_cc_p_ubarsea_bbar;
    nucc_channels[15] = kHEDIS_v_cc_p_cbarsea_dbar;
    nucc_channels[16] = kHEDIS_v_cc_p_cbarsea_sbar;
    nucc_channels[17] = kHEDIS_v_cc_p_cbarsea_bbar;
    nucc_channels[18] = kHEDIS_v_cc_n_dval_u;
    nucc_channels[19] = kHEDIS_v_cc_n_dval_c;
    nucc_channels[20] = kHEDIS_v_cc_n_dval_t;
    nucc_channels[21] = kHEDIS_v_cc_n_dsea_u;
    nucc_channels[22] = kHEDIS_v_cc_n_dsea_c;
    nucc_channels[23] = kHEDIS_v_cc_n_dsea_t;
    nucc_channels[24] = kHEDIS_v_cc_n_ssea_u;
    nucc_channels[25] = kHEDIS_v_cc_n_ssea_c;
    nucc_channels[26] = kHEDIS_v_cc_n_ssea_t;
    nucc_channels[27] = kHEDIS_v_cc_n_bsea_u;
    nucc_channels[28] = kHEDIS_v_cc_n_bsea_c;
    nucc_channels[29] = kHEDIS_v_cc_n_bsea_t;
    nucc_channels[30] = kHEDIS_v_cc_n_ubarsea_dbar;
    nucc_channels[31] = kHEDIS_v_cc_n_ubarsea_sbar;
    nucc_channels[32] = kHEDIS_v_cc_n_ubarsea_bbar;
    nucc_channels[33] = kHEDIS_v_cc_n_cbarsea_dbar;
    nucc_channels[34] = kHEDIS_v_cc_n_cbarsea_sbar;
    nucc_channels[35] = kHEDIS_v_cc_n_cbarsea_bbar;
    nunc_channels[0]  = kHEDIS_v_nc_p_dval_d;
    nunc_channels[1]  = kHEDIS_v_nc_p_uval_u;
    nunc_channels[2]  = kHEDIS_v_nc_p_dsea_d;
    nunc_channels[3]  = kHEDIS_v_nc_p_usea_u;
    nunc_channels[4]  = kHEDIS_v_nc_p_ssea_s;
    nunc_channels[5]  = kHEDIS_v_nc_p_csea_c;
    nunc_channels[6]  = kHEDIS_v_nc_p_bsea_b;
    nunc_channels[7]  = kHEDIS_v_nc_p_dbarsea_dbar;
    nunc_channels[8]  = kHEDIS_v_nc_p_ubarsea_ubar;
    nunc_channels[9]  = kHEDIS_v_nc_p_sbarsea_sbar;
    nunc_channels[10] = kHEDIS_v_nc_p_cbarsea_cbar;
    nunc_channels[11] = kHEDIS_v_nc_p_bbarsea_bbar;
    nunc_channels[12] = kHEDIS_v_nc_n_dval_d;
    nunc_channels[13] = kHEDIS_v_nc_n_uval_u;
    nunc_channels[14] = kHEDIS_v_nc_n_dsea_d;
    nunc_channels[15] = kHEDIS_v_nc_n_usea_u;
    nunc_channels[16] = kHEDIS_v_nc_n_ssea_s;
    nunc_channels[17] = kHEDIS_v_nc_n_csea_c;
    nunc_channels[18] = kHEDIS_v_nc_n_bsea_b;
    nunc_channels[19] = kHEDIS_v_nc_n_dbarsea_dbar;
    nunc_channels[20] = kHEDIS_v_nc_n_ubarsea_ubar;
    nunc_channels[21] = kHEDIS_v_nc_n_sbarsea_sbar;
    nunc_channels[22] = kHEDIS_v_nc_n_cbarsea_cbar;
    nunc_channels[23] = kHEDIS_v_nc_n_bbarsea_bbar;
  } 
  else if ( pdg::IsAntiNeutrino(nupdg) ) {
    nucc_channels[0]  = kHEDIS_vbar_cc_p_uval_d;
    nucc_channels[1]  = kHEDIS_vbar_cc_p_uval_s;
    nucc_channels[2]  = kHEDIS_vbar_cc_p_uval_b;
    nucc_channels[3]  = kHEDIS_vbar_cc_p_usea_d;
    nucc_channels[4]  = kHEDIS_vbar_cc_p_usea_s;
    nucc_channels[5]  = kHEDIS_vbar_cc_p_usea_b;
    nucc_channels[6]  = kHEDIS_vbar_cc_p_csea_d;
    nucc_channels[7]  = kHEDIS_vbar_cc_p_csea_s;
    nucc_channels[8]  = kHEDIS_vbar_cc_p_csea_b;
    nucc_channels[9]  = kHEDIS_vbar_cc_p_dbarsea_ubar;
    nucc_channels[10] = kHEDIS_vbar_cc_p_dbarsea_cbar;
    nucc_channels[11] = kHEDIS_vbar_cc_p_dbarsea_tbar;
    nucc_channels[12] = kHEDIS_vbar_cc_p_sbarsea_ubar;
    nucc_channels[13] = kHEDIS_vbar_cc_p_sbarsea_cbar;
    nucc_channels[14] = kHEDIS_vbar_cc_p_sbarsea_tbar;
    nucc_channels[15] = kHEDIS_vbar_cc_p_bbarsea_ubar;
    nucc_channels[16] = kHEDIS_vbar_cc_p_bbarsea_cbar;
    nucc_channels[17] = kHEDIS_vbar_cc_p_bbarsea_tbar;
    nucc_channels[18] = kHEDIS_vbar_cc_n_uval_d;
    nucc_channels[19] = kHEDIS_vbar_cc_n_uval_s;
    nucc_channels[20] = kHEDIS_vbar_cc_n_uval_b;
    nucc_channels[21] = kHEDIS_vbar_cc_n_usea_d;
    nucc_channels[22] = kHEDIS_vbar_cc_n_usea_s;
    nucc_channels[23] = kHEDIS_vbar_cc_n_usea_b;
    nucc_channels[24] = kHEDIS_vbar_cc_n_csea_d;
    nucc_channels[25] = kHEDIS_vbar_cc_n_csea_s;
    nucc_channels[26] = kHEDIS_vbar_cc_n_csea_b;
    nucc_channels[27] = kHEDIS_vbar_cc_n_dbarsea_ubar;
    nucc_channels[28] = kHEDIS_vbar_cc_n_dbarsea_cbar;
    nucc_channels[29] = kHEDIS_vbar_cc_n_dbarsea_tbar;
    nucc_channels[30] = kHEDIS_vbar_cc_n_sbarsea_ubar;
    nucc_channels[31] = kHEDIS_vbar_cc_n_sbarsea_cbar;
    nucc_channels[32] = kHEDIS_vbar_cc_n_sbarsea_tbar;
    nucc_channels[33] = kHEDIS_vbar_cc_n_bbarsea_ubar;
    nucc_channels[34] = kHEDIS_vbar_cc_n_bbarsea_cbar;
    nucc_channels[35] = kHEDIS_vbar_cc_n_bbarsea_tbar;
    nunc_channels[0]  = kHEDIS_vbar_nc_p_dval_d;
    nunc_channels[1]  = kHEDIS_vbar_nc_p_uval_u;
    nunc_channels[2]  = kHEDIS_vbar_nc_p_dsea_d;
    nunc_channels[3]  = kHEDIS_vbar_nc_p_usea_u;
    nunc_channels[4]  = kHEDIS_vbar_nc_p_ssea_s;
    nunc_channels[5]  = kHEDIS_vbar_nc_p_csea_c;
    nunc_channels[6]  = kHEDIS_vbar_nc_p_bsea_b;
    nunc_channels[7]  = kHEDIS_vbar_nc_p_dbarsea_dbar;
    nunc_channels[8]  = kHEDIS_vbar_nc_p_ubarsea_ubar;
    nunc_channels[9]  = kHEDIS_vbar_nc_p_sbarsea_sbar;
    nunc_channels[10] = kHEDIS_vbar_nc_p_cbarsea_cbar;
    nunc_channels[11] = kHEDIS_vbar_nc_p_bbarsea_bbar;
    nunc_channels[12] = kHEDIS_vbar_nc_n_dval_d;
    nunc_channels[13] = kHEDIS_vbar_nc_n_uval_u;
    nunc_channels[14] = kHEDIS_vbar_nc_n_dsea_d;
    nunc_channels[15] = kHEDIS_vbar_nc_n_usea_u;
    nunc_channels[16] = kHEDIS_vbar_nc_n_ssea_s;
    nunc_channels[17] = kHEDIS_vbar_nc_n_csea_c;
    nunc_channels[18] = kHEDIS_vbar_nc_n_bsea_b;
    nunc_channels[19] = kHEDIS_vbar_nc_n_dbarsea_dbar;
    nunc_channels[20] = kHEDIS_vbar_nc_n_ubarsea_ubar;
    nunc_channels[21] = kHEDIS_vbar_nc_n_sbarsea_sbar;
    nunc_channels[22] = kHEDIS_vbar_nc_n_cbarsea_cbar;
    nunc_channels[23] = kHEDIS_vbar_nc_n_bbarsea_bbar;
  } 
  else {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  bool hasP = (init_state.Tgt().Z() > 0);
  bool hasN = (init_state.Tgt().N() > 0);

  InteractionList * intlist = new InteractionList;

  if (fIsCC) {
    for(int i=0; i<n_nucc_channels; i++) {
      int struck_nucleon = HEDISChannel::HitNuclPdg(nucc_channels[i]);
      if( (struck_nucleon == kPdgProton && hasP) || (struck_nucleon == kPdgNeutron && hasN) ) {
        ProcessInfo proc_info(kScHEDIS, kIntWeakCC);
        Interaction * interaction = new Interaction(init_state, proc_info);
        Target * target = interaction->InitStatePtr()->TgtPtr();
        target->SetHitNucPdg(struck_nucleon);
        this->AddFinalStateInfo(interaction, nucc_channels[i]);
        intlist->push_back(interaction);
      }
    }
  }
  else if (fIsNC) {
    // NC
    for(int i=0; i<n_nunc_channels; i++) {
      int struck_nucleon = HEDISChannel::HitNuclPdg(nunc_channels[i]);
      if( (struck_nucleon == kPdgProton && hasP) || (struck_nucleon == kPdgNeutron && hasN) ) {
        ProcessInfo proc_info(kScHEDIS, kIntWeakNC);
        Interaction * interaction = new Interaction(init_state, proc_info);
        interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(struck_nucleon);
        this->AddFinalStateInfo(interaction, nunc_channels[i]);
        intlist->push_back(interaction);
      }
    }//nc channels   
  }

  if(intlist->size() == 0) {
     LOG("IntLst", pERROR)
         << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;

}
//___________________________________________________________________________
void HEDISInteractionListGenerator::AddFinalStateInfo(
                       Interaction * interaction, HEDISChannel_t hedischan) const
{

  bool iq_sea = HEDISChannel::HitQuarkSea(hedischan);
  int  iq_pdg = HEDISChannel::HitQuarkPdg(hedischan);
  int  fq_pdg = HEDISChannel::FnlQuarkPdg(hedischan);

  interaction->InitStatePtr()->TgtPtr()->SetHitSeaQrk(iq_sea);
  interaction->InitStatePtr()->TgtPtr()->SetHitQrkPdg(iq_pdg);

  XclsTag exclusive_tag;
  exclusive_tag.SetFinalQuark (fq_pdg);
  exclusive_tag.SetHEDISChannel (hedischan);
  interaction->SetExclTag(exclusive_tag);

}
//___________________________________________________________________________
void HEDISInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void HEDISInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void HEDISInteractionListGenerator::LoadConfigData(void)
{

  GetParamDef("is-CC", fIsCC, false ) ;
  GetParamDef("is-NC", fIsNC, false ) ;

}

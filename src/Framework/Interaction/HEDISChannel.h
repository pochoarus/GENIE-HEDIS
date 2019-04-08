#ifndef _HEDIS_CHANNEL_H_
#define _HEDIS_CHANNEL_H_

#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Interaction/InteractionType.h"

namespace genie {

  typedef enum EHEDISNuclChannel {

    kHEDISNuclNull = 0,

    kHEDISNucl_v_cc_p,
    kHEDISNucl_v_cc_n,
    kHEDISNucl_vbar_cc_p,
    kHEDISNucl_vbar_cc_n,
    kHEDISNucl_v_nc_p,
    kHEDISNucl_v_nc_n,
    kHEDISNucl_vbar_nc_p,
    kHEDISNucl_vbar_nc_n,
    kHEDISNucl_numofchannels

  } HEDISNuclChannel_t;

  typedef enum EHEDISChannel {

    kHEDISNull = 0,

    kHEDIS_v_cc_p_dval_u,
    kHEDIS_v_cc_p_dval_c,
    kHEDIS_v_cc_p_dval_t,
    kHEDIS_v_cc_p_dsea_u,
    kHEDIS_v_cc_p_dsea_c,
    kHEDIS_v_cc_p_dsea_t,
    kHEDIS_v_cc_p_ssea_u,
    kHEDIS_v_cc_p_ssea_c,
    kHEDIS_v_cc_p_ssea_t,
    kHEDIS_v_cc_p_bsea_u,
    kHEDIS_v_cc_p_bsea_c,
    kHEDIS_v_cc_p_bsea_t,
    kHEDIS_v_cc_p_ubarsea_dbar,
    kHEDIS_v_cc_p_ubarsea_sbar,
    kHEDIS_v_cc_p_ubarsea_bbar,
    kHEDIS_v_cc_p_cbarsea_dbar,
    kHEDIS_v_cc_p_cbarsea_sbar,
    kHEDIS_v_cc_p_cbarsea_bbar,

    kHEDIS_v_cc_n_dval_u,
    kHEDIS_v_cc_n_dval_c,
    kHEDIS_v_cc_n_dval_t,
    kHEDIS_v_cc_n_dsea_u,
    kHEDIS_v_cc_n_dsea_c,
    kHEDIS_v_cc_n_dsea_t,
    kHEDIS_v_cc_n_ssea_u,
    kHEDIS_v_cc_n_ssea_c,
    kHEDIS_v_cc_n_ssea_t,
    kHEDIS_v_cc_n_bsea_u,
    kHEDIS_v_cc_n_bsea_c,
    kHEDIS_v_cc_n_bsea_t,
    kHEDIS_v_cc_n_ubarsea_dbar,
    kHEDIS_v_cc_n_ubarsea_sbar,
    kHEDIS_v_cc_n_ubarsea_bbar,
    kHEDIS_v_cc_n_cbarsea_dbar,
    kHEDIS_v_cc_n_cbarsea_sbar,
    kHEDIS_v_cc_n_cbarsea_bbar,

    kHEDIS_vbar_cc_p_uval_d,
    kHEDIS_vbar_cc_p_uval_s,
    kHEDIS_vbar_cc_p_uval_b,
    kHEDIS_vbar_cc_p_usea_d,
    kHEDIS_vbar_cc_p_usea_s,
    kHEDIS_vbar_cc_p_usea_b,
    kHEDIS_vbar_cc_p_csea_d,
    kHEDIS_vbar_cc_p_csea_s,
    kHEDIS_vbar_cc_p_csea_b,
    kHEDIS_vbar_cc_p_dbarsea_ubar,
    kHEDIS_vbar_cc_p_dbarsea_cbar,
    kHEDIS_vbar_cc_p_dbarsea_tbar,
    kHEDIS_vbar_cc_p_sbarsea_ubar,
    kHEDIS_vbar_cc_p_sbarsea_cbar,
    kHEDIS_vbar_cc_p_sbarsea_tbar,
    kHEDIS_vbar_cc_p_bbarsea_ubar,
    kHEDIS_vbar_cc_p_bbarsea_cbar,
    kHEDIS_vbar_cc_p_bbarsea_tbar,

    kHEDIS_vbar_cc_n_uval_d,
    kHEDIS_vbar_cc_n_uval_s,
    kHEDIS_vbar_cc_n_uval_b,
    kHEDIS_vbar_cc_n_usea_d,
    kHEDIS_vbar_cc_n_usea_s,
    kHEDIS_vbar_cc_n_usea_b,
    kHEDIS_vbar_cc_n_csea_d,
    kHEDIS_vbar_cc_n_csea_s,
    kHEDIS_vbar_cc_n_csea_b,
    kHEDIS_vbar_cc_n_dbarsea_ubar,
    kHEDIS_vbar_cc_n_dbarsea_cbar,
    kHEDIS_vbar_cc_n_dbarsea_tbar,
    kHEDIS_vbar_cc_n_sbarsea_ubar,
    kHEDIS_vbar_cc_n_sbarsea_cbar,
    kHEDIS_vbar_cc_n_sbarsea_tbar,
    kHEDIS_vbar_cc_n_bbarsea_ubar,
    kHEDIS_vbar_cc_n_bbarsea_cbar,
    kHEDIS_vbar_cc_n_bbarsea_tbar,

    kHEDIS_v_nc_p_dval_d,
    kHEDIS_v_nc_p_uval_u,
    kHEDIS_v_nc_p_dsea_d,
    kHEDIS_v_nc_p_usea_u,
    kHEDIS_v_nc_p_ssea_s,
    kHEDIS_v_nc_p_csea_c,
    kHEDIS_v_nc_p_bsea_b,
    kHEDIS_v_nc_p_dbarsea_dbar,
    kHEDIS_v_nc_p_ubarsea_ubar,
    kHEDIS_v_nc_p_sbarsea_sbar,
    kHEDIS_v_nc_p_cbarsea_cbar,
    kHEDIS_v_nc_p_bbarsea_bbar,

    kHEDIS_v_nc_n_dval_d,
    kHEDIS_v_nc_n_uval_u,
    kHEDIS_v_nc_n_dsea_d,
    kHEDIS_v_nc_n_usea_u,
    kHEDIS_v_nc_n_ssea_s,
    kHEDIS_v_nc_n_csea_c,
    kHEDIS_v_nc_n_bsea_b,
    kHEDIS_v_nc_n_dbarsea_dbar,
    kHEDIS_v_nc_n_ubarsea_ubar,
    kHEDIS_v_nc_n_sbarsea_sbar,
    kHEDIS_v_nc_n_cbarsea_cbar,
    kHEDIS_v_nc_n_bbarsea_bbar,

    kHEDIS_vbar_nc_p_dval_d,
    kHEDIS_vbar_nc_p_uval_u,
    kHEDIS_vbar_nc_p_dsea_d,
    kHEDIS_vbar_nc_p_usea_u,
    kHEDIS_vbar_nc_p_ssea_s,
    kHEDIS_vbar_nc_p_csea_c,
    kHEDIS_vbar_nc_p_bsea_b,
    kHEDIS_vbar_nc_p_dbarsea_dbar,
    kHEDIS_vbar_nc_p_ubarsea_ubar,
    kHEDIS_vbar_nc_p_sbarsea_sbar,
    kHEDIS_vbar_nc_p_cbarsea_cbar,
    kHEDIS_vbar_nc_p_bbarsea_bbar,

    kHEDIS_vbar_nc_n_dval_d,
    kHEDIS_vbar_nc_n_uval_u,
    kHEDIS_vbar_nc_n_dsea_d,
    kHEDIS_vbar_nc_n_usea_u,
    kHEDIS_vbar_nc_n_ssea_s,
    kHEDIS_vbar_nc_n_csea_c,
    kHEDIS_vbar_nc_n_bsea_b,
    kHEDIS_vbar_nc_n_dbarsea_dbar,
    kHEDIS_vbar_nc_n_ubarsea_ubar,
    kHEDIS_vbar_nc_n_sbarsea_sbar,
    kHEDIS_vbar_nc_n_cbarsea_cbar,
    kHEDIS_vbar_nc_n_bbarsea_bbar,
    kHEDIS_numofchannels

  } HEDISChannel_t;


  class HEDISChannel
  {
    public:


      //__________________________________________________________________________
      static string AsString(HEDISNuclChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNucl_v_cc_p)       : return "nu_CC_p";      break;
          case (kHEDISNucl_v_cc_n)       : return "nu_CC_n";      break;
          case (kHEDISNucl_vbar_cc_p)    : return "nubar_CC_p";   break;
          case (kHEDISNucl_vbar_cc_n)    : return "nubar_CC_n";   break;
          case (kHEDISNucl_v_nc_p)       : return "nu_NC_p";      break;
          case (kHEDISNucl_v_nc_n)       : return "nu_NC_n";      break;
          case (kHEDISNucl_vbar_nc_p)    : return "nubar_NC_p";   break;
          case (kHEDISNucl_vbar_nc_n)    : return "nubar_NC_n";   break;

          default : return "Unknown";  break;
        }
        return "Unknown";
      }
      //__________________________________________________________________________
      static HEDISChannel_t GetFirstHEDISChannel(HEDISNuclChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNucl_v_cc_p)       : return kHEDIS_v_cc_p_dval_u;      break;
          case (kHEDISNucl_v_cc_n)       : return kHEDIS_v_cc_n_dval_u;      break;
          case (kHEDISNucl_vbar_cc_p)    : return kHEDIS_vbar_cc_p_uval_d;   break;
          case (kHEDISNucl_vbar_cc_n)    : return kHEDIS_vbar_cc_n_uval_d;   break;
          case (kHEDISNucl_v_nc_p)       : return kHEDIS_v_nc_p_dval_d;      break;
          case (kHEDISNucl_v_nc_n)       : return kHEDIS_v_nc_n_dval_d;      break;
          case (kHEDISNucl_vbar_nc_p)    : return kHEDIS_vbar_nc_p_dval_d;   break;
          case (kHEDISNucl_vbar_nc_n)    : return kHEDIS_vbar_nc_n_dval_d;   break;

          default : return kHEDISNull;  break;
        }
        return kHEDISNull;
      }
      //__________________________________________________________________________
      static HEDISChannel_t GetLastHEDISChannel(HEDISNuclChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNucl_v_cc_p)       : return kHEDIS_v_cc_p_cbarsea_bbar;      break;
          case (kHEDISNucl_v_cc_n)       : return kHEDIS_v_cc_n_cbarsea_bbar;      break;
          case (kHEDISNucl_vbar_cc_p)    : return kHEDIS_vbar_cc_p_bbarsea_tbar;   break;
          case (kHEDISNucl_vbar_cc_n)    : return kHEDIS_vbar_cc_n_bbarsea_tbar;   break;
          case (kHEDISNucl_v_nc_p)       : return kHEDIS_v_nc_p_bbarsea_bbar;      break;
          case (kHEDISNucl_v_nc_n)       : return kHEDIS_v_nc_n_bbarsea_bbar;      break;
          case (kHEDISNucl_vbar_nc_p)    : return kHEDIS_vbar_nc_p_bbarsea_bbar;   break;
          case (kHEDISNucl_vbar_nc_n)    : return kHEDIS_vbar_nc_n_bbarsea_bbar;   break;

          default : return kHEDISNull;  break;
        }
        return kHEDISNull;
      }
      static InteractionType_t InteractionType(HEDISNuclChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNucl_v_nc_p)       : return kIntWeakNC;   break;
          case (kHEDISNucl_v_nc_n)       : return kIntWeakNC;   break;
          case (kHEDISNucl_vbar_nc_p)    : return kIntWeakNC;   break;
          case (kHEDISNucl_vbar_nc_n)    : return kIntWeakNC;   break;

          default : return kIntWeakCC;  break;
        }
        return kIntNull;
      }
      //__________________________________________________________________________
      static bool IsNu(HEDISNuclChannel_t channel)
      {
        switch (channel) {
          case (kHEDISNucl_vbar_cc_p)    : return false;   break;
          case (kHEDISNucl_vbar_cc_n)    : return false;   break;
          case (kHEDISNucl_vbar_nc_p)    : return false;   break;
          case (kHEDISNucl_vbar_nc_n)    : return false;   break;

          default : return true;  break;
        }
        return true;
      }
      //__________________________________________________________________________
      static int HitNuclPdg(HEDISNuclChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNucl_v_cc_n)       : return kPdgNeutron;      break;
          case (kHEDISNucl_vbar_cc_n)    : return kPdgNeutron;   break;
          case (kHEDISNucl_v_nc_n)       : return kPdgNeutron;      break;
          case (kHEDISNucl_vbar_nc_n)    : return kPdgNeutron;   break;

          default : return kPdgProton;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      static string AsString(HEDISChannel_t channel)
      {

        switch (channel) {
          case (kHEDIS_v_cc_p_dval_u)       : return "nu_CC_p_dval_u";   break;
          case (kHEDIS_v_cc_p_dval_c)       : return "nu_CC_p_dval_c";   break;
          case (kHEDIS_v_cc_p_dval_t)       : return "nu_CC_p_dval_t";   break;
          case (kHEDIS_v_cc_p_dsea_u)       : return "nu_CC_p_dsea_u";   break;
          case (kHEDIS_v_cc_p_dsea_c)       : return "nu_CC_p_dsea_c";   break;
          case (kHEDIS_v_cc_p_dsea_t)       : return "nu_CC_p_dsea_t";   break;
          case (kHEDIS_v_cc_p_ssea_u)       : return "nu_CC_p_ssea_u";   break;
          case (kHEDIS_v_cc_p_ssea_c)       : return "nu_CC_p_ssea_c";   break;
          case (kHEDIS_v_cc_p_ssea_t)       : return "nu_CC_p_ssea_t";   break;
          case (kHEDIS_v_cc_p_bsea_u)       : return "nu_CC_p_bsea_u";   break;
          case (kHEDIS_v_cc_p_bsea_c)       : return "nu_CC_p_bsea_c";   break;
          case (kHEDIS_v_cc_p_bsea_t)       : return "nu_CC_p_bsea_t";   break;
          case (kHEDIS_v_cc_p_ubarsea_dbar) : return "nu_CC_p_ubarsea_dbar";   break;
          case (kHEDIS_v_cc_p_ubarsea_sbar) : return "nu_CC_p_ubarsea_sbar";   break;
          case (kHEDIS_v_cc_p_ubarsea_bbar) : return "nu_CC_p_ubarsea_bbar";   break;
          case (kHEDIS_v_cc_p_cbarsea_dbar) : return "nu_CC_p_cbarsea_dbar";   break;
          case (kHEDIS_v_cc_p_cbarsea_sbar) : return "nu_CC_p_cbarsea_sbar";   break;
          case (kHEDIS_v_cc_p_cbarsea_bbar) : return "nu_CC_p_cbarsea_bbar";   break;

          case (kHEDIS_v_cc_n_dval_u)       : return "nu_CC_n_dval_u";   break;
          case (kHEDIS_v_cc_n_dval_c)       : return "nu_CC_n_dval_c";   break;
          case (kHEDIS_v_cc_n_dval_t)       : return "nu_CC_n_dval_t";   break;
          case (kHEDIS_v_cc_n_dsea_u)       : return "nu_CC_n_dsea_u";   break;
          case (kHEDIS_v_cc_n_dsea_c)       : return "nu_CC_n_dsea_c";   break;
          case (kHEDIS_v_cc_n_dsea_t)       : return "nu_CC_n_dsea_t";   break;
          case (kHEDIS_v_cc_n_ssea_u)       : return "nu_CC_n_ssea_u";   break;
          case (kHEDIS_v_cc_n_ssea_c)       : return "nu_CC_n_ssea_c";   break;
          case (kHEDIS_v_cc_n_ssea_t)       : return "nu_CC_n_ssea_t";   break;
          case (kHEDIS_v_cc_n_bsea_u)       : return "nu_CC_n_bsea_u";   break;
          case (kHEDIS_v_cc_n_bsea_c)       : return "nu_CC_n_bsea_c";   break;
          case (kHEDIS_v_cc_n_bsea_t)       : return "nu_CC_n_bsea_t";   break;
          case (kHEDIS_v_cc_n_ubarsea_dbar) : return "nu_CC_n_ubarsea_dbar";   break;
          case (kHEDIS_v_cc_n_ubarsea_sbar) : return "nu_CC_n_ubarsea_sbar";   break;
          case (kHEDIS_v_cc_n_ubarsea_bbar) : return "nu_CC_n_ubarsea_bbar";   break;
          case (kHEDIS_v_cc_n_cbarsea_dbar) : return "nu_CC_n_cbarsea_dbar";   break;
          case (kHEDIS_v_cc_n_cbarsea_sbar) : return "nu_CC_n_cbarsea_sbar";   break;
          case (kHEDIS_v_cc_n_cbarsea_bbar) : return "nu_CC_n_cbarsea_bbar";   break;

          case (kHEDIS_vbar_cc_p_uval_d)       : return "nubar_CC_p_uval_d";   break;
          case (kHEDIS_vbar_cc_p_uval_s)       : return "nubar_CC_p_uval_s";   break;
          case (kHEDIS_vbar_cc_p_uval_b)       : return "nubar_CC_p_uval_b";   break;
          case (kHEDIS_vbar_cc_p_usea_d)       : return "nubar_CC_p_usea_d";   break;
          case (kHEDIS_vbar_cc_p_usea_s)       : return "nubar_CC_p_usea_s";   break;
          case (kHEDIS_vbar_cc_p_usea_b)       : return "nubar_CC_p_usea_b";   break;
          case (kHEDIS_vbar_cc_p_csea_d)       : return "nubar_CC_p_csea_d";   break;
          case (kHEDIS_vbar_cc_p_csea_s)       : return "nubar_CC_p_csea_s";   break;
          case (kHEDIS_vbar_cc_p_csea_b)       : return "nubar_CC_p_csea_b";   break;
          case (kHEDIS_vbar_cc_p_dbarsea_ubar) : return "nubar_CC_p_dbarsea_ubar";   break;
          case (kHEDIS_vbar_cc_p_dbarsea_cbar) : return "nubar_CC_p_dbarsea_cbar";   break;
          case (kHEDIS_vbar_cc_p_dbarsea_tbar) : return "nubar_CC_p_dbarsea_tbar";   break;
          case (kHEDIS_vbar_cc_p_sbarsea_ubar) : return "nubar_CC_p_sbarsea_ubar";   break;
          case (kHEDIS_vbar_cc_p_sbarsea_cbar) : return "nubar_CC_p_sbarsea_cbar";   break;
          case (kHEDIS_vbar_cc_p_sbarsea_tbar) : return "nubar_CC_p_sbarsea_tbar";   break;
          case (kHEDIS_vbar_cc_p_bbarsea_ubar) : return "nubar_CC_p_bbarsea_ubar";   break;
          case (kHEDIS_vbar_cc_p_bbarsea_cbar) : return "nubar_CC_p_bbarsea_cbar";   break;
          case (kHEDIS_vbar_cc_p_bbarsea_tbar) : return "nubar_CC_p_bbarsea_tbar";   break;

          case (kHEDIS_vbar_cc_n_uval_d)       : return "nubar_CC_n_uval_d";   break;
          case (kHEDIS_vbar_cc_n_uval_s)       : return "nubar_CC_n_uval_s";   break;
          case (kHEDIS_vbar_cc_n_uval_b)       : return "nubar_CC_n_uval_b";   break;
          case (kHEDIS_vbar_cc_n_usea_d)       : return "nubar_CC_n_usea_d";   break;
          case (kHEDIS_vbar_cc_n_usea_s)       : return "nubar_CC_n_usea_s";   break;
          case (kHEDIS_vbar_cc_n_usea_b)       : return "nubar_CC_n_usea_b";   break;
          case (kHEDIS_vbar_cc_n_csea_d)       : return "nubar_CC_n_csea_d";   break;
          case (kHEDIS_vbar_cc_n_csea_s)       : return "nubar_CC_n_csea_s";   break;
          case (kHEDIS_vbar_cc_n_csea_b)       : return "nubar_CC_n_csea_b";   break;
          case (kHEDIS_vbar_cc_n_dbarsea_ubar) : return "nubar_CC_n_dbarsea_ubar";   break;
          case (kHEDIS_vbar_cc_n_dbarsea_cbar) : return "nubar_CC_n_dbarsea_cbar";   break;
          case (kHEDIS_vbar_cc_n_dbarsea_tbar) : return "nubar_CC_n_dbarsea_tbar";   break;
          case (kHEDIS_vbar_cc_n_sbarsea_ubar) : return "nubar_CC_n_sbarsea_ubar";   break;
          case (kHEDIS_vbar_cc_n_sbarsea_cbar) : return "nubar_CC_n_sbarsea_cbar";   break;
          case (kHEDIS_vbar_cc_n_sbarsea_tbar) : return "nubar_CC_n_sbarsea_tbar";   break;
          case (kHEDIS_vbar_cc_n_bbarsea_ubar) : return "nubar_CC_n_bbarsea_ubar";   break;
          case (kHEDIS_vbar_cc_n_bbarsea_cbar) : return "nubar_CC_n_bbarsea_cbar";   break;
          case (kHEDIS_vbar_cc_n_bbarsea_tbar) : return "nubar_CC_n_bbarsea_tbar";   break;

          case (kHEDIS_v_nc_p_dval_d)       : return "nu_NC_p_dval_d";   break;
          case (kHEDIS_v_nc_p_uval_u)       : return "nu_NC_p_uval_u";   break;
          case (kHEDIS_v_nc_p_dsea_d)       : return "nu_NC_p_dsea_d";   break;
          case (kHEDIS_v_nc_p_usea_u)       : return "nu_NC_p_usea_u";   break;
          case (kHEDIS_v_nc_p_ssea_s)       : return "nu_NC_p_ssea_s";   break;
          case (kHEDIS_v_nc_p_csea_c)       : return "nu_NC_p_csea_c";   break;
          case (kHEDIS_v_nc_p_bsea_b)       : return "nu_NC_p_bsea_b";   break;
          case (kHEDIS_v_nc_p_dbarsea_dbar) : return "nu_NC_p_dbarsea_dbar";   break;
          case (kHEDIS_v_nc_p_ubarsea_ubar) : return "nu_NC_p_ubarsea_ubar";   break;
          case (kHEDIS_v_nc_p_sbarsea_sbar) : return "nu_NC_p_sbarsea_sbar";   break;
          case (kHEDIS_v_nc_p_cbarsea_cbar) : return "nu_NC_p_cbarsea_cbar";   break;
          case (kHEDIS_v_nc_p_bbarsea_bbar) : return "nu_NC_p_bbarsea_bbar";   break;

          case (kHEDIS_v_nc_n_dval_d)       : return "nu_NC_n_dval_d";   break;
          case (kHEDIS_v_nc_n_uval_u)       : return "nu_NC_n_uval_u";   break;
          case (kHEDIS_v_nc_n_dsea_d)       : return "nu_NC_n_dsea_d";   break;
          case (kHEDIS_v_nc_n_usea_u)       : return "nu_NC_n_usea_u";   break;
          case (kHEDIS_v_nc_n_ssea_s)       : return "nu_NC_n_ssea_s";   break;
          case (kHEDIS_v_nc_n_csea_c)       : return "nu_NC_n_csea_c";   break;
          case (kHEDIS_v_nc_n_bsea_b)       : return "nu_NC_n_bsea_b";   break;
          case (kHEDIS_v_nc_n_dbarsea_dbar) : return "nu_NC_n_dbarsea_dbar";   break;
          case (kHEDIS_v_nc_n_ubarsea_ubar) : return "nu_NC_n_ubarsea_ubar";   break;
          case (kHEDIS_v_nc_n_sbarsea_sbar) : return "nu_NC_n_sbarsea_sbar";   break;
          case (kHEDIS_v_nc_n_cbarsea_cbar) : return "nu_NC_n_cbarsea_cbar";   break;
          case (kHEDIS_v_nc_n_bbarsea_bbar) : return "nu_NC_n_bbarsea_bbar";   break;

          case (kHEDIS_vbar_nc_p_dval_d)       : return "nubar_NC_p_dval_d";   break;
          case (kHEDIS_vbar_nc_p_uval_u)       : return "nubar_NC_p_uval_u";   break;
          case (kHEDIS_vbar_nc_p_dsea_d)       : return "nubar_NC_p_dsea_d";   break;
          case (kHEDIS_vbar_nc_p_usea_u)       : return "nubar_NC_p_usea_u";   break;
          case (kHEDIS_vbar_nc_p_ssea_s)       : return "nubar_NC_p_ssea_s";   break;
          case (kHEDIS_vbar_nc_p_csea_c)       : return "nubar_NC_p_csea_c";   break;
          case (kHEDIS_vbar_nc_p_bsea_b)       : return "nubar_NC_p_bsea_b";   break;
          case (kHEDIS_vbar_nc_p_dbarsea_dbar) : return "nubar_NC_p_dbarsea_dbar";   break;
          case (kHEDIS_vbar_nc_p_ubarsea_ubar) : return "nubar_NC_p_ubarsea_ubar";   break;
          case (kHEDIS_vbar_nc_p_sbarsea_sbar) : return "nubar_NC_p_sbarsea_sbar";   break;
          case (kHEDIS_vbar_nc_p_cbarsea_cbar) : return "nubar_NC_p_cbarsea_cbar";   break;
          case (kHEDIS_vbar_nc_p_bbarsea_bbar) : return "nubar_NC_p_bbarsea_bbar";   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return "nubar_NC_n_dval_d";   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return "nubar_NC_n_uval_u";   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return "nubar_NC_n_dsea_d";   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return "nubar_NC_n_usea_u";   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return "nubar_NC_n_ssea_s";   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return "nubar_NC_n_csea_c";   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return "nubar_NC_n_bsea_b";   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return "nubar_NC_n_dbarsea_dbar";   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return "nubar_NC_n_ubarsea_ubar";   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return "nubar_NC_n_sbarsea_sbar";   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return "nubar_NC_n_cbarsea_cbar";   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return "nubar_NC_n_bbarsea_bbar";   break;


          default : return "Unknown";  break;
        }
        return "Unknown";
      }
      //__________________________________________________________________________
      static InteractionType_t InteractionType(HEDISChannel_t channel)
      {

        switch (channel) {
          case (kHEDIS_v_nc_p_dval_d)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_uval_u)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_usea_u)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_csea_c)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_p_bbarsea_bbar) : return kIntWeakNC;   break;

          case (kHEDIS_v_nc_n_dval_d)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_uval_u)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_usea_u)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_csea_c)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDIS_v_nc_n_bbarsea_bbar) : return kIntWeakNC;   break;

          case (kHEDIS_vbar_nc_p_dval_d)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_uval_u)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_usea_u)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_csea_c)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_p_bbarsea_bbar) : return kIntWeakNC;   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return kIntWeakNC;   break;

          default : return kIntWeakCC;  break;
        }
        return kIntNull;
      }
      //__________________________________________________________________________
      static bool IsNu(HEDISChannel_t channel)
      {
        switch (channel) {
          case (kHEDIS_vbar_cc_p_uval_d)       : return false;   break;
          case (kHEDIS_vbar_cc_p_uval_s)       : return false;   break;
          case (kHEDIS_vbar_cc_p_uval_b)       : return false;   break;
          case (kHEDIS_vbar_cc_p_usea_d)       : return false;   break;
          case (kHEDIS_vbar_cc_p_usea_s)       : return false;   break;
          case (kHEDIS_vbar_cc_p_usea_b)       : return false;   break;
          case (kHEDIS_vbar_cc_p_csea_d)       : return false;   break;
          case (kHEDIS_vbar_cc_p_csea_s)       : return false;   break;
          case (kHEDIS_vbar_cc_p_csea_b)       : return false;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_tbar) : return false;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_tbar) : return false;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_tbar) : return false;   break;

          case (kHEDIS_vbar_cc_n_uval_d)       : return false;   break;
          case (kHEDIS_vbar_cc_n_uval_s)       : return false;   break;
          case (kHEDIS_vbar_cc_n_uval_b)       : return false;   break;
          case (kHEDIS_vbar_cc_n_usea_d)       : return false;   break;
          case (kHEDIS_vbar_cc_n_usea_s)       : return false;   break;
          case (kHEDIS_vbar_cc_n_usea_b)       : return false;   break;
          case (kHEDIS_vbar_cc_n_csea_d)       : return false;   break;
          case (kHEDIS_vbar_cc_n_csea_s)       : return false;   break;
          case (kHEDIS_vbar_cc_n_csea_b)       : return false;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_tbar) : return false;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_tbar) : return false;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_tbar) : return false;   break;

          case (kHEDIS_vbar_nc_p_dval_d)       : return false;   break;
          case (kHEDIS_vbar_nc_p_uval_u)       : return false;   break;
          case (kHEDIS_vbar_nc_p_dsea_d)       : return false;   break;
          case (kHEDIS_vbar_nc_p_usea_u)       : return false;   break;
          case (kHEDIS_vbar_nc_p_ssea_s)       : return false;   break;
          case (kHEDIS_vbar_nc_p_csea_c)       : return false;   break;
          case (kHEDIS_vbar_nc_p_bsea_b)       : return false;   break;
          case (kHEDIS_vbar_nc_p_dbarsea_dbar) : return false;   break;
          case (kHEDIS_vbar_nc_p_ubarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_nc_p_sbarsea_sbar) : return false;   break;
          case (kHEDIS_vbar_nc_p_cbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_nc_p_bbarsea_bbar) : return false;   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return false;   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return false;   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return false;   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return false;   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return false;   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return false;   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return false;   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return false;   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return false;   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return false;   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return false;   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return false;   break;

          default : return true;  break;
        }
        return true;
      }

      //__________________________________________________________________________
      static int HitNuclPdg(HEDISChannel_t channel)
      {

        switch (channel) {

          case (kHEDIS_v_cc_n_dval_u)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_dval_c)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_dval_t)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_dsea_u)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_dsea_c)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_dsea_t)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_ssea_u)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_ssea_c)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_ssea_t)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_bsea_u)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_bsea_c)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_bsea_t)       : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_ubarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_ubarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_ubarsea_bbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_cbarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_cbarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_cc_n_cbarsea_bbar) : return kPdgNeutron;   break;

          case (kHEDIS_vbar_cc_n_uval_d)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_uval_s)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_uval_b)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_usea_d)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_usea_s)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_usea_b)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_csea_d)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_csea_s)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_csea_b)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_tbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_tbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_tbar) : return kPdgNeutron;   break;

          case (kHEDIS_v_nc_n_dval_d)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_uval_u)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_dsea_d)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_usea_u)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_ssea_s)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_csea_c)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_bsea_b)       : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_dbarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_ubarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_sbarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_cbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDIS_v_nc_n_bbarsea_bbar) : return kPdgNeutron;   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return kPdgNeutron;   break;

          default : return kPdgProton;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      static bool HitQuarkSea(HEDISChannel_t channel)
      {
        switch (channel) {
          case (kHEDIS_v_cc_p_dval_u)       : return false;  break;
          case (kHEDIS_v_cc_p_dval_c)       : return false;  break;
          case (kHEDIS_v_cc_p_dval_t)       : return false;  break;
          case (kHEDIS_v_cc_n_dval_u)       : return false;  break;
          case (kHEDIS_v_cc_n_dval_c)       : return false;  break;
          case (kHEDIS_v_cc_n_dval_t)       : return false;  break;
          case (kHEDIS_vbar_cc_p_uval_d)    : return false;  break;
          case (kHEDIS_vbar_cc_p_uval_s)    : return false;  break;
          case (kHEDIS_vbar_cc_p_uval_b)    : return false;  break;
          case (kHEDIS_vbar_cc_n_uval_d)    : return false;  break;
          case (kHEDIS_vbar_cc_n_uval_s)    : return false;  break;
          case (kHEDIS_vbar_cc_n_uval_b)    : return false;  break;
          case (kHEDIS_v_nc_p_dval_d)       : return false;  break;
          case (kHEDIS_v_nc_p_uval_u)       : return false;  break;
          case (kHEDIS_v_nc_n_dval_d)       : return false;  break;
          case (kHEDIS_v_nc_n_uval_u)       : return false;  break;
          case (kHEDIS_vbar_nc_p_dval_d)    : return false;  break;
          case (kHEDIS_vbar_nc_p_uval_u)    : return false;  break;
          case (kHEDIS_vbar_nc_n_dval_d)    : return false;  break;
          case (kHEDIS_vbar_nc_n_uval_u)    : return false;  break;

          default : return true;  break;
        }
        return true;
      }
      //__________________________________________________________________________
      static int HitQuarkPdg(HEDISChannel_t channel)
      {
        switch (channel) {
          case (kHEDIS_v_cc_p_dval_u)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_p_dval_c)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_p_dval_t)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_p_dsea_u)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_p_dsea_c)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_p_dsea_t)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_p_ssea_u)       : return kPdgSQuark;   break;
          case (kHEDIS_v_cc_p_ssea_c)       : return kPdgSQuark;   break;
          case (kHEDIS_v_cc_p_ssea_t)       : return kPdgSQuark;   break;
          case (kHEDIS_v_cc_p_bsea_u)       : return kPdgBQuark;   break;
          case (kHEDIS_v_cc_p_bsea_c)       : return kPdgBQuark;   break;
          case (kHEDIS_v_cc_p_bsea_t)       : return kPdgBQuark;   break;
          case (kHEDIS_v_cc_p_ubarsea_dbar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_cc_p_ubarsea_sbar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_cc_p_ubarsea_bbar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_cc_p_cbarsea_dbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_cc_p_cbarsea_sbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_cc_p_cbarsea_bbar) : return kPdgAntiCQuark;   break;

          case (kHEDIS_v_cc_n_dval_u)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_n_dval_c)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_n_dval_t)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_n_dsea_u)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_n_dsea_c)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_n_dsea_t)       : return kPdgDQuark;   break;
          case (kHEDIS_v_cc_n_ssea_u)       : return kPdgSQuark;   break;
          case (kHEDIS_v_cc_n_ssea_c)       : return kPdgSQuark;   break;
          case (kHEDIS_v_cc_n_ssea_t)       : return kPdgSQuark;   break;
          case (kHEDIS_v_cc_n_bsea_u)       : return kPdgBQuark;   break;
          case (kHEDIS_v_cc_n_bsea_c)       : return kPdgBQuark;   break;
          case (kHEDIS_v_cc_n_bsea_t)       : return kPdgBQuark;   break;
          case (kHEDIS_v_cc_n_ubarsea_dbar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_cc_n_ubarsea_sbar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_cc_n_ubarsea_bbar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_cc_n_cbarsea_dbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_cc_n_cbarsea_sbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_cc_n_cbarsea_bbar) : return kPdgAntiCQuark;   break;

          case (kHEDIS_vbar_cc_p_uval_d)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_p_uval_s)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_p_uval_b)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_p_usea_d)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_p_usea_s)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_p_usea_b)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_p_csea_d)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_cc_p_csea_s)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_cc_p_csea_b)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_ubar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_cbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_tbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_ubar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_cbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_tbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_ubar) : return kPdgAntiBQuark;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_cbar) : return kPdgAntiBQuark;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_tbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_vbar_cc_n_uval_d)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_n_uval_s)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_n_uval_b)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_n_usea_d)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_n_usea_s)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_n_usea_b)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_cc_n_csea_d)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_cc_n_csea_s)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_cc_n_csea_b)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_ubar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_cbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_tbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_ubar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_cbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_tbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_ubar) : return kPdgAntiBQuark;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_cbar) : return kPdgAntiBQuark;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_tbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_v_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_v_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_v_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_v_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_v_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_v_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_vbar_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          default : return 0;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      static int FnlQuarkPdg(HEDISChannel_t channel)
      {
        switch (channel) {
          case (kHEDIS_v_cc_p_dval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_p_dval_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_p_dval_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_p_dsea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_p_dsea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_p_dsea_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_p_ssea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_p_ssea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_p_ssea_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_p_bsea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_p_bsea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_p_bsea_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_p_ubarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_cc_p_ubarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_cc_p_ubarsea_bbar) : return kPdgAntiBQuark;   break;
          case (kHEDIS_v_cc_p_cbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_cc_p_cbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_cc_p_cbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_v_cc_n_dval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_n_dval_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_n_dval_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_n_dsea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_n_dsea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_n_dsea_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_n_ssea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_n_ssea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_n_ssea_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_n_bsea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_cc_n_bsea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_cc_n_bsea_t)       : return kPdgTQuark;   break;
          case (kHEDIS_v_cc_n_ubarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_cc_n_ubarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_cc_n_ubarsea_bbar) : return kPdgAntiBQuark;   break;
          case (kHEDIS_v_cc_n_cbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_cc_n_cbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_cc_n_cbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_vbar_cc_p_uval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_cc_p_uval_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_cc_p_uval_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_cc_p_usea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_cc_p_usea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_cc_p_usea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_cc_p_csea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_cc_p_csea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_cc_p_csea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_tbar) : return kPdgAntiTQuark;   break;

          case (kHEDIS_vbar_cc_n_uval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_cc_n_uval_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_cc_n_uval_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_cc_n_usea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_cc_n_usea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_cc_n_usea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_cc_n_csea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_cc_n_csea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_cc_n_csea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_tbar) : return kPdgAntiTQuark;   break;

          case (kHEDIS_v_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_v_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_v_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_v_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_v_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_v_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_v_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_v_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_v_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_v_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_v_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_v_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_v_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_vbar_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          default : return 0;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      //__________________________________________________________________________
      static HEDISNuclChannel_t HEDISNuclChannel(HEDISChannel_t channel)
      {
        switch (channel) {
          case (kHEDIS_v_cc_p_dval_u)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_dval_c)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_dval_t)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_dsea_u)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_dsea_c)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_dsea_t)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_ssea_u)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_ssea_c)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_ssea_t)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_bsea_u)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_bsea_c)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_bsea_t)       : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_ubarsea_dbar) : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_ubarsea_sbar) : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_ubarsea_bbar) : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_cbarsea_dbar) : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_cbarsea_sbar) : return kHEDISNucl_v_cc_p;   break;
          case (kHEDIS_v_cc_p_cbarsea_bbar) : return kHEDISNucl_v_cc_p;   break;

          case (kHEDIS_v_cc_n_dval_u)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_dval_c)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_dval_t)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_dsea_u)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_dsea_c)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_dsea_t)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_ssea_u)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_ssea_c)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_ssea_t)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_bsea_u)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_bsea_c)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_bsea_t)       : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_ubarsea_dbar) : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_ubarsea_sbar) : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_ubarsea_bbar) : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_cbarsea_dbar) : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_cbarsea_sbar) : return kHEDISNucl_v_cc_n;   break;
          case (kHEDIS_v_cc_n_cbarsea_bbar) : return kHEDISNucl_v_cc_n;   break;

          case (kHEDIS_vbar_cc_p_uval_d)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_uval_s)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_uval_b)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_usea_d)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_usea_s)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_usea_b)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_csea_d)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_csea_s)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_csea_b)       : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_ubar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_cbar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_dbarsea_tbar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_ubar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_cbar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_sbarsea_tbar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_ubar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_cbar) : return kHEDISNucl_vbar_cc_p;   break;
          case (kHEDIS_vbar_cc_p_bbarsea_tbar) : return kHEDISNucl_vbar_cc_p;   break;

          case (kHEDIS_vbar_cc_n_uval_d)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_uval_s)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_uval_b)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_usea_d)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_usea_s)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_usea_b)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_csea_d)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_csea_s)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_csea_b)       : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_ubar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_cbar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_dbarsea_tbar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_ubar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_cbar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_sbarsea_tbar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_ubar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_cbar) : return kHEDISNucl_vbar_cc_n;   break;
          case (kHEDIS_vbar_cc_n_bbarsea_tbar) : return kHEDISNucl_vbar_cc_n;   break;

          case (kHEDIS_v_nc_p_dval_d)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_uval_u)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_dsea_d)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_usea_u)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_ssea_s)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_csea_c)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_bsea_b)       : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_dbarsea_dbar) : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_ubarsea_ubar) : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_sbarsea_sbar) : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_cbarsea_cbar) : return kHEDISNucl_v_nc_p;   break;
          case (kHEDIS_v_nc_p_bbarsea_bbar) : return kHEDISNucl_v_nc_p;   break;

          case (kHEDIS_v_nc_n_dval_d)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_uval_u)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_dsea_d)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_usea_u)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_ssea_s)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_csea_c)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_bsea_b)       : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_dbarsea_dbar) : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_ubarsea_ubar) : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_sbarsea_sbar) : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_cbarsea_cbar) : return kHEDISNucl_v_nc_n;   break;
          case (kHEDIS_v_nc_n_bbarsea_bbar) : return kHEDISNucl_v_nc_n;   break;

          case (kHEDIS_vbar_nc_p_dval_d)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_uval_u)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_dsea_d)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_usea_u)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_ssea_s)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_csea_c)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_bsea_b)       : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_dbarsea_dbar) : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_ubarsea_ubar) : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_sbarsea_sbar) : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_cbarsea_cbar) : return kHEDISNucl_vbar_nc_p;   break;
          case (kHEDIS_vbar_nc_p_bbarsea_bbar) : return kHEDISNucl_vbar_nc_p;   break;

          case (kHEDIS_vbar_nc_n_dval_d)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_uval_u)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_dsea_d)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_usea_u)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_ssea_s)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_csea_c)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_bsea_b)       : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_dbarsea_dbar) : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_ubarsea_ubar) : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_sbarsea_sbar) : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_cbarsea_cbar) : return kHEDISNucl_vbar_nc_n;   break;
          case (kHEDIS_vbar_nc_n_bbarsea_bbar) : return kHEDISNucl_vbar_nc_n;   break;

          default : return kHEDISNuclNull;  break;
        }
        return kHEDISNuclNull;
      }
      //__________________________________________________________________________





  };

}      // genie namespace

#endif // _HEDIS_CHANNEL_H_

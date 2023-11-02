// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \brief Table definitions for Heavy flavor Tagging jets
///
/// \file   JetTagging.h
/// \author Hanseo Park

#ifndef PWGJE_DATAMODEL_JETTAGGING_H
#define PWGJE_DATAMODEL_JETTAGGING_H

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"

namespace o2::aod
{
// selection of heavy flavor tagging
namespace jet_hf_tag_sel
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(IslfJet, islfjet, int);
DECLARE_SOA_COLUMN(IscJet, iscjet, int);
DECLARE_SOA_COLUMN(IsbJet, isbjet, int);
} // namespace jet_hf_tag_sel

// general properties of tagged constituents 
namespace jet_hf_tag_constituents
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(SignVertex, signVertex, int);
DECLARE_SOA_COLUMN(ImpXY, impXY, float);
DECLARE_SOA_DYNAMIC_COLUMN(SignImpXY, signImpXY,
                          [](int sign, float impXY) -> float { return sign * impXY;});
DECLARE_SOA_COLUMN(ImpXYSig, impXYSig, float);
DECLARE_SOA_DYNAMIC_COLUMN(SignImpXYSig, signImpXYSig,
                          [](int sign, float impXYSig) -> float { return sign * impXYSig;});
DECLARE_SOA_COLUMN(ImpXYZ, impXYZ, float);
DECLARE_SOA_DYNAMIC_COLUMN(SignImpXYZ, signImpXYZ,
                          [](int sign, float impXYZ) -> float { return sign * impXYZ;});
DECLARE_SOA_COLUMN(ImpXYZSig, impXYZSig, float);
DECLARE_SOA_DYNAMIC_COLUMN(SignImpXYZSig, signImpXYZSig,
                          [](int sign, float impXYZSig) -> float { return sign * impXYZSig;});

} // namespace jet_hf_tag_constituents
DECLARE_SOA_TABLE(TagJetConstBase, "AOD", "TAGJETCONSTBASE",
                  o2::soa::Index<>,
                  jet_hf_tag_constituents::CollisionId,
                  jet_hf_tag_constituents::SignVertex,
                  jet_hf_tag_constituents::ImpXY,
                  jet_hf_tag_constituents::ImpXYSig,
                  jet_hf_tag_constituents::SignImpXY<jet_hf_tag_constituents::SignVertex, jet_hf_tag_constituents::ImpXY>,
                  jet_hf_tag_constituents::SignImpXYSig<jet_hf_tag_constituents::SignVertex, jet_hf_tag_constituents::ImpXYSig>);

// track and jet probability of tagged constituents
namespace jet_hf_tag_constituents_ext_jp
{
DECLARE_SOA_COLUMN(TP, tp, float); // track probability
DECLARE_SOA_COLUMN(JP, jp, float); // Jet probability
} // namespace jet_hf_tag_constituents_ext_jp

DECLARE_SOA_TABLE(TagJetConstExtJP, "AOD", "TAGJETCONSTJP",
                  jet_hf_tag_constituents_ext_jp::TP,
                  jet_hf_tag_constituents_ext_jp::JP);

using TagJetConst = soa::Join<TagJetConstBase, TagJetConstExtJP>;

    //DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(TagJet, tagjet, int32_t, TagJetConstBase, "_tagtrack"); \
                    //_name_##constituents::TagJetConstBase##Ids, \

// Defines the tagged jet table
#define CONSTITUENTS_TABLE_DEF(_jet_type_, _name_, _Description_, _track_type_)    \
  namespace _name_##constituents \
  {   \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                                                  \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_track_type_, trackstemp); \
  } \
  DECLARE_SOA_TABLE(_jet_type_##Constituents, "AOD", _Description_ "CONSTS", \
                    _name_##constituents::_jet_type_##Id, \
                    _name_##constituents::_track_type_##Ids);

//#define CONSTITUENTS_BASIC_TABLE_DEF(_jet_type_, _name_, _Description_) \
//  DECLARE_SOA_TABLE(_jet_type_##ConstituentsBasic, "AOD", _Description_ "CONSTTAGBAS") \
//                    _name_##constituents::_jet_type_##Id, \
//                    jet_hf_tag_constituents::SignVertex, \
//                    jet_hf_tag_constituents::ImpXY, \
//                    jet_hf_tag_constituents::ImpXYSig, \
//                    jet_hf_tag_constituents::SignImpXY<jet_hf_tag_constituents::SignVertex, jet_hf_tag_constituents::ImpXY>, \
//                    jet_hf_tag_constituents::SignImpXYSig<jet_hf_tag_constituents::SignVertex, jet_hf_tag_constituents::ImpXYSig>);
//
////#define CONSTITUENTS_EXT_TABLE_DEF(_jet_type_, _name_, _Description_ )

#define JETTAGGING_TABLE_DEF(_collision_name_, _jet_type_, _track_type_, _description_) \
  DECLARE_SOA_TABLE(_jet_type_##s, "AOD", _description_, \
                    o2::soa::Index<>, \
                    jet_hf_tag_sel::IslfJet, \
                    jet_hf_tag_sel::IscJet, \
                    jet_hf_tag_sel::IsbJet); \

#define JETTAGGING_TABLES_DEF(_collision_name_, _jet_type_, _track_type_, _description_) \
  JETTAGGING_TABLE_DEF(_collision_name_, _jet_type_##Jet, _track_type_, _description_);       \
  using _jet_type_##Jet = _jet_type_##Jet##s::iterator; \
  CONSTITUENTS_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_, _track_type_); \
  using _jet_type_##Jet##Constituent = _jet_type_##Jet##Constituents::iterator; \

JETTAGGING_TABLES_DEF(Collision, TaggedMCDetectorLevel, Track, "TAGMCD");
//JETTAGGING_TABLES_DEF(McCollision, TaggedMCParticleLevel, McParticle, "TAGMCP");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H

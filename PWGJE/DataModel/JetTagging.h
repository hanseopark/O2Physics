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
/// \brief Table definitions for hf jet tagging
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETTAGGING_H_
#define PWGJE_DATAMODEL_JETTAGGING_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2::analysis;

namespace o2::aod
{

// Defines tagger track extention
namespace jtracktagext
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_COLUMN(DCAX, dcaX, float);
DECLARE_SOA_COLUMN(DCAY, dcaY, float);
} // namespace jtracktagext

DECLARE_SOA_TABLE(JTracksTagExt, "AOD", "JTracksTagExt",
                  o2::soa::Index<>,
                  jtracktagext::CollisionId,
                  jtracktagext::DCAX,
                  jtracktagext::DCAY);

using JTrackTagExt = JTracksTagExt::iterator;

DECLARE_SOA_TABLE(JTrackTagExtPIs, "AOD", "JTrackTagExtPIs",
                  jtracktagext::TrackId);

// Defines tagger trackIU extention
namespace jtrackiutagext
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_COLUMN(DCAX, dcaX, float);
DECLARE_SOA_COLUMN(DCAY, dcaY, float);
} // namespace jtrackiutagext

DECLARE_SOA_TABLE(JTracksIUTagExt, "AOD", "JTracksIUTagExt",
                  o2::soa::Index<>,
                  jtrackiutagext::CollisionId,
                  jtrackiutagext::DCAX,
                  jtrackiutagext::DCAY);

using JTrackIUTagExt = JTracksIUTagExt::iterator;

DECLARE_SOA_TABLE(JTrackIUTagExtPIs, "AOD", "JTraIUTagExtPIs",
                  jtrackiutagext::TrackId);

// Defines the jet substrcuture table definition
#define JETTAGGING_TABLE_DEF(_jet_type_, _name_, _description_) \
  namespace _name_##tagging                                     \
  {                                                             \
    DECLARE_SOA_COLUMN(Origin, origin, int);                    \
    DECLARE_SOA_COLUMN(Algorithm1, algorithm1, int);            \
    DECLARE_SOA_COLUMN(Algorithm2, algorithm2, int);            \
    DECLARE_SOA_COLUMN(Algorithm3, algorithm3, int);            \
  }                                                             \
  DECLARE_SOA_TABLE(_jet_type_##Tags, "AOD", _description_ "Tags", _name_##tagging::Origin, _name_##tagging::Algorithm1, _name_##tagging::Algorithm2, _name_##tagging::Algorithm3);

#define JETTAGGING_TABLES_DEF(_jet_type_, _description_)                                                    \
  JETTAGGING_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_)                                     \
  JETTAGGING_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD") \
  JETTAGGING_TABLE_DEF(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _description_ "MCP")

JETTAGGING_TABLES_DEF(Charged, "C");
JETTAGGING_TABLES_DEF(Full, "F");
JETTAGGING_TABLES_DEF(Neutral, "N");
JETTAGGING_TABLES_DEF(D0Charged, "D0");
JETTAGGING_TABLES_DEF(LcCharged, "Lc");
JETTAGGING_TABLES_DEF(BplusCharged, "BPL");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H_

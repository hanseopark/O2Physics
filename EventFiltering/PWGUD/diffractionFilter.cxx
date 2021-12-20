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
// O2 includes
///
/// \brief A filter task for diffractive events
///        requires: EvSels, o2-analysis-event-selection
///                  TrackSelection, o2-analysis-trackselection
///                  TracksExtended, o2-analysis-trackextension
///                  pidTOF*, o2-analysis-pid-tof
///                  pidTPC*, o2-analysis-pid-tpc
///                  Timestamps, o2-analysis-timestamp
///        usage: o2-analysis-timestamp --aod-file AO2D.root |
///               o2-analysis-event-selection --processRun3 1 --isMC 1 |
///               o2-analysis-trackextension |
///               o2-analysis-trackselection --isRun3 |
///               o2-analysis-pid-tof |
///               o2-analysis-pid-tpc |
///               o2-analysis-diffraction-filter
/// \author P. Buehler , paul.buehler@oeaw.ac.at
/// \since June 1, 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/LHCConstants.h"

#include "cutHolder.h"
#include "diffractionSelectors.h"
#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Run 3
struct DGFilterRun3 {

  // Productions
  Produces<aod::DiffractionFilters> filterTable;

  // configurable cutHolder
  MutableConfigurable<cutHolder> diffCuts{"DiffCuts", {}, "Diffractive events cuts"};

  // DG selector
  DGSelector dgSelector;

  HistogramRegistry registry{
    "registry",
    {
      {"aftercut", "#aftercut", {HistType::kTH1F, {{7, -0.5, 6.5}}}},
    }};

  // some general Collisions and Tracks filter

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::TOFSignal, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  // using MFs = aod::MFTTracks;
  using FWs = aod::FwdTracks;

  void process(CC const& collision,
               BCs const& bcs,
               TCs& tracks,
               // MFs& mfttracks,
               FWs& fwdtracks,
               aod::Zdcs& zdcs,
               aod::FT0s& ft0s,
               aod::FV0As& fv0as,
               aod::FV0Cs& fv0cs,
               aod::FDDs& fdds)
  {
    // nominal BC
    auto bc = collision.bc_as<BCs>();

    // obtain slice of compatible BCs
    auto bcRange = compatibleBCs(collision, diffCuts->NDtcoll(), bcs);
    LOGF(info, "  Number of compatible BCs in +- %i dtcoll: %i", diffCuts->NDtcoll(), bcRange.size());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bc, bcRange, tracks, fwdtracks);

    // fill filterTable
    if (isDGEvent == 0) {
      LOGF(info, "This collision is a DG candidate!");
    }
    filterTable(isDGEvent == 0);

    // update histogram
    registry.get<TH1>(HIST("aftercut"))->Fill(isDGEvent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGFilterRun3>(cfgc, TaskName{"DGfilterRun3"}),
  };
}

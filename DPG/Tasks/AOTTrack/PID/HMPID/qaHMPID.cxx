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
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <TTree.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct pidHmpidQa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsP{"nBinsP", 500, "Number of momentum bins"};
  Configurable<float> minP{"minP", 0.01f, "Minimum momentum plotted (GeV/c)"};
  Configurable<float> maxP{"maxP", 10.f, "Maximum momentum plotted (GeV/c)"};
  Configurable<float> maxDCA{"maxDCA", 3.f, "Maximum DCA xy use for the plot (cm)"};
  Configurable<float> maxDistance{"maxDistance", 5.f, "Maximum HMPID distance between the track and the cluster (cm)"};
  Configurable<float> minCharge{"minCharge", 120.f, "Minimum HMPID charge collected in the cluster"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec momAxis{nBinsP, minP, maxP};

    histos.add("hmpidSignal", "hmpidSignal", kTH1F, {{1000, 0, 1}});
    histos.add("hmpidMomvsTrackMom", "hmpidMomvsTrackMom", kTH2F, {{1200, 0, 30, "Track #it{p} (GeV/#it{c})"}, {1200, 0, 30, "HMP #it{p} (GeV/#it{c})"}});
    histos.add("PhivsEta", "PhivsEta", kTH2F, {{550, -0.55, 0.55, "#eta"}, {550, 0, 1.1, "#phi (rad)"}});
    histos.add("hmpidCkovvsMom", "hmpidCkovvsMom", kTH2F, {{500, 0, 10, "#it{p} (GeV/#it{c})"}, {1000, 0, 1, "Cherenkov angle (rad)"}});
    histos.add("hmpidXTrack", "hmpidXTrack", kTH1F, {{280, 0, 140, "X track (cm)"}});
    histos.add("hmpidYTrack", "hmpidYTrack", kTH1F, {{280, 0, 140, "Y track (cm)"}});
    histos.add("hmpidXMip", "hmpidXMip", kTH1F, {{280, 0, 140, "X mip (cm)"}});
    histos.add("hmpidYMip", "hmpidYMip", kTH1F, {{280, 0, 140, "X mip (cm)"}});
    histos.add("hmpidXResiduals", "hmpidXResiduals", kTH1F, {{400, -20, 20, "X Residuals (cm)"}});
    histos.add("hmpidYResiduals", "hmpidYResiduals", kTH1F, {{400, -20, 20, "Y Residuals (cm)"}});
    histos.add("hmpidNPhotons", "hmpidNPhotons", kTH1F, {{50, 0, 50, "Number of photons"}});
    histos.add("hmpidQMip", "hmpidQMip", kTH1F, {{2000, 200, 2200, "Charge (ADCD)"}});
    histos.add("hmpidClusSize", "hmpidClusSize", kTH1F, {{15, 0, 15, "MIP Cluster size"}});
    histos.add("TrackMom", "TrackMom", kTH1F, {{1200, -30, 30, "#it{p} (GeV/#it{c})"}});
    histos.add("hmpidMom", "hmpidMom", kTH1F, {{1200, -30, 30, "#it{p} (GeV/#it{c})"}});
    histos.add("hmpidPhotsCharge", "hmpidPhotsCharge", kTH1F, {{300, 0, 300}});
    for (int iCh = 0; iCh < 7; iCh++) {
      histos.add(Form("hmpidXTrack%i", iCh), Form("hmpidXTrack%i", iCh), kTH1F, {{280, 0, 140, "X track (cm)"}});
      histos.add(Form("hmpidYTrack%i", iCh), Form("hmpidYTrack%i", iCh), kTH1F, {{280, 0, 140, "Y track (cm)"}});
      histos.add(Form("hmpidXMip%i", iCh), Form("hmpidXMip%i", iCh), kTH1F, {{280, 0, 140, "X mip (cm)"}});
      histos.add(Form("hmpidYMip%i", iCh), Form("hmpidYMip%i", iCh), kTH1F, {{280, 0, 140, "X mip (cm)"}});
      histos.add(Form("hmpidXResiduals%i", iCh), Form("hmpidXResiduals%i", iCh), kTH1F, {{400, -20, 20, "X Residuals (cm)"}});
      histos.add(Form("hmpidYResiduals%i", iCh), Form("hmpidYResiduals%i", iCh), kTH1F, {{400, -20, 20, "Y Residuals (cm)"}});
      histos.add(Form("hmpidNPhotons%i", iCh), Form("hmpidNPhotons%i", iCh), kTH1F, {{50, 0, 50, "Number of photons"}});
      histos.add(Form("hmpidQMip%i", iCh), Form("hmpidQMip%i", iCh), kTH1F, {{2000, 200, 2200, "Charge (ADCD)"}});
      histos.add(Form("hmpidClusSize%i", iCh), Form("hmpidClusSize%i", iCh), kTH1F, {{15, 0, 15, "MIP Cluster size"}});
      histos.add(Form("TrackMom%i", iCh), Form("TrackMom%i", iCh), kTH1F, {{1200, -30, 30, "#it{p} (GeV/#it{c})"}});
      histos.add(Form("hmpidMom%i", iCh), Form("hmpidMom%i", iCh), kTH1F, {{1200, -30, 30, "#it{p} (GeV/#it{c})"}});
      histos.add(Form("hmpidPhotsCharge%i", iCh), Form("hmpidPhotsCharge%i", iCh), kTH1F, {{300, 0, 300}});
    }
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  void process(const aod::HMPIDs& hmpids,
               const TrackCandidates& /*tracks*/,
               const aod::Collisions& /*colls*/)
  {

    for (const auto& t : hmpids) {
      if (t.track_as<TrackCandidates>().isGlobalTrack() != (uint8_t) true) {
        continue;
      }
      if (abs(t.track_as<TrackCandidates>().dcaXY()) > maxDCA) {
        continue;
      }

      histos.fill(HIST("hmpidSignal"), t.hmpidSignal());
      histos.fill(HIST("PhivsEta"), t.track_as<TrackCandidates>().eta(), t.track_as<TrackCandidates>().phi());
      histos.fill(HIST("hmpidMomvsTrackMom"), t.track_as<TrackCandidates>().p(), abs(t.hmpidMom()));
      histos.fill(HIST("hmpidCkovvsMom"), abs(t.hmpidMom()), t.hmpidSignal());
      histos.fill(HIST("hmpidXTrack"), t.hmpidYTrack());
      histos.fill(HIST("hmpidYTrack"), t.hmpidYTrack());
      histos.fill(HIST("hmpidXMip"), t.hmpidXMip());
      histos.fill(HIST("hmpidYMip"), t.hmpidYMip());
      if (t.track_as<TrackCandidates>().p() > 1.5) {
        histos.fill(HIST("hmpidXResiduals"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals"), t.hmpidYMip() - t.hmpidYTrack());
      }
      histos.fill(HIST("hmpidNPhotons"), t.hmpidNPhotons());
      histos.fill(HIST("hmpidQMip"), t.hmpidQMip());
      histos.fill(HIST("hmpidClusSize"), (t.hmpidClusSize() % 1000000) / 1000);
      histos.fill(HIST("TrackMom"), t.track_as<TrackCandidates>().p());
      histos.fill(HIST("hmpidMom"), abs(t.hmpidMom()));
      for (int i = 0; i < 10; i++) {
        if (t.hmpidPhotsCharge()[i] > 0)
          histos.fill(HIST("hmpidPhotsCharge"), t.hmpidPhotsCharge()[i]);
      }

      if (t.hmpidClusSize() / 1000000 == 0) {
        histos.fill(HIST("hmpidXTrack0"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack0"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip0"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip0"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals0"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals0"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons0"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip0"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize0"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom0"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom0"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge0"), t.hmpidPhotsCharge()[i]);
        }
      }

      if (t.hmpidClusSize() / 1000000 == 1) {
        histos.fill(HIST("hmpidXTrack1"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack1"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip1"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip1"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals1"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals1"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons1"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip1"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize1"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom1"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom1"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge1"), t.hmpidPhotsCharge()[i]);
        }
      }

      if (t.hmpidClusSize() / 1000000 == 2) {
        histos.fill(HIST("hmpidXTrack2"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack2"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip2"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip2"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals2"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals2"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons2"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip2"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize2"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom2"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom2"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge2"), t.hmpidPhotsCharge()[i]);
        }
      }

      if (t.hmpidClusSize() / 1000000 == 3) {
        histos.fill(HIST("hmpidXTrack3"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack3"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip3"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip3"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals3"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals3"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons3"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip3"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize3"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom3"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom3"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge3"), t.hmpidPhotsCharge()[i]);
        }
      }

      if (t.hmpidClusSize() / 1000000 == 4) {
        histos.fill(HIST("hmpidXTrack4"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack4"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip4"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip4"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals4"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals4"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons4"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip4"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize4"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom4"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom4"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge4"), t.hmpidPhotsCharge()[i]);
        }
      }

      if (t.hmpidClusSize() / 1000000 == 5) {
        histos.fill(HIST("hmpidXTrack5"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack5"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip5"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip5"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals5"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals5"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons5"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip5"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize5"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom5"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom5"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge5"), t.hmpidPhotsCharge()[i]);
        }
      }

      if (t.hmpidClusSize() / 1000000 == 6) {
        histos.fill(HIST("hmpidXTrack6"), t.hmpidXTrack());
        histos.fill(HIST("hmpidYTrack6"), t.hmpidYTrack());
        histos.fill(HIST("hmpidXMip6"), t.hmpidXMip());
        histos.fill(HIST("hmpidYMip6"), t.hmpidYMip());
        histos.fill(HIST("hmpidXResiduals6"), t.hmpidXMip() - t.hmpidXTrack());
        histos.fill(HIST("hmpidYResiduals6"), t.hmpidYMip() - t.hmpidYTrack());
        histos.fill(HIST("hmpidNPhotons6"), t.hmpidNPhotons());
        histos.fill(HIST("hmpidQMip6"), t.hmpidQMip());
        histos.fill(HIST("hmpidClusSize6"), (t.hmpidClusSize() % 1000000) / 1000);
        histos.fill(HIST("TrackMom6"), t.track_as<TrackCandidates>().p());
        histos.fill(HIST("hmpidMom6"), abs(t.hmpidMom()));
        for (int i = 0; i < 10; i++) {
          if (t.hmpidPhotsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge6"), t.hmpidPhotsCharge()[i]);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<pidHmpidQa>(cfg)}; }

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

// jet heavy flavor tagging task
//
// Authors: Hanseo Park

#include "TF1.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
//#include "PWGJE/TableProducer/jetHfTagging.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

template <typename T>
int getJetFlavor(T const& jet)
{
  if (jet.isbjet() > 0) return 2; // bjet
  else if (jet.iscjet() > 0) return 1; // cjet
  else return 0; // lfjet
}

namespace hf_tagging_jet_cut
{
  static constexpr int nBinsJetPt = 18;
  constexpr double binsJetPt[nBinsJetPt + 1] = {
    0,
    2.0,
    4.0,
    10.0,
    15.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
    45.0,
    50.0,
    55.0,
    60.0,
    65.0,
    70.0,
    80.0,
    90.0,
    100.0};
  auto vecBinsJetPt = std::vector<double>{binsJetPt, binsJetPt + nBinsJetPt + 1};

} // namespace hf_tagging_jet_cut

struct JetHfTaggingTask {

  static constexpr std::string_view charJetFlavor[] = {"inc_jet", "lfjet", "cjet", "bjet"};

  // Binning
  ConfigurableAxis binPdgCode{"binPdgCode", {5000, -2500.f, 2500.f}, ""};
  ConfigurableAxis binStatus{"binStatusCode", {200, -99.5f, 100.5f}, ""};
  ConfigurableAxis binJetPt{"binJetPt", {200, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, -0.5, 99.5}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {1000, -0.4f, 0.4f}, ""};
  ConfigurableAxis binImpactParameterSignificanceXY{"binImpactParameterSignificanceXY", {1000, -40.f, 40.f}, ""};
  ConfigurableAxis binImpactParameterXYZ{"binImpactParameterXYZ", {1000, -0.4f, 0.4f}, ""};
  //ConfigurableAxis binImpactParameterSignificanceXYZ{"binImpactParameterSignificanceXYZ", {1000, -40.f, 40.f}, ""};
  ConfigurableAxis binImpactParameterSignificanceXYZ{"binImpactParameterSignificanceXYZ", {1000, -0.1f, 0.1f}, ""};
  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbability{"binJetProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbabilityLog{"binJetProbabilityLog", {100, 0.f, 10.f}, ""};
  ConfigurableAxis binJetEntries{"binEntries", {3, 0.f, 3.f}, ""};
  Configurable<std::vector<double>> binsJetPt{"binsJetPt", std::vector<double>{hf_tagging_jet_cut::vecBinsJetPt}, "pT bin limits"};

  // Axis
  AxisSpec jetPtRangeAxis{(std::vector<double>)binsJetPt, "#it{p}_{T}^{Jet} (GeV/#it{c})"};
  AxisSpec pdgCodeAxis = {binPdgCode, "PDG code"};
  AxisSpec statusAxis = {binStatus, "Status code"};
  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
  AxisSpec etaAxis = {binEta, "#eta"};
  AxisSpec phiAxis = {binPhi, "#phi"};
  AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
  AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
  AxisSpec ImpactParameterXYAxis = {binImpactParameterXY, "IP_{XY} [cm]"};
  AxisSpec ImpactParameterSignificanceXYAxis = {binImpactParameterSignificanceXY, "IPs_{XY} [cm]"};
  AxisSpec ImpactParameterXYZAxis = {binImpactParameterXYZ, "IP_{XYZ} [cm]"};
  AxisSpec ImpactParameterSignificanceXYZAxis = {binImpactParameterSignificanceXYZ, "IPs_{XYZ} [cm]"};
  AxisSpec TrackProbabilityAxis = {binTrackProbability, "track probability"};
  AxisSpec JetProbabilityAxis = {binJetProbability, "JP"};
  AxisSpec JetProbabilityLogAxis = {binJetProbabilityLog, "-Log(JP)"};
  AxisSpec JetEntries = {binJetEntries, "lf=1, c=2, b=3"};

  TF1* fResoFunccjet = new TF1("fResoFunccjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1* fResoFuncbjet = new TF1("fResoFuncbjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1* fResoFunclfjet = new TF1("fResoFunclfjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1* fResoFuncincjet = new TF1("fResoFuncincjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    for (int i = 0; i < 4; i++) {
      // common
      registry.add(Form("h_%s_pt_jet_eta", charJetFlavor[i].data()), Form("%s #it{p}_{T}", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, etaAxis}}, true);
      registry.add(Form("h_%s_pt", charJetFlavor[i].data()), Form("%s #it{p}_{T}", charJetFlavor[i].data()), {HistType::kTH1F, {jetPtAxis}}, true);
      registry.add(Form("h_%s_eta", charJetFlavor[i].data()), Form("%s #eta", charJetFlavor[i].data()), {HistType::kTH1F, {etaAxis}}, true);
      registry.add(Form("h_%s_phi", charJetFlavor[i].data()), Form("%s #phi", charJetFlavor[i].data()), {HistType::kTH1F, {phiAxis}}, true);
      registry.add(Form("h_%s_ntracks", charJetFlavor[i].data()), Form("%s N tracks", charJetFlavor[i].data()), {HistType::kTH1F, {ntracksAxis}}, true);
      registry.add(Form("h_%s_track_pt", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH1F, {trackPtAxis}}, true);
      registry.add(Form("h_%s_track_eta", charJetFlavor[i].data()), Form("#eta of track in %s", charJetFlavor[i].data()), {HistType::kTH1F, {etaAxis}}, true);
      registry.add(Form("h_%s_track_phi", charJetFlavor[i].data()), Form("#phi of track in %s", charJetFlavor[i].data()), {HistType::kTH1F, {phiAxis}}, true);
      registry.add(Form("h2_%s_track_pt", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtAxis, trackPtAxis}}, true);
      registry.add(Form("h2_%s_track_eta", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtAxis, etaAxis}}, true);
      registry.add(Form("h2_%s_track_phi", charJetFlavor[i].data()), Form("#it{p}_{T} of track in %s", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtAxis, phiAxis}}, true);
      registry.add(Form("h_%s_particle_pdgcode", charJetFlavor[i].data()), Form("particles' pdg code in %s", charJetFlavor[i].data()), {HistType::kTH1F, {pdgCodeAxis}});
      registry.add(Form("h_%s_particle_statuscode", charJetFlavor[i].data()), Form("particles' status code in %s", charJetFlavor[i].data()), {HistType::kTH1F, {statusAxis}});
      registry.add(Form("h_%s_impact_parameter_xy", charJetFlavor[i].data()), Form("%s impact parameter dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xy", charJetFlavor[i].data()), Form("%s sign impact parameter dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYAxis}});
      registry.add(Form("h_%s_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_pt_sign_impact_parameter_xy_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xy}", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYAxis}});
      registry.add(Form("h_%s_impact_parameter_xyz", charJetFlavor[i].data()), Form("%s impact parameter dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYZAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xyz", charJetFlavor[i].data()), Form("%s sign impact parameter dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterXYZAxis}});
      registry.add(Form("h_%s_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
      registry.add(Form("h_%s_sign_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
      registry.add(Form("h_%s_pt_sign_impact_parameter_xyz_significance", charJetFlavor[i].data()), Form("%s sign impact parameter significance dca_{xyz}", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYZAxis}});

      for (int j = 1; j < 4; j++) { // Track counting
        registry.add(Form("h_%s_impact_parameter_xy_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xy_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYAxis}});
        registry.add(Form("h_%s_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_pt_sign_impact_parameter_xy_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} N=%d", charJetFlavor[i].data(), j), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYAxis}});
        registry.add(Form("h_%s_impact_parameter_xyz_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYZAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xyz_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterXYZAxis}});
        registry.add(Form("h_%s_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
        registry.add(Form("h_%s_sign_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
        registry.add(Form("h_%s_pt_sign_impact_parameter_xyz_significance_N%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xyz} N=%d", charJetFlavor[i].data(), j), {HistType::kTH2F, {jetPtRangeAxis, ImpactParameterSignificanceXYZAxis}});
      }
      // track and jet probability
      registry.add(Form("h_%s_pos_track_probability", charJetFlavor[i].data()), Form("%s positive sign track probability", charJetFlavor[i].data()), {HistType::kTH1F, {TrackProbabilityAxis}});
      registry.add(Form("h_%s_neg_track_probability", charJetFlavor[i].data()), Form("%s negative sign track probability", charJetFlavor[i].data()), {HistType::kTH1F, {TrackProbabilityAxis}});
      registry.add(Form("h_%s_pt_pos_track_probability", charJetFlavor[i].data()), Form("%s positive sign track probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, TrackProbabilityAxis}});
      registry.add(Form("h_%s_pt_neg_track_probability", charJetFlavor[i].data()), Form("%s negative sign track probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, TrackProbabilityAxis}});
      registry.add(Form("h_%s_JP", charJetFlavor[i].data()), Form("%s jet probability", charJetFlavor[i].data()), {HistType::kTH1F, {JetProbabilityAxis}});
      registry.add(Form("h_%s_JP_Log", charJetFlavor[i].data()), Form("Log %s jet probability", charJetFlavor[i].data()), {HistType::kTH1F, {JetProbabilityLogAxis}});
      registry.add(Form("h_%s_pt_JP", charJetFlavor[i].data()), Form("%s jet probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, JetProbabilityAxis}});
      registry.add(Form("h_%s_pt_JP_Log", charJetFlavor[i].data()), Form("Log %s jet probability", charJetFlavor[i].data()), {HistType::kTH2F, {jetPtRangeAxis, JetProbabilityLogAxis}});
    }

    registry.add("h_jet_entries", "Entries of jets (lf=1, c=2, b=3)", {HistType::kTH1F, {JetEntries}});
    registry.add("h_jet_track_chi2", "track ch2", {HistType::kTH1F, {{100, -10., 10.}}});

    // for Resolution function

    for (int i = 0; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        registry.add(Form("h_%s_sign_impact_parameter_xy_significance_class%d", charJetFlavor[i].data(), j), Form("%s sign impact parameter significance dca_{xy} (class%d)", charJetFlavor[i].data(), j), {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
      }
    }
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>;
  //using TagJetTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels, aod::TagJetConstBase>>;
  using TagJetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels, aod::TagJetConstBase>;

  template <typename T, typename U, typename V>
  void fillHistogramsTC(T const& collision, U const& jets, V const& mcParticles)
  {
    for (const auto& jet : jets) {
      const int jetflavor = getJetFlavor(jet);
      if (jetflavor == 2) LOGF(info, "b jet");
      //for (const auto& track : jet.template tracks_as<TagJetTracks>()) {
//      for (const auto& track : jet.template tracks_as<JetTracks>()) {
//      //for (const auto& track : jet.template tracks_as<aod::TagJetConstBase>()) {
//        //LOGF(info, Form("track sign vertex: %f", track.impXY()));
//        LOGF(info, Form("track sign vertex: %f", track.pt()));
//      }
      //for (const auto& track: jet.template trackstemp_as<aod:::TagJetConstBase>()) {
      //for (const auto& track: jet.template trackstemp_as<JetTracks>()) {
      for (const auto& track : jet.template trackstemp_as<TagJetTracks>()) {
        LOGF(info, Form("track sign vertex: %f", track.pt()));
      }
//      for (const auto& track : jet.template tracks_as<aod::TagJetConstBase>()) {
//        LOGF(info, "in tracks");
//        //LOGF(info, Form("track sign vertex: %f", track.pt()));
//        LOGF(info, Form("track sign vertex: %d", track.signVertex()));
////        if (!track.has_mcParticle()) {
////          LOGF(warning, "No MC particle for track, skip...");
////          continue;
////        }
////        if (!track.isPrimaryTrack()) {
////          LOGF(warning, "No primary track, skip...");
////          continue;
////        }
////        //LOGF(info, track.ImpXY());
//      }
    }
  }

  void processDummy(o2::aod::TagJetConstBase const& tracks)
  {
    LOGF(info, "PROCESS DUMMY");
  }
  PROCESS_SWITCH(JetHfTaggingTask, processDummy, "Dummy process function turned on by default", false);

  void processTCHfJetsMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::TaggedMCDetectorLevelJets> const& jets,
                         //soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::TaggedMCDetectorLevelJets, aod::TaggedMCDetectorLevelJetConstituents> const& jets,
                         soa::Join<aod::TaggedMCDetectorLevelJets, aod::TaggedMCDetectorLevelJetConstituents> const& tagjets,
                         JetTracks const& tracks,
                         aod::McCollisions const& mcCollisions,
                         aod::McParticles const& mcParticles)
  {
    LOGF(info, "PROCESS TC");
    //fillHistogramsTC(collision, jets, mcParticles);
    fillHistogramsTC(collision, tagjets, mcParticles);
  }
  PROCESS_SWITCH(JetHfTaggingTask, processTCHfJetsMCD, "Task of impact parameter for track counting (MCD)", false);

  //void processTEMP(soa::Join<aod::Tracks, aod::TagJetConstBase> const& tracks)
  void processTEMP(aod::TagJetConstBase const& tracks)
  {
    for (auto track: tracks) {
      LOGF(info, Form("track sign vertex: %f", track.impXY()));
    }
  }
  PROCESS_SWITCH(JetHfTaggingTask, processTEMP, "Task of impact parameter for track counting (MCD)", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetHfTaggingTask>(cfgc, TaskName{"jet-hf-tagging-task"})}; }


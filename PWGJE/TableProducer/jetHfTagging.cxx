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

// jet hf tagging
//
// Authors: Hanseo Park

#include <string>
#include "DCAFitter/DCAFitterN.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "TF1.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/TableProducer/jetHfTagging.h"

#include "EventFiltering/filterTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetHfTagging {

  //HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Produces<o2::aod::TaggedMCDetectorLevelJets> taggedJets;
  Produces<o2::aod::TaggedMCDetectorLevelJetConstituents> taggedJetConstituents;
  //Produces<o2::aod::TaggedMCDetectorLevelJetConstituentsBasic> taggedJetConstituentsBasic;
  //Produces<o2::aod::TaggedMCParticleLevelJets> taggedPartJets;
  Produces<o2::aod::TagJetConstBase> taggedJetConstituent;

  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>;

//  OutputObj<TH1F> hLfJetPtWoHfShower{"h_lfjet_pt_wo_hf_shower"};
//  OutputObj<TH1F> hCJetPtWoHfShower{"h_bjet_pt_wo_hf_shower"};
//  OutputObj<TH1F> hBJetPtWoHfShower{"h_cjet_pt_wo_hf_shower"};
//  OutputObj<TH1F> hLfJetPtWithHfShower{"h_lfjet_pt_with_hf_shower"};
//  OutputObj<TH1F> hCJetPtWithHfShower{"h_bjet_pt_with_hf_shower"};
//  OutputObj<TH1F> hBJetPtWithHfShower{"h_cjet_pt_with_hf_shower"};

  template <typename T, typename U, typename V, typename P>
  void fillHistogramMCDQA(T const& collision, U const& jets, V const& tracks, P const& mcParticles)
  {
    for (auto part : mcParticles) {
      auto partIndex = JetTaggingUtilities::getOriginalMotherIndex<P>(part);
      auto hfPartIndex = JetTaggingUtilities::getOriginalHFMotherIndex<P>(part);
      if (partIndex == hfPartIndex) LOGF(info, Form("heavy flavor particle index: %d", hfPartIndex));
    }
    int jetFlavor = 0;
    int jetFlavorWJetAxis = 0;
    typename V::iterator hftrack;
    for (auto& jet : jets) {
      jetFlavor = JetTaggingUtilities::jetTrackFromHFShower(jet, mcParticles, tracks, hftrack);
      jetFlavorWJetAxis = JetTaggingUtilities::mcdJetFromHFShower(jet, mcParticles, tracks);
      LOGF(info, Form("jetWoHfShower color: %d", jetFlavor));
      LOGF(info, Form("jet color: %d", jetFlavorWJetAxis));
    }
  }

  template <typename T, typename U, typename V, typename P>
  void fillHistogramMCPQA(T const& collision, U const& jets, V const& tracks, P const& mcParticles)
  {
//    for (auto part : mcParticles) {
//      auto partIndex = JetTaggingUtilities::getOriginalMotherIndex<P>(part);
//      auto hfPartIndex = JetTaggingUtilities::getOriginalHFMotherIndex<P>(part);
//      if (partIndex == hfPartIndex) LOGF(info, Form("heavy flavor particle index: %d", hfPartIndex));
//    }
    int jetFlavor = 0;
    //int jetFlavorWJetAxis = 0;
    typename P::iterator hfpart;
    for (auto& jet : jets) {
      jetFlavor = JetTaggingUtilities::jetParticleFromHFShower(jet, mcParticles, hfpart);
      //jetFlavorWJetAxis = JetTaggingUtilities::mcpJetFromHFShower(jet, mcParticles);
      if (jetFlavor>0) LOGF(info, Form("jetWoHfShower color: %d", jetFlavor));
      //LOGF(info, Form("jet color: %d", jetFlavorWJetAxis));
    }

  }

  template <typename T, typename U, typename V, typename P>
  void createTaggingMCD(T const& collision, U const& jets, V const& tracks, P const& mcParticles)
  {
    typename V::iterator hftrack;
    for (const auto& jet : jets) {
      int jetflavor = 0;
      jetflavor = JetTaggingUtilities::jetTrackFromHFShower(jet, mcParticles, tracks, hftrack);
      int islfjet = 0;
      int iscjet = 0;
      int isbjet = 0;
      if (jetflavor == 0) islfjet=1;
      if (jetflavor == 1) iscjet =1;
      if (jetflavor == 2) isbjet =1;
      //taggedJets(collision.globalIndex(), islfjet, iscjet, isbjet);
      taggedJets(islfjet, iscjet, isbjet);

      std::vector<int> tracktemp;
      for (auto& track : jet.template tracks_as<JetTracks>()) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }

        if (!track.isPrimaryTrack()) {
          LOGF(warning, "No primary track, skip...");
          continue;
        }
        double sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        SetSignImpactParameterSignificance(collision, jet, track, sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig);
        taggedJetConstituent(track.globalIndex(), sign, varImpXY, varImpXYSig);
        tracktemp.push_back(track.globalIndex());
        //taggedJetConstituentsBasic(taggedJets.lastIndex(), sign, varImpXY, varImpXYSig);
      }
      taggedJetConstituents(jet.globalIndex(), tracktemp);
    }
  }

  template <typename T, typename U, typename V, typename P>
  void createTaggingMCP(T const& collision, U const& jets, V const& tracks, P const& mcParticles)
  {
    typename V::iterator hftrack;
    for (const auto& jet : jets) {
      int jetflavor = 0;
      jetflavor = JetTaggingUtilities::jetTrackFromHFShower(jet, mcParticles, tracks, hftrack);
      int islfjet = 0;
      int iscjet = 0;
      int isbjet = 0;
      if (jetflavor == 0) islfjet=1;
      if (jetflavor == 1) iscjet =1;
      if (jetflavor == 2) isbjet =1;
      //taggedPartJets(collision.globalIndex(), islfjet, iscjet, isbjet);

      for (auto& track : jet.template tracks_as<JetTracks>()) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }

        if (!track.isPrimaryTrack()) {
          LOGF(warning, "No primary track, skip...");
          continue;
        }
        double sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        SetSignImpactParameterSignificance(collision, jet, track, sign, varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig);
        taggedJetConstituent(track.globalIndex(), sign, varImpXY, varImpXYSig);
      }
    }
  }

  void processDummy(aod::Tracks const& track)
  {
    LOGF(info, "PROCESS DUMMY");
  }
  PROCESS_SWITCH(JetHfTagging, processDummy, "Process dummy", false);

  void processTaggingMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         JetTracks const& tracks,
                         aod::McParticles const& mcParticles)
  {
    createTaggingMCD(collision, jets, tracks, mcParticles);
  }
  PROCESS_SWITCH(JetHfTagging, processTaggingMCD, "Task of impact parameter for heavy flavor tagging flavor (MCD)", false);

  void processTaggingMCDQA(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         JetTracks const& tracks,
                         aod::McParticles const& mcParticles)
  {
    fillHistogramMCDQA(collision, jets, tracks, mcParticles);
  }
  PROCESS_SWITCH(JetHfTagging, processTaggingMCDQA, "Task of impact parameter for heavy flavor tagging flavor QA (MCD)", false);

  void processTaggingMCP(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                         JetTracks const& tracks,
                         aod::McParticles const& mcParticles)
  {
    //createTaggingMCP(mcCollision, jets, tracks, mcParticles);
    createTaggingMCP(collision, jets, tracks, mcParticles);
  }
  PROCESS_SWITCH(JetHfTagging, processTaggingMCP, "Task of impact parameter for heavy flavor tagging flavor (MCP)", false);

  void processTaggingMCPQA(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                         JetTracks const& tracks,
                         aod::McParticles const& mcParticles)
  {
    //fillHistogramMCPQA(mcCollision, jets, tracks, mcParticles);
    fillHistogramMCPQA(collision, jets, tracks, mcParticles);
  }
  PROCESS_SWITCH(JetHfTagging, processTaggingMCPQA, "Task of impact parameter for heavy flavor tagging flavor QA (MCP)", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetHfTagging>(cfgc, TaskName{"jet-hf-tagging"})}; }

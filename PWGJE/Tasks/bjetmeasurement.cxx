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

// \file jetfragmentationhf.cxx
// \note Extended from jetfinderhfQA.cxx and jetfragmentation.cxx

// HF jets fragmentation function task
//
// Authors: hanseo.park@cern.ch

#include "DCAFitter/DCAFitterN.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include <string>
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

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGHF/HFJ/Utils/HFTaggingUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

#include "Framework/runDataProcessing.h"

struct bJetMeasurementTask {

  TrackSelection customTrackCuts;
  //static constexpr std::string_view charSelection[8] = {"1", "2", "3", "4", "5", "6", "7", "8"};
  static constexpr std::string_view charJetFlavor[] = {"undefined", "lfjet", "cjet", "bjet"};
  std::vector<int> pdgVector = {1,2,3,4,5,21};
  std::vector<std::string> partonVector = {"d", "u", "s", "c", "b", "g"};
  std::array<double, 3> FitParlfjet;
  std::array<double, 3> FitParbjet;
  std::array<double, 3> FitParcjet;
  std::vector<double> N1Tracks;
  std::vector<double> N2Tracks;
  std::vector<double> N3Tracks;

  Configurable<double> yCandMax{"yCandMax", -1, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<bool> doFitResoFunc{"doFitResoFunc", false, "Do Fit resolution function"};

  // Binning
  ConfigurableAxis binPdgCode{"binPdgCode", {5000, -2500.f, 2500.f}, ""};
  ConfigurableAxis binStatus{"binStatusCode", {200, -99.5f, 100.5f}, ""};
  ConfigurableAxis binJetPt{"binJetPt", {200, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, -0.5, 99.5}, ""};
  ConfigurableAxis binZ{"binZ", {100, -5e-3f, 1.f + 5e-3f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binD0CandidatePt{"binD0CandidatePt", {360, 0.f, 36.f}, ""};
  ConfigurableAxis binD0CandidateMass{"binD0CandidatMass", {500, 0.f, 5.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {1000, -0.4f, 0.4f}, ""};
  ConfigurableAxis binImpactParameterSignificanceXY{"binImpactParameterSignificanceXY", {1000, -40.f, 40.f}, ""};
  ConfigurableAxis binImpactParameterXYZ{"binImpactParameterXYZ", {1000, -0.4f, 0.4f}, ""};
  ConfigurableAxis binImpactParameterSignificanceXYZ{"binImpactParameterSignificanceXYZ", {1000, -40.f, 40.f}, ""};
  ConfigurableAxis binJetProbability{"binJetProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbabilityLog{"binJetProbabilityLog", {100, 0.f, 10.f}, ""};
  ConfigurableAxis binJetEntries{"binEntries", {3, 0.f, 3.f}, ""};

  // Axis
  AxisSpec pdgCodeAxis = {binPdgCode, "PDG code"};
  AxisSpec statusAxis = {binStatus, "Status code"};
  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
  AxisSpec etaAxis = {binEta, "#eta"};
  AxisSpec phiAxis = {binPhi, "#phi"};
  AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
  AxisSpec zAxis = {binZ, "#it{z}"};
  AxisSpec jetRAxis = {binJetR, "#it{R}"};
  AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
  AxisSpec d0candidatePtAxis = {binD0CandidatePt, "#it{p}_{T}^{D0 candidate}"};
  AxisSpec d0candidateMassAxis = {binD0CandidateMass, "inv. mass(#pi K)"};
  AxisSpec ImpactParameterXYAxis = {binImpactParameterXY, "IP_{XY} [cm]"};
  AxisSpec ImpactParameterSignificanceXYAxis = {binImpactParameterSignificanceXY, "IPs_{XY} [cm]"};
  AxisSpec ImpactParameterXYZAxis = {binImpactParameterXYZ, "IP_{XYZ} [cm]"};
  AxisSpec ImpactParameterSignificanceXYZAxis = {binImpactParameterSignificanceXYZ, "IPs_{XYZ} [cm]"};
  AxisSpec JetProbabilityAxis = {binJetProbability, "JP"};
  AxisSpec JetProbabilityLogAxis = {binJetProbabilityLog, "-Log(JP)"};
  AxisSpec JetEntries = {binJetEntries, "lf=1, c=2, b=3"};
  //const AxisSpec axisSelections{10, 0.5, 10.5f, "Selection"};

  TF1 *fResoFunccjet = new TF1("fResoFunccjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1 *fResoFuncbjet = new TF1("fResoFuncbjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  TF1 *fResoFunclfjet = new TF1("fResoFunclfjet", "expo(0)+expo(2)+expo(4)+gaus(6)", -40, 0);
  HistogramRegistry registry{"registry",
                             {{"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}}, true},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_lfjet_ntracks", "lf-jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_cjet_ntracks", "c-jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_bjet_ntracks", "b-jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_lfjet_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_cjet_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_bjet_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_undefined_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_track_phi", "track #phi;#phi_{track};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h2_jet_pt_track_pt", ";#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}}, true},
                              {"h2_jet_pt_track_eta", ";#it{p}_{T,jet} (GeV/#it{c}); #eta_{track}", {HistType::kTH2F, {jetPtAxis, etaAxis}}, true},
                              {"h2_jet_pt_track_phi", ";#it{p}_{T,jet} (GeV/#it{c}); #phi_{track}", {HistType::kTH2F, {jetPtAxis, phiAxis}}, true},
                              {"h_part_jet_pt", "Particle Level jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}}, true},
                              {"h_part_jet_eta", "Particle Level jet #eta;#eta_{jet};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_part_jet_phi", "Particle Level jet #phi;#phi_{jet};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h_part_jet_ntracks", "Particle Level jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_part_track_pt", "Particle Level track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_part_track_eta", "Particle Level track #eta;#eta_{track};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_part_track_phi", "Particle Level track #phi;#phi_{track};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h2_part_jet_pt_part_track_pt", "Particle Level;#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}}, true},
                              {"h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {{100, 0.0, 10.0}}}, true}}};

  void init(o2::framework::InitContext&)
  {


    auto vbins = (std::vector<double>)binsPt;
    // Data and MC (Detector level, Rec)
    registry.add("h_d0candidate_pt", "D0 candidates #it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_d0candidate_eta", "D0 candidates #eta;#eta_{candidate};entries", {HistType::kTH1F, {etaAxis}});
    registry.add("h_d0candidate_phi", "D0 candidates #phi;#phi_{candidate};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("h2_jet_pt_d0candidate_pt", "#it{p}_{T, jet} vs #it{t}_{T, D0 candidates};#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, d0candidatePtAxis}});
    registry.add("h2_jet_pt_d0candidate_eta", "#it{p}_{T,jet} vs #eta_{D0 candidates};#it{p}_{T,jet} (GeV/#it{c}); #eta_{candidate}", {HistType::kTH2F, {jetPtAxis, etaAxis}});
    registry.add("h2_jet_pt_d0candidate_phi", "#it{p}_{T,jet} vs #phi_{D0 candidates};#it{p}_{T,jet} (GeV/#it{c}); #phi_{candidate}", {HistType::kTH2F, {jetPtAxis, phiAxis}});
    registry.add("h_d0candidate_mass", "Inv. mass distribution of D0 candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {d0candidateMassAxis}});
    registry.add("h2_d0candidate_mass_d0candidate_pt", "Inv. mass of D0 candidates vs #it{p}_{T, D0 candidates};inv. mass (#pi K) (GeV/#it{c}^{2});entries ", {HistType::kTH2F, {d0candidateMassAxis, {vbins, "#it{p}_{T,candidate} (GeV/#it{c})"}}});
    registry.add("h2_d0candidate_mass_d0candidate_phi", "Inv. mass of D0 candidates vs #phi_{D0 candidates};inv. mass (#pi K) (GeV/#it{c}^{2});entries ", {HistType::kTH2F, {d0candidateMassAxis, phiAxis}});
    registry.add("h_jet_chargefrag", "#it{p}_{T,track}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});
    registry.add("h_jet_d0frag", "#it{p}_{T,D0 candidates}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});
    // prong
    registry.add("h_d0candidate_prong0_pt", "D0 candidates prong0#it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_d0candidate_prong1_pt", "D0 candidates prong1#it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_decay_length", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} GeV/#it{c}"}}});
    registry.add("h_decay_length_xy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} GeV/#it{c}"}}});
    registry.add("h_decay_length_error", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_decay_length_xy_error", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong01D", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_impact_parameter_d0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong11D", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_impact_parameter_d0Prong0_error", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong1_error", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_jet_d0Prong0frag", "#it{p}_{T,D0 candidates Prong0}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});
    registry.add("h_jet_d0Prong1frag", "#it{p}_{T,D0 candidates Prong1}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});

    // MC (Partcle level, Gen)
    registry.add("h_part_d0candidate_pt", "Particle Level D0 candidates #it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c}) ();entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_part_d0candidate_eta", "Particle Level D0 candidates #eta;#eta_{candidate};entries", {HistType::kTH1F, {etaAxis}});
    registry.add("h_part_d0candidate_phi", "Particle Level D0 candidates #phi;#phi_{candidate};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("h2_part_jet_pt_part_d0candidate_pt", "#it{p}_{T, jet}^{Part} vs #it{t}_{T, D0 candidates}^{Part};#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, d0candidatePtAxis}});
    registry.add("h2_part_jet_pt_part_d0candidate_eta", "#it{p}_{T,jet}^{Part} vs #eta_{D0 candidates}^{Part};#it{p}_{T,jet} (GeV/#it{c}); #eta_{candidate}", {HistType::kTH2F, {jetPtAxis, etaAxis}});
    registry.add("h2_part_jet_pt_part_d0candidate_phi", "#it{p}_{T,jet}^{Part} vs #phi_{D0 candidates}^{Part};#it{p}_{T,jet} (GeV/#it{c}); #phi_{candidate}", {HistType::kTH2F, {jetPtAxis, phiAxis}});
    registry.add("h_part_jet_part_chargefrag", "#it{p}_{T,track}^{Part}/#it{p}_{T,jet}^{Part};entries", {HistType::kTH1F, {zAxis}});
    registry.add("h_part_jet_part_d0frag", "#it{p}_{T,D0 candidates}^{Part}/#it{p}_{T,jet}^{Part};entries", {HistType::kTH1F, {zAxis}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);D^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);D^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    // Matched
    registry.add("h2_jet_pt_part_jet_pt", "#it{p}_{T,jet}^{part} vs #it{p}_{T,jet}^{Dec} ;#it{p}_{T,jet}^{DeC} (GeV/#it{c}); #it{p}_{T,jet}^{part}", {HistType::kTH2F, {binJetPt, binJetPt}});

    // Comp. hf jets and inclusive jets

    // Tagging method
    registry.add("h_mc_pdgcode", "pdgcode", {HistType::kTH1F, {{1000, 0., 1000}}});
    //registry.add("h_jet_mother_track_pdgcode", "pdg code", {HistType::kTH1F, {{30, -0.5, 29.5}}});
    registry.add("h_no_mother_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_no_mother_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_no_mother_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_jet_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_jet_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_jet_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_jet_mother_1_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_jet_mother_1_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_jet_mother_1_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_jet_mother_2_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_jet_mother_2_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_jet_mother_2_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_jet_mother_3_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_jet_mother_3_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_jet_mother_3_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_jet_mother_4_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_jet_mother_4_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_jet_mother_4_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_jet_mother_5_particle_pdgcode", "pdg code", {HistType::kTH1F, {{5000, -2500, 2500}}});
    registry.add("h_jet_mother_5_particle_statuscode", "status code", {HistType::kTH1F, {{200, -100.5, 99.5}}});
    registry.add("h_jet_mother_5_particle_pdgcode_statuscode", "status code", {HistType::kTH2F, {{5000, -2500, 2500}, {200, -100.5, 99.5}}});
    registry.add("h_bjet_particle_pdgcode", "particles' pdg code in b jet", {HistType::kTH1F, {pdgCodeAxis}});
    registry.add("h_cjet_particle_pdgcode", "particles' pdg code in c jet", {HistType::kTH1F, {pdgCodeAxis}});
    registry.add("h_lfjet_particle_pdgcode", "particles' pdg code in lf jet", {HistType::kTH1F, {pdgCodeAxis}});
    registry.add("h_undefined_particle_pdgcode", "undefined particles' pdg code", {HistType::kTH1F, {pdgCodeAxis}});
    registry.add("h_bjet_particle_statuscode", "particles' status code in b jet", {HistType::kTH1F, {statusAxis}});
    registry.add("h_cjet_particle_statuscode", "particles' status code in c jet", {HistType::kTH1F, {statusAxis}});
    registry.add("h_lfjet_particle_statuscode", "particles' status code in lf jet", {HistType::kTH1F, {statusAxis}});
    registry.add("h_undefined_particle_statuscode", "undefined particles' status code", {HistType::kTH1F, {statusAxis}});
    //registry.add("h_bjet_particle_pdgcode", "particles' pdg code, status code in b jet", {HistType::kTH1F, {pdgCodeAxis}});

    registry.add("h_inc_jet_impact_parameter_xy", "inclusive jet impact parameter dca_{xy}", {HistType::kTH1F, {{1000, -0.1, 0.1}}});
    registry.add("h_inc_jet_impact_parameter_xy_significance", "inclusive jet impact parameter significance dca_{xy}", {HistType::kTH1F, {{1000, -0.1, 0.1}}});
    registry.add("h_inc_jet_impact_parameter_xyz", "inclusive jet impact parameter dca_{xyz}", {HistType::kTH1F, {{1000, -0.1, 0.1}}});
    registry.add("h_inc_jet_impact_parameter_xyz_significance", "inclusive jet impact parameter significance dca_{xy}", {HistType::kTH1F, {{1000, -0.1, 0.1}}});
    registry.add("h_bjet_impact_parameter_xy", "bjet impact parameter dca_{xy}", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_cjet_impact_parameter_xy", "cjet impact parameter dca_{xy}", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy", "lfjet impact paramete dca_{xy}", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance", "bjet impact parameter significance dca_{xy}", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance", "cjet impact parameter significance dca_{xy}", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_significance", "lfjet impact parameter significance dca_{xy}", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_N1", "bjet impact parameter dca_{xy} N=1", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_N2", "bjet impact parameter dca_{xy} N=2", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_N3", "bjet impact parameter dca_{xy} N=3", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_N1", "cjet impact parameter dca_{xy} N=1", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_N2", "cjet impact parameter dca_{xy} N=2", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_N3", "cjet impact parameter dca_{xy} N=3", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_N1", "lfjet impact parameter dca_{xy} N=1", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_N2", "lfjet impact parameter dca_{xy} N=2", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_N3", "lfjet impact parameter dca_{xy} N=3", {HistType::kTH1F, {ImpactParameterXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_N1", "bjet impact parameter dca_{xy} N=1", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_N2", "bjet impact parameter dca_{xy} N=2", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_N3", "bjet impact parameter dca_{xy} N=3", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_N1", "cjet impact parameter dca_{xy} N=1", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_N2", "cjet impact parameter dca_{xy} N=2", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_N3", "cjet impact parameter dca_{xy} N=3", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_significance_N1", "lfjet impact parameter dca_{xy} N=1", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_significance_N2", "lfjet impact parameter dca_{xy} N=2", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_significance_N3", "lfjet impact parameter dca_{xy} N=3", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

    registry.add("h_bjet_impact_parameter_xyz", "bjet impact parameter dca_{xyz}", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz", "cjet impact parameter dca_{xyz}", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz", "lfjet impact paramete dca_{xyz}", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_significance", "bjet impact parameter significance dca_{xyz}", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_significance", "cjet impact parameter significance dca_{xyz}", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_significance", "lfjet impact parameter significance dca_{xyz}", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_N1", "bjet impact parameter dca_{xyz} N=1", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_N2", "bjet impact parameter dca_{xyz} N=2", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_N3", "bjet impact parameter dca_{xyz} N=3", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_N1", "cjet impact parameter dca_{xyz} N=1", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_N2", "cjet impact parameter dca_{xyz} N=2", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_N3", "cjet impact parameter dca_{xyz} N=3", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_N1", "lfjet impact parameter dca_{xyz} N=1", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_N2", "lfjet impact parameter dca_{xyz} N=2", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_N3", "lfjet impact parameter dca_{xyz} N=3", {HistType::kTH1F, {ImpactParameterXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_significance_N1", "bjet impact parameter dca_{xyz} N=1", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_significance_N2", "bjet impact parameter dca_{xyz} N=2", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_bjet_impact_parameter_xyz_significance_N3", "bjet impact parameter dca_{xyz} N=3", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_significance_N1", "cjet impact parameter dca_{xyz} N=1", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_significance_N2", "cjet impact parameter dca_{xyz} N=2", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_cjet_impact_parameter_xyz_significance_N3", "cjet impact parameter dca_{xyz} N=3", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_significance_N1", "lfjet impact parameter dca_{xyz} N=1", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_significance_N2", "lfjet impact parameter dca_{xyz} N=2", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_lfjet_impact_parameter_xyz_significance_N3", "lfjet impact parameter dca_{xyz} N=3", {HistType::kTH1F, {ImpactParameterSignificanceXYZAxis}});
    registry.add("h_bjet_decay_length", "bjet decay length", {HistType::kTH1F, {{1000, -0.1, 0.1}}});
    registry.add("h_cjet_decay_length", "cjet dacay length", {HistType::kTH1F, {{1000, -0.1, 0.1}}});
    registry.add("h_lfjet_decay_length", "lfjet decay length", {HistType::kTH1F, {{1000, -0.1, 0.1}}});

    registry.add("h_bjet_JP", "bjet jet probability", {HistType::kTH1F, {JetProbabilityAxis}});
    registry.add("h_cjet_JP", "cjet jet probability", {HistType::kTH1F, {JetProbabilityAxis}});
    registry.add("h_lfjet_JP", "lfjet jet probability", {HistType::kTH1F, {JetProbabilityAxis}});
    registry.add("h_bjet_JP_Log", "Log bjet jet probability", {HistType::kTH1F, {JetProbabilityLogAxis}});
    registry.add("h_cjet_JP_Log", "Log cjet jet probability", {HistType::kTH1F, {JetProbabilityLogAxis}});
    registry.add("h_lfjet_JP_Log", "Log lfjet jet probability", {HistType::kTH1F, {JetProbabilityLogAxis}});
    registry.add("h_jet_entries", "Entries of jets (lf=1, c=2, b=3)", {HistType::kTH1F, {JetEntries}});

    registry.add("h_jet_track_chi2", "track ch2",{HistType::kTH1F, {{100, -10., 10.}}});
    registry.add("h_track_dcaXY", "dcaXY",{HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_track_dcaZ", "dcaZ",{HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_particle_pt_temp", "particle",{HistType::kTH1F, {{10000, 0., 100.}}});
    registry.add("h_track_pt_temp", "track",{HistType::kTH1F, {{10000, 0., 100.}}});
    registry.add("h_track_pt_jet_temp", "track in jet",{HistType::kTH1F, {{10000, 0., 100.}}});
    registry.add("h_hftrack_pt_jet_temp", "hftrack in jet",{HistType::kTH1F, {{10000, 0., 100.}}});
    registry.add("h_constituent_pt", "constituent",{HistType::kTH1F, {{10000, 0., 100.}}});

    // for Resolution function
    registry.add("h_lfjet_impact_parameter_xy_significance_class1", "lfjet impact parameter significance dca_{xy} (class1)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class1", "cjet impact parameter significance dca_{xy} (class1)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class1", "bjet impact parameter significance dca_{xy} (class1)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_lfjet_impact_parameter_xy_significance_class2", "lfjet impact parameter significance dca_{xy} (class2)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class2", "cjet impact parameter significance dca_{xy} (class2)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class2", "bjet impact parameter significance dca_{xy} (class2)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

    registry.add("h_lfjet_impact_parameter_xy_significance_class3", "lfjet impact parameter significance dca_{xy} (class3)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class3", "cjet impact parameter significance dca_{xy} (class3)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class3", "bjet impact parameter significance dca_{xy} (class3)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

    registry.add("h_lfjet_impact_parameter_xy_significance_class4", "lfjet impact parameter significance dca_{xy} (class4)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class4", "cjet impact parameter significance dca_{xy} (class4)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class4", "bjet impact parameter significance dca_{xy} (class4)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

    registry.add("h_lfjet_impact_parameter_xy_significance_class5", "lfjet impact parameter significance dca_{xy} (class5)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class5", "cjet impact parameter significance dca_{xy} (class5)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class5", "bjet impact parameter significance dca_{xy} (class5)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

    registry.add("h_lfjet_impact_parameter_xy_significance_class6", "lfjet impact parameter significance dca_{xy} (class6)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class6", "cjet impact parameter significance dca_{xy} (class6)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class6", "bjet impact parameter significance dca_{xy} (class6)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

    registry.add("h_lfjet_impact_parameter_xy_significance_class7", "lfjet impact parameter significance dca_{xy} (class7)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_cjet_impact_parameter_xy_significance_class7", "cjet impact parameter significance dca_{xy} (class7)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});
    registry.add("h_bjet_impact_parameter_xy_significance_class7", "bjet impact parameter significance dca_{xy} (class7)", {HistType::kTH1F, {ImpactParameterSignificanceXYAxis}});

  }

  //using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::McTrackLabels>;
  using CandidateD0Data = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  //using CandidateD0MC = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec, aod::HfCand2ProngMcGen>;
  using CandidateD0MC = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>;
  using JetParticles2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

//  template <typename T>
//  int GetJetFlavor(T const& particle){
//	int pdgcode=TMath::Abs(particle.pdgCode());
//    if (pdgcode==1 || pdgcode==2 || pdgcode==3 || pdgcode==21) return 1; // light flavor
//	else if (pdgcode==4) return 2; // charm
//    else if (pdgcode==5) return 3; // beauty
//    else return 0; // undefined
//  }

//  template <typename T, typename U>
//  void SearchingMothers(T const& particle, U &mothersid){
//	if (particle.has_mothers()){
//	  //auto pdg = TMath::Abs(particle.pdgCode());
//	  auto status = TMath::Abs(particle.getGenStatusCode());
//	  if (status==23) mothersid.push_back(particle.globalIndex());
//	  for (auto& mother : particle.template mothers_as<aod::McParticles>()) {
//	    SearchingMothers(mother, mothersid);
//      }
//	}
//	//if ((mothersid.size()%1000)==0) LOGF(info, "mothers id size: %d", mothersid.size());
//  } 

//  template <typename T>
//  bool CheckGluonSpliting(T const& track){
//	if (!track.has_mcParticle()) return false;
//	auto particle=track.mcParticle();
//
//	if (!particle.has_mothers()) return false;
//
//	if (!particle.isPhysicalPrimary()) return false;
//
//	//for (auto& mothers : particle.template mothers_as<aod::McParticles>()){
//	  auto mothers = particle.template mothers_first_as<aod::McParticles>();
//	  int motherid = -999;
//	  if (mothers.size()==1){
//	    motherid = mothers.mothersIds()[0];
//	  } else {
//		motherid = mothers.mothersIds()[0];
//	  }
//	  while(motherid > -1) {
//		auto mp = mothers.iteratorAt(motherid);
//		auto pdg = TMath::Abs(mp.pdgCode());
//		if (pdg == 4 || pdg == 5) {
//		  if (mp.has_mothers()){
//			LOGF(info, "mother's id: %d, pdg: %d", motherid, pdg);
//			motherid = mp.mothersIds()[0];
//		  } else {
//			motherid = -999;
//		  }
//		}
//      }
//
//    //}
//	return true;
//	
//  }



//  template <typename T>
//  int GetJetFlavor(T const& jet){
//	const int arraySize=99;
//	int candCharmLfCode[arraySize], count=0;
//	for (auto& track : jet.template tracks_as<JetTracks>()) {
//	  if (!track.has_mcParticle()){
//		  LOGF(warning, "No MC particle for track, skip...");
//		  continue;
//	  } 
//	  auto particle=track.mcParticle();
//	  if (!particle.has_mothers()){
//		LOGF(warning, "No mother particle for particle, skip...");
//		continue;
//      }
//	  for (auto& mother : particle.template mothers_as<aod::McParticles>()){
//		//auto m=mother.template mothers_first_as<aod::McParticles>();
//		int pdgcode=TMath::Abs(mother.pdgCode());
//		if (pdgcode==21 || (pdgcode>0 && pdgcode<6)){
//		  if(pdgcode==5) {
//		    return 3; // beauty
//		  }else {
//			candCharmLfCode[count]=pdgcode;
//			count++;
//		  }
//		}
//	  }
//
//	}
//	if (count>0){
//	  for (int i=0; i<count; i++){
//		if(candCharmLfCode[i]==4) return 2; // charm
//	  }
//	} else {
//	  return 0;
//	}
//    if (count>0) return 1; // light flavor
//	return 0;
//  }

  template <typename T>
  bool trackSelectionJP(T const& jet){
	if (jet.tracks().size()<2) return 0;

	return 1;
  }
  
  void SetFitResoFunc(){
	if (doFitResoFunc){
	  fResoFunclfjet->SetParameters(FitParlfjet[0], FitParlfjet[1], FitParlfjet[2]); // refer to IPs distribution
	  fResoFunccjet->SetParameters(FitParcjet[0], FitParcjet[1], FitParcjet[2]); // refer to IPs distribution
	  fResoFuncbjet->SetParameters(FitParbjet[0], FitParbjet[1], FitParbjet[2]); // refer to IPs distribution
    } else {
	  fResoFunclfjet->SetParameters(2.57403e+00,1.20758e-01, 1.44472e+00, 1.02265e+00, 7.57231e+00, 1.03262e+00, 6.23147e+03, 1.26972e-02, 9.06710e-01); // refer to IPs distribution
	  fResoFunccjet->SetParameters(1.39055e+00, 6.24134e-02, 1.20702e+00, 1.05035e+00, 5.65658e+00, 8.66233e-01, 1.33617e+03, 7.10938e-02, 0.32481e-01); // refer to IPs distribution
	  fResoFuncbjet->SetParameters(-8.06871e-03,4.05482e-02, 4.89640e-01, 4.05976e-02, 5.51132e+00, 1.24328e+00, 2.12108e+02,-3.72092e-01, 7.18397e-01); // refer to IPs distribution
	  }
  }

  //template <typename T, typename U, typename V = int, typename D = double>
  //bool CalculateJP(T const& jet, U const& mcParticles, V &jetflavor ,D &JP){
  template <typename T, typename U>
  double CalculateJP(T const& jet, U const& mcParticles){
    double JP = -1.;
    if (!trackSelectionJP(jet)) return -1.;
    SetFitResoFunc();
    std::vector<double> jetTracksPt;

    double trackjetProb = 1.;
    const int jetflavor = getJetFlavorAll(jet);
    if (jetflavor==0) return -1;
    for (auto& track : jet.template tracks_as<JetTracks>()) { 
      if (!track.has_mcParticle()){
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      } 
      auto particle=track.mcParticle();
      if (!particle.has_mothers()){
        LOGF(warning, "No mother particle for particle, skip...");
        continue;
      }
      if (track.sign()<0) continue; // only take positive track

      double ProbTrack=0.;

      double dcaXYSig = GetDCAXYSig(track); //track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2());
      if (TMath::Abs(dcaXYSig)>100) dcaXYSig = 99.9; // Limit to function definition range
      if (jetflavor == 1){ // lf
        ProbTrack=fResoFunclfjet->Integral(-100, -1*TMath::Abs(dcaXYSig)) / fResoFunclfjet->Integral(-100,0);
      }	
      if (jetflavor == 2) { // charm
        ProbTrack=fResoFunccjet->Integral(-100, -1*TMath::Abs(dcaXYSig)) / fResoFunccjet->Integral(-100,0);
      }
      if (jetflavor == 3) { // beauty
        ProbTrack=fResoFuncbjet->Integral(-100, -1*TMath::Abs(dcaXYSig)) / fResoFuncbjet->Integral(-100,0);
      }
      //LOGF(info, "track ");
      trackjetProb = trackjetProb*TMath::Abs(ProbTrack);
      jetTracksPt.push_back(track.pt());
    }
    double sumjetProb = 0.;

    if(jetTracksPt.size()<2) return -1;

    for (int i=0; i<jetTracksPt.size(); i++){ // JP
      sumjetProb = sumjetProb+(TMath::Power(-1*TMath::Log(trackjetProb), i) / TMath::Factorial(i));
    }

    JP = trackjetProb*sumjetProb;

    return JP;
  }

  template <typename T>
  double GetDCAXY(T const& track){
    return track.dcaXY() * track.sign();
  }
  template <typename T>
  double GetDCAXYSig(T const& track){
	  return track.dcaXY() * track.sign() / TMath::Sqrt(track.sigmaDcaXY2());
  }

  template <typename T>
  double GetDCAXYZ(T const& track){
	  //return TMath::Sqrt(track.dcaXY()*track.dcaXY()+track.dcaZ()*track.dcaZ()); // track.getSigmaZY()
    auto temp = getTrackParCov(track);
	  return temp.getSigmaZY();
  }

  template <typename T>
  double GetDCAXYZSig(T const& track){
	  return GetDCAXYZ(track) / TMath::Sqrt(TMath::Sqrt(track.sigmaDcaXY2()*track.sigmaDcaXY2()+track.sigmaDcaZ2()*track.sigmaDcaZ2())); // TODO: cal cov
  }

  template <typename T>
  int TrackCounting(T const& jet){
	if (!jet.tracks.size()>2) return 0;
	for (auto& track : jet.template tracks_as<JetTracks>()) {
	  return 0;
	}
	return 0;
  }

  template <typename T>
  int getQualityTrack(T const& track) {
      
    if (!track.hasITS()) return 0;
    //LOGF(info, "its Chi2NCl: %f, pt: %f", track.itsChi2NCl(), track.pt());
    if (track.itsChi2NCl()>2) return 1;
      
    if (track.itsChi2NCl()<2 && track.pt()<2) return 2;
    if (track.itsChi2NCl()<2 && track.pt()>2) return 3;

    return 0;
  }

//  template <typename T>
//  bool GetXYZAt(T const* track, double r, double b, double *xyz)
//  {
//    double dx=r-track.getX();
//    //if(TMath::Abs(dx)<=0.000004) return 1; // return track.getX(), track.getY(), track.getZ()
//    
//    vector<double> fP;
//    fP[0] = track.y();
//    fP[1] = track.z();
//    fP[2] = TMath::Sin(track.phi());
//    fP[3] = track.pz()/track.pt();
//    fP[4] = TMath::Sign(1/track.pt(),track.sign());
//    
//    double crv=1; //GetC(b)
//    double x2r=crv*dx;
//    double f1=track.pz(), f2=f1+dx*crv;
//
//    //if (TMath::Abs(f1)>=0.000004) return kFALSE;
//    //if (TMath::Abs(f2)>=0.000004) return kFALSE;
//
//    double r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
//    double dy2dx = (f1+f2)/(r1+r2);
//    xyz[0] = r;
//    xyz[1] = track.px()+dx*dy2dx;
//    xyz[2] = track.py()+dx*(r2+f2*dy2dx)*
//
//  }

  template <typename T>
  void fillDataHistograms(T const& jets, float weight = 1.0)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);
      for (auto& track : jet.template tracks_as<JetTracks>()) {
        auto chargeFrag = track.pt() / jet.pt();
        registry.fill(HIST("h_jet_chargefrag"), chargeFrag, weight);
        registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), track.pt(), weight);
        registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), track.eta(), weight);
        registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), track.phi(), weight);
        registry.fill(HIST("h_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_track_phi"), track.phi(), weight);
        registry.fill(HIST("h_track_dcaXY"), track.dcaXY(), weight);
        registry.fill(HIST("h_track_dcaZ"), track.dcaZ(), weight);
      }
      for (auto& d0candidate : jet.template hfcandidates_as<CandidateD0Data>()) {
        auto massD0 = invMassD0ToPiK(d0candidate);
        auto D0Frag = d0candidate.pt() / jet.pt();
        auto D0Prong0Frag = d0candidate.ptProng0() / jet.pt();
        auto D0Prong1Frag = d0candidate.ptProng1() / jet.pt();
        registry.fill(HIST("h_jet_d0frag"), D0Frag, weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_pt"), jet.pt(), d0candidate.pt(), weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_eta"), jet.pt(), d0candidate.eta(), weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_phi"), jet.pt(), d0candidate.phi(), weight);
        registry.fill(HIST("h_d0candidate_pt"), d0candidate.pt(), weight);
        registry.fill(HIST("h_d0candidate_eta"), d0candidate.eta(), weight);
        registry.fill(HIST("h_d0candidate_phi"), d0candidate.phi(), weight);

        registry.fill(HIST("h_jet_d0Prong0frag"), D0Prong0Frag, weight);
        registry.fill(HIST("h_jet_d0Prong1frag"), D0Prong1Frag, weight);
        registry.fill(HIST("h_d0candidate_prong0_pt"), d0candidate.ptProng0(), weight);
        registry.fill(HIST("h_d0candidate_prong1_pt"), d0candidate.ptProng1(), weight);
        registry.fill(HIST("h_decay_length"), d0candidate.decayLength(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_xy"), d0candidate.decayLengthXY(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_error"), d0candidate.errorDecayLength(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_xy_error"), d0candidate.errorDecayLengthXY(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong0"), d0candidate.impactParameter0(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong01D"), d0candidate.impactParameter0(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong1"), d0candidate.impactParameter1(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong11D"), d0candidate.impactParameter1(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong0_error"), d0candidate.errorImpactParameter0(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong1_error"), d0candidate.errorImpactParameter1(), d0candidate.pt(), weight);

        if (d0candidate.isSelD0() >= 1 || d0candidate.isSelD0bar() >= 1) {
          registry.fill(HIST("h_d0candidate_mass"), massD0);
          registry.fill(HIST("h2_d0candidate_mass_d0candidate_pt"), massD0, d0candidate.pt());
          registry.fill(HIST("h2_d0candidate_mass_d0candidate_phi"), massD0, d0candidate.phi());
        }
      }
    }
  }

  template <typename T, typename U>
  void fillCommonMCPHistograms(T const& tracks, U const& mcparticles)
  {
    if (1) {
	for (const auto& track : tracks) {
	  if (!track.has_mcParticle()){
		LOGF(warning, "No MC particle for track, skip...");
		continue;
	  } 
	  auto particle=track.mcParticle();
	  if (!particle.has_mothers()){
		  LOGF(warning, "No mother particle for particle, skip...");
		  continue;
      }

	  registry.fill(HIST("h_jet_mother_1_particle_pdgcode"), particle.pdgCode()); 
	  registry.fill(HIST("h_jet_mother_1_particle_statuscode"), particle.getGenStatusCode()); 
	  auto mm = particle.template mothers_first_as<aod::McParticles>();
	  if  (!mm.has_mothers()) continue;
	  auto m = mm.template mothers_first_as<aod::McParticles>();

	  auto status = TMath::Abs(m.getGenStatusCode());
	  if (status== 23 || status ==33){
	    registry.fill(HIST("h_jet_mother_1_particle_pdgcode"), m.pdgCode()); 
	    registry.fill(HIST("h_jet_mother_1_particle_statuscode"), m.getGenStatusCode()); 
	    registry.fill(HIST("h_jet_mother_1_particle_pdgcode_statuscode"), m.pdgCode(), m.getGenStatusCode()); 
	  }
    }
    } else {
      for(const auto& particle : mcparticles){
		registry.fill(HIST("h_jet_particle_pdgcode"), particle.pdgCode()); 
		registry.fill(HIST("h_jet_particle_statuscode"), particle.getGenStatusCode()); 
		registry.fill(HIST("h_jet_particle_pdgcode_statuscode"), particle.pdgCode(), particle.getGenStatusCode()); 
      }
    }
  }


//  template <typename T>
//  void fillD0JetMCPHistograms(T const& d0jets, float weight = 1.0)
//  {
//    for (const auto& d0jet : d0jets) {
//      const int jetflavor = getJetFlavorHadron(d0jet);
//      if (jetflavor==1) {
//
//      }
//      if (jetflavor==2) {
//
//      }
//      if (jetflavor==3) {
//
//      }
//      registry.fill(HIST("h_jet_pt"), d0jet.pt(), weight);
//
//    }
//  }

  template <typename T>
  void fillResoFuncQualityClass(T const& jets)
  {
    for (const auto& jet : jets){
      const int jetflavor = getJetFlavorAll(jet);
      for (auto& track : jet.template tracks_as<JetTracks>()) {
        int QualityClass = 0;
        if (getQualityTrack(track)) QualityClass = getQualityTrack(track);
		    auto dcaXY=track.dcaXY();
		    auto dcaXYSig=dcaXY / TMath::Sqrt(track.sigmaDcaXY2());
        switch (QualityClass) {
          case 1:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class1"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class1"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class1"), dcaXYSig); 
            }
            continue;
          case 2:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class2"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class2"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class2"), dcaXYSig); 
            }
            continue;
          case 3:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class3"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class3"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class3"), dcaXYSig); 
            }
            continue;
          case 4:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class4"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class4"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class4"), dcaXYSig); 
            }
            continue;
          case 5:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class5"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class5"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class5"), dcaXYSig); 
            }
            continue;
          case 6:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class6"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class6"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class6"), dcaXYSig); 
            }
            continue;
          case 7:
            if (jetflavor==1){
		          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_class7"), dcaXYSig); 
            }
            if (jetflavor==2){
		          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_class7"), dcaXYSig); 
            }
            if (jetflavor==3){
		          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_class7"), dcaXYSig); 
            }
            continue;
        }
      }
    }
  }

  template <typename T, typename U>
  void fillIPsMCPHistograms(T const& jets, U const& mcParticles, float weight = 1.0)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      std::vector<double> lfjetTracksImpXY, cjetTracksImpXY, bjetTracksImpXY;
      std::vector<double> lfjetTracksImpXYSig, cjetTracksImpXYSig, bjetTracksImpXYSig;
      std::vector<double> lfjetTracksImpXYZ, cjetTracksImpXYZ, bjetTracksImpXYZ;
      std::vector<double> lfjetTracksImpXYZSig, cjetTracksImpXYZSig, bjetTracksImpXYZSig;

      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);

      //const int jetflavor = GetJetFlavor(jet, mcParticles);
      //const int jetflavor = getJetFlavorAll(jet);
      for (auto& particle: jet.template tracks_as<aod::McParticles>()) {
        registry.fill(HIST("h_no_mother_particle_pdgcode"), particle.pdgCode()); 
        registry.fill(HIST("h_no_mother_particle_statuscode"), particle.getGenStatusCode()); 
//        if (!track.has_mcParticle()){
//          LOGF(warning, "No MC particle for track, skip...");
//          continue;
//        } 

//        if(!track.isPrimaryTrack()) {
//          LOGF(warning, "No primary track, skip...");
//          continue;
//        }
//
//        auto particle=track.mcParticle();
//        if (!particle.has_mothers()){
//          LOGF(warning, "No mother particle for particle, skip...");
//          registry.fill(HIST("h_no_mother_particle_pdgcode"), particle.pdgCode()); 
//          registry.fill(HIST("h_no_mother_particle_statuscode"), particle.getGenStatusCode()); 
//          continue;
//        }
//        // IP and IPs with XY and XYZ
//        Preslice<aod::Tracks> perCol = aod::track::dcaXY;
//        auto track = particle.sliceBy(perCol, particle.globalIndex());
//        auto dcaXY=track.dcaXY();
//        auto dcaXYSig=dcaXY / TMath::Sqrt(track.sigmaDcaXY2());
//        auto dcaXYZ=GetDCAXYZ(track);
//        auto dcaXYZSig=dcaXYZ/TMath::Sqrt(TMath::Sqrt(track.sigmaDcaXY2()*track.sigmaDcaXY2()+track.sigmaDcaZ2()*track.sigmaDcaZ2())); // TODO: cal cov to add sigmaDcaZY()
//        registry.fill(HIST("h_inc_jet_impact_parameter_xy"), dcaXY); 
//        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance"), dcaXYSig); 
//        registry.fill(HIST("h_inc_jet_impact_parameter_xyz"), dcaXYZ); 
//        registry.fill(HIST("h_inc_jet_impact_parameter_xyz_significance"), dcaXYZSig); 
//        if(jetflavor==1){ // lf
//          registry.fill(HIST("h_lfjet_track_pt"), track.pt()); 
//          registry.fill(HIST("h_lfjet_ntracks"), jet.tracks().size(), weight);
//          registry.fill(HIST("h_lfjet_impact_parameter_xy"), dcaXY); 
//          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance"), dcaXYSig); 
//          registry.fill(HIST("h_lfjet_impact_parameter_xyz"), dcaXYZ); 
//          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance"), dcaXYZSig); 
//          lfjetTracksImpXY.push_back(dcaXY);
//          lfjetTracksImpXYSig.push_back(dcaXYSig);
//          lfjetTracksImpXYZ.push_back(dcaXYZ);
//          lfjetTracksImpXYZSig.push_back(dcaXYZSig);
//        }
//        if(jetflavor==2){ // charm
//          registry.fill(HIST("h_cjet_track_pt"), track.pt()); 
//          registry.fill(HIST("h_cjet_ntracks"), jet.tracks().size(), weight);
//          registry.fill(HIST("h_cjet_impact_parameter_xy"), dcaXY); 
//          registry.fill(HIST("h_cjet_impact_parameter_xy_significance"), dcaXYSig); 
//          registry.fill(HIST("h_cjet_impact_parameter_xyz"), dcaXYZ); 
//          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance"), dcaXYZSig); 
//          cjetTracksImpXY.push_back(dcaXY);
//          cjetTracksImpXYSig.push_back(dcaXYSig);
//          cjetTracksImpXYZ.push_back(dcaXYZ);
//          cjetTracksImpXYZSig.push_back(dcaXYZSig);
//        }
//        if(jetflavor==3){ // beauty
//          registry.fill(HIST("h_bjet_track_pt"), track.pt()); 
//          registry.fill(HIST("h_bjet_ntracks"), jet.tracks().size(), weight);
//          registry.fill(HIST("h_bjet_impact_parameter_xy"), dcaXY); 
//          registry.fill(HIST("h_bjet_impact_parameter_xy_significance"), dcaXYSig); 
//          registry.fill(HIST("h_bjet_impact_parameter_xyz"), dcaXYZ); 
//          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance"), dcaXYZSig); 
//          bjetTracksImpXY.push_back(dcaXY);
//          bjetTracksImpXYSig.push_back(dcaXYSig);
//          bjetTracksImpXYZ.push_back(dcaXYZ);
//          bjetTracksImpXYZSig.push_back(dcaXYZSig);
//        }
      }
    }
  }

  template <typename T, typename U>
  void fillIPsMCDHistograms(T const& jets, U const& mcParticles, float weight = 1.0)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      std::vector<double> lfjetTracksImpXY, cjetTracksImpXY, bjetTracksImpXY;
      std::vector<double> lfjetTracksImpXYSig, cjetTracksImpXYSig, bjetTracksImpXYSig;
      std::vector<double> lfjetTracksImpXYZ, cjetTracksImpXYZ, bjetTracksImpXYZ;
      std::vector<double> lfjetTracksImpXYZSig, cjetTracksImpXYZSig, bjetTracksImpXYZSig;

      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);


      //const int jetflavor = GetJetFlavor(jet, mcParticles);
      const int jetflavor = getJetFlavorAll(jet);
      for (auto& track : jet.template tracks_as<JetTracks>()) {
        if (!track.has_mcParticle()){
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        } 

        if(!track.isPrimaryTrack()) {
          LOGF(warning, "No primary track, skip...");
          continue;
        }

        auto particle=track.mcParticle();
        if (!particle.has_mothers()){
          LOGF(warning, "No mother particle for particle, skip...");
          registry.fill(HIST("h_no_mother_particle_pdgcode"), particle.pdgCode()); 
          registry.fill(HIST("h_no_mother_particle_statuscode"), particle.getGenStatusCode()); 
          continue;
        }

        ////if (!CheckGluonSpliting(track)) continue;
        //std::vector<int> mothersid;
        //SearchingMothers(particle, mothersid);
        //LOGF(info, "track size: %d, mothers size: %d", jet.tracks().size(), mothersid.size());
        //		for (auto &motherid : mothersid) {
        //		  auto mp = mcParticles.iteratorAt(motherid);
        //		  registry.fill(HIST("h_jet_particle_pdgcode"), mp.pdgCode()); 
        //		}
        registry.fill(HIST("h_jet_particle_pdgcode"), particle.pdgCode()); 
        registry.fill(HIST("h_jet_particle_statuscode"), particle.getGenStatusCode()); 
        registry.fill(HIST("h_jet_particle_pdgcode_statuscode"), particle.pdgCode(), particle.getGenStatusCode()); 

        // common
        auto chargeFrag = track.pt() / jet.pt();
        registry.fill(HIST("h_jet_chargefrag"), chargeFrag, weight);
        registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), track.pt(), weight);
        registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), track.eta(), weight);
        registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), track.phi(), weight);
        registry.fill(HIST("h_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_track_phi"), track.phi(), weight);
        registry.fill(HIST("h_track_dcaXY"), track.dcaXY(), weight);
        registry.fill(HIST("h_track_dcaZ"), track.dcaZ(), weight);

        //	    for (auto& mother : particle.template mothers_as<aod::McParticles>()){
        //		  //auto m=mother.template mothers_first_as<aod::McParticles>();
        ////		  auto pdgcode = TMath::Abs(m.getGenStatusCode());
        ////		  if (pdgcode==23 && pdgcode==33) {
        ////			registry.fill(HIST("h_jet_mother_particle_pdgcode"), m.pdgCode()); 
        ////		  }
        //		  
        //		  registry.fill(HIST("h_jet_mother_particle_pdgcode"), mother.pdgCode()); 
        //		  registry.fill(HIST("h_jet_mother_particle_statuscode"), mother.getGenStatusCode()); 
        //		  registry.fill(HIST("h_jet_mother_particle_pdgcode_statuscode"), mother.pdgCode(), mother.getGenStatusCode()); 
        //		}

        //		LOGF(info, "------------------------");
        //		for (auto &m : particle.template mothers_as<aod::McParticles>()) {
        //		  LOGF(info, "pdg: %d, status: %d", m.pdgCode(), m.getGenStatusCode());
        //		}
        //		LOGF(info, "------------------------");
        auto m1 = particle.template mothers_first_as<aod::McParticles>();
        //		LOGF(info, "first mother pdg: %d, status: %d", m1.pdgCode(), m1.getGenStatusCode());
        if (m1.isPhysicalPrimary()){
          registry.fill(HIST("h_jet_mother_1_particle_pdgcode"), m1.pdgCode()); 
          registry.fill(HIST("h_jet_mother_1_particle_statuscode"), m1.getGenStatusCode()); 
          registry.fill(HIST("h_jet_mother_1_particle_pdgcode_statuscode"), m1.pdgCode(), m1.getGenStatusCode()); 
        }
        if (!m1.has_mothers()) continue;
        auto m2 = m1.template mothers_first_as<aod::McParticles>();
        registry.fill(HIST("h_jet_mother_2_particle_pdgcode"), m2.pdgCode()); 
        registry.fill(HIST("h_jet_mother_2_particle_statuscode"), m2.getGenStatusCode()); 
        registry.fill(HIST("h_jet_mother_2_particle_pdgcode_statuscode"), m2.pdgCode(), m2.getGenStatusCode()); 
        if (!m2.has_mothers()) continue;
        auto m3 = m2.template mothers_first_as<aod::McParticles>();
        registry.fill(HIST("h_jet_mother_3_particle_pdgcode"), m3.pdgCode()); 
        registry.fill(HIST("h_jet_mother_3_particle_statuscode"), m3.getGenStatusCode()); 
        registry.fill(HIST("h_jet_mother_3_particle_pdgcode_statuscode"), m3.pdgCode(), m3.getGenStatusCode()); 
        if (!m3.has_mothers()) continue;
        auto m4 = m3.template mothers_first_as<aod::McParticles>();
        registry.fill(HIST("h_jet_mother_4_particle_pdgcode"), m4.pdgCode()); 
        registry.fill(HIST("h_jet_mother_4_particle_statuscode"), m4.getGenStatusCode()); 
        registry.fill(HIST("h_jet_mother_4_particle_pdgcode_statuscode"), m4.pdgCode(), m4.getGenStatusCode()); 
        if (!m4.has_mothers()) continue;
        auto m5 = m4.template mothers_first_as<aod::McParticles>();
        registry.fill(HIST("h_jet_mother_5_particle_pdgcode"), m5.pdgCode()); 
        registry.fill(HIST("h_jet_mother_5_particle_statuscode"), m5.getGenStatusCode()); 
        registry.fill(HIST("h_jet_mother_5_particle_pdgcode_statuscode"), m5.pdgCode(), m5.getGenStatusCode()); 
        //LOGF(info, "PDG m1: %d, m2: %d, m3: %d, m4: %d, m5: %d", m1.pdgCode(), m2.pdgCode(), m3.pdgCode(), m4.pdgCode(), m5.pdgCode());
        //LOGF(info, "Status m1: %d, m2: %d, m3: %d, m4: %d, m5: %d", m1.getGenStatusCode(), m2.getGenStatusCode(), m3.getGenStatusCode(), m4.getGenStatusCode(), m5.getGenStatusCode());

        //		int res =0;
        //		while(1){
        //		  auto m = particle.template mothers_first_as<aod::McParticles>();
        //		  if (!m.has_mothers()) continue;
        //		  res++;
        //		  auto statusCode = TMath::Abs(m.getGenStatusCode());
        //		  if (statusCode==23) break;
        //		  particle = m.template mothers_first_as<aod::McParticles>();
        //		  if (res>100) break;
        //		}
        //	    LOGF(info, "MOTHER SIZE: %d", res);


        //	    for (auto& mother : particle.template mothers_as<aod::mcparticles>()){
        //		  auto m=mother.template mothers_first_as<aod::mcparticles>();
        //		  registry.fill(hist("h_jet_first_mother_particle_pdgcode"), m.pdgcode()); 
        //		  registry.fill(hist("h_jet_first_mother_particle_statuscode"), m.getgenstatuscode()); 
        //		  registry.fill(hist("h_jet_first_mother_particle_pdgcode_statuscode"), m.pdgcode(), m.getgenstatuscode()); 
        //
        //		}

        //		for (auto& mother : particle.template mothers_as<aod::McParticles>()){
        //		  auto m=mother.template mothers_first_as<aod::McParticles>();
        //		  for (auto& mmm : m.template mothers_as<aod::McParticles>()) {
        //			//auto mmm = mm.template mothers_first_as<aod::McParticles>();
        //		    registry.fill(HIST("h_jet_first_mother_particle_pdgcode"), mmm.pdgCode()); 
        //		    registry.fill(HIST("h_jet_first_mother_particle_statuscode"), mmm.getGenStatusCode()); 
        //		    registry.fill(HIST("h_jet_first_mother_particle_pdgcode_statuscode"), mmm.pdgCode(), mmm.getGenStatusCode()); 
        //		  }
        //		}

        //		for (auto& mother : particle.template mothers_as<aod::McParticles>()){
        //		  auto m=mother.template mothers_first_as<aod::McParticles>();
        //		  for (auto& mm : m.template mothers_as<aod::McParticles>()) {
        //			auto mmm = mm.template mothers_first_as<aod::McParticles>();
        //			for (auto& mmmm : mmm.template mothers_as<aod::McParticles>()){
        //		      registry.fill(HIST("h_jet_first_mother_particle_pdgcode"), mmmm.pdgCode()); 
        //		      registry.fill(HIST("h_jet_first_mother_particle_statuscode"), mmmm.getGenStatusCode()); 
        //		      registry.fill(HIST("h_jet_first_mother_particle_pdgcode_statuscode"), mmmm.pdgCode(), mmmm.getGenStatusCode()); 
        //			}
        //		  }
        //		}

        // IP and IPs with XY and XYZ
        auto dcaXY=track.dcaXY()*track.sign();
        auto dcaXYSig=dcaXY / TMath::Sqrt(track.sigmaDcaXY2());
        auto dcaXYZ=GetDCAXYZ(track);
        auto dcaXYZSig=dcaXYZ/TMath::Sqrt(TMath::Sqrt(track.sigmaDcaXY2()*track.sigmaDcaXY2()+track.sigmaDcaZ2()*track.sigmaDcaZ2())); // TODO: cal cov to add sigmaDcaZY()
        registry.fill(HIST("h_inc_jet_impact_parameter_xy"), dcaXY); 
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance"), dcaXYSig); 
        registry.fill(HIST("h_inc_jet_impact_parameter_xy"), dcaXY); 
        registry.fill(HIST("h_inc_jet_impact_parameter_xy_significance"), dcaXYSig); 
        if(jetflavor==1){ // lf
                          //registry.fill(HIST("h_") + HIST(charJetFlavor[jetflavor]) + HIST("_track_pt"), track.pt()); 
          registry.fill(HIST("h_lfjet_track_pt"), track.pt()); 
          registry.fill(HIST("h_lfjet_ntracks"), jet.tracks().size(), weight);
          registry.fill(HIST("h_lfjet_impact_parameter_xy"), dcaXY); 
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance"), dcaXYSig); 
          registry.fill(HIST("h_lfjet_impact_parameter_xyz"), dcaXYZ); 
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance"), dcaXYZSig); 
          registry.fill(HIST("h_lfjet_particle_pdgcode"), track.mcParticle().pdgCode()); 
          registry.fill(HIST("h_lfjet_particle_statuscode"), track.mcParticle().getGenStatusCode()); 
          for (auto& mother : particle.template mothers_as<aod::McParticles>()){
            auto m=mother.template mothers_first_as<aod::McParticles>();
            registry.fill(HIST("h_lfjet_particle_pdgcode"), m.pdgCode()); 
            registry.fill(HIST("h_lfjet_particle_statuscode"), m.getGenStatusCode()); 
          }
          lfjetTracksImpXY.push_back(dcaXY);
          lfjetTracksImpXYSig.push_back(dcaXYSig);
          lfjetTracksImpXYZ.push_back(dcaXYZ);
          lfjetTracksImpXYZSig.push_back(dcaXYZSig);
        }
        if(jetflavor==2){ // charm
                          //registry.fill(HIST(subDir[0]), track.pt()); 
          registry.fill(HIST("h_cjet_track_pt"), track.pt()); 
          registry.fill(HIST("h_cjet_ntracks"), jet.tracks().size(), weight);
          registry.fill(HIST("h_cjet_impact_parameter_xy"), dcaXY); 
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance"), dcaXYSig); 
          registry.fill(HIST("h_cjet_impact_parameter_xyz"), dcaXYZ); 
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance"), dcaXYZSig); 
          for (auto& mother : particle.template mothers_as<aod::McParticles>()){
            auto m=mother.template mothers_first_as<aod::McParticles>();
            registry.fill(HIST("h_cjet_particle_pdgcode"), m.pdgCode()); 
            registry.fill(HIST("h_cjet_particle_statuscode"), m.getGenStatusCode()); 
          }
          cjetTracksImpXY.push_back(dcaXY);
          cjetTracksImpXYSig.push_back(dcaXYSig);
          cjetTracksImpXYZ.push_back(dcaXYZ);
          cjetTracksImpXYZSig.push_back(dcaXYZSig);
        }
        if(jetflavor==3){ // beauty
          registry.fill(HIST("h_bjet_track_pt"), track.pt()); 
          registry.fill(HIST("h_bjet_ntracks"), jet.tracks().size(), weight);
          registry.fill(HIST("h_bjet_impact_parameter_xy"), dcaXY); 
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance"), dcaXYSig); 
          registry.fill(HIST("h_bjet_impact_parameter_xyz"), dcaXYZ); 
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance"), dcaXYZSig); 
          for (auto& mother : particle.template mothers_as<aod::McParticles>()){
            auto m=mother.template mothers_first_as<aod::McParticles>();
            registry.fill(HIST("h_bjet_particle_pdgcode"), m.pdgCode()); 
            registry.fill(HIST("h_bjet_particle_statuscode"), m.getGenStatusCode()); 
          }
          bjetTracksImpXY.push_back(dcaXY);
          bjetTracksImpXYSig.push_back(dcaXYSig);
          bjetTracksImpXYZ.push_back(dcaXYZ);
          bjetTracksImpXYZSig.push_back(dcaXYZSig);
        }
        if(jetflavor==0){ // undefiend
          for (auto& mother : particle.template mothers_as<aod::McParticles>()){
            auto m=mother.template mothers_first_as<aod::McParticles>();
            registry.fill(HIST("h_undefined_particle_pdgcode"), m.pdgCode()); 
            registry.fill(HIST("h_undefined_particle_statuscode"), m.getGenStatusCode()); 
          }
        }
      }
      // Track Counting 
      if(jetflavor==1 && !(lfjetTracksImpXY.empty())){ // lfjet
        sort(lfjetTracksImpXY.begin(), lfjetTracksImpXY.end(), std::greater<double>());
        sort(lfjetTracksImpXYSig.begin(), lfjetTracksImpXYSig.end(), std::greater<double>());
        sort(lfjetTracksImpXYZ.begin(), lfjetTracksImpXYZ.end(), std::greater<double>());
        sort(lfjetTracksImpXYZSig.begin(), lfjetTracksImpXYZSig.end(), std::greater<double>());
        if (lfjetTracksImpXY.size()>0) { // N1
          registry.fill(HIST("h_lfjet_impact_parameter_xy_N1"), lfjetTracksImpXY[0]);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N1"), lfjetTracksImpXYSig[0]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_N1"), lfjetTracksImpXYZ[0]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N1"), lfjetTracksImpXYZSig[0]);
        }
        if (lfjetTracksImpXY.size()>1) { // N2
          registry.fill(HIST("h_lfjet_impact_parameter_xy_N2"), lfjetTracksImpXY[1]);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N2"), lfjetTracksImpXYSig[1]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_N2"), lfjetTracksImpXYZ[1]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N2"), lfjetTracksImpXYZSig[1]);
        }
        if (lfjetTracksImpXY.size()>2) { // N3
          registry.fill(HIST("h_lfjet_impact_parameter_xy_N3"), lfjetTracksImpXY[2]);
          registry.fill(HIST("h_lfjet_impact_parameter_xy_significance_N3"), lfjetTracksImpXYSig[2]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_N3"), lfjetTracksImpXYZ[2]);
          registry.fill(HIST("h_lfjet_impact_parameter_xyz_significance_N3"), lfjetTracksImpXYZSig[2]);
        }
      }
      if(jetflavor==2 && !(cjetTracksImpXY.empty())) { // cjet
        sort(cjetTracksImpXY.begin(), cjetTracksImpXY.end(), std::greater<double>());
        sort(cjetTracksImpXYSig.begin(), cjetTracksImpXYSig.end(), std::greater<double>());
        sort(cjetTracksImpXYZ.begin(), cjetTracksImpXYZ.end(), std::greater<double>());
        sort(cjetTracksImpXYZSig.begin(), cjetTracksImpXYZSig.end(), std::greater<double>());
        if (cjetTracksImpXY.size()>0) { // N1
          registry.fill(HIST("h_cjet_impact_parameter_xy_N1"), cjetTracksImpXY[0]);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N1"), cjetTracksImpXYSig[0]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_N1"), cjetTracksImpXYZ[0]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N1"), cjetTracksImpXYZSig[0]);
        }
        if (cjetTracksImpXY.size()>1) { // N2
          registry.fill(HIST("h_cjet_impact_parameter_xy_N2"), cjetTracksImpXY[1]);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N2"), cjetTracksImpXYSig[1]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_N2"), cjetTracksImpXYZ[1]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N2"), cjetTracksImpXYZSig[1]);
        }
        if (cjetTracksImpXY.size()>2) { // N3
          registry.fill(HIST("h_cjet_impact_parameter_xy_N3"), cjetTracksImpXY[2]);
          registry.fill(HIST("h_cjet_impact_parameter_xy_significance_N3"), cjetTracksImpXYSig[2]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_N3"), cjetTracksImpXYZ[2]);
          registry.fill(HIST("h_cjet_impact_parameter_xyz_significance_N3"), cjetTracksImpXYZSig[2]);
        }
      }
      if(jetflavor==3 && !(bjetTracksImpXY.empty())){ // bjet
        sort(bjetTracksImpXY.begin(), bjetTracksImpXY.end(), std::greater<double>());
        sort(bjetTracksImpXYSig.begin(), bjetTracksImpXYSig.end(), std::greater<double>());
        sort(bjetTracksImpXYZ.begin(), bjetTracksImpXYZ.end(), std::greater<double>());
        sort(bjetTracksImpXYZSig.begin(), bjetTracksImpXYZSig.end(), std::greater<double>());
        if (bjetTracksImpXY.size()>0) { // N1
          registry.fill(HIST("h_bjet_impact_parameter_xy_N1"), bjetTracksImpXY[0]);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N1"), bjetTracksImpXYSig[0]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_N1"), bjetTracksImpXYZ[0]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N1"), bjetTracksImpXYZSig[0]);
        }
        if (bjetTracksImpXY.size()>1) { // N2
          registry.fill(HIST("h_bjet_impact_parameter_xy_N2"), bjetTracksImpXY[1]);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N2"), bjetTracksImpXYSig[1]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_N2"), bjetTracksImpXYZ[1]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N2"), bjetTracksImpXYZSig[1]);
        }
        if (bjetTracksImpXY.size()>2) { // N3
          registry.fill(HIST("h_bjet_impact_parameter_xy_N3"), bjetTracksImpXY[2]);
          registry.fill(HIST("h_bjet_impact_parameter_xy_significance_N3"), bjetTracksImpXYSig[2]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_N3"), bjetTracksImpXYZ[2]);
          registry.fill(HIST("h_bjet_impact_parameter_xyz_significance_N3"), bjetTracksImpXYZSig[2]);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillJetProbMCDHistograms(T const& jets, U const& mcParticles, float weight = 1.0)
  {
    for (const auto& jet : jets) {

      //int jetflavor = GetJetFlavor(jet, mcParticles);
      const int jetflavor = getJetFlavorAll(jet);
      double JP = -1.;
      JP = CalculateJP(jet, mcParticles);
      if (JP < 0) {
        //LOGF(warning, "No JP due to track size, skip...");
        continue;
      }
      //LOGF(info, "flavor: %d, JP: %f", jetflavor, JP);

      if(jetflavor==1) { // lf
        registry.fill(HIST("h_lfjet_JP"), JP); 
        registry.fill(HIST("h_lfjet_JP_Log"), -1*TMath::Log(JP)); 
      }
      if(jetflavor==2) { // charm 
        registry.fill(HIST("h_cjet_JP"), JP); 
        registry.fill(HIST("h_cjet_JP_Log"), -1*TMath::Log(JP)); 
      }
      if(jetflavor==3) { // beauty
        registry.fill(HIST("h_bjet_JP"), JP); 
        registry.fill(HIST("h_bjet_JP_Log"), -1*TMath::Log(JP)); 
      }
    }
  }

  template <typename T, typename U>
  void fillMCDHistogramsHf(T const& recjets, U const& genjets, float weight = 1.0)
  {
    for (const auto& jet: genjets){
      registry.fill(HIST("h_part_jet_pt"), jet.pt(), weight);
      for (auto& constituent : jet.template tracks_as<JetParticles2Prong>()) {
	    LOGF(info, "TEST: %d", constituent.originMcGen());
        registry.fill(HIST("h_constituent_pt"), constituent.pt(), weight);
	  }
    }
    for (const auto& jet : recjets) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);
//      for (auto& track : jet.template tracks_as<JetTracks>()){
//		LOGF(info, "track pt : %f", track.pt());
//		if (!track.has_mcParticle()){
//		  LOGF(warning, "No MC particle for track, skip...");
//		  continue;
//		} 
//		auto particle =track.mcParticle();
//		LOGF(info, "part: pt : %f", particle.pt());
//	  }
      for (auto& d0candidate : jet.template hfcandidates_as<CandidateD0MC>()) {
      //for (auto& d0candidate : jet.template tracks_as<JetParticles2Prong>()) {
		//LOGF(info, "TEST: %d", d0candidate.originMcGen());
        registry.fill(HIST("h2_jet_pt_d0candidate_pt"), jet.pt(), d0candidate.pt(), weight);
	  }
    }
  }

  template <typename T, typename U, typename V, typename K>
  void fillMCPHistograms(T const& collision, U const& partjets, V const& particles, K const& tracks, float weight = 1.0)
  {
	for (auto& particle : particles){
      registry.fill(HIST("h_particle_pt_temp"), particle.pt(), weight);
    }
	for (auto& track : tracks){
      registry.fill(HIST("h_track_pt_temp"), track.pt(), weight);
	}
    for (const auto& jet : partjets) {
      //registry.fill(HIST("h_part_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_part_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_part_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_part_jet_ntracks"), jet.tracks().size(), weight);
	  o2::vertexing::DCAFitterN<2> df;
      for (auto& track : jet.template tracks_as<JetParticles2Prong>()) {
        //registry.fill(HIST("h_part_jet_pt"), track.pt(), weight);
		registry.fill(HIST("h_mc_pdgcode"), track.pdgCode(), weight);
		registry.fill(HIST("h_track_pt_jet_temp"), track.pt(), weight);
      }
      for (auto& d0candidate : jet.template hfcandidates_as<JetParticles2Prong>()) {
		registry.fill(HIST("h_mc_pdgcode"), d0candidate.pdgCode(), weight);
		registry.fill(HIST("h_hftrack_pt_jet_temp"), d0candidate.pt(), weight);
      }
    }
  }

  template <typename T>
  void fillMCMatchedHistograms(T const& mcdjet, float weight = 1.0)
  {
    if (mcdjet.has_matchedJetCand() && mcdjet.matchedJetCandId() >= 0) {
      const auto& mcpjet = mcdjet.template matchedJetCand_as<soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>>();
      registry.fill(HIST("h2_jet_pt_part_jet_pt"), mcpjet.pt(), mcdjet.pt(), weight);
    }
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(bJetMeasurementTask, processDummy, "Dummy process function turned on by default", true);

  void processFitResoFunc(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
                      soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                      //CandidateD0MC const& candidates,
                      JetTracks const& tracks,
					  aod::McCollisions const& mcCollisions,
					  aod::McParticles const& mcParticles)
  {
    TH1D *hlfJetIPsXY = new TH1D("hlfJetIPsXY", "", 1000,-40,40);
    TH1D *hcJetIPsXY = new TH1D("hcJetIPsXY", "", 1000,-40,40);
    TH1D *hbJetIPsXY = new TH1D("hbJetIPsXY", "", 1000,-40,40);
    TF1 *fResoFitFunclfjet = new TF1("fResoFitFunclfjet", "gaus(0)", -40, 40);
    TF1 *fResoFitFunccjet = new TF1("fResoFitFunccjet", "gaus(0)", -40, 40);
    TF1 *fResoFitFuncbjet = new TF1("fResoFitFuncbjet", "gaus(0)", -40, 40);

	for (auto& collision : collisions) {
	  if(!collision.posZ()) continue;
      for (auto& jet : jets){
	    //int jetflavor = GetJetFlavor(jet, mcParticles);
	    int jetflavor = getJetFlavorAll(jet);
      for (auto& track : jet.template tracks_as<JetTracks>()) {
	      if (!track.has_mcParticle()){
		    LOGF(warning, "No MC particle for track, skip...");
		    continue;
	    } 
	    auto particle=track.mcParticle();
	    if (!particle.has_mothers()){
		    LOGF(warning, "No mother particle for particle, skip...");
		    continue;
      }
	    //if (track.sign()<0) continue; // only take positive track

		  double dcaXYSig = track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2());
          //if (TMath::Abs(dcaXYSig)>100) dcaXYSig = 99.9; // Limit to function definition range
		  if(jetflavor==1){ // lf
		    hlfJetIPsXY->Fill(dcaXYSig);
		  }
		  if(jetflavor==2){ // charm
		    hcJetIPsXY->Fill(dcaXYSig);
		  }
		  if(jetflavor==3){ // beauty
		    hbJetIPsXY->Fill(dcaXYSig);
		  }
        }
      }
    }

    // TODO:: set fitting function
    hlfJetIPsXY->Fit("fResoFitFunclfjet");
    hcJetIPsXY->Fit("fResoFitFunccjet");
    hbJetIPsXY->Fit("fResoFitFuncbjet");
	for (int i=0; i<3; i++){
	  FitParlfjet[i] = fResoFitFunclfjet->GetParameter(i);
	  FitParcjet[i] = fResoFitFunccjet->GetParameter(i);
	  FitParbjet[i] = fResoFitFuncbjet->GetParameter(i);
    }
  }
  PROCESS_SWITCH(bJetMeasurementTask, processFitResoFunc, "Task of jet fragmentation for heavy flavor (MCD)", false);

  void processJetsData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                       soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                       CandidateD0Data const& candidates,
                       JetTracks const& tracks)
  {
    fillDataHistograms(jets);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processJetsData, "Task of jet fragmentation for heavy flavor (Data)", false);

  void processFitResoFuncMCP(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                      soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                      //CandidateD0MC const& candidates,
                      JetTracks const& tracks,
					  aod::McCollisions const& mcCollisions,
					  aod::McParticles const& mcParticles)
  {
    fillIPsMCDHistograms(jets, mcParticles);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processFitResoFuncMCP, "Task of jet fragmentation for heavy flavor (MCD)", false);

//  void processIPsHfJetsMCP(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
//                      soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents> const& jets,
//                      //CandidateD0MC const& candidates,
//                      JetTracks const& tracks,
//					            aod::McCollisions const& mcCollisions,
//					            aod::McParticles const& mcParticles)
//  {
//	  //LOGF(info, "vtx-z (data) = %f | vtx-z (MC) =%f", collision.posZ(), collision.mcCollision().posZ());
//    fillD0JetMCPHistograms(jets);
//	  //fillCommonMCPHistograms(tracks, mcParticles);
//  }
//  PROCESS_SWITCH(bJetMeasurementTask, processIPsHfJetsMCP, "Task of jet fragmentation for heavy flavor (MCD)", false);

  void processIPsJetsMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                      soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                      //CandidateD0MC const& candidates,
                      JetTracks const& tracks,
					  aod::McCollisions const& mcCollisions,
					  aod::McParticles const& mcParticles)
  {
	//LOGF(info, "vtx-z (data) = %f | vtx-z (MC) =%f", collision.posZ(), collision.mcCollision().posZ());
    fillIPsMCDHistograms(jets, mcParticles);
	//fillCommonMCPHistograms(tracks, mcParticles);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processIPsJetsMCD, "Task of jet fragmentation for heavy flavor (MCD)", false);

  void processIPsJetsMCP(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                      soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                      //CandidateD0MC const& candidates,
                      JetTracks const& tracks,
					  aod::McCollisions const& mcCollisions,
					  aod::McParticles const& mcParticles)
  {
	//LOGF(info, "vtx-z (data) = %f | vtx-z (MC) =%f", collision.posZ(), collision.mcCollision().posZ());
    fillIPsMCPHistograms(jets, mcParticles);
	//fillCommonMCPHistograms(tracks, mcParticles);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processIPsJetsMCP, "Task of jet fragmentation for heavy flavor (MCD)", false);



  void processCreateResoFunc(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                      soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                      //CandidateD0MC const& candidates,
                      JetTracks const& tracks,
					  aod::McCollisions const& mcCollisions,
					  aod::McParticles const& mcParticles)
  {
	//LOGF(info, "vtx-z (data) = %f | vtx-z (MC) =%f", collision.posZ(), collision.mcCollision().posZ());
    fillResoFuncQualityClass(jets);
	//fillCommonMCPHistograms(tracks, mcParticles);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processCreateResoFunc, "Task of jet fragmentation for heavy flavor (MCD)", false);



  void processJetProbMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                      soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                      //CandidateD0MC const& candidates,
                      JetTracks const& tracks,
					  aod::McCollisions const& mcCollisions,
					  aod::McParticles const& mcParticles)
  {
	//LOGF(info, "vtx-z (data) = %f | vtx-z (MC) =%f", collision.posZ(), collision.mcCollision().posZ());
    fillJetProbMCDHistograms(jets, mcParticles);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processJetProbMCD, "Task of jet fragmentation for heavy flavor (MCD)", false);


  void processHfJetsMCD(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                      soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents> const& hfjets,
                      soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents> const& genjets,
                      CandidateD0MC const& candidates,
                      JetParticles2Prong const& particles,
                      JetTracks const& tracks)
  {
	//LOGF(info, "vtx-z (data) = %f | vtx-z (MC) =%f", collision.posZ(), collision.mcCollision().posZ());
    fillMCDHistogramsHf(hfjets, genjets);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processHfJetsMCD, "Task of jet fragmentation for heavy flavor (MCD)", false);


  void processJetsMCP(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
//                      soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents> const& hfjets,
                      soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& hfjets,
                      JetParticles2Prong const& particles,
                      JetTracks const& tracks)
  {
    fillMCPHistograms(collision, hfjets, particles, tracks);
  }
  PROCESS_SWITCH(bJetMeasurementTask, processJetsMCP, "Task of jet fragmentaion for heavy flavor (MCP)", false);

  void processJetsMCPMCDMatched(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                                soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets> const& mcdjets,
                                soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets> const& mcpjets,
                                JetParticles2Prong const& particles,
                                JetTracks const& tracks)
  {
    for (const auto& mcdjet : mcdjets) {
      fillMCMatchedHistograms(mcdjet);
    }
  }
  PROCESS_SWITCH(bJetMeasurementTask, processJetsMCPMCDMatched, "Matching of detector level jets and particle level jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<bJetMeasurementTask>(cfgc, TaskName{"bjet-measurement"})}; }

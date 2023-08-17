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

/// \commonly used for HF tagging analysis.
/// \author hanseo.park@cern.ch

#ifndef PWGJE_UTILS_HFTAGGINGUTILITIES_H_
#define PWGJE_UTILS_HFTAGGINGUTILITIES_H_

#include "Framework/AnalysisTask.h"

using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::McTrackLabels>;

template <typename TTrack, typename Vec>
void SearchingMothers(TTrack const& particle, Vec &mothersid){
  if (particle.has_mothers()){
    //auto pdg = TMath::Abs(particle.pdgCode());
    auto status = TMath::Abs(particle.getGenStatusCode());
    if (status==23) {
      if (std::find(mothersid.begin(), mothersid.end(), particle.globalIndex()) == mothersid.end()) mothersid.push_back(particle.globalIndex());
      //LOGF(info, "pdg: %d", pdg);
    }
    for (auto& mother : particle.template mothers_as<o2::aod::McParticles>()) {
      //if (!mother.isPhysicalPrimary()) continue;
      SearchingMothers(mother, mothersid);
    }
  }
} 

template <typename T, typename U>
int GetJetFlavor(T const& jet, U const& mcparticles){
  const int arraySize=99;
  int candCharmLfCode[arraySize], count=0;
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

    std::vector<int> mothersid;
    SearchingMothers(particle, mothersid);
    if (mothersid.size()==0) return 0;
    //LOGF(info, "mothers size: %d", mothersid.size());

    for (auto &motherid : mothersid) {
      auto mp = mcparticles.iteratorAt(motherid);
      int pdgcode=TMath::Abs(mp.pdgCode());
      if (pdgcode==21 || (pdgcode>0 && pdgcode<6)){
        if(pdgcode==5) {
          return 3; // beauty
        }else {
          candCharmLfCode[count]=pdgcode;
          count++;
        }
      }
    }
  }
  if (count>0){
    for (int i=0; i<count; i++){
      if(candCharmLfCode[i]==4) return 2; // charm
    }
  } else {
    return 0;
  }
  if (count>0) return 1; // light flavor
  return 0;
}

template <typename T>
int getJetFlavorAll(T const& jet){
  const int arraySize=99;
  int candCharmLfCode[arraySize], count=0;
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

    if (!particle.isPhysicalPrimary()) continue;
    for (auto& mother : particle.template mothers_as<aod::McParticles>()){
    //auto mother = particle.template mothers_first_as<aod::McParticles>();
      //if (!mother.isPhysicalPrimary()) continue;
      int pdgcode=TMath::Abs(mother.pdgCode());
      if (pdgcode==21 || (pdgcode>0 && pdgcode<6)){
        if(pdgcode==5) {
          return 3; // beauty
        }else {
          candCharmLfCode[count]=pdgcode;
          count++;
        }
      }
    }
  }
  if (count>0){
    for (int i=0; i<count; i++){
      if(candCharmLfCode[i]==4) return 2; // charm
    }
  } else {
    return 0;
  }
  if (count>0) return 1; // light flavor
  return 0;
}

template <typename T>
int getJetFlavorParton(T const& jet){
  const int arraySize=99;
  int candCharmLfCode[arraySize], count=0;
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

    int pdgcode=TMath::Abs(particle.pdgCode());
    if (pdgcode==21 || (pdgcode>0 && pdgcode<6)){
      if(pdgcode==5) {
        return 3; // beauty
      }else {
        candCharmLfCode[count]=pdgcode;
        count++;
      }
    }
  }
  if (count>0){
    for (int i=0; i<count; i++){
      if(candCharmLfCode[i]==4) return 2; // charm
    }
  } else {
    return 0;
  }
  if (count>0) return 1; // light flavor
  return 0;
}

template <typename T>
int getJetFlavorHadron(T const& jet) {
  const int arraySize=99;
  int candCharmLfCode[arraySize], count=0;
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

    int pdgcode=TMath::Abs(particle.pdgCode());
    if (pdgcode==21 || (pdgcode>0 && pdgcode<6)){
      if(pdgcode==5) {
        return 3; // beauty
      }else {
        candCharmLfCode[count]=pdgcode;
        count++;
      }
    }
  }
  if (count>0){
    for (int i=0; i<count; i++){
      if(candCharmLfCode[i]==4) return 2; // charm
    }
  } else {
    return 0;
  }
  if (count>0) return 1; // light flavor
  return 0;
}

#endif // PWGEM_UTILS_HFTAGGINGUTILTIES_H_

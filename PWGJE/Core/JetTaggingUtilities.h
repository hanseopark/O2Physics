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

/// \file JetTaggingUtilities.h
/// \brief Jet tagging related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETTAGGINGUTILITIES_H_
#define PWGJE_CORE_JETTAGGINGUTILITIES_H_

#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

#include "Framework/Logger.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/Core/JetUtilities.h"

namespace JetTaggingUtilities
{
/**
 * returns the globalIndex of the earliest mother of a particle in the shower. returns -1 if a suitable mother is not found
 *
 * @param particle MCParticle whose mother is to be found
 */

template <typename T>
int getOriginalMotherIndex(const typename T::iterator& particle)
{

  auto mother = particle;

  while (mother.has_mothers()) {

    mother = mother.template mothers_first_as<T>();

    int motherStatusCode = std::abs(mother.getGenStatusCode());

    if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63) {
      return mother.globalIndex();
    }
  }
  return -1.0;
}

/**
 * returns the globalIndex of the earliest HF mother of a particle in the shower. returns -1 if a suitable mother is not found. Should be used only on already identified HF particles
 *
 * @param hfparticle MCParticle whose mother is to be found
 */

//template <typename T>
//int getOriginalHFMotherIndex(T &particles, const typename T::iterator& hfparticle)
template <typename T>
int getOriginalHFMotherIndex(const typename T::iterator& hfparticle)
{

  auto mother = hfparticle;

  while (mother.has_mothers()) {

    mother = mother.template mothers_first_as<T>();

    int motherStatusCode = std::abs(mother.getGenStatusCode());

    if (motherStatusCode == 23 || motherStatusCode == 33 || motherStatusCode == 43 || motherStatusCode == 63 || (motherStatusCode == 51 && mother.template mothers_first_as<T>().pdgCode() == 21)) {
      return mother.globalIndex();
    }
  }
  return -1.0;
}

/**
 * checks if atrack in a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The first track originating from an HF shower can be extracted by reference
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param hftrack track passed as reference which is then replaced by the first track that originated from an HF shower
 */

template <typename T, typename U, typename V>
int jetTrackFromHFShower(T const& jet, U const& particles, V const& tracks, typename V::iterator& hftrack)
{#
  std::vector<std::pair<int, typename V::iterator>> candCharmLf;
  int origin = -1;
  for (auto& track : jet.template tracks_as<V>()) {
    if (!track.has_mcParticle()) {
      continue;
    }
      auto const& particle = track.template mcParticle_as<U>();
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, true);

    hftrack = track;
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hftrack = track;
      return origin; 
    }

    if (origin == 2) {
      hftrack = track;
      return origin; // 2=beauty
    }
    candCharmLf.push_back(std::make_pair(origin, track));
  }
  for (const auto& cand : candCharmLf) {
    if(cand.first==1) {
      hftrack = cand.second;
      return 1; // 1-charm
    }
  }
  return 0;
}

/**
 * checks if a particle in a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The first particle originating from an HF shower can be extracted by reference
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param hfparticle particle passed as reference which is then replaced by the first track that originated from an HF shower
 */
template <typename T, typename U>
int jetParticleFromHFShower(T const& jet, U const& particles, typename U::iterator& hfparticle)
{

  int origin = -1;
  for (auto& particle : jet.template tracks_as<U>()) {
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, true);
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hfparticle = particle;
      return origin;
    }
  }
  return 0;
}

/**
 * returns if a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U, typename V>
int mcdJetFromHFShower(T const& jet, U const& particles, V const& tracks, float dRMax = 0.25)
{
  typename V::iterator hftrack;
  int origin = jetTrackFromHFShower(jet, particles, tracks, hftrack);
  if (origin == 0) return 0;
  if (!hftrack.has_mcParticle()) {
    return 0;
  }
  auto const& hfparticle = hftrack.template mcParticle_as<U>();
  if (origin == 1 || origin == 2) {

    int originalHFMotherIndex = getOriginalHFMotherIndex<U>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (JetUtilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return 0;
      }

    } else {
      return 0;
    }
  } else {

    return 0;
  }
}

/**
 * checks if a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

//template <typename T, typename U, typename V>
template <typename T, typename U>
int mcpJetFromHFShower(T const& jet, U const& particles, float dRMax = 0.25)
{

  typename U::iterator hfparticle;
  int origin = jetParticleFromHFShower(jet, particles, hfparticle);
  if (origin == 1 || origin == 2) {

    //int originalHFMotherIndex = getOriginalHFMotherIndex(particles, hfparticle);
    int originalHFMotherIndex = getOriginalHFMotherIndex<U>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (JetUtilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return 0;
      }

    } else {
      return 0;
    }

  } else {

    return 0;
  }
}

/**
 * returns the pdg code of the original scattered parton closest to the jet axis, with the restriction that the parton and jet axis must be within dRMax in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U>
int jetOrigin(T const& jet, U const& particles, float dRMax = 0.25)
{
  bool firstPartonFound = false;
  typename U::iterator parton1;
  typename U::iterator parton2;
  for (auto& particle : particles) {
    if (std::abs(particle.getGenStatusCode() == 23)) {
      if (!firstPartonFound) {
        parton1 = particle;
        firstPartonFound = true;
      } else {
        parton2 = particle;
      }
    }
  }

  float dR1 = JetUtilities::deltaR(jet, parton1);
  float dR2 = JetUtilities::deltaR(jet, parton2);

  if (dR1 <= dR2 && dR1 < dRMax) {

    return parton1.pdgCode();
  }
  if (dR2 <= dR1 && dR2 < dRMax) {

    return parton2.pdgCode();
  }

  return 0;
}

}; // namespace JetTaggingUtilities

#endif // PWGJE_CORE_JETTAGGINGUTILITIES_H_

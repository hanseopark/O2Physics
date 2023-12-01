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
/// \author Hanseo Park <hanseo.park@cern.ch>

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

enum JetTaggingSpecies {
  none = 0,
  charm = 1,
  beauty = 2,
  lightflavour = 3,
  lightquark = 4,
  gluon = 5
};

namespace JetTaggingBinCut
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
}; // namespace JetTaggingBinCut

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

  // if (!particle) {
  //   return -1.0;
  // }
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

template <typename T>
int getOriginalHFMotherIndex(const typename T::iterator& hfparticle)
{

  // if (!hfparticle) {
  //   return -1.0;
  // }
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
int jetTrackFromHFShower(T const& jet, U const& tracks, V const& particles, typename U::iterator& hftrack)
{

  bool hasMcParticle = false;
  int origin = -1;
  for (auto& track : jet.template tracks_as<U>()) {
    if (!track.has_mcParticle()) {
      continue;
    }
    hasMcParticle = true;
    auto const& particle = track.template mcParticle_as<V>();
    origin = RecoDecay::getCharmHadronOrigin(particles, particle, true);
    if (origin == 1 || origin == 2) { // 1=charm , 2=beauty
      hftrack = track;
      if (origin == 1) {
        return JetTaggingSpecies::charm;
      }
      if (origin == 2) {
        return JetTaggingSpecies::beauty;
      }
    }
  }
  if (hasMcParticle) {
    return JetTaggingSpecies::lightflavour;
  } else {
    return JetTaggingSpecies::none;
  }
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
      if (origin == 1) {
        return JetTaggingSpecies::charm;
      }
      if (origin == 2) {
        return JetTaggingSpecies::beauty;
      }
    }
  }
  return JetTaggingSpecies::lightflavour;
}

/**
 * returns if a reco level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U, typename V>
int mcdJetFromHFShower(T const& jet, U const& tracks, V const& particles, float dRMax = 0.25)
{

  typename U::iterator hftrack;
  int origin = jetTrackFromHFShower(jet, tracks, particles, hftrack);
  if (origin == JetTaggingSpecies::charm || origin == JetTaggingSpecies::beauty) {
    if (!hftrack.has_mcParticle()) {
      return JetTaggingSpecies::none;
    }
    auto const& hfparticle = hftrack.template mcParticle_as<V>();

    int originalHFMotherIndex = getOriginalHFMotherIndex<V>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      return origin;
      if (JetUtilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return JetTaggingSpecies::none;
      }

    } else {
      return JetTaggingSpecies::none;
    }

  } else {

    return JetTaggingSpecies::lightflavour;
  }
}

/**
 * checks if a generator level jet originates from a HF shower. 0:no HF shower, 1:charm shower, 2:beauty shower. The requirement is that the jet contains a particle from an HF shower and that the original HF quark is within dRMax of the jet axis in eta-phi
 *
 * @param jet
 * @param particles table of generator level particles to be searched through
 * @param dRMax maximum distance in eta-phi of initiating heavy-flavour quark from the jet axis
 */

template <typename T, typename U, typename V>
int mcpJetFromHFShower(T const& jet, U const& particles, float dRMax = 0.25)
{

  typename U::iterator hfparticle;
  int origin = jetParticleFromHFShower(jet, particles, hfparticle);
  if (origin == JetTaggingSpecies::charm || origin == JetTaggingSpecies::beauty) {

    int originalHFMotherIndex = getOriginalHFMotherIndex<U>(hfparticle);
    if (originalHFMotherIndex > -1.0) {

      if (JetUtilities::deltaR(jet, particles.iteratorAt(originalHFMotherIndex)) < dRMax) {

        return origin;

      } else {
        return JetTaggingSpecies::none;
      }

    } else {
      return JetTaggingSpecies::none;
    }

  } else {

    return JetTaggingSpecies::lightflavour;
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

template <typename T, typename U, typename V>
int getGeoSign(T const& collision, U const& jet, V const& track)
{
  auto sgn = TMath::Sign(1, (track.dcaX() - collision.posX()) * jet.px() + (track.dcaY() - collision.posY()) * jet.py() + (track.dcaZ() - collision.posZ())*jet.pz());
//  std::cout << "dcax: " << track.dcaX() << " collx: " << collision.posX() << " jetpx: " << jet.px() << "\n";
//  std::cout << "dcay: " << track.dcaY() << " collx: " << collision.posY() << " jetpx: " << jet.py() << "\n";
//  std::cout << "dcaz: " << track.dcaZ() << " collx: " << collision.posZ() << " jetpx: " << jet.pz() << "\n";
//  std::cout << "det: " << (track.dcaX() - collision.posX()) * jet.px() + (track.dcaY() - collision.posY()) * jet.py() + (track.dcaZ() - collision.posZ())*jet.pz()
//  std::cout << "sign: " << sgn << "\n";
  return sgn;  
}

template <typename AnyCollision, typename AnyJet, typename Constituent>
int getSignIP(AnyCollision const& collision, AnyJet const& jet, Constituent const& constit) /// Determine Impact paramter sign
{
  // Calculate Sign
  double posdcatrack[3] = {0., 0., 0.};
  posdcatrack[0] = constit.x();
  posdcatrack[1] = constit.y();
  posdcatrack[2] = constit.z();

  // Performs local->global transformation of the track position.
  double cs = TMath::Cos(constit.alpha()), sn = TMath::Sin(constit.alpha()), x = posdcatrack[0];
  posdcatrack[0] = x * cs - posdcatrack[1] * sn;
  posdcatrack[1] = x * sn + posdcatrack[1] * cs;

  // auto collision = jet.template collision_as<aod::Collisions>();

  double ipvector3[3] = {posdcatrack[0] - collision.posX(), posdcatrack[1] - collision.posY(), posdcatrack[2] - collision.posZ()};
  int sign = TMath::Sign(1., ipvector3[0] * jet.px() + ipvector3[1] * jet.py() + ipvector3[2] * jet.pz());

  return sign;
}

template <typename AnyCollision, typename AnyJet, typename Constituent>
float getdcaXGlo(AnyCollision const& collision, AnyJet const& jet, Constituent const& constit) /// Determine Impact paramter sign
{
  float posdcatrack[3] = {0., 0., 0.};
  posdcatrack[0] = constit.x();
  posdcatrack[1] = constit.y();
  posdcatrack[2] = constit.z();

  // Performs local->global transformation of the track position.
  float cs = TMath::Cos(constit.alpha()), sn = TMath::Sin(constit.alpha()), x = posdcatrack[0];
  posdcatrack[0] = x * cs - posdcatrack[1] * sn;
  posdcatrack[1] = x * sn + posdcatrack[1] * cs;

  float ipvector3[3] = {posdcatrack[0] - collision.posX(), posdcatrack[1] - collision.posY(), posdcatrack[2] - collision.posZ()};
  return ipvector3[0];

}

//template <typename T, typename U, typename V>
//void SetSgnImpactParameterSignificance(T const& collision, U const& jet, V const& track, int& sgn, double& IPxy, double& SgnIPxy, double& IPxySig, double& SgnIPxySig, double& IPxyz, double& SgnIPxyz, double& IPxyzSig, double& SgnIPxyzSig)
template <typename T, typename U, typename V>
void SetSgnImpactParameterSignificance(T const& collision, U const& jet, V const& track, int& sgn, float& IPxy, float& SgnIPxy, float& IPxySig, float& SgnIPxySig, float& IPxyz, float& SgnIPxyz, float& IPxyzSig, float& SgnIPxyzSig)
{
  //LOGF(info, Form("coll x: %f, track x: %f", collision.posX(), track.x()));
  //LOGF(info, Form("coll y: %f, track y: %f", collision.posY(), track.y()));
  //LOGF(info, Form("track z: %f", track.z()));
  //sgn = TMath::Sign(1, (track.x() - collision.posX()) * jet.px() + (track.y() - collision.posY()) * jet.py()); // geometric sign +1: from SV, +1,-1: from PV
  //sgn = TMath::Sign(1, (track.x() - collision.posX()) * jet.px() + (track.y() - collision.posY()) * jet.py() + (track.z()-collision.posZ()) * jet.pz()); // geometric sign +1: from SV, +1,-1: from PV
  //sgn = TMath::Sign(1, (track.x() - collision.posX()) * jet.px() + (track.y() - collision.posY()) * jet.py() + (track.z()-collision.posZ()) * jet.pz()); // geometric sign +1: from SV, +1,-1: from PV
  //std::cout << "track x: " << track.x() << " collision x: " << collision.posX() << " jet px: " << jet.px() << "\n";
  sgn = getGeoSign(collision, jet, track);
  //sgn = getSignIP(collision, jet, track);
  // IPxy
  IPxy = track.dcaXY();
  SgnIPxy = sgn * TMath::Abs(IPxy);
  IPxySig = IPxy / TMath::Sqrt(track.sigmaDcaXY2());
  SgnIPxySig = sgn * TMath::Abs(IPxySig);

  // IPxyz
  IPxyz = TMath::Sqrt(track.dcaXY() * track.dcaXY() + track.dcaZ() * track.dcaZ());
  IPxyzSig = IPxyz;
  SgnIPxyz = sgn * TMath::Abs(IPxyz);
  float dFdxy = 2 * track.dcaXY() / IPxyz;
  float dFdz = 2 * track.dcaZ() / IPxyz;
  IPxyzSig /= TMath::Sqrt(track.cYY() * dFdxy * dFdxy + track.cZZ() * dFdz * dFdz + 2 * track.cZY() * dFdxy * dFdz);
  SgnIPxyzSig = sgn * TMath::Abs(IPxyzSig);
}

}; // namespace JetTaggingUtilities

#endif // PWGJE_CORE_JETTAGGINGUTILITIES_H_

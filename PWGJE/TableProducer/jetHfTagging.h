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

// jet hf tagging header file
//
// Authors: Hanseo Park

#ifndef PWGJE_TABLEPRODUCER_JETHFTAGGING_H_
#define PWGJE_TABLEPRODUCER_JETHFTAGGING_H_

#include <array>
#include <vector>
#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;

//using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels>;

// functions for tagging jets
template <typename T, typename U, typename V>
void SetSignImpactParameterSignificance(T const& collision, U const& jet, V const& track, double& sign, double& IPxy, double& SignIPxy, double& IPxySig, double& SignIPxySig, double& IPxyz, double& SignIPxyz, double& IPxyzSig, double& SignIPxyzSig)
{
  sign = TMath::Sign(1, (track.x() - collision.posX()) * jet.px() + (track.y() - collision.posY()) * jet.py() + (track.z() - collision.posZ()) * jet.pz());
  // IPxy
  IPxy = track.dcaXY();
  SignIPxy = sign * TMath::Abs(IPxy);
  IPxySig = IPxy / TMath::Sqrt(track.sigmaDcaXY2());
  SignIPxySig = sign * TMath::Abs(IPxySig);

  // IPxyz
  IPxyz = TMath::Sqrt(track.dcaXY() * track.dcaXY() + track.dcaZ() * track.dcaZ());
  SignIPxyz = sign * TMath::Abs(IPxyz);
  double dFdxy = 2*track.dcaXY()/IPxyz;
  double dFdz = 2*track.dcaZ()/IPxyz;
  IPxyzSig /= TMath::Sqrt(track.cYY()*dFdxy*dFdxy+track.cZZ()*dFdz*dFdz+2*track.cZY()*dFdxy*dFdz);
  //double det = track.cYY() * track.cZZ() - track.cZY() * track.cZY();
  //double sigmaDcaXYZ2 = (track.dcaXY()*track.dcaXY()*track.cZZ() + track.dcaZ()*track.dcaZ()*track.cYY() - 2*track.dcaXY()*track.dcaZ()*track.cZY())/det;
  //IPxyzSig =  IPxyz / TMath::Sqrt(TMath::Abs(sigmaDcaXYZ2));
  SignIPxyzSig = sign * TMath::Abs(IPxyzSig);
}

template <typename T>
int mcdHfJetTagging(T const& jets)
{
  for (auto &jet : jets){}

}

#endif // PWGJE_TABLEPRODUCER_JETHFTAGGING_H_


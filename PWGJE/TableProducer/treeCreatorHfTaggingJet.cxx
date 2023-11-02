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

// \file treeCreatorHfTaggingJet.cxx

// tree creator for HF trgging jet 
//
// Authors: hanseo.park@cern.ch

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

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

#include "EventFiltering/filterTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace hfjet
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace full

DECLARE_SOA_TABLE(hfjetTrack, "AOD", "HFJETTRACK",
                  hfjet::Sign,
                  hfjet::ImpactParameterXY)

}

struct treeCreatorHfTaggingJet {

  Produces<aod::jet> train;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {

  }

  template <typename T>
  void fillJet(const T& collision)
  {

  }
  
  template <typename T>
  void fillJetTrack(const T& collision)
  {

  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(treeCreatorHfTaggingJet, processDummy, "Dummy process function turned on by default", true);

  void processMC(aod::Collisions const& collisions, aod::Tracks const& track)
  {
    fillJet(collisions);
    fillJetTrack(collisions);
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<treeCreatorHfTaggingJet>(cfgc, TaskName{"tree-creator-hf-tagging-jet"})}; }

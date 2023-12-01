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

// Task to produce a table joinable to the jet tables for hf jet tagging
//
/// copy from Common/TableProducer/trackPropagation.cxx on 23.Nov.23
/// \author Hanseo Park <hanseo.park@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "DetectorsBase/Propagator.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetTaggerHFTrackTask {

  Produces<aod::JTracksTagExt> jTracksTagExtTable;
  Produces<aod::JTrackTagExtPIs> jTracksParentIndexTable;
  Produces<aod::JTracksIUTagExt> jTracksIUTagExtTable;
  Produces<aod::JTrackIUTagExtPIs> jTracksIUParentIndexTable;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  const o2::dataformats::MeanVertexObject* mVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current() << " A for run " << bc.runNumber() << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }


  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTrackTask, processDummy, "Dummy process", true);

  //void processTracks(aod::Collisions const&, soa::Join<aod::Tracks, aod::TracksDCA> const& tracks, aod::BCsWithTimestamps const& bcs)
  void processTracks(aod::Collisions const&, aod::Tracks const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());
    gpu::gpustd::array<float, 2> dcaInfo;

    for (auto& track : tracks) {
      dcaInfo[0] = 999;
      dcaInfo[1] = 999;
      //aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      auto trackPar = getTrackPar(track);
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        } else {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({mVtx->getX(), mVtx->getY(), mVtx->getZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        }
        //trackType = aod::track::Track;
      }
      auto xyz = trackPar.getXYZGlo();
      float dcaX =999.;
      float dcaY= 999.;
      if (track.has_collision()){
        auto const& collision = track.collision();
        dcaX = xyz.X()-collision.posX();
        dcaY = xyz.Y()-collision.posY();
      } else {
        dcaX = xyz.X()-mVtx->getX();
        dcaY = xyz.Y()-mVtx->getY();
      }
      jTracksTagExtTable(track.collisionId(), dcaX, dcaY);
      jTracksParentIndexTable(track.globalIndex());
    }
  }
  PROCESS_SWITCH(JetTaggerHFTrackTask, processTracks, "Fill tagging decision for mcd jets", false);


  void processTracksIU(aod::Collisions const&, aod::StoredTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());
    gpu::gpustd::array<float, 2> dcaInfo;

    for (auto& track : tracks) {
      dcaInfo[0] = 999;
      dcaInfo[1] = 999;
      //aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      auto trackPar = getTrackPar(track);
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        } else {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({mVtx->getX(), mVtx->getY(), mVtx->getZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        }
        //trackType = aod::track::Track;
      }
      auto xyz = trackPar.getXYZGlo();
      float dcaX =999.;
      float dcaY= 999.;
      if (track.has_collision()){
        auto const& collision = track.collision();
        dcaX = xyz.X()-collision.posX();
        dcaY = xyz.Y()-collision.posY();
      } else {
        dcaX = xyz.X()-mVtx->getX();
        dcaY = xyz.Y()-mVtx->getY();
      }
      jTracksIUTagExtTable(track.collisionId(), dcaX, dcaY);
      jTracksIUParentIndexTable(track.globalIndex());
      //std::cout << "dca x: " << dcaX << " dca y: " << dcaY << "\n";
      //std::cout << "dca xy: " << track.dcaXY() << " dca xy from par: " << dcaInfo[0] << "\n";
      //std::cout << " dca xy from par: " << dcaInfo[0] << "\n";
    }
  }
  PROCESS_SWITCH(JetTaggerHFTrackTask, processTracksIU, "Fill tagging decision for mcd jets", false);
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTaggerHFTrackTask>(cfgc, TaskName{"jet-taggerhf-track"})};
}


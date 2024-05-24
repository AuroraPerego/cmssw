#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HGCalReco/interface/MtdHostCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

using namespace edm;

class MTDSoAProducer : public edm::stream::EDProducer<> {
public:
  MTDSoAProducer(const ParameterSet& pset);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event& ev, const edm::EventSetup& es) final;

private:
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> betaToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> MVAQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> momentumWithMTDToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> timePiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> timeKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> timePToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmaTimePiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmaTimeKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmaTimePToken_;
};

MTDSoAProducer::MTDSoAProducer(const ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
      trackAssocToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"))),
      t0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"))),
      sigmat0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"))),
      tmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtdSrc"))),
      sigmatmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtdSrc"))),
      betaToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("betamtd"))),
      pathToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathmtd"))),
      MVAQualityToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mvaquality"))),
      momentumWithMTDToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("momentum"))),
      probPiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"))),
      probKToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"))),
      probPToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"))),
      btlMatchChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2"))),
      btlMatchTimeChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2"))),
      etlMatchChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchChi2"))),
      etlMatchTimeChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2"))),
      timePiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("generalTrackTofPi"))),
      timeKToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("generalTrackTofK"))),
      timePToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("generalTrackTofP"))),
      sigmaTimePiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("generalTrackSigmaTofPi"))),
      sigmaTimeKToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("generalTrackSigmaTofK"))),
      sigmaTimePToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("generalTrackSigmaTofP"))) {
  produces<MtdHostCollection>();
}

// Configuration descriptions
void MTDSoAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("betamtd", edm::InputTag("trackExtenderWithMTD:generalTrackBeta"));
  desc.add<edm::InputTag>("pathmtd", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("mvaquality", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("posmtd", edm::InputTag("trackExtenderWithMTD:generalTrackmtdpos"));
  desc.add<edm::InputTag>("momentum", edm::InputTag("trackExtenderWithMTD:generalTrackp"));
  desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));
  desc.add<edm::InputTag>("btlMatchChi2", edm::InputTag("trackExtenderWithMTD:btlMatchChi2"));
  desc.add<edm::InputTag>("btlMatchTimeChi2", edm::InputTag("trackExtenderWithMTD:btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchChi2", edm::InputTag("trackExtenderWithMTD:etlMatchChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2", edm::InputTag("trackExtenderWithMTD:etlMatchTimeChi2"));
  desc.add<edm::InputTag>("generalTrackTofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("generalTrackTofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("generalTrackTofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("generalTrackSigmaTofPi", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("generalTrackSigmaTofK", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("generalTrackSigmaTofP", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));

  descriptions.add("mtdSoAProducer", desc);
}

void MTDSoAProducer::produce(edm::Event& ev, const edm::EventSetup& es) {
  edm::Handle<reco::TrackCollection> tracksH;
  ev.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  const auto& trackAssoc = ev.get(trackAssocToken_);

  const auto& t0 = ev.get(t0Token_);
  const auto& sigmat0 = ev.get(sigmat0Token_);

  const auto& tmtd = ev.get(tmtdToken_);
  const auto& sigmatmtd = ev.get(sigmatmtdToken_);

  const auto& beta = ev.get(betaToken_);
  const auto& path = ev.get(pathToken_);
  const auto& MVAquality = ev.get(MVAQualityToken_);
  const auto& momentum = ev.get(momentumWithMTDToken_);
  const auto& probPi = ev.get(probPiToken_);
  const auto& probK = ev.get(probKToken_);
  const auto& probP = ev.get(probPToken_);

  const auto& btlMatchChi2 = ev.get(btlMatchChi2Token_);
  const auto& btlMatchTimeChi2 = ev.get(btlMatchTimeChi2Token_);
  const auto& etlMatchChi2 = ev.get(etlMatchChi2Token_);
  const auto& etlMatchTimeChi2 = ev.get(etlMatchTimeChi2Token_);
  const auto& timePi = ev.get(timePiToken_);
  const auto& timeK = ev.get(timeKToken_);
  const auto& timeP = ev.get(timePToken_);
  const auto& sigmaTimePi = ev.get(sigmaTimePiToken_);
  const auto& sigmaTimeK = ev.get(sigmaTimeKToken_);
  const auto& sigmaTimeP = ev.get(sigmaTimePToken_);

  auto MtdInfo = std::make_unique<MtdHostCollection>(tracks.size(), cms::alpakatools::host());

  auto& MtdInfoView = MtdInfo->view();
  for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack) {
    const reco::TrackRef trackref(tracksH, iTrack);

    if (trackAssoc[trackref] == -1) {
      MtdInfoView.trackAsocMTD()[iTrack] = -1;
      MtdInfoView.time0()[iTrack] = 0.f;
      MtdInfoView.time0Err()[iTrack] = -1.f;
      MtdInfoView.time()[iTrack] = 0.f;
      MtdInfoView.timeErr()[iTrack] = -1.f;
      MtdInfoView.MVAquality()[iTrack] = 0.f;
      MtdInfoView.pathLength()[iTrack] = 0.f;
      MtdInfoView.beta()[iTrack] = 0.f;
      MtdInfoView.momentumWithMTD()[iTrack] = 0.f;
      MtdInfoView.probPi()[iTrack] = 0.f;
      MtdInfoView.probK()[iTrack] = 0.f;
      MtdInfoView.probP()[iTrack] = 0.f;
      MtdInfoView.btlMatchChi2()[iTrack] = -1.f;
      MtdInfoView.btlMatchTimeChi2()[iTrack] = -1.f;
      MtdInfoView.etlMatchChi2()[iTrack] = -1.f;
      MtdInfoView.etlMatchTimeChi2()[iTrack] = -1.f;
      MtdInfoView.timePi()[iTrack] = 0.f;
      MtdInfoView.timeK()[iTrack] = 0.f;
      MtdInfoView.timeP()[iTrack] = 0.f;
      MtdInfoView.sigmaTimePi()[iTrack] = -1.f;
      MtdInfoView.sigmaTimeK()[iTrack] = -1.f;
      MtdInfoView.sigmaTimeP()[iTrack] = -1.f;
      continue;
    }

    MtdInfoView.trackAsocMTD()[iTrack] = trackAssoc[trackref];
    MtdInfoView.time0()[iTrack] = t0[trackref];
    MtdInfoView.time0Err()[iTrack] = sigmat0[trackref];
    MtdInfoView.time()[iTrack] = tmtd[trackref];
    MtdInfoView.timeErr()[iTrack] = sigmatmtd[trackref];
    MtdInfoView.MVAquality()[iTrack] = MVAquality[trackref];
    MtdInfoView.pathLength()[iTrack] = path[trackref];
    MtdInfoView.beta()[iTrack] = beta[trackref];
    MtdInfoView.momentumWithMTD()[iTrack] = momentum[trackref];
    MtdInfoView.probPi()[iTrack] = probPi[trackref];
    MtdInfoView.probK()[iTrack] = probK[trackref];
    MtdInfoView.probP()[iTrack] = probP[trackref];
    MtdInfoView.btlMatchChi2()[iTrack] = btlMatchChi2[trackref];
    MtdInfoView.btlMatchTimeChi2()[iTrack] = btlMatchTimeChi2[trackref];
    MtdInfoView.etlMatchChi2()[iTrack] = etlMatchChi2[trackref];
    MtdInfoView.etlMatchTimeChi2()[iTrack] = etlMatchTimeChi2[trackref];
    MtdInfoView.timePi()[iTrack] = timePi[trackref];
    MtdInfoView.timeK()[iTrack] = timeK[trackref];
    MtdInfoView.timeP()[iTrack] = timeP[trackref];
    MtdInfoView.sigmaTimePi()[iTrack] = sigmaTimePi[trackref];
    MtdInfoView.sigmaTimeK()[iTrack] = sigmaTimeK[trackref];
    MtdInfoView.sigmaTimeP()[iTrack] = sigmaTimeP[trackref];
  }

  ev.put(std::move(MtdInfo));
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MTDSoAProducer);

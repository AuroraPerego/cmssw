#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>  // unique_ptr
#include <string>
#include <tuple>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Associations/interface/TrackToGenParticleAssociator.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // Adding header files for electrons
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"  // Adding header files for electrons
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class MtdEleIsoNtupler : public edm::one::EDAnalyzer<
                             edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit MtdEleIsoNtupler(const edm::ParameterSet&);
  ~MtdEleIsoNtupler() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef math::XYZVector Vector;
  typedef std::vector<double> Vec;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  // ------------ member data ------------

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<reco::GsfElectronCollection>
      GsfElectronToken_EB_;  // Adding token for electron collection
  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_EE_;
  edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;

  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;

  // cuts
  const auto min_pt_cut_EB = 0.7;  // GeV
  const auto min_pt_cut_EE = 0.4;  // GeV
 
  // Output tree
  TTree* tree_;

  void clearVariables();

  unsigned int event_index;

  // Variables for branches
  unsigned int ev_event_;

  // electrons
  unsigned int nEle_;
  unsigned int nPromptEle_;
  std::vector<float> ele_energy_;
  std::vector<float> ele_pt_;
  std::vector<float> ele_sim_pt_;
  std::vector<float> ele_eta_;
  std::vector<float> ele_phi_;
  std::vector<float> ele_dz_;
  std::vector<float> ele_dxy_;
  std::vector<int> ele_TrkRef_;
  std::vector<bool> ele_barrel_;
  std::vector<bool> ele_Promt_;
  std::vector<double> ele_time_;
  std::vector<double> ele_timeErr_;
  std::vector<double> ele_mva_;
  std::vector<double> ele_sim_time_;
  std::vector<float> sum_pT_;
  std::vector<int> nTracks_;

  // tracks
  std::vector<std::vector<float>> track_pt_;
  std::vector<std::vector<double>> track_sim_pt_;
  std::vector<std::vector<float>> track_dtEle_;
  std::vector<std::vector<float>> track_dtVtx_;
  std::vector<std::vector<float>> track_dzEle_;
  std::vector<std::vector<float>> track_weight_;
  std::vector<std::vector<double>> track_time_;
  std::vector<std::vector<double>> track_timeErr_;
  std::vector<std::vector<double>> track_mva_;
  std::vector<std::vector<double>> track_sim_time_;
  std::vector<std::vector<int>> track_genMatched_;

  // vertex
  float vertex_time_;
  float vertex_timeErr_;
  int number_of_tracks_;

  TTree* electrons_tree_;
  TTree* tracks_tree_;
  TTree* vertices_tree_;
};

void MtdEleIsoNtupler::clearVariables() {
  ev_event_ = 0;
  nPromptEle_ = 0;
  nEle_ = 0;

  ele_energy_.clear();
  ele_pt_.clear();
  ele_sim_pt_.clear();
  ele_eta_.clear();
  ele_phi_.clear();
  ele_dz_.clear();
  ele_dxy_.clear();
  ele_TrkRef_.clear();
  ele_barrel_.clear();
  ele_Promt_.clear();
  ele_time_.clear();
  ele_timeErr_.clear();
  ele_mva_.clear();
  ele_sim_time_.clear();
  sum_pT_.clear();
  nTracks_.clear();

  track_pt_.clear();
  track_sim_pt_.clear();
  track_dtEle_.clear();
  track_weight_.clear();
  track_dtVtx_.clear();
  track_dzEle_.clear();
  track_time_.clear();
  track_timeErr_.clear();
  track_mva_.clear();
  track_sim_time_.clear();
  track_genMatched_.clear();

  vertex_time_ = 0;
  vertex_timeErr_ = 0;
  number_of_tracks_ = 0;
};

// ------------ constructor and destructor --------------
MtdEleIsoNtupler::MtdEleIsoNtupler(const edm::ParameterSet& iConfig) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(
      iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(
      iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ =
      consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>(
          "inputTag_vtx"));  // Vtx 4D collection

  GsfElectronToken_EB_ =
      consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>(
          "inputEle_EB"));  // Barrel electron collection input/token
  GsfElectronToken_EE_ =
      consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>(
          "inputEle_EE"));  // Endcap electron collection input/token
  GenParticleToken_ = consumes<reco::GenParticleCollection>(
      iConfig.getParameter<edm::InputTag>("inputGenP"));

  HepMCProductToken_ = consumes<edm::HepMCProduct>(
      iConfig.getParameter<edm::InputTag>("inputTagH"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(
      iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(
      iConfig.getParameter<edm::InputTag>("trackMVAQual"));

  trackingParticleCollectionToken_ = consumes<TrackingParticleCollection>(
      iConfig.getParameter<edm::InputTag>("SimTag"));
  recoToSimAssociationToken_ = consumes<reco::RecoToSimCollection>(
      iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
};

MtdEleIsoNtupler::~MtdEleIsoNtupler() {
  clearVariables();
};

void MtdEleIsoNtupler::beginRun(edm::Run const&, edm::EventSetup const& es) {}

void MtdEleIsoNtupler::beginJob() {
  edm::Service<TFileService> fs;
  electrons_tree_ = fs->make<TTree>("electrons", "electrons info");
  tracks_tree_ = fs->make<TTree>("tracks", "general tracks info");
  vertices_tree_ = fs->make<TTree>("vertices", "vertex info");

  electrons_tree_->Branch("nEle", &nEle_);
  electrons_tree_->Branch("ele_energy", &ele_energy_);
  electrons_tree_->Branch("ele_pt", &ele_pt_);
  electrons_tree_->Branch("ele_sim_pt", &ele_sim_pt_);
  electrons_tree_->Branch("ele_eta", &ele_eta_);
  electrons_tree_->Branch("ele_phi", &ele_phi_);
  electrons_tree_->Branch("ele_dz", &ele_dz_);
  electrons_tree_->Branch("ele_dxy", &ele_dxy_);
  electrons_tree_->Branch("ele_track", &ele_TrkRef_);
  electrons_tree_->Branch("ele_barrel", &ele_barrel_);
  electrons_tree_->Branch("ele_prompt", &ele_Promt_);
  electrons_tree_->Branch("ele_time", &ele_time_);
  electrons_tree_->Branch("ele_timeErr", &ele_timeErr_);
  electrons_tree_->Branch("ele_mva", &ele_mva_);
  electrons_tree_->Branch("ele_sim_time", &ele_sim_time_);
  electrons_tree_->Branch("sum_pT", &sum_pT_);
  electrons_tree_->Branch("nTracks", &nTracks_);

  tracks_tree_->Branch("track_pt", &track_pt_);
  tracks_tree_->Branch("track_sim_pt", &track_sim_pt_);
  tracks_tree_->Branch("track_dt_ele", &track_dtEle_);
  tracks_tree_->Branch("track_PVweight", &track_weight_);
  tracks_tree_->Branch("track_dt_vtx", &track_dtVtx_);
  tracks_tree_->Branch("track_dz_ele", &track_dzEle_);
  tracks_tree_->Branch("track_time", &track_time_);
  tracks_tree_->Branch("track_timeErr", &track_timeErr_);
  tracks_tree_->Branch("track_mva", &track_mva_);
  tracks_tree_->Branch("track_sim_time", &track_sim_time_);
  tracks_tree_->Branch("track_gen_matched", &track_genMatched_);

  vertices_tree_->Branch("vertex_time", &vertex_time_);
  vertices_tree_->Branch("vertex_timeErr", &vertex_timeErr_);
  vertices_tree_->Branch("vertex_NumberofTracks", &number_of_tracks_);
  vertices_tree_->Branch("event", &ev_event_);

  event_index = 0;
}

// ------------ method called for each event  ------------
void MtdEleIsoNtupler::analyze(const edm::Event& iEvent,
                               const edm::EventSetup& iSetup) {
  event_index++;
  clearVariables();
  ev_event_ = event_index;

  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = iEvent.getHandle(GenRecTrackToken_);
  auto VertexHandle = iEvent.getHandle(RecVertexToken_);
  std::vector<reco::Vertex> vertices = *VertexHandle;

  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  // const auto& tMtd = iEvent.get(tmtdToken_);
  // const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  // const auto& t0Src = iEvent.get(t0SrcToken_);
  // const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  // const auto& t0Safe = iEvent.get(t0SafePidToken_);
  // const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);

  auto eleHandle_EB = makeValid(iEvent.getHandle(GsfElectronToken_EB_));
  reco::GsfElectronCollection eleColl_EB = *(eleHandle_EB.product());

  auto eleHandle_EE = makeValid(iEvent.getHandle(GsfElectronToken_EE_));
  reco::GsfElectronCollection eleColl_EE = *(eleHandle_EE.product());

  edm::Handle<reco::GenParticleCollection> GenPartHandle;
  iEvent.getByToken(GenParticleToken_, GenPartHandle);
  const reco::GenParticleCollection GenPartColl = *(GenPartHandle.product());

  auto tpHandle = makeValid(iEvent.getHandle(trackingParticleCollectionToken_));
  TrackingParticleCollection tpColl = *(tpHandle.product());

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  const reco::RecoToSimCollection* r2s_ = recoToSimH.product();

//  std::cout << "genParticles list:" << std::endl;
//  int nelet = 0;
//  for (const auto& gp : GenPartColl){
//    if (std::abs(gp.pdgId())==11){
//	std::cout << "  pdg " << gp.pdgId() << " pt " << gp.pt() << " eta " << gp.eta() << " is prompt " << gp.isPromptFinalState() << " from tau " << gp.isDirectPromptTauDecayProductFinalState() << " decayed " << gp.isPromptDecayed() << " last copy " << gp.isLastCopy() << " status " << gp.status() ;
//	if (gp.numberOfMothers() > 0)
//	std::cout << " mother " << gp.mother()->pdgId();
//        if (gp.numberOfDaughters() > 0) 
//	std::cout << " daughter " << gp.daughter(0)->pdgId();
//        std::cout << std::endl;
////    if (std::abs(gp.pdgId())==11){
//        nelet++;
//    }
// }
//std::cout << "there are " << GenPartColl.size() << " genParticles, of which " << nelet << " electrons" << std::endl;
//std::cout << "there are " << (eleColl_EB.size() + eleColl_EE.size()) << " reco electrons, " << GenRecTrackHandle->size() << " recotracks and " << tpColl.size() << " tracking particles" << std::endl; 

  // Creating combined electron collection
  std::vector<reco::GsfElectron> localEleCollection;
  localEleCollection.reserve(eleColl_EB.size() + eleColl_EE.size());
  for (const auto& ele_EB : eleColl_EB) {
    if (ele_EB.isEB()) {
      localEleCollection.emplace_back(ele_EB);
    }
  }
  for (const auto& ele_EE : eleColl_EE) {
    if (ele_EE.isEE()) {
      localEleCollection.emplace_back(ele_EE);
    }
  }
  localEleCollection.shrink_to_fit();
  nEle_ = localEleCollection.size();
  // Selecting the PV from 3D and 4D vertex collections
  reco::Vertex Vtx_chosen;

  // This part has to be included, because in ~1% of the events, the "good"
  // vertex is the 1st one not the 0th one in the collection
  for (int iVtx = 0; iVtx < (int)vertices.size(); iVtx++) {
    const reco::Vertex& vertex = vertices.at(iVtx);

    if (!vertex.isFake() && vertex.ndof() >= 4) {
      Vtx_chosen = vertex;
      break;
    }
  }
  vertex_time_ = Vtx_chosen.t();
  vertex_timeErr_ = Vtx_chosen.tError();
  number_of_tracks_ = Vtx_chosen.tracksSize();
  // Vertex selection ends

  auto isPrompt = [](int pdg) {
    pdg = std::abs(pdg);
    return (pdg == 23 or pdg == 24 or pdg==15 or pdg==11);
  };

  for (const auto& ele : localEleCollection) {
    bool ele_Promt = false;

    float ele_track_source_dz = fabs(ele.gsfTrack()->dz(Vtx_chosen.position()));
    float ele_track_source_dxy =
        fabs(ele.gsfTrack()->dxy(Vtx_chosen.position()));

    if (ele.pt() < 10 || fabs(ele.eta()) > 2.8)
      continue;

    const reco::TrackRef ele_TrkRef = ele.core()->ctfTrack();
    const reco::TrackBaseRef trkrefb(ele_TrkRef);
    
    double tsim = -1.;
    double ele_sim_pt = -1.;
    auto found = r2s_->find(trkrefb);
    if (found != r2s_->end()) {
      const auto& tp = (found->val)[0];
      tsim = (tp.first)->parentVertex()->position().t() * 1e9;
      ele_sim_pt = (tp.first)->pt();
      // check that the genParticle vector is not empty
      if (((found->val)[0]).first->status() != -99) {
        const auto genParticle = *(tp.first->genParticles()[0]);
        // check if prompt (not from hadron, muon, or tau decay) and final state
        // or if is a direct decay product of a prompt tau and is final state
        if ((genParticle.isPromptFinalState() or genParticle.isDirectPromptTauDecayProductFinalState()) and isPrompt(genParticle.mother()->pdgId())) {
          ele_Promt = true;
        }
      }
    }
    std::cout << std::endl;
    // const reco::GsfTrackRef ele_gsfTrkRef = ele.gsfTrack();

    math::XYZVector EleSigTrackMomentumAtVtx = ele.gsfTrack()->momentum();
    double EleSigTrackEtaAtVtx = ele.gsfTrack()->eta();

    // if we found a track match, we add MTD timing information for this matched
    // track and do the isolation check
    if (ele_TrkRef.isNonnull()) {
      ele_dz_.push_back(ele_track_source_dz);
      ele_dxy_.push_back(ele_track_source_dxy);
      ele_energy_.push_back(ele.caloEnergy());
      ele_Promt_.push_back(ele_Promt);
      ele_TrkRef_.push_back(ele_TrkRef.key());

      // for track pT/dz cuts (Different for EB and EE in TDR)
      bool Barrel_ele = ele.isEB();  
      ele_barrel_.push_back(Barrel_ele);
      float min_pt_cut = Barrel_ele ? min_pt_cut_EB : min_pt_cut_EE;

      double ele_sigTrkTime = t0Pid[ele_TrkRef]; 
      double ele_sigTrkTimeErr = Sigmat0Pid[ele_TrkRef];
      double ele_sigTrkMtdMva = mtdQualMVA[ele_TrkRef];

      ele_time_.push_back(ele_sigTrkTime);
      ele_timeErr_.push_back(ele_sigTrkTimeErr);
      ele_mva_.push_back(ele_sigTrkMtdMva);
      ele_sim_time_.push_back(tsim);
      ele_pt_.push_back(ele.pt());
      ele_sim_pt_.push_back(ele_sim_pt);
      ele_eta_.push_back(ele.eta());
      ele_phi_.push_back(ele.phi());

      int N_tracks_noMTD = 0;  // values for no MTD case
      double pT_sum_noMTD = 0;

      std::vector<float> track_pt;
      std::vector<double> track_sim_pt;
      std::vector<float> track_dtEle;
      std::vector<float> track_dtVtx;
      std::vector<float> track_dzEle;
      std::vector<float> track_weight;
      std::vector<double> track_time;
      std::vector<double> track_timeErr;
      std::vector<double> track_mva;
      std::vector<double> track_sim_time;
      std::vector<int> track_genMatched;

      int general_index = 0;
      for (const auto& trackGen : *GenRecTrackHandle) {
        const reco::TrackRef trackref_general(GenRecTrackHandle, general_index);
        general_index++;

        if (trackref_general == ele_TrkRef)
          continue;

        if (trackGen.pt() < min_pt_cut) {  // track pT cut
          continue;
        }

        if (fabs(trackGen.vz() - ele.gsfTrack()->vz()) >
            1) {  // general track vs signal track dz cut
          continue;
        }

        double dr_check =
            reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx);
        double deta = fabs(trackGen.eta() - EleSigTrackEtaAtVtx);

        const auto min_dR_cut = 0.01;
        const auto max_dR_cut = 0.3;
        const auto min_strip_cut = 0.01;
        if (dr_check < min_dR_cut || dr_check > max_dR_cut ||
            deta < min_strip_cut) {
          continue;
        }

        const reco::TrackBaseRef trkrefBase(trackref_general);

        auto TPmatched = r2s_->find(trkrefBase);
        double tsim_trk = -1.;
        double trk_ptSim = -1.;
	int genMatched = 0;
        if (TPmatched != r2s_->end()) {
          // reco track matched to a TP
          const auto& tp = (TPmatched->val)[0];
          tsim_trk = (tp.first)->parentVertex()->position().t() * 1e9;
          trk_ptSim = (tp.first)->pt();
          // check that the genParticle vector is not empty
            if (((TPmatched->val)[0]).first->status() != -99) {
		genMatched = 1;
	    }
        }
	track_genMatched.push_back(genMatched);
        track_dzEle.push_back(fabs(trackGen.vz() - ele.gsfTrack()->vz()));
        track_weight.push_back(Vtx_chosen.trackWeight(trackref_general));

        track_pt.push_back(trackGen.pt());
        track_sim_pt.push_back(trk_ptSim);
        ++N_tracks_noMTD;
        pT_sum_noMTD += trackGen.pt();

        // checking the MTD timing cuts
        double TrkMTDTime = t0Pid[trackref_general];
        double TrkMTDTimeErr = Sigmat0Pid[trackref_general];
        double TrkMTDMva = mtdQualMVA[trackref_general];

        track_time.push_back(TrkMTDTime);
        track_timeErr.push_back(TrkMTDTimeErr);
        track_mva.push_back(TrkMTDMva);
        track_sim_time.push_back(tsim_trk);

        track_dtEle.push_back(fabs(TrkMTDTime - ele_sigTrkTime));
        track_dtVtx.push_back(fabs(TrkMTDTime - Vtx_chosen.t()));
      }
      nTracks_.push_back(N_tracks_noMTD);
      sum_pT_.push_back(pT_sum_noMTD);

      track_pt_.push_back(track_pt);
      track_sim_pt_.push_back(track_sim_pt);
      track_dtEle_.push_back(track_dtEle);
      track_dtVtx_.push_back(track_dtVtx);
      track_dzEle_.push_back(track_dzEle);
      track_genMatched_.push_back(track_genMatched);
      track_weight_.push_back(track_weight);
      track_time_.push_back(track_time);
      track_timeErr_.push_back(track_timeErr);
      track_mva_.push_back(track_mva);
      track_sim_time_.push_back(track_sim_time);

    }  // electron matched to a track
  }    // electron collection inside single event

  electrons_tree_->Fill();
  tracks_tree_->Fill();
  vertices_tree_->Fill();
}

void MtdEleIsoNtupler::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MtdEleIsoNtupler::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>(
      "inputTag_vtx",
      edm::InputTag(
          "offlinePrimaryVertices4D"));  //  "offlinePrimaryVertices4D" or
                                         // "offlinePrimaryVertices" (3D case)
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));

  desc.add<edm::InputTag>("inputEle_EB",
                          edm::InputTag("gedGsfElectrons"));  // barrel ecal and
                                                              // track driven
                                                              // electrons
  desc.add<edm::InputTag>(
      "inputEle_EE",
      edm::InputTag("ecalDrivenGsfElectronsHGC"));  // only endcap electrons
  desc.add<edm::InputTag>("inputGenP", edm::InputTag("genParticles"));

  desc.add<edm::InputTag>(
      "tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>(
      "sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src",
                          edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>(
      "sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>(
           "trackAssocSrc",
           edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>(
      "pathLengthSrc",
      edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID",
                          edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("trackMVAQual",
                          edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>(
      "TPtoRecoTrackAssoc",
      edm::InputTag("trackingParticleRecoTrackAsssociation"));
  // desc.add<edm::InputTag>("gsfTracks", edm::InputTag("electronGsfTracks"));

  desc.add<double>("trackMinimumPt", 1.0);  // [GeV]

  descriptions.add("mtdEleIsoNtupler", desc);
}

DEFINE_FWK_MODULE(MtdEleIsoNtupler);

#include <string>
#include <tuple>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // Adding header files for electrons
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"         // Adding header files for electrons
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class MtdEleIsoAnalyzer : public DQMEDAnalyzer {
public:
  explicit MtdEleIsoAnalyzer(const edm::ParameterSet&);
  ~MtdEleIsoAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const bool mvaGenSel(const HepMC::GenParticle&, const float&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);
  const bool mvaGenRecMatch(const HepMC::GenParticle&, const double&, const reco::TrackBase&);

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMinEta_;
  const float trackMaxEta_;
  const double rel_iso_cut_;
  // const std::vector<double> pTcuts_;

  bool track_match_PV_;
  bool dt_sig_vtx_;
  bool dt_sig_track_;

  static constexpr float min_dR_cut = 0.01;
  static constexpr float max_dR_cut = 0.3;
  static constexpr float min_pt_cut_EB = 0.7;  // 0.4 for endcap
  static constexpr float min_pt_cut_EE = 0.4;
  static constexpr float max_dz_cut_EB =
      0.2;  // will check 0.15/0.20/0.25/0.30 for tracks in iso cone // 0.1 for barrel, 0.2 for endcap
  static constexpr float max_dz_cut_EE = 0.2;
  static constexpr float max_dz_vtx_cut = 0.5;
  static constexpr float max_dxy_vtx_cut = 0.2;
  const std::vector<double> max_dt_vtx_cut{
      0.30,
      0.27,
      0.24,
      0.21,
      0.18,
      0.15,
      0.12};  // has to be 7 values here!!! Code created for 7 dt values!! Change lists and histogram count if number of dt values change.
  const std::vector<double> max_dt_track_cut{0.30, 0.27, 0.24, 0.21, 0.18, 0.15, 0.12};
  const int number_of_cuts = max_dt_track_cut.size();
  static constexpr float min_strip_cut = 0.01;
  static constexpr float min_track_mtd_mva_cut = 0.5;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_EB_;  // Adding token for electron collection
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

  // edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  // Signal histograms
  MonitorElement* meEleISO_Ntracks_Sig_EB_;  // Adding histograms for barrel electrons (isolation stuff)
  MonitorElement* meEleISO_chIso_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_1_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_1_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_2_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_2_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_3_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_3_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_4_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_4_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_5_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_5_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_6_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_6_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_7_Sig_EB_;
  MonitorElement* meEleISO_chIso_MTD_7_Sig_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_Sig_EB_;

  MonitorElement* meEle_pt_tot_Sig_EB_;
  MonitorElement* meEle_eta_tot_Sig_EB_;
  MonitorElement* meEle_phi_tot_Sig_EB_;

  MonitorElement* meEle_pt_MTD_1_Sig_EB_;
  MonitorElement* meEle_eta_MTD_1_Sig_EB_;
  MonitorElement* meEle_phi_MTD_1_Sig_EB_;

  MonitorElement* meEle_pt_MTD_2_Sig_EB_;
  MonitorElement* meEle_eta_MTD_2_Sig_EB_;
  MonitorElement* meEle_phi_MTD_2_Sig_EB_;

  MonitorElement* meEle_pt_MTD_3_Sig_EB_;
  MonitorElement* meEle_eta_MTD_3_Sig_EB_;
  MonitorElement* meEle_phi_MTD_3_Sig_EB_;

  MonitorElement* meEle_pt_MTD_4_Sig_EB_;
  MonitorElement* meEle_eta_MTD_4_Sig_EB_;
  MonitorElement* meEle_phi_MTD_4_Sig_EB_;

  MonitorElement* meEle_pt_MTD_5_Sig_EB_;
  MonitorElement* meEle_eta_MTD_5_Sig_EB_;
  MonitorElement* meEle_phi_MTD_5_Sig_EB_;

  MonitorElement* meEle_pt_MTD_6_Sig_EB_;
  MonitorElement* meEle_eta_MTD_6_Sig_EB_;
  MonitorElement* meEle_phi_MTD_6_Sig_EB_;

  MonitorElement* meEle_pt_MTD_7_Sig_EB_;
  MonitorElement* meEle_eta_MTD_7_Sig_EB_;
  MonitorElement* meEle_phi_MTD_7_Sig_EB_;

  MonitorElement* meEle_pt_noMTD_Sig_EB_;
  MonitorElement* meEle_eta_noMTD_Sig_EB_;
  MonitorElement* meEle_phi_noMTD_Sig_EB_;

  MonitorElement* meEleISO_Ntracks_Sig_EE_;  // Adding histograms for endcap electrons (isolation stuff)
  MonitorElement* meEleISO_chIso_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_1_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_1_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_2_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_2_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_3_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_3_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_4_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_4_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_5_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_5_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_6_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_6_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_Sig_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_7_Sig_EE_;
  MonitorElement* meEleISO_chIso_MTD_7_Sig_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_Sig_EE_;

  MonitorElement* meEle_pt_tot_Sig_EE_;
  MonitorElement* meEle_eta_tot_Sig_EE_;
  MonitorElement* meEle_phi_tot_Sig_EE_;

  MonitorElement* meEle_pt_MTD_1_Sig_EE_;
  MonitorElement* meEle_eta_MTD_1_Sig_EE_;
  MonitorElement* meEle_phi_MTD_1_Sig_EE_;

  MonitorElement* meEle_pt_MTD_2_Sig_EE_;
  MonitorElement* meEle_eta_MTD_2_Sig_EE_;
  MonitorElement* meEle_phi_MTD_2_Sig_EE_;

  MonitorElement* meEle_pt_MTD_3_Sig_EE_;
  MonitorElement* meEle_eta_MTD_3_Sig_EE_;
  MonitorElement* meEle_phi_MTD_3_Sig_EE_;

  MonitorElement* meEle_pt_MTD_4_Sig_EE_;
  MonitorElement* meEle_eta_MTD_4_Sig_EE_;
  MonitorElement* meEle_phi_MTD_4_Sig_EE_;

  MonitorElement* meEle_pt_MTD_5_Sig_EE_;
  MonitorElement* meEle_eta_MTD_5_Sig_EE_;
  MonitorElement* meEle_phi_MTD_5_Sig_EE_;

  MonitorElement* meEle_pt_MTD_6_Sig_EE_;
  MonitorElement* meEle_eta_MTD_6_Sig_EE_;
  MonitorElement* meEle_phi_MTD_6_Sig_EE_;

  MonitorElement* meEle_pt_MTD_7_Sig_EE_;
  MonitorElement* meEle_eta_MTD_7_Sig_EE_;
  MonitorElement* meEle_phi_MTD_7_Sig_EE_;

  MonitorElement* meEle_pt_noMTD_Sig_EE_;
  MonitorElement* meEle_eta_noMTD_Sig_EE_;
  MonitorElement* meEle_phi_noMTD_Sig_EE_;
  // Signal histograms end

  // Background histograms
  MonitorElement* meEleISO_Ntracks_Bkg_EB_;  // Adding histograms for barrel electrons (isolation stuff)
  MonitorElement* meEleISO_chIso_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_1_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_1_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_2_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_2_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_3_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_3_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_4_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_4_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_5_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_5_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_6_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_6_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_MTD_7_Bkg_EB_;
  MonitorElement* meEleISO_chIso_MTD_7_Bkg_EB_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_Bkg_EB_;

  MonitorElement* meEle_pt_tot_Bkg_EB_;
  MonitorElement* meEle_eta_tot_Bkg_EB_;
  MonitorElement* meEle_phi_tot_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_1_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_1_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_1_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_2_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_2_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_2_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_3_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_3_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_3_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_4_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_4_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_4_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_5_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_5_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_5_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_6_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_6_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_6_Bkg_EB_;

  MonitorElement* meEle_pt_MTD_7_Bkg_EB_;
  MonitorElement* meEle_eta_MTD_7_Bkg_EB_;
  MonitorElement* meEle_phi_MTD_7_Bkg_EB_;

  MonitorElement* meEle_pt_noMTD_Bkg_EB_;
  MonitorElement* meEle_eta_noMTD_Bkg_EB_;
  MonitorElement* meEle_phi_noMTD_Bkg_EB_;

  MonitorElement* meEleISO_Ntracks_Bkg_EE_;  // Adding histograms for endcap electrons (isolation stuff)
  MonitorElement* meEleISO_chIso_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_1_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_1_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_1_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_2_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_2_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_2_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_3_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_3_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_3_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_4_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_4_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_4_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_5_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_5_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_5_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_6_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_6_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_6_Bkg_EE_;

  MonitorElement* meEleISO_Ntracks_MTD_7_Bkg_EE_;
  MonitorElement* meEleISO_chIso_MTD_7_Bkg_EE_;
  MonitorElement* meEleISO_rel_chIso_MTD_7_Bkg_EE_;

  MonitorElement* meEle_pt_tot_Bkg_EE_;
  MonitorElement* meEle_eta_tot_Bkg_EE_;
  MonitorElement* meEle_phi_tot_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_1_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_1_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_1_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_2_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_2_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_2_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_3_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_3_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_3_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_4_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_4_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_4_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_5_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_5_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_5_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_6_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_6_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_6_Bkg_EE_;

  MonitorElement* meEle_pt_MTD_7_Bkg_EE_;
  MonitorElement* meEle_eta_MTD_7_Bkg_EE_;
  MonitorElement* meEle_phi_MTD_7_Bkg_EE_;

  MonitorElement* meEle_pt_noMTD_Bkg_EE_;
  MonitorElement* meEle_eta_noMTD_Bkg_EE_;
  MonitorElement* meEle_phi_noMTD_Bkg_EE_;
  // Background histograms end

  // promt part for histogram vectors
  std::vector<MonitorElement*> Ntracks_EB_list_Sig;
  std::vector<MonitorElement*> ch_iso_EB_list_Sig;
  std::vector<MonitorElement*> rel_ch_iso_EB_list_Sig;

  std::vector<MonitorElement*> Ntracks_EE_list_Sig;
  std::vector<MonitorElement*> ch_iso_EE_list_Sig;
  std::vector<MonitorElement*> rel_ch_iso_EE_list_Sig;

  std::vector<MonitorElement*> Ele_pT_MTD_EB_list_Sig;
  std::vector<MonitorElement*> Ele_eta_MTD_EB_list_Sig;
  std::vector<MonitorElement*> Ele_phi_MTD_EB_list_Sig;

  std::vector<MonitorElement*> Ele_pT_MTD_EE_list_Sig;
  std::vector<MonitorElement*> Ele_eta_MTD_EE_list_Sig;
  std::vector<MonitorElement*> Ele_phi_MTD_EE_list_Sig;

  // Non-promt part for histogram vectors
  std::vector<MonitorElement*> Ntracks_EB_list_Bkg;
  std::vector<MonitorElement*> ch_iso_EB_list_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_EB_list_Bkg;

  std::vector<MonitorElement*> Ntracks_EE_list_Bkg;
  std::vector<MonitorElement*> ch_iso_EE_list_Bkg;
  std::vector<MonitorElement*> rel_ch_iso_EE_list_Bkg;

  std::vector<MonitorElement*> Ele_pT_MTD_EB_list_Bkg;
  std::vector<MonitorElement*> Ele_eta_MTD_EB_list_Bkg;
  std::vector<MonitorElement*> Ele_phi_MTD_EB_list_Bkg;

  std::vector<MonitorElement*> Ele_pT_MTD_EE_list_Bkg;
  std::vector<MonitorElement*> Ele_eta_MTD_EE_list_Bkg;
  std::vector<MonitorElement*> Ele_phi_MTD_EE_list_Bkg;
};

// ------------ constructor and destructor --------------
MtdEleIsoAnalyzer::MtdEleIsoAnalyzer(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMinEta_(iConfig.getParameter<double>("trackMinimumEta")),
      trackMaxEta_(iConfig.getParameter<double>("trackMaximumEta")),
      rel_iso_cut_(iConfig.getParameter<double>("rel_iso_cut")),
      //pTcuts_(iConfig.getParameter<std::vector<double>>(
      //    "pt_cuts")),  // example that works for double, but not for MonitorElement*
      track_match_PV_(iConfig.getUntrackedParameter<bool>("optionTrackMatchToPV")),
      dt_sig_vtx_(iConfig.getUntrackedParameter<bool>("option_dtToPV")),
      dt_sig_track_(iConfig.getUntrackedParameter<bool>("option_dtToTrack")) {
  //Ntracks_EB_list_Sig(iConfig.getParameter<std::vector<MonitorElement*>>("Ntracks_EB_list_Sig_test")) { // Example that does not work, but does work for double type

  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ =
      consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTag_vtx"));  // Vtx 4D collection

  GsfElectronToken_EB_ = consumes<reco::GsfElectronCollection>(
      iConfig.getParameter<edm::InputTag>("inputEle_EB"));  // Barrel electron collection input/token
  GsfElectronToken_EE_ = consumes<reco::GsfElectronCollection>(
      iConfig.getParameter<edm::InputTag>("inputEle_EE"));  // Endcap electron collection input/token
  GenParticleToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("inputGenP"));

  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  //particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();
}

MtdEleIsoAnalyzer::~MtdEleIsoAnalyzer() {}

// ------------ method called for each event  ------------
void MtdEleIsoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = iEvent.getHandle(GenRecTrackToken_);
  auto VertexHandle = iEvent.getHandle(RecVertexToken_);
  std::vector<reco::Vertex> vertices = *VertexHandle;

  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  //const auto& tMtd = iEvent.get(tmtdToken_);
  //const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  //const auto& t0Src = iEvent.get(t0SrcToken_);
  //const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  //const auto& t0Safe = iEvent.get(t0SafePidToken_);
  //const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  
  auto eleHandle_EB = makeValid(iEvent.getHandle(GsfElectronToken_EB_));
  reco::GsfElectronCollection eleColl_EB = *(eleHandle_EB.product());

  auto eleHandle_EE = makeValid(iEvent.getHandle(GsfElectronToken_EE_));
  reco::GsfElectronCollection eleColl_EE = *(eleHandle_EE.product());

  auto GenPartHandle = makeValid(iEvent.getHandle(GenParticleToken_));
  reco::GenParticleCollection GenPartColl = *(GenPartHandle.product());
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
  // localEleCollection.shrink_to_fit();

  // Selecting the PV from 3D and 4D vertex collections
  reco::Vertex Vtx_chosen;

  // This part has to be included, because in ~1% of the events, the "good" vertex is the 1st one not the 0th one in the collection
  for (int iVtx = 0; iVtx < (int)vertices.size(); iVtx++) {
    const reco::Vertex& vertex = vertices.at(iVtx);

    if (!vertex.isFake() && vertex.ndof() >= 4) {
      Vtx_chosen = vertex;
      break;
    }
  }
  // Vertex selection ends

  for (const auto& ele : localEleCollection) {
    //bool ele_Matched = false;
    bool ele_Promt = false;

    float ele_track_source_dz =  fabs(ele.gsfTrack()->dz(Vtx_chosen.position()));
    float ele_track_source_dxy = fabs(ele.gsfTrack()->dxy(Vtx_chosen.position()));

    if (ele.pt() > 10 && fabs(ele.eta()) < 2.4 && ele_track_source_dz < max_dz_vtx_cut &&
        ele_track_source_dxy < max_dxy_vtx_cut) {  // selecting "good" RECO electrons

      math::XYZVector EleMomentum = ele.momentum();

      for (const auto& genParticle : GenPartColl) {
        math::XYZVector GenPartMomentum = genParticle.momentum();
        double dr_match = reco::deltaR(GenPartMomentum, EleMomentum);
        if (((genParticle.pdgId() == 11 && ele.charge() == -1) || (genParticle.pdgId() == -11 && ele.charge() == 1)) &&
            dr_match < 0.10) {  // dR match check and charge match check

          //ele_Matched = true; // GEN-RECO electrons only pass here

          // check if ele mother is Z,W or tau (W->tau->ele decays, as these still should be quite isolated)
          if (genParticle.mother()->pdgId() == 23 || fabs(genParticle.mother()->pdgId()) == 24 ||
              fabs(genParticle.mother()->pdgId()) == 15) { 
            ele_Promt = true;  // ele is defined as promt
          }
          break;
        } 
      }
    } else {
      continue;
    }

    const reco::TrackRef ele_TrkRef = ele.core()->ctfTrack();
    math::XYZVector EleSigTrackMomentumAtVtx = ele.gsfTrack()->momentum(); // not sure if correct, but works
    double EleSigTrackEtaAtVtx = ele.gsfTrack()->eta();

    // Electron timing information check (whether electron signal track can be found in general track collection, for 99% we can.)
    // if we found a track match, we add MTD timing information for this matched track and do the isolation check
    if (ele_TrkRef.isNonnull()) { 

      bool Barrel_ele = ele.isEB();  // for track pT/dz cuts (Different for EB and EE in TDR)

      float min_pt_cut = Barrel_ele ? min_pt_cut_EB : min_pt_cut_EE;
      float max_dz_cut = Barrel_ele ? max_dz_cut_EB : max_dz_cut_EE;

      double ele_sigTrkTime = t0Pid[ele_TrkRef]; // electron signal track MTD information
      double ele_sigTrkTimeErr = Sigmat0Pid[ele_TrkRef];
      double ele_sigTrkMtdMva = mtdQualMVA[ele_TrkRef];
      ele_sigTrkTimeErr =
          (ele_sigTrkMtdMva > min_track_mtd_mva_cut) ? ele_sigTrkTimeErr : -1;  // track MTD MVA cut/check

      if (ele_Promt) {
        // For signal (promt)
        if (Barrel_ele) {
          meEle_pt_tot_Sig_EB_->Fill(ele.pt());  // All selected electron information for efficiency plots later
          meEle_eta_tot_Sig_EB_->Fill(ele.eta());
          meEle_phi_tot_Sig_EB_->Fill(ele.phi());
        } else {
          meEle_pt_tot_Sig_EE_->Fill(ele.pt());  // All selected electron information for efficiency plots later
          meEle_eta_tot_Sig_EE_->Fill(ele.eta());
          meEle_phi_tot_Sig_EE_->Fill(ele.phi());
        }
      } else {
        // For background (non-promt)
        if (Barrel_ele) {
          meEle_pt_tot_Bkg_EB_->Fill(ele.pt());
          meEle_eta_tot_Bkg_EB_->Fill(ele.eta());
          meEle_phi_tot_Bkg_EB_->Fill(ele.phi());
        } else {
          meEle_pt_tot_Bkg_EE_->Fill(ele.pt());
          meEle_eta_tot_Bkg_EE_->Fill(ele.eta());
          meEle_phi_tot_Bkg_EE_->Fill(ele.phi());
        }
      }

      int N_tracks_noMTD = 0;  // values for no MTD case
      double pT_sum_noMTD = 0;
      double rel_pT_sum_noMTD = 0;
      std::vector<int> N_tracks_MTD(number_of_cuts, 0);       
      std::vector<double> pT_sum_MTD(number_of_cuts, 0.);     
      std::vector<double> rel_pT_sum_MTD(number_of_cuts, 0.); 

      int general_index = 0;
      for (const auto& trackGen : *GenRecTrackHandle) {
        const reco::TrackRef trackref_general(GenRecTrackHandle, general_index);
        general_index++;

        if (trackGen.pt() < min_pt_cut) {  // track pT cut
          continue;
        }

        if (fabs(trackGen.vz() - ele.gsfTrack()->vz()) > max_dz_cut) {  // general track vs signal track dz cut
          continue;
        }

        if (track_match_PV_) {
          if (Vtx_chosen.trackWeight(trackref_general) < 0.5) {  // cut for general track matching to PV
            continue;
          }
        }

        double dr_check = reco::deltaR(trackGen.momentum(), EleSigTrackMomentumAtVtx);
        double deta = fabs(trackGen.eta() - EleSigTrackEtaAtVtx);

        if (dr_check > min_dR_cut && dr_check < max_dR_cut &&
            deta > min_strip_cut) {  // checking if the track is inside isolation cone

          ++N_tracks_noMTD;
          pT_sum_noMTD += trackGen.pt();
        } else {
          continue;
        }

        // checking the MTD timing cuts

        double TrkMTDTime = t0Pid[trackref_general];  // MTD timing info for particular track fron general tracks
        double TrkMTDTimeErr = Sigmat0Pid[trackref_general];
        double TrkMTDMva = mtdQualMVA[trackref_general];
        TrkMTDTimeErr = (TrkMTDMva > min_track_mtd_mva_cut) ? TrkMTDTimeErr : -1;  // track MTD MVA cut/check

        if (dt_sig_track_) {
          double dt_sigTrk = 0;  // dt regular track vs signal track (absolute value)

          if (TrkMTDTimeErr > 0 && ele_sigTrkTimeErr > 0) {
            dt_sigTrk = fabs(TrkMTDTime - ele_sigTrkTime);
            //dt_sigTrk_signif = dt_sigTrk/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + ele_sigTrkTimeErr*ele_sigTrkTimeErr);
    std::cout << "dt_sigTrk "  << dt_sigTrk << " err trk " << TrkMTDTimeErr << " err ele " << ele_sigTrkTimeErr << " tot err " << std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + ele_sigTrkTimeErr*ele_sigTrkTimeErr)*3 << std::endl;

            for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
              if (dt_sigTrk < max_dt_track_cut[i]) {
                N_tracks_MTD[i] += 1;
                pT_sum_MTD[i]   += trackGen.pt();
              }
            }
          } else {
            // if there is no error for MTD information, we count the MTD isolation case same as noMTD
            for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
              N_tracks_MTD[i] += 1;
              pT_sum_MTD[i]   += trackGen.pt();
            }
          }
        }

        if (dt_sig_vtx_) {
          double dt_vtx = 0;  // dt regular track vs vtx

          if (TrkMTDTimeErr > 0 && Vtx_chosen.tError() > 0) {
            dt_vtx = TrkMTDTime - Vtx_chosen.t();
            //dt_vtx_signif = dt_vtx/std::sqrt(TrkMTDTimeErr*TrkMTDTimeErr + Vtx_chosen.tError()*Vtx_chosen.tError());
    std::cout << "dt_vtx "  << dt_vtx << " err trk " << TrkMTDTimeErr << " err vtx " << Vtx_chosen.tError() << " tot err " << TrkMTDTimeErr*3 << std::endl;

            for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
              if (dt_vtx < max_dt_vtx_cut[i]) {
                N_tracks_MTD[i] += 1;
                pT_sum_MTD[i]   += trackGen.pt();
              }
            }
          } else {
            // if there is no error for MTD information, we count the MTD isolation case same as noMTD
            for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
              N_tracks_MTD[i] += 1;
              pT_sum_MTD[i]   += trackGen.pt(); 
            }
          }
        }
      }

      rel_pT_sum_noMTD = pT_sum_noMTD / ele.gsfTrack()->pt();  // rel_ch_iso calculation
 std::cout << "no mtd: pt sum " << pT_sum_noMTD << " rel " << rel_pT_sum_noMTD << " ele pt " << ele.gsfTrack()->pt() << " ntracks " << N_tracks_noMTD << std::endl;
      for (long unsigned int i = 0; i < N_tracks_MTD.size(); i++) {
        rel_pT_sum_MTD[i] = pT_sum_MTD[i] / ele.gsfTrack()->pt();
 std::cout << "mtd:" << i << " pt sum " << pT_sum_MTD[i] << " rel " << rel_pT_sum_MTD[i] << " ele pt " << ele.gsfTrack()->pt() << " ntracks " << N_tracks_MTD[i] << std::endl;
      }

      // defining vectors for more efficient hist filling
      // Promt part
      Ntracks_EB_list_Sig = {meEleISO_Ntracks_MTD_1_Sig_EB_,
                             meEleISO_Ntracks_MTD_2_Sig_EB_,
                             meEleISO_Ntracks_MTD_3_Sig_EB_,
                             meEleISO_Ntracks_MTD_4_Sig_EB_,
                             meEleISO_Ntracks_MTD_5_Sig_EB_,
                             meEleISO_Ntracks_MTD_6_Sig_EB_,
                             meEleISO_Ntracks_MTD_7_Sig_EB_};
      ch_iso_EB_list_Sig = {meEleISO_chIso_MTD_1_Sig_EB_,
                            meEleISO_chIso_MTD_2_Sig_EB_,
                            meEleISO_chIso_MTD_3_Sig_EB_,
                            meEleISO_chIso_MTD_4_Sig_EB_,
                            meEleISO_chIso_MTD_5_Sig_EB_,
                            meEleISO_chIso_MTD_6_Sig_EB_,
                            meEleISO_chIso_MTD_7_Sig_EB_};
      rel_ch_iso_EB_list_Sig = {meEleISO_rel_chIso_MTD_1_Sig_EB_,
                                meEleISO_rel_chIso_MTD_2_Sig_EB_,
                                meEleISO_rel_chIso_MTD_3_Sig_EB_,
                                meEleISO_rel_chIso_MTD_4_Sig_EB_,
                                meEleISO_rel_chIso_MTD_5_Sig_EB_,
                                meEleISO_rel_chIso_MTD_6_Sig_EB_,
                                meEleISO_rel_chIso_MTD_7_Sig_EB_};

      Ntracks_EE_list_Sig = {meEleISO_Ntracks_MTD_1_Sig_EE_,
                             meEleISO_Ntracks_MTD_2_Sig_EE_,
                             meEleISO_Ntracks_MTD_3_Sig_EE_,
                             meEleISO_Ntracks_MTD_4_Sig_EE_,
                             meEleISO_Ntracks_MTD_5_Sig_EE_,
                             meEleISO_Ntracks_MTD_6_Sig_EE_,
                             meEleISO_Ntracks_MTD_7_Sig_EE_};
      ch_iso_EE_list_Sig = {meEleISO_chIso_MTD_1_Sig_EE_,
                            meEleISO_chIso_MTD_2_Sig_EE_,
                            meEleISO_chIso_MTD_3_Sig_EE_,
                            meEleISO_chIso_MTD_4_Sig_EE_,
                            meEleISO_chIso_MTD_5_Sig_EE_,
                            meEleISO_chIso_MTD_6_Sig_EE_,
                            meEleISO_chIso_MTD_7_Sig_EE_};
      rel_ch_iso_EE_list_Sig = {meEleISO_rel_chIso_MTD_1_Sig_EE_,
                                meEleISO_rel_chIso_MTD_2_Sig_EE_,
                                meEleISO_rel_chIso_MTD_3_Sig_EE_,
                                meEleISO_rel_chIso_MTD_4_Sig_EE_,
                                meEleISO_rel_chIso_MTD_5_Sig_EE_,
                                meEleISO_rel_chIso_MTD_6_Sig_EE_,
                                meEleISO_rel_chIso_MTD_7_Sig_EE_};

      Ele_pT_MTD_EB_list_Sig = {meEle_pt_MTD_1_Sig_EB_,
                                meEle_pt_MTD_2_Sig_EB_,
                                meEle_pt_MTD_3_Sig_EB_,
                                meEle_pt_MTD_4_Sig_EB_,
                                meEle_pt_MTD_5_Sig_EB_,
                                meEle_pt_MTD_6_Sig_EB_,
                                meEle_pt_MTD_7_Sig_EB_};
      Ele_eta_MTD_EB_list_Sig = {meEle_eta_MTD_1_Sig_EB_,
                                 meEle_eta_MTD_2_Sig_EB_,
                                 meEle_eta_MTD_3_Sig_EB_,
                                 meEle_eta_MTD_4_Sig_EB_,
                                 meEle_eta_MTD_5_Sig_EB_,
                                 meEle_eta_MTD_6_Sig_EB_,
                                 meEle_eta_MTD_7_Sig_EB_};
      Ele_phi_MTD_EB_list_Sig = {meEle_phi_MTD_1_Sig_EB_,
                                 meEle_phi_MTD_2_Sig_EB_,
                                 meEle_phi_MTD_3_Sig_EB_,
                                 meEle_phi_MTD_4_Sig_EB_,
                                 meEle_phi_MTD_5_Sig_EB_,
                                 meEle_phi_MTD_6_Sig_EB_,
                                 meEle_phi_MTD_7_Sig_EB_};

      Ele_pT_MTD_EE_list_Sig = {meEle_pt_MTD_1_Sig_EE_,
                                meEle_pt_MTD_2_Sig_EE_,
                                meEle_pt_MTD_3_Sig_EE_,
                                meEle_pt_MTD_4_Sig_EE_,
                                meEle_pt_MTD_5_Sig_EE_,
                                meEle_pt_MTD_6_Sig_EE_,
                                meEle_pt_MTD_7_Sig_EE_};
      Ele_eta_MTD_EE_list_Sig = {meEle_eta_MTD_1_Sig_EE_,
                                 meEle_eta_MTD_2_Sig_EE_,
                                 meEle_eta_MTD_3_Sig_EE_,
                                 meEle_eta_MTD_4_Sig_EE_,
                                 meEle_eta_MTD_5_Sig_EE_,
                                 meEle_eta_MTD_6_Sig_EE_,
                                 meEle_eta_MTD_7_Sig_EE_};
      Ele_phi_MTD_EE_list_Sig = {meEle_phi_MTD_1_Sig_EE_,
                                 meEle_phi_MTD_2_Sig_EE_,
                                 meEle_phi_MTD_3_Sig_EE_,
                                 meEle_phi_MTD_4_Sig_EE_,
                                 meEle_phi_MTD_5_Sig_EE_,
                                 meEle_phi_MTD_6_Sig_EE_,
                                 meEle_phi_MTD_7_Sig_EE_};

      // Non-promt part
      Ntracks_EB_list_Bkg = {meEleISO_Ntracks_MTD_1_Bkg_EB_,
                             meEleISO_Ntracks_MTD_2_Bkg_EB_,
                             meEleISO_Ntracks_MTD_3_Bkg_EB_,
                             meEleISO_Ntracks_MTD_4_Bkg_EB_,
                             meEleISO_Ntracks_MTD_5_Bkg_EB_,
                             meEleISO_Ntracks_MTD_6_Bkg_EB_,
                             meEleISO_Ntracks_MTD_7_Bkg_EB_};
      ch_iso_EB_list_Bkg = {meEleISO_chIso_MTD_1_Bkg_EB_,
                            meEleISO_chIso_MTD_2_Bkg_EB_,
                            meEleISO_chIso_MTD_3_Bkg_EB_,
                            meEleISO_chIso_MTD_4_Bkg_EB_,
                            meEleISO_chIso_MTD_5_Bkg_EB_,
                            meEleISO_chIso_MTD_6_Bkg_EB_,
                            meEleISO_chIso_MTD_7_Bkg_EB_};
      rel_ch_iso_EB_list_Bkg = {meEleISO_rel_chIso_MTD_1_Bkg_EB_,
                                meEleISO_rel_chIso_MTD_2_Bkg_EB_,
                                meEleISO_rel_chIso_MTD_3_Bkg_EB_,
                                meEleISO_rel_chIso_MTD_4_Bkg_EB_,
                                meEleISO_rel_chIso_MTD_5_Bkg_EB_,
                                meEleISO_rel_chIso_MTD_6_Bkg_EB_,
                                meEleISO_rel_chIso_MTD_7_Bkg_EB_};

      Ntracks_EE_list_Bkg = {meEleISO_Ntracks_MTD_1_Bkg_EE_,
                             meEleISO_Ntracks_MTD_2_Bkg_EE_,
                             meEleISO_Ntracks_MTD_3_Bkg_EE_,
                             meEleISO_Ntracks_MTD_4_Bkg_EE_,
                             meEleISO_Ntracks_MTD_5_Bkg_EE_,
                             meEleISO_Ntracks_MTD_6_Bkg_EE_,
                             meEleISO_Ntracks_MTD_7_Bkg_EE_};
      ch_iso_EE_list_Bkg = {meEleISO_chIso_MTD_1_Bkg_EE_,
                            meEleISO_chIso_MTD_2_Bkg_EE_,
                            meEleISO_chIso_MTD_3_Bkg_EE_,
                            meEleISO_chIso_MTD_4_Bkg_EE_,
                            meEleISO_chIso_MTD_5_Bkg_EE_,
                            meEleISO_chIso_MTD_6_Bkg_EE_,
                            meEleISO_chIso_MTD_7_Bkg_EE_};
      rel_ch_iso_EE_list_Bkg = {meEleISO_rel_chIso_MTD_1_Bkg_EE_,
                                meEleISO_rel_chIso_MTD_2_Bkg_EE_,
                                meEleISO_rel_chIso_MTD_3_Bkg_EE_,
                                meEleISO_rel_chIso_MTD_4_Bkg_EE_,
                                meEleISO_rel_chIso_MTD_5_Bkg_EE_,
                                meEleISO_rel_chIso_MTD_6_Bkg_EE_,
                                meEleISO_rel_chIso_MTD_7_Bkg_EE_};

      Ele_pT_MTD_EB_list_Bkg = {meEle_pt_MTD_1_Bkg_EB_,
                                meEle_pt_MTD_2_Bkg_EB_,
                                meEle_pt_MTD_3_Bkg_EB_,
                                meEle_pt_MTD_4_Bkg_EB_,
                                meEle_pt_MTD_5_Bkg_EB_,
                                meEle_pt_MTD_6_Bkg_EB_,
                                meEle_pt_MTD_7_Bkg_EB_};
      Ele_eta_MTD_EB_list_Bkg = {meEle_eta_MTD_1_Bkg_EB_,
                                 meEle_eta_MTD_2_Bkg_EB_,
                                 meEle_eta_MTD_3_Bkg_EB_,
                                 meEle_eta_MTD_4_Bkg_EB_,
                                 meEle_eta_MTD_5_Bkg_EB_,
                                 meEle_eta_MTD_6_Bkg_EB_,
                                 meEle_eta_MTD_7_Bkg_EB_};
      Ele_phi_MTD_EB_list_Bkg = {meEle_phi_MTD_1_Bkg_EB_,
                                 meEle_phi_MTD_2_Bkg_EB_,
                                 meEle_phi_MTD_3_Bkg_EB_,
                                 meEle_phi_MTD_4_Bkg_EB_,
                                 meEle_phi_MTD_5_Bkg_EB_,
                                 meEle_phi_MTD_6_Bkg_EB_,
                                 meEle_phi_MTD_7_Bkg_EB_};

      Ele_pT_MTD_EE_list_Bkg = {meEle_pt_MTD_1_Bkg_EE_,
                                meEle_pt_MTD_2_Bkg_EE_,
                                meEle_pt_MTD_3_Bkg_EE_,
                                meEle_pt_MTD_4_Bkg_EE_,
                                meEle_pt_MTD_5_Bkg_EE_,
                                meEle_pt_MTD_6_Bkg_EE_,
                                meEle_pt_MTD_7_Bkg_EE_};
      Ele_eta_MTD_EE_list_Bkg = {meEle_eta_MTD_1_Bkg_EE_,
                                 meEle_eta_MTD_2_Bkg_EE_,
                                 meEle_eta_MTD_3_Bkg_EE_,
                                 meEle_eta_MTD_4_Bkg_EE_,
                                 meEle_eta_MTD_5_Bkg_EE_,
                                 meEle_eta_MTD_6_Bkg_EE_,
                                 meEle_eta_MTD_7_Bkg_EE_};
      Ele_phi_MTD_EE_list_Bkg = {meEle_phi_MTD_1_Bkg_EE_,
                                 meEle_phi_MTD_2_Bkg_EE_,
                                 meEle_phi_MTD_3_Bkg_EE_,
                                 meEle_phi_MTD_4_Bkg_EE_,
                                 meEle_phi_MTD_5_Bkg_EE_,
                                 meEle_phi_MTD_6_Bkg_EE_,
                                 meEle_phi_MTD_7_Bkg_EE_};

      if (ele_Promt) {  // promt part
        if (Barrel_ele) {
          meEleISO_Ntracks_Sig_EB_->Fill(N_tracks_noMTD);  // Filling hists for Ntraks and chIso sums for noMTD case //
          meEleISO_chIso_Sig_EB_->Fill(pT_sum_noMTD);
          meEleISO_rel_chIso_Sig_EB_->Fill(rel_pT_sum_noMTD);

          for (long unsigned int j = 0; j < Ntracks_EB_list_Sig.size(); j++) {
            Ntracks_EB_list_Sig[j]->Fill(N_tracks_MTD[j]);
            ch_iso_EB_list_Sig[j]->Fill(pT_sum_MTD[j]);
            rel_ch_iso_EB_list_Sig[j]->Fill(rel_pT_sum_MTD[j]);
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meEle_pt_noMTD_Sig_EB_->Fill(ele.pt());
            meEle_eta_noMTD_Sig_EB_->Fill(ele.eta());
            meEle_phi_noMTD_Sig_EB_->Fill(ele.phi());
          }

          for (long unsigned int k = 0; k < Ntracks_EB_list_Sig.size(); k++) {
            if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
              Ele_pT_MTD_EB_list_Sig[k]->Fill(ele.pt());
              Ele_eta_MTD_EB_list_Sig[k]->Fill(ele.eta());
              Ele_phi_MTD_EB_list_Sig[k]->Fill(ele.phi());
            }
          }

        } else {  // for endcap

          meEleISO_Ntracks_Sig_EE_->Fill(N_tracks_noMTD);  // Filling hists for Ntraks and chIso sums for noMTD case //
          meEleISO_chIso_Sig_EE_->Fill(pT_sum_noMTD);
          meEleISO_rel_chIso_Sig_EE_->Fill(rel_pT_sum_noMTD);

          for (long unsigned int j = 0; j < Ntracks_EE_list_Sig.size(); j++) {
            Ntracks_EE_list_Sig[j]->Fill(N_tracks_MTD[j]);
            ch_iso_EE_list_Sig[j]->Fill(pT_sum_MTD[j]);
            rel_ch_iso_EE_list_Sig[j]->Fill(rel_pT_sum_MTD[j]);
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meEle_pt_noMTD_Sig_EE_->Fill(ele.pt());
            meEle_eta_noMTD_Sig_EE_->Fill(ele.eta());
            meEle_phi_noMTD_Sig_EE_->Fill(ele.phi());
          }

          for (long unsigned int k = 0; k < Ntracks_EE_list_Sig.size(); k++) {
            if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
              Ele_pT_MTD_EE_list_Sig[k]->Fill(ele.pt());
              Ele_eta_MTD_EE_list_Sig[k]->Fill(ele.eta());
              Ele_phi_MTD_EE_list_Sig[k]->Fill(ele.phi());
            }
          }
        }
      } else {  // non-promt part
        if (Barrel_ele) {
          meEleISO_Ntracks_Bkg_EB_->Fill(N_tracks_noMTD);  // Filling hists for Ntraks and chIso sums for noMTD case //
          meEleISO_chIso_Bkg_EB_->Fill(pT_sum_noMTD);
          meEleISO_rel_chIso_Bkg_EB_->Fill(rel_pT_sum_noMTD);

          for (long unsigned int j = 0; j < Ntracks_EB_list_Bkg.size(); j++) {
            Ntracks_EB_list_Bkg[j]->Fill(N_tracks_MTD[j]);
            ch_iso_EB_list_Bkg[j]->Fill(pT_sum_MTD[j]);
            rel_ch_iso_EB_list_Bkg[j]->Fill(rel_pT_sum_MTD[j]);
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meEle_pt_noMTD_Bkg_EB_->Fill(ele.pt());
            meEle_eta_noMTD_Bkg_EB_->Fill(ele.eta());
            meEle_phi_noMTD_Bkg_EB_->Fill(ele.phi());
          }

          for (long unsigned int k = 0; k < Ntracks_EB_list_Bkg.size(); k++) {
            if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
              Ele_pT_MTD_EB_list_Bkg[k]->Fill(ele.pt());
              Ele_eta_MTD_EB_list_Bkg[k]->Fill(ele.eta());
              Ele_phi_MTD_EB_list_Bkg[k]->Fill(ele.phi());
            }
          }

        } else {                                           // for endcap
          meEleISO_Ntracks_Bkg_EE_->Fill(N_tracks_noMTD);  // Filling hists for Ntraks and chIso sums for noMTD case //
          meEleISO_chIso_Bkg_EE_->Fill(pT_sum_noMTD);
          meEleISO_rel_chIso_Bkg_EE_->Fill(rel_pT_sum_noMTD);

          for (long unsigned int j = 0; j < Ntracks_EE_list_Bkg.size(); j++) {
            Ntracks_EE_list_Bkg[j]->Fill(N_tracks_MTD[j]);
            ch_iso_EE_list_Bkg[j]->Fill(pT_sum_MTD[j]);
            rel_ch_iso_EE_list_Bkg[j]->Fill(rel_pT_sum_MTD[j]);
          }

          if (rel_pT_sum_noMTD < rel_iso_cut_) {  // filling hists for iso efficiency calculations
            meEle_pt_noMTD_Bkg_EE_->Fill(ele.pt());
            meEle_eta_noMTD_Bkg_EE_->Fill(ele.eta());
            meEle_phi_noMTD_Bkg_EE_->Fill(ele.phi());
          }

          for (long unsigned int k = 0; k < Ntracks_EE_list_Bkg.size(); k++) {
            if (rel_pT_sum_MTD[k] < rel_iso_cut_) {
              Ele_pT_MTD_EE_list_Bkg[k]->Fill(ele.pt());
              Ele_eta_MTD_EE_list_Bkg[k]->Fill(ele.eta());
              Ele_phi_MTD_EE_list_Bkg[k]->Fill(ele.phi());
            }
          }
        }
      }
    }  // electron matched to a track
  }    // electron collection inside single event
}

// ------------ method for histogram booking ------------
void MtdEleIsoAnalyzer::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking

  // signal
  meEleISO_Ntracks_Sig_EB_ = ibook.book1D("Ele_Iso_Ntracks_Sig_EB",
                                          "Tracks in isolation cone around electron track after basic cuts",
                                          20,
                                          0,
                                          20);  // hists for electrons
  meEleISO_chIso_Sig_EB_ = ibook.book1D(
      "Ele_chIso_sum_Sig_EB", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_1_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_1_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_1_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_1_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_1_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_1_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_2_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_2_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_2_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_2_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_2_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_2_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_3_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_3_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_3_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_3_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_3_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_3_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_4_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_4_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_4_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_4_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_4_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_4_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_5_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_5_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_5_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_5_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_5_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_5_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_6_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_6_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_6_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_6_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_6_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_6_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_7_Sig_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_7_Sig_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_7_Sig_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_7_Sig_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_7_Sig_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_7_Sig_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEle_pt_tot_Sig_EB_ =
      ibook.book1D("Ele_pT_tot_Sig_EB", "Electron pT tot", 30, 10, 100);  // hists for ele isto stuff start
  meEle_pt_noMTD_Sig_EB_ = ibook.book1D("Ele_pT_noMTD_Sig_EB", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_Sig_EB_ = ibook.book1D("Ele_eta_tot_Sig_EB", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_Sig_EB_ = ibook.book1D("Ele_eta_noMTD_Sig_EB", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_Sig_EB_ = ibook.book1D("Ele_phi_tot_Sig_EB", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_Sig_EB_ = ibook.book1D("Ele_phi_noMTD_Sig_EB", "Electron phi noMTD", 128, -3.2, 3.2);

  meEle_pt_MTD_1_Sig_EB_ = ibook.book1D("Ele_pT_MTD_1_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_Sig_EB_ = ibook.book1D("Ele_eta_MTD_1_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_Sig_EB_ = ibook.book1D("Ele_phi_MTD_1_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_Sig_EB_ = ibook.book1D("Ele_pT_MTD_2_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_Sig_EB_ = ibook.book1D("Ele_eta_MTD_2_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_Sig_EB_ = ibook.book1D("Ele_phi_MTD_2_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_Sig_EB_ = ibook.book1D("Ele_pT_MTD_3_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_Sig_EB_ = ibook.book1D("Ele_eta_MTD_3_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_Sig_EB_ = ibook.book1D("Ele_phi_MTD_3_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_Sig_EB_ = ibook.book1D("Ele_pT_MTD_4_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_Sig_EB_ = ibook.book1D("Ele_eta_MTD_4_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_Sig_EB_ = ibook.book1D("Ele_phi_MTD_4_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_Sig_EB_ = ibook.book1D("Ele_pT_MTD_5_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_Sig_EB_ = ibook.book1D("Ele_eta_MTD_5_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_Sig_EB_ = ibook.book1D("Ele_phi_MTD_5_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_Sig_EB_ = ibook.book1D("Ele_pT_MTD_6_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_Sig_EB_ = ibook.book1D("Ele_eta_MTD_6_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_Sig_EB_ = ibook.book1D("Ele_phi_MTD_6_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_Sig_EB_ = ibook.book1D("Ele_pT_MTD_7_Sig_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_Sig_EB_ = ibook.book1D("Ele_eta_MTD_7_Sig_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_Sig_EB_ = ibook.book1D("Ele_phi_MTD_7_Sig_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEleISO_Ntracks_Sig_EE_ = ibook.book1D("Ele_Iso_Ntracks_Sig_EE",
                                          "Tracks in isolation cone around electron track after basic cuts",
                                          20,
                                          0,
                                          20);  // hists for electrons
  meEleISO_chIso_Sig_EE_ = ibook.book1D(
      "Ele_chIso_sum_Sig_EE", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_1_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_1_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_1_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_1_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_1_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_1_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_2_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_2_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_2_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_2_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_2_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_2_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_3_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_3_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_3_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_3_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_3_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_3_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_4_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_4_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_4_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_4_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_4_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_4_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_5_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_5_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_5_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_5_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_5_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_5_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_6_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_6_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_6_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_6_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_6_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_6_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_7_Sig_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_7_Sig_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_7_Sig_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_7_Sig_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_7_Sig_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_7_Sig_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEle_pt_tot_Sig_EE_ =
      ibook.book1D("Ele_pT_tot_Sig_EE", "Electron pT tot", 30, 10, 100);  // hists for ele isto stuff start
  meEle_pt_noMTD_Sig_EE_ = ibook.book1D("Ele_pT_noMTD_Sig_EE", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_Sig_EE_ = ibook.book1D("Ele_eta_tot_Sig_EE", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_Sig_EE_ = ibook.book1D("Ele_eta_noMTD_Sig_EE", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_Sig_EE_ = ibook.book1D("Ele_phi_tot_Sig_EE", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_Sig_EE_ = ibook.book1D("Ele_phi_noMTD_Sig_EE", "Electron phi noMTD", 128, -3.2, 3.2);

  meEle_pt_MTD_1_Sig_EE_ = ibook.book1D("Ele_pT_MTD_1_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_Sig_EE_ = ibook.book1D("Ele_eta_MTD_1_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_Sig_EE_ = ibook.book1D("Ele_phi_MTD_1_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_Sig_EE_ = ibook.book1D("Ele_pT_MTD_2_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_Sig_EE_ = ibook.book1D("Ele_eta_MTD_2_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_Sig_EE_ = ibook.book1D("Ele_phi_MTD_2_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_Sig_EE_ = ibook.book1D("Ele_pT_MTD_3_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_Sig_EE_ = ibook.book1D("Ele_eta_MTD_3_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_Sig_EE_ = ibook.book1D("Ele_phi_MTD_3_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_Sig_EE_ = ibook.book1D("Ele_pT_MTD_4_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_Sig_EE_ = ibook.book1D("Ele_eta_MTD_4_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_Sig_EE_ = ibook.book1D("Ele_phi_MTD_4_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_Sig_EE_ = ibook.book1D("Ele_pT_MTD_5_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_Sig_EE_ = ibook.book1D("Ele_eta_MTD_5_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_Sig_EE_ = ibook.book1D("Ele_phi_MTD_5_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_Sig_EE_ = ibook.book1D("Ele_pT_MTD_6_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_Sig_EE_ = ibook.book1D("Ele_eta_MTD_6_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_Sig_EE_ = ibook.book1D("Ele_phi_MTD_6_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_Sig_EE_ = ibook.book1D("Ele_pT_MTD_7_Sig_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_Sig_EE_ = ibook.book1D("Ele_eta_MTD_7_Sig_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_Sig_EE_ = ibook.book1D("Ele_phi_MTD_7_Sig_EE", "Electron phi MTD", 128, -3.2, 3.2);

  // background
  meEleISO_Ntracks_Bkg_EB_ = ibook.book1D("Ele_Iso_Ntracks_Bkg_EB",
                                          "Tracks in isolation cone around electron track after basic cuts",
                                          20,
                                          0,
                                          20);  // hists for electrons
  meEleISO_chIso_Bkg_EB_ = ibook.book1D(
      "Ele_chIso_sum_Bkg_EB", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_1_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_1_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_1_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_1_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_1_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_1_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_2_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_2_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_2_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_2_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_2_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_2_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_3_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_3_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_3_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_3_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_3_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_3_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_4_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_4_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_4_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_4_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_4_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_4_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_5_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_5_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_5_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_5_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_5_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_5_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_6_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_6_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_6_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_6_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_6_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_6_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_7_Bkg_EB_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_7_Bkg_EB",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_7_Bkg_EB_ =
      ibook.book1D("Ele_chIso_sum_MTD_7_Bkg_EB",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_7_Bkg_EB_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_7_Bkg_EB",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEle_pt_tot_Bkg_EB_ =
      ibook.book1D("Ele_pT_tot_Bkg_EB", "Electron pT tot", 30, 10, 100);  // hists for ele isto stuff start
  meEle_pt_noMTD_Bkg_EB_ = ibook.book1D("Ele_pT_noMTD_Bkg_EB", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_Bkg_EB_ = ibook.book1D("Ele_eta_tot_Bkg_EB", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_Bkg_EB_ = ibook.book1D("Ele_eta_noMTD_Bkg_EB", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_Bkg_EB_ = ibook.book1D("Ele_phi_tot_Bkg_EB", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_Bkg_EB_ = ibook.book1D("Ele_phi_noMTD_Bkg_EB", "Electron phi noMTD", 128, -3.2, 3.2);

  meEle_pt_MTD_1_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_1_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_1_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_1_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_2_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_2_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_2_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_3_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_3_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_3_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_4_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_4_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_4_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_5_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_5_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_5_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_6_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_6_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_6_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_Bkg_EB_ = ibook.book1D("Ele_pT_MTD_7_Bkg_EB", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_Bkg_EB_ = ibook.book1D("Ele_eta_MTD_7_Bkg_EB", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_Bkg_EB_ = ibook.book1D("Ele_phi_MTD_7_Bkg_EB", "Electron phi MTD", 128, -3.2, 3.2);

  meEleISO_Ntracks_Bkg_EE_ = ibook.book1D("Ele_Iso_Ntracks_Bkg_EE",
                                          "Tracks in isolation cone around electron track after basic cuts",
                                          20,
                                          0,
                                          20);  // hists for electrons
  meEleISO_chIso_Bkg_EE_ = ibook.book1D(
      "Ele_chIso_sum_Bkg_EE", "Track pT sum in isolation cone around electron track after basic cuts", 400, 0, 20);
  meEleISO_rel_chIso_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_1_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_1_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_1_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_1_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_1_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_1_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_2_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_2_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_2_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_2_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_2_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_2_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_3_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_3_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_3_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_3_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_3_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_3_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_4_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_4_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_4_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_4_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_4_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_4_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_5_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_5_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_5_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_5_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_5_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_5_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_6_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_6_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_6_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_6_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_6_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_6_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEleISO_Ntracks_MTD_7_Bkg_EE_ =
      ibook.book1D("Ele_Iso_Ntracks_MTD_7_Bkg_EE",
                   "Tracks in isolation cone around electron track after basic cuts with MTD",
                   20,
                   0,
                   20);  // hists for electrons
  meEleISO_chIso_MTD_7_Bkg_EE_ =
      ibook.book1D("Ele_chIso_sum_MTD_7_Bkg_EE",
                   "Track pT sum in isolation cone around electron track after basic cuts with MTD",
                   400,
                   0,
                   20);
  meEleISO_rel_chIso_MTD_7_Bkg_EE_ =
      ibook.book1D("Ele_rel_chIso_sum_MTD_7_Bkg_EE",
                   "Track relative pT sum in isolation cone around electron track after basic cuts with MTD",
                   200,
                   0,
                   4);

  meEle_pt_tot_Bkg_EE_ =
      ibook.book1D("Ele_pT_tot_Bkg_EE", "Electron pT tot", 30, 10, 100);  // hists for ele isto stuff start
  meEle_pt_noMTD_Bkg_EE_ = ibook.book1D("Ele_pT_noMTD_Bkg_EE", "Electron pT noMTD", 30, 10, 100);

  meEle_eta_tot_Bkg_EE_ = ibook.book1D("Ele_eta_tot_Bkg_EE", "Electron eta tot", 128, -3.2, 3.2);
  meEle_eta_noMTD_Bkg_EE_ = ibook.book1D("Ele_eta_noMTD_Bkg_EE", "Electron eta noMTD", 128, -3.2, 3.2);

  meEle_phi_tot_Bkg_EE_ = ibook.book1D("Ele_phi_tot_Bkg_EE", "Electron phi tot", 128, -3.2, 3.2);
  meEle_phi_noMTD_Bkg_EE_ = ibook.book1D("Ele_phi_noMTD_Bkg_EE", "Electron phi noMTD", 128, -3.2, 3.2);

  meEle_pt_MTD_1_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_1_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_1_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_1_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_1_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_1_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_2_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_2_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_2_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_2_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_2_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_2_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_3_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_3_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_3_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_3_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_3_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_3_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_4_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_4_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_4_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_4_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_4_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_4_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_5_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_5_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_5_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_5_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_5_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_5_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_6_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_6_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_6_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_6_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_6_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_6_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);

  meEle_pt_MTD_7_Bkg_EE_ = ibook.book1D("Ele_pT_MTD_7_Bkg_EE", "Electron pT MTD", 30, 10, 100);
  meEle_eta_MTD_7_Bkg_EE_ = ibook.book1D("Ele_eta_MTD_7_Bkg_EE", "Electron eta MTD", 128, -3.2, 3.2);
  meEle_phi_MTD_7_Bkg_EE_ = ibook.book1D("Ele_phi_MTD_7_Bkg_EE", "Electron phi MTD", 128, -3.2, 3.2);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdEleIsoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/ElectronIso");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>(
      "inputTag_vtx",
      edm::InputTag("offlinePrimaryVertices4D"));  //  "offlinePrimaryVertices4D" or "offlinePrimaryVertices" (3D case)
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));

  desc.add<edm::InputTag>(
      "inputEle_EB",
      edm::InputTag(
          "gedGsfElectrons"));  // Adding the electron collection name ! (At least I think I'm doing that) // barrel ecal and track driven electrons
  //desc.add<edm::InputTag>("inputEle", edm::InputTag("ecalDrivenGsfElectrons")); // barrel + endcap, but without track seeded electrons
  desc.add<edm::InputTag>("inputEle_EE", edm::InputTag("ecalDrivenGsfElectronsHGC"));  // only endcap electrons
  desc.add<edm::InputTag>("inputGenP", edm::InputTag("genParticles"));

  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("trackMinimumPt", 1.0);  // [GeV]
  desc.add<double>("trackMinimumEta", 1.5);
  desc.add<double>("trackMaximumEta", 3.2);
  desc.add<double>("rel_iso_cut", 0.08);
//  desc.add<std::vector<double>>(
//      "pt_cuts",
//      {0.08,
//       0.10,
//       0.12});  // example for vector initialization, works for int, does not work for predifined pointers, like in the next line
//  //desc.add<std::vector<MonitorElement*>>("Ntracks_EB_list_Sig_test", {meEleISO_Ntracks_MTD_1_Sig_EB_,meEleISO_Ntracks_MTD_2_Sig_EB_,meEleISO_Ntracks_MTD_3_Sig_EB_,meEleISO_Ntracks_MTD_4_Sig_EB_,meEleISO_Ntracks_MTD_5_Sig_EB_,meEleISO_Ntracks_MTD_6_Sig_EB_,meEleISO_Ntracks_MTD_7_Sig_EB_}); // example that does not work...
  desc.addUntracked<bool>("optionTrackMatchToPV", false);
  desc.addUntracked<bool>("option_dtToPV", true);
  desc.addUntracked<bool>("option_dtToTrack", false);

  descriptions.add("mtdEleIsoAnalyzer", desc);
}

DEFINE_FWK_MODULE(MtdEleIsoAnalyzer);

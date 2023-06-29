#ifndef DataFormats_HGCalReco_TICLCandidate_h
#define DataFormats_HGCalReco_TICLCandidate_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// A TICLCandidate is a lightweight physics object made from one or multiple Tracksters.

class TICLCandidate : public reco::LeafCandidate {
public:
  typedef ticl::Trackster::ParticleType ParticleType;

  TICLCandidate(Charge q, const LorentzVector& p4)
      : LeafCandidate(q, p4), CALOtime_(0.f), CALOtimeError_(-1.f), t0Mtd_(0.f), t0MtdError_(-1.f), tMtd_(0.f), tMtdError_(-1.f), speed_(0.f), tMtdSim_(-1.f), tMtdPos_{0.f, 0.f, 0.f}, rawEnergy_(0.f), idProbabilities_{} {}

  TICLCandidate() : LeafCandidate(), CALOtime_(0.f), CALOtimeError_(-1.f), t0Mtd_(0.f), t0MtdError_(-1.f), tMtd_(0.f), tMtdError_(-1.f), speed_(0.f), tMtdSim_(-1.f), tMtdPos_{0.f, 0.f, 0.f}, rawEnergy_(0.f), idProbabilities_{} {}

  TICLCandidate(const edm::Ptr<ticl::Trackster>& trackster)
      : LeafCandidate(),
        CALOtime_(trackster->time()),
        CALOtimeError_(trackster->timeError()),
        t0Mtd_(trackster->t0Mtd()),
        t0MtdError_(trackster->t0MtdError()),
        tMtd_(trackster->tMtd()),
        tMtdError_(trackster->tMtdError()),
        speed_(trackster->speed()),
        tMtdSim_(trackster->MTDSimTime()),
        tMtdPos_(trackster->tMtdPos()),
        rawEnergy_(0.f),
        tracksters_({trackster}),
        idProbabilities_{} {}

  inline const float time() const { return CALOtime_; }
  inline const float timeError() const { return CALOtimeError_; }
  inline const float t0Mtd() const { return t0Mtd_; }
  inline const float t0MtdError() const { return t0MtdError_; }
  inline const float tMtd() const { return tMtd_; }
  inline const float tMtdError() const { return tMtdError_; }
  inline const float speed() const { return speed_; }
  inline const float MTDSimTime() const { return tMtdSim_; }
  inline const GlobalPoint tMtdPos() const { return tMtdPos_; }

  void setTime(float time) { CALOtime_ = time; };
  void setTimeError(float timeError) { CALOtimeError_ = timeError; }

  inline void sett0MtdTimeAndError(float t, float tError) {
    t0Mtd_ = t;
    t0MtdError_ = tError;
  }
  inline void settMtdTimeAndError(float t, float tError) {
    tMtd_ = t;
    tMtdError_ = tError;
  }
  inline void setSpeed(float v) {
    speed_ = v;
  }
  inline void setMTDSimTime(float t) {
    tMtdSim_ = t;
  }
  inline void settMtdPos(GlobalPoint pos) {
    tMtdPos_ = pos;
  }

  inline const edm::Ptr<reco::Track> trackPtr() const { return trackPtr_; }
  void setTrackPtr(const edm::Ptr<reco::Track>& trackPtr) { trackPtr_ = trackPtr; }

  inline float rawEnergy() const { return rawEnergy_; }
  void setRawEnergy(float rawEnergy) { rawEnergy_ = rawEnergy; }

  inline const std::vector<edm::Ptr<ticl::Trackster> > tracksters() const { return tracksters_; };

  void setTracksters(const std::vector<edm::Ptr<ticl::Trackster> >& tracksters) { tracksters_ = tracksters; }
  void addTrackster(const edm::Ptr<ticl::Trackster>& trackster) { tracksters_.push_back(trackster); }
  // convenience method to return the ID probability for a certain particle type
  inline float id_probability(ParticleType type) const {
    // probabilities are stored in the same order as defined in the ParticleType enum
    return idProbabilities_[(int)type];
  }

  void zeroProbabilities() {
    for (auto& p : idProbabilities_) {
      p = 0.f;
    }
  }

  void setIdProbabilities(const std::array<float, 8>& idProbs) { idProbabilities_ = idProbs; }
  inline void setIdProbability(ParticleType type, float value) { idProbabilities_[int(type)] = value; }

private:
  float CALOtime_;
  float CALOtimeError_;
  // MTD time
  float t0Mtd_;
  float t0MtdError_;
  float tMtd_;
  float tMtdError_;
  float speed_;

  // float simt0Mtd_; maybe?
  float tMtdSim_;
  GlobalPoint tMtdPos_;  

  edm::Ptr<reco::Track> trackPtr_;

  float rawEnergy_;

  // vector of Ptr so Tracksters can come from different collections
  // and there can be derived classes
  std::vector<edm::Ptr<ticl::Trackster> > tracksters_;

  // Since it contains multiple tracksters, duplicate the probability interface
  std::array<float, 8> idProbabilities_;
};
#endif

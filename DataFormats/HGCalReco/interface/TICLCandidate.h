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
      : LeafCandidate(q, p4),
        VertexTime_{},
        BoundaryTime_{},
        CALOtime_{},
        t0Mtd_{},
        tMtd_{},
        betaMtd_(0.f),
        pathMtd_(0.f),
        tMtdPos_{},
        rawEnergy_(0.f),
        idProbabilities_{} {}

  TICLCandidate()
      : LeafCandidate(),
        VertexTime_{},
        BoundaryTime_{},
        CALOtime_{},
        t0Mtd_{},
        tMtd_{},
        betaMtd_(0.f),
        pathMtd_(0.f),
        tMtdPos_{},
        rawEnergy_(0.f),
        idProbabilities_{} {}

  TICLCandidate(const edm::Ptr<ticl::Trackster>& trackster)
      : LeafCandidate(),
        VertexTime_(trackster->VertexTimeAndErr()),
        BoundaryTime_(trackster->BoundaryTimeAndErr()),
        CALOtime_(trackster->timeAndErr()),
        t0Mtd_(trackster->t0MtdAndErr()),
        tMtd_(trackster->tMtdAndErr()),
        betaMtd_(trackster->betaMtd()),
        pathMtd_(trackster->pathMtd()),
        tMtdPos_(trackster->tMtdPos()),
        rawEnergy_(0.f),
        tracksters_({trackster}),
        idProbabilities_{} {}

  inline const float BoundaryTime() const { return BoundaryTime_.time; }
  inline const float BoundaryTimeError() const { return BoundaryTime_.timeErr; }
  inline const Time BoundaryTimeAndErr() const { return BoundaryTime_; }
  inline const float VertexTime() const { return VertexTime_.time; }
  inline const float VertexTimeError() const { return VertexTime_.timeErr; }
  inline const Time VertexTimeAndErr() const { return VertexTime_; }
  inline const float time() const { return CALOtime_.time; }
  inline const float timeError() const { return CALOtime_.timeErr; }
  inline const Time timeAndErr() const { return CALOtime_; }
  inline const float t0Mtd() const { return t0Mtd_.time; }
  inline const float t0MtdError() const { return t0Mtd_.timeErr; }
  inline const Time t0MtdAndErr() const { return t0Mtd_; }
  inline const float tMtd() const { return tMtd_.time; }
  inline const float tMtdError() const { return tMtd_.timeErr; }
  inline const Time tMtdAndErr() const { return tMtd_; }
  inline const float betaMtd() const { return betaMtd_; }
  inline const float pathMtd() const { return pathMtd_; }
  inline const GlobalPoint tMtdPos() const { return tMtdPos_; }

  void setTime(float time) { CALOtime_.time = time; };
  void setTimeError(float timeError) { CALOtime_.timeErr = timeError; }

  inline void setBoundaryTimeAndError(float t, float tError) {
    BoundaryTime_ = Time{t,tError};
  }
  inline void setVertexTimeAndError(float t, float tError) {
    VertexTime_ = Time{t,tError};
  }
  inline void setTimeAndError(float t, float tError) {
    CALOtime_ = Time{t,tError};
  }
  inline void sett0MtdTimeAndError(float t, float tError) {
    t0Mtd_ = Time{t,tError};
  }
  inline void settMtdTimeAndError(float t, float tError) {
    tMtd_ = Time{t, tError};
  }
  inline void setBetaMtd(float b) { betaMtd_ = b; }
  inline void setPathMtd(float b) { pathMtd_ = b; }
  inline void settMtdPos(GlobalPoint pos) { tMtdPos_ = pos; }

  inline const edm::Ptr<reco::Track> trackPtr() const { return trackPtr_; }
  void setTrackPtr(const edm::Ptr<reco::Track>& trackPtr) { trackPtr_ = trackPtr; }

  inline float rawEnergy() const { return rawEnergy_; }
  void setRawEnergy(float rawEnergy) { rawEnergy_ = rawEnergy; }

  inline const std::vector<edm::Ptr<ticl::Trackster> > tracksters() const { return tracksters_; };

  void setTracksters(const std::vector<edm::Ptr<ticl::Trackster> >& tracksters) { tracksters_ = tracksters; }
  void addTrackster(const edm::Ptr<ticl::Trackster>& trackster) {
    tracksters_.push_back(trackster);
 // not needed now, bc also for candidates with 1 trackster we have to do the propagation
 //   CALOtime_ = trackster->time();
 //   CALOtimeError_ = trackster->timeError();
  }

  // convenience method to return the ID probability for a certain particle type
  inline float id_probability(ParticleType type) const {
    // probabilities are stored in the same order as defined in the ParticleType enum
    return idProbabilities_[(int)type];
  }

  inline const std::array<float, 8>& idProbabilities() const { return idProbabilities_; }

  void zeroProbabilities() {
    for (auto& p : idProbabilities_) {
      p = 0.f;
    }
  }

  void setIdProbabilities(const std::array<float, 8>& idProbs) { idProbabilities_ = idProbs; }
  inline void setIdProbability(ParticleType type, float value) { idProbabilities_[int(type)] = value; }

private:
  // sim
  Time VertexTime_;
  Time BoundaryTime_;
  // reco
  Time CALOtime_;
  // mixed
  Time t0Mtd_;
  Time tMtd_;
  float betaMtd_;
  float pathMtd_;
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

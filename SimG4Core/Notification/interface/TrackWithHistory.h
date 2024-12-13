#ifndef SimG4Core_TrackWithHistory_H
#define SimG4Core_TrackWithHistory_H

#include "G4Track.hh"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "G4Allocator.hh"

class G4PrimaryParticle;
/** The part of the information about a SimTrack that we need from
 *  a G4Track
 */

class TrackWithHistory {
public:
  /** The constructor is called at time,
     *  when some of the information may not available yet.
     */
  TrackWithHistory(const G4Track *g4track, int pID);
  TrackWithHistory(const G4PrimaryParticle *, int trackID, const math::XYZVectorD &pos, const double time);
  ~TrackWithHistory() = default;

  inline void *operator new(std::size_t);
  inline void operator delete(void *TrackWithHistory);

  int trackID() const { return trackID_; }
  int particleID() const { return pdgID_; }
  int parentID() const { return parentID_; }
  int genParticleID() const { return isPrimary() ? genParticleID_ : -1; }
  int vertexID() const { return vertexID_; }
  int processType() const { return procType_; }
  int getIDAtBoundary() const { return idAtBoundary_; }

  void setTrackID(int i) { trackID_ = i; }
  void setParentID(int i) { parentID_ = i; }
  void setVertexID(int i) { vertexID_ = i; }
  void setGenParticleID(int i) { genParticleID_ = i; }

  double totalEnergy() const { return totalEnergy_; }
  double time() const { return time_; }
  double weight() const { return weight_; }
  void setToBeSaved() { trackInfo_ |= 1 << 4;; }
  bool storeTrack() const { return (trackInfo_ >> 3) & 1; }
  bool saved() const { return (trackInfo_ >> 4) & 1; }
  bool crossedBoundary() const { return (trackInfo_ >> 2) & 1; }

  const math::XYZVectorD &momentum() const { return momentum_; }
  const math::XYZVectorD &vertexPosition() const { return vertexPosition_; }

  // Boundary crossing variables
  void setCrossedBoundaryPosMom(int id,
                                const math::XYZTLorentzVectorF &position,
                                const math::XYZTLorentzVectorF &momentum) {
    trackInfo_ |= 1 << 2;
    idAtBoundary_ = id;
    positionAtBoundary_ = position;
    momentumAtBoundary_ = momentum;
  }
  const math::XYZTLorentzVectorF &getPositionAtBoundary() const { return positionAtBoundary_; }
  const math::XYZTLorentzVectorF &getMomentumAtBoundary() const { return momentumAtBoundary_; }

  // tracker surface
  const math::XYZVectorD &trackerSurfacePosition() const { return tkSurfacePosition_; }
  const math::XYZTLorentzVectorD &trackerSurfaceMomentum() const { return tkSurfaceMomentum_; }
  void setSurfacePosMom(const math::XYZVectorD &pos, const math::XYZTLorentzVectorD &mom) {
    tkSurfacePosition_ = pos;
    tkSurfaceMomentum_ = mom;
  }
  bool isFromBackScattering() const { return (trackInfo_ >> 0) & 1; }
  void setFromBackScattering() { trackInfo_ |= 1 << 0; }

  bool isPrimary() const { return (trackInfo_ >> 1) & 1; }
  void setIsPrimary() { trackInfo_ |= 1 << 1; }
  int getPrimaryID() const { return genParticleID_; }

private:
  int trackID_;
  int pdgID_;
  int parentID_;
  int genParticleID_{-1};
  int vertexID_{-1};
  int idAtBoundary_{-1};
  int procType_{0};
  double totalEnergy_;
  double time_;  // lab system
  double weight_;
  math::XYZVectorD momentum_;
  math::XYZVectorD vertexPosition_;
  math::XYZTLorentzVectorF positionAtBoundary_{math::XYZTLorentzVectorF(0.f, 0.f, 0.f, 0.f)};
  math::XYZTLorentzVectorF momentumAtBoundary_{math::XYZTLorentzVectorF(0.f, 0.f, 0.f, 0.f)};
  math::XYZVectorD tkSurfacePosition_{math::XYZVectorD(0., 0., 0.)};
  math::XYZTLorentzVectorD tkSurfaceMomentum_{math::XYZTLorentzVectorD(0., 0., 0., 0.)};
  uint8_t trackInfo_; // 0 = isFromBackScattering, 1 = isPrimary, 2 = crossedBoundary, 3 = storeTrack, 4 = saved
};

extern G4ThreadLocal G4Allocator<TrackWithHistory> *fpTrackWithHistoryAllocator;

inline void *TrackWithHistory::operator new(std::size_t) {
  if (!fpTrackWithHistoryAllocator)
    fpTrackWithHistoryAllocator = new G4Allocator<TrackWithHistory>;
  return (void *)fpTrackWithHistoryAllocator->MallocSingle();
}

inline void TrackWithHistory::operator delete(void *aTwH) {
  fpTrackWithHistoryAllocator->FreeSingle((TrackWithHistory *)aTwH);
}

#endif

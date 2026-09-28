#ifndef EgammaReco_SuperCluster_h
#define EgammaReco_SuperCluster_h
/** \class reco::SuperCluster SuperCluster.h DataFormats/EgammaReco/interface/SuperCluster.h
 *
 * A SuperCluster reconstructed in the Electromagnetic Calorimeter
 * contains references to seed and constituent BasicClusters
 *
 * \author Luca Lista, INFN
 *
 *
 */
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

namespace reco {
  namespace io_v1 {
    class SuperCluster : public CaloCluster {
    public:
      typedef math::XYZPoint Point;

      /// default constructor
      SuperCluster()
          : CaloCluster(0., Point(0., 0., 0.)),
            rawEnergy_(-1.),
            phiWidth_(0),
            etaWidth_(0) {}

      /// constructor defined by CaloCluster - will have to use setSeed and add() separately
      SuperCluster(double energy, const Point& position);

      SuperCluster(double energy,
                   const Point& position,
                   const CaloClusterPtr& seed,
                   const CaloClusterPtrVector& clusters,
                   double phiWidth = 0.,
                   double etaWidth = 0.);

      /// raw uncorrected energy (sum of energies of component BasicClusters)
      double rawEnergy() const { return rawEnergy_; }

      /// obtain phi and eta width of the Super Cluster
      double phiWidth() const { return phiWidth_; }
      double etaWidth() const { return etaWidth_; }

      void setPhiWidth(double pw) { phiWidth_ = pw; }
      void setEtaWidth(double ew) { etaWidth_ = ew; }

      /// seed BasicCluster
      const CaloClusterPtr& seed() const { return seed_; }

      /// const access to the cluster list itself
      const CaloClusterPtrVector& clusters() const { return clusters_; }

      /// fist iterator over BasicCluster constituents
      CaloCluster_iterator clustersBegin() const { return clusters_.begin(); }

      /// last iterator over BasicCluster constituents
      CaloCluster_iterator clustersEnd() const { return clusters_.end(); }

      /// number of BasicCluster constituents
      size_t clustersSize() const { return clusters_.size(); }

      /// list of used xtals by DetId // now inherited by CaloCluster
      //std::vector<DetId> getHitsByDetId() const { return usedHits_; }

      /// set reference to seed BasicCluster
      void setSeed(const CaloClusterPtr& r) { seed_ = r; }

      //(re)-set clusters
      void setClusters(const CaloClusterPtrVector& clusters) {
        clusters_ = clusters;
        computeRawEnergy();
      }

      //clear hits and fractions vector (for slimming)
      void clearHitsAndFractions() { hitsAndFractions_.clear(); }

      /// add reference to constituent BasicCluster
      void addCluster(const CaloClusterPtr& r) {
        clusters_.push_back(r);
        computeRawEnergy();
      }

      const int seedCrysIEtaOrIx() const {
        auto detid = seed_->seed();
        int ietaorix = 0;
        if (detid.subdetId() == EcalBarrel) {
          EBDetId ebdetid(detid);
          ietaorix = ebdetid.ieta();
        } else if (detid.subdetId() == EcalEndcap) {
          EEDetId eedetid(detid);
          ietaorix = eedetid.ix();
        }
        return ietaorix;
      }

      const int seedCrysIPhiOrIy() const {
        auto detid = seed_->seed();
        int iphioriy = 0;
        if (detid.subdetId() == EcalBarrel) {
          EBDetId ebdetid(detid);
          iphioriy = ebdetid.iphi();
        } else if (detid.subdetId() == EcalEndcap) {
          EEDetId eedetid(detid);
          iphioriy = eedetid.iy();
        }
        return iphioriy;
      }
      double preshowerEnergy() const { return preshowerEnergy_; }

      void setPreshowerEnergy(double preshowerEnergy) { preshowerEnergy_ = preshowerEnergy; };


    private:
      void computeRawEnergy();

      /// reference to BasicCluster seed
      CaloClusterPtr seed_;

      /// references to BasicCluster constitunets
      CaloClusterPtrVector clusters_;

      /// used hits by detId - retrieved from BC constituents -- now inherited from CaloCluster
      //std::vector<DetId> usedHits_;

      double preshowerEnergy_ = 0.;

      double rawEnergy_;

      double phiWidth_;
      double etaWidth_;
    };

  }  // namespace io_v1
  using SuperCluster = io_v1::SuperCluster;
}  // namespace reco
#endif

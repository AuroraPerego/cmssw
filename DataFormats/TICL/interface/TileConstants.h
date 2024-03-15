#ifndef DataFormats_TICL_TileConstants_h
#define DataFormats_TICL_TileConstants_h

#include <vector>
#include <array>
#include <cstdint>
#include <cmath>

namespace ticl {
  struct TileConstantsGlobal_EtaPhi {
  static constexpr float tileSize = 0.15f;
  static constexpr float minDim1 = -3.f;
  static constexpr float maxDim1 = 3.f;
  static constexpr float minDim2 = -M_PI;
  static constexpr float maxDim2 = M_PI;
  static constexpr bool wrapped = true;
  };

  struct TileConstantsEndcapNeg_EtaPhi {
    static constexpr float tileSize = 0.15f;
    static constexpr float minDim1 = -3.f;
    static constexpr float maxDim1 = -1.5f;
    static constexpr float minDim2 = -M_PI;
    static constexpr float maxDim2 = M_PI;
    static constexpr bool wrapped = true;

  };

    struct TileConstantsEndcapPos_EtaPhi {
    static constexpr float tileSize = 0.15f;
    static constexpr float minDim1 = 1.5f;
    static constexpr float maxDim1 = 3.f;
    static constexpr float minDim2 = -M_PI;
    static constexpr float maxDim2 = M_PI;
    static constexpr bool wrapped = true;

  };

  struct TileConstantsBarrel_EtaPhi {
    static constexpr float tileSize = 3*0.087f;
    static constexpr float minDim1 = -1.5f;
    static constexpr float maxDim1 = 1.5f;
    static constexpr float minDim2 = -M_PI;
    static constexpr float maxDim2 = M_PI;
    static constexpr bool wrapped = true;
  };


  struct TileConstantsEndcap_XY {
    static constexpr float tileSize = 5.f;
    static constexpr float minDim1 = -285.f;
    static constexpr float maxDim1 = 285.f;
    static constexpr float minDim2 = -285.f;
    static constexpr float maxDim2 = 285.f;
    static constexpr bool wrapped = false;

  };

  struct TileConstantsHFNose_XY {
    static constexpr float tileSize = 5.f;
    static constexpr float minDim1 = -110.f;
    static constexpr float maxDim1 = 110.f;
    static constexpr float minDim2 = -110.f;
    static constexpr float maxDim2 = 110.f;
    static constexpr bool wrapped = false;
  };

  struct TileConstants {
    static constexpr float minEta = 1.5f;
    static constexpr float maxEta = 3.2f;
    static constexpr int nEtaBins = 34;
    static constexpr int nPhiBins = 126;
    static constexpr int nLayers = 104;
    static constexpr int iterations = 4;
    static constexpr int nBins = nEtaBins * nPhiBins;
  };

  struct TileConstantsHFNose {
    static constexpr float minEta = 3.0f;
    static constexpr float maxEta = 4.2f;
    static constexpr int nEtaBins = 24;
    static constexpr int nPhiBins = 126;
    static constexpr int nLayers = 16;  // 8x2
    static constexpr int iterations = 4;
    static constexpr int nBins = nEtaBins * nPhiBins;
  };

}  // namespace ticl

namespace ticl {
  typedef std::vector<std::pair<unsigned int, float> > TICLClusterFilterMask;
}  // namespace ticl

#endif

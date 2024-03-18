// Authors: Marco Rovere, Felice Pantaleo - marco.rovere@cern.ch, felice.pantaleo@cern.ch
// Date: 05/2019

#ifndef DataFormats_HGCalReco_TICLLayerTile_h
#define DataFormats_HGCalReco_TICLLayerTile_h

#include "DataFormats/TICL/interface/TICLTile.h"

namespace ticl {
  using TICLLayerTile = Tile<TileConstants>;
  using Tiles = std::array<TICLLayerTile, TileConstants::nLayers>;
  using TracksterTiles = std::array<TICLLayerTile, TileConstants::iterations>;

  using TICLLayerTileHFNose = Tile<TileConstantsHFNose>;
  using TilesHFNose = std::array<TICLLayerTileHFNose, TileConstantsHFNose::nLayers>;
  using TracksterTilesHFNose = std::array<TICLLayerTileHFNose, TileConstantsHFNose::iterations>;

}  // namespace ticl

template <typename T>
class TICLGenericTile {
public:
  // value_type_t is the type of the type of the array used by the incoming <T> type.
  using constants_type_t = typename T::value_type::type;
  using tile_type_t = typename T::value_type;
  // This class represents a generic collection of Tiles. The additional index
  // numbering is not handled internally. It is the user's responsibility to
  // properly use and consistently access it here.
  const auto& operator[](int index) const { return tiles_[index]; }
  void fill(int index, float eta, float phi, unsigned int objectId) { tiles_[index].fill(eta, phi, objectId); }

private:
  T tiles_;
};

using TICLLayerTiles = TICLGenericTile<ticl::Tiles>;
using TICLTracksterTiles = TICLGenericTile<ticl::TracksterTiles>;
using TICLLayerTilesHFNose = TICLGenericTile<ticl::TilesHFNose>;
using TICLTracksterTilesHFNose = TICLGenericTile<ticl::TracksterTilesHFNose>;

#endif

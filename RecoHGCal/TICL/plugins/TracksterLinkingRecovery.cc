#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingRecovery.h"

using namespace ticl;

void TracksterLinkingRecovery::linkTracksters(
    const Inputs& input,
    std::vector<Trackster>& resultTracksters,
    std::vector<unsigned int>& maskTracksters,
    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {
  resultTracksters.reserve(input.tracksters.size());
  linkedTracksterIdToInputTracksterId.resize(input.tracksters.size());
  // Merge all trackster collections into a single collection
  for (size_t i = 0; i < input.tracksters.size(); ++i) {
    resultTracksters.push_back(input.tracksters[i]);
    linkedTracksterIdToInputTracksterId[i].push_back(resultTracksters.size() - 1);
    maskTracksters[i] = 0;
  }
}

#ifndef RecoHGCal_TICL_MinCostFlow_h
#define RecoHGCal_TICL_MinCostFlow_h

#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace ticl {

  class MinCostFlow {
  public:
    enum Status { OPTIMAL, INFEASIBLE };

    explicit MinCostFlow(int n)
        : graph_(n), potential_(n, 0), n_(n) {}

    // Add arc u→v with given capacity and unit cost. Returns arc id.
    int addArc(int u, int v, int cap, int64_t cost) {
      assert(u >= 0 && u < n_ && v >= 0 && v < n_);
      int arcId = static_cast<int>(arcMeta_.size());
      arcMeta_.push_back({u, static_cast<int>(graph_[u].size()), cap});
      graph_[u].push_back({v, static_cast<int>(graph_[v].size()), cap, cost});
      graph_[v].push_back({u, static_cast<int>(graph_[u].size()) - 1, 0, -cost});
      return arcId;
    }

    // Set supply (positive) or demand (negative) on a node.
    void setNodeSupply(int v, int supply) {
      assert(v >= 0 && v < n_);
      supply_.resize(n_, 0);
      supply_[v] = supply;
    }

    // Solve by converting supplies/demands to a source/sink formulation,
    // then running Successive Shortest Paths with Dijkstra + Johnson potentials.
    Status solve() {
      supply_.resize(n_, 0);

      // Create super-source (n_) and super-sink (n_+1)
      int S = n_, T = n_ + 1;
      graph_.push_back({});  // node S
      graph_.push_back({});  // node T
      potential_.push_back(0);
      potential_.push_back(0);
      int totalSupply = 0;

      for (int v = 0; v < n_; ++v) {
        if (supply_[v] > 0) {
          // S → v
          graph_[S].push_back({v, static_cast<int>(graph_[v].size()), supply_[v], 0});
          graph_[v].push_back({S, static_cast<int>(graph_[S].size()) - 1, 0, 0});
          totalSupply += supply_[v];
        } else if (supply_[v] < 0) {
          // v → T
          graph_[v].push_back({T, static_cast<int>(graph_[T].size()), -supply_[v], 0});
          graph_[T].push_back({v, static_cast<int>(graph_[v].size()) - 1, 0, 0});
        }
      }

      int N = static_cast<int>(graph_.size());  // n_ + 2

      // Initialize potentials via Bellman-Ford from S
      potential_.assign(N, INF);
      potential_[S] = 0;
      for (int iter = 0; iter < N - 1; ++iter)
        for (int u = 0; u < N; ++u)
          if (potential_[u] < INF)
            for (const auto& e : graph_[u])
              if (e.cap > 0 && potential_[u] + e.cost < potential_[e.to])
                potential_[e.to] = potential_[u] + e.cost;

      // Successive Shortest Paths
      int pushed = 0;
      while (pushed < totalSupply) {
        // Dijkstra with Johnson reweighting
        std::vector<int64_t> dist(N, INF);
        std::vector<int> prevv(N, -1), preve(N, -1);
        dist[S] = 0;
        using P = std::pair<int64_t, int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        pq.push({0, S});

        while (!pq.empty()) {
          auto [d, u] = pq.top(); pq.pop();
          if (d > dist[u]) continue;
          for (int i = 0; i < static_cast<int>(graph_[u].size()); ++i) {
            const auto& e = graph_[u][i];
            if (e.cap <= 0) continue;
            // Reweighted cost: always non-negative after Bellman-Ford init
            int64_t newCost = e.cost + potential_[u] - potential_[e.to];
            if (newCost < 0) newCost = 0;  // numerical safety
            int64_t nd = dist[u] + newCost;
            if (nd < dist[e.to]) {
              dist[e.to] = nd;
              prevv[e.to] = u;
              preve[e.to] = i;
              pq.push({nd, e.to});
            }
          }
        }

        if (dist[T] == INF) return INFEASIBLE;

        // Update potentials
        for (int v = 0; v < N; ++v)
          if (dist[v] < INF)
            potential_[v] += dist[v];

        // Bottleneck capacity along shortest path
        int aug = totalSupply - pushed;
        for (int v = T; v != S; v = prevv[v])
          aug = std::min(aug, graph_[prevv[v]][preve[v]].cap);

        // Augment
        for (int v = T; v != S; v = prevv[v]) {
          auto& e = graph_[prevv[v]][preve[v]];
          e.cap -= aug;
          graph_[v][e.rev].cap += aug;
        }
        pushed += aug;
      }
      return OPTIMAL;
    }

    // Flow on a user arc (original cap - remaining cap)
    int flow(int arcId) const {
      const auto& m = arcMeta_[arcId];
      return m.origCap - graph_[m.node][m.edgeIdx].cap;
    }

    int tail(int arcId) const { return arcMeta_[arcId].node; }

    int head(int arcId) const {
      const auto& m = arcMeta_[arcId];
      return graph_[m.node][m.edgeIdx].to;
    }

    int numArcs() const { return static_cast<int>(arcMeta_.size()); }

  private:
    static constexpr int64_t INF = std::numeric_limits<int64_t>::max() / 2;

    struct Edge {
      int to, rev, cap;
      int64_t cost;
    };
    struct ArcMeta {
      int node, edgeIdx, origCap;
    };

    std::vector<std::vector<Edge>> graph_;
    std::vector<int64_t> potential_;
    std::vector<int> supply_;
    std::vector<ArcMeta> arcMeta_;
    int n_;
  };

}  // namespace ticl

#endif

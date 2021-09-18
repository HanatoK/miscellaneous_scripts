#ifndef GRAPH_H
#define GRAPH_H

#include <cstddef>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <functional>
#include <queue>
#include <utility>
#include <initializer_list>
#include <limits>
#include <chrono>

class Graph {
public:
  enum class FindPathMode {
    SumOfEdges,
    // MaximumEdges,
    // MFEPMode,
  };
  struct Node {
    size_t mIndex;
    double mWeight;
    // Image mImage;
  };
  struct Edge {
    size_t mSource;
    size_t mDestination;
    double mWeight;
  };
  struct FindPathResult {
    size_t mNumLoops;
    std::vector<bool> mVisitedNodes;
    std::vector<size_t> mPathNodes;
    std::vector<double> mDistances;
    void dump() const;
  };
  Graph(bool directed = false);
  Graph(const size_t num_nodes, bool directed = false);
  // special constructor for SMwST images
  // Graph(const std::vector<Image>& images);
  bool setEdge(size_t source, size_t destination, double weight = 1.0);
  bool setEdges(const std::vector<Edge> &edges);
  bool getEdge(size_t source, size_t destination, double &weight) const;
  void printGraph(std::ostream &os) const;
  void summary(std::ostream &os) const;
  size_t totalEdges() const;
  void DFS(size_t start, std::function<void(const Node &)> func) const;
  void DFS(size_t start, size_t end, std::function<void(const Node &)> func) const;
  FindPathResult Dijkstra(size_t start, size_t end, FindPathMode mode) const;
  template <typename DistanceType>
  FindPathResult Dijkstra(
      size_t start, size_t end, const DistanceType &dist_start,
      const DistanceType &dist_infinity,
      std::function<DistanceType(DistanceType, double)> calc_new_dist) const;
protected:
  size_t mNumNodes;
  bool mIsDirected;
  std::vector<std::deque<Node>> mHead;
  bool setEdgeHelper(size_t source, size_t destination, double weight = 1.0);
  void DFSHelper(size_t i, std::vector<bool> &visited,
                 std::function<void(const Node &)> func) const;
  void DFSHelper(size_t i, size_t j, std::vector<bool> &visited,
                 std::function<void(const Node &)> func) const;
};

template <typename DistanceType>
Graph::FindPathResult Graph::Dijkstra(
    size_t start, size_t end, const DistanceType &dist_start,
    const DistanceType &dist_infinity,
    std::function<DistanceType(DistanceType, double)> calc_new_dist) const {
  using std::make_pair;
  using std::priority_queue;
  typedef std::pair<DistanceType, size_t> DistNodePair;
  std::vector<bool> visited(mNumNodes, false);
  std::vector<size_t> previous(mNumNodes);
  priority_queue<DistNodePair, std::vector<DistNodePair>,
                 std::greater<DistNodePair>>
      pq;
  std::vector<DistanceType> distances(mNumNodes);
  for (size_t i = 0; i < mNumNodes; ++i) {
    distances[i] = (i == start) ? dist_start : dist_infinity;
    previous[i] = mNumNodes;
  }
  pq.push(make_pair(dist_start, start));
  size_t loop = 0;
  const auto start_time = std::chrono::high_resolution_clock::now();
  while (!pq.empty()) {
#ifdef DEBUG_DIJKSTRA
    std::cout << "Loop " << loop << " ============= \n";
    debug_priority_queue(pq, "Current priority search queue:");
    std::cout << "Distance from " << start << " : " << distances << "\n";
    std::cout << "Previous visited vertices array: " << previous << "\n";
    std::cout << "Visited vertices: " << visited << "\n";
#endif
    const size_t to_visit = pq.top().second;
#ifdef DEBUG_DIJKSTRA
    std::cout << "Visiting neighbor vertices of vertex " << to_visit << " :\n";
#endif
    pq.pop();
    auto neighbor_node = std::next(mHead[to_visit].cbegin(), 1);
    while (neighbor_node != mHead[to_visit].cend()) {
      const size_t neighbor_index = neighbor_node->mIndex;
      if (visited[neighbor_index] == false) {
#ifdef DEBUG_DIJKSTRA
        std::cout << "Neighbor vertex " << neighbor_index << " is not visited\n";
#endif
        const DistanceType new_distance =
            calc_new_dist(distances[to_visit], neighbor_node->mWeight);
#ifdef DEBUG_DIJKSTRA
        std::cout << "Current distance of vertex " << neighbor_index << " is "
                  << distances[neighbor_index] << "\n";
#endif
        if (new_distance < distances[neighbor_index]) {
#ifdef DEBUG_DIJKSTRA
          std::cout << "Distance is updated to " << new_distance << "\n";
#endif
          distances[neighbor_index] = new_distance;
          previous[neighbor_index] = to_visit;
          pq.push(make_pair(distances[neighbor_index], neighbor_index));
        }
      } else {
#ifdef DEBUG_DIJKSTRA
        std::cout << "Neighbor vertex " << neighbor_index
                  << " is already visited. Skip it\n";
#endif
      }
      std::advance(neighbor_node, 1);
    }
    if (to_visit == end) {
      break;
    }
    visited[to_visit] = true;
    ++loop;
#ifdef DEBUG_DIJKSTRA
    qDebug() << "===============================================";
#endif
  }
  const auto end_time = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed_time = end_time - start_time;
  std::cout << "Dijkstra's algorithm takes " << elapsed_time.count()
            << " milliseconds; total number of loops: " << loop << "\n";
  std::vector<size_t> path;
  size_t target = end;
  while (previous[target] != mNumNodes) {
    path.push_back(target);
    target = previous[target];
  }
  if (path.back() != end)
    path.push_back(start);
  std::reverse(path.begin(), path.end());
  std::vector<double> res_distance(distances.size());
  for (size_t i = 0; i < distances.size(); ++i) {
    res_distance[i] = static_cast<double>(distances[i]);
  }
  FindPathResult result{loop, visited, path, res_distance};
  return result;
}

#endif // GRAPH_H

#include "Graph.h"


Graph::Graph(bool directed) : mNumNodes(0), mIsDirected(directed), mHead(0) {}

Graph::Graph(const size_t numNodes, bool directed)
    : mNumNodes(numNodes), mIsDirected(directed), mHead(mNumNodes) {
  for (size_t i = 0; i < mNumNodes; ++i) {
    mHead[i].push_back(Node{i, 0.0});
  }
}

// Graph::Graph(const std::vector<Image>& images)
//   : mNumNodes(images.size()), mIsDirected(false), mHead(images.size()) {
//   for (size_t i = 0; i < mNumNodes; ++i) {
//     mHead[i].push_back(Node{i, 0.0, images[i]});
//   }
// }

bool Graph::setEdge(size_t source, size_t destination, double weight) {
  bool success = setEdgeHelper(source, destination, weight);
  if (!success)
    return success;
  if (mIsDirected == false) {
    success = success && setEdgeHelper(destination, source, weight);
  }
  return success;
}

bool Graph::setEdges(const std::vector<Edge> &edges) {
  for (size_t i = 0; i < edges.size(); ++i) {
    if (!setEdge(edges[i].mSource, edges[i].mDestination, edges[i].mWeight)) {
      return false;
    }
  }
  return true;
}

bool Graph::getEdge(size_t source, size_t destination, double &weight) const {
  if (source >= mNumNodes || destination >= mNumNodes || source == destination)
    return false;
  const auto &headList = mHead[source];
  auto it = headList.cbegin();
  while (it != headList.cend()) {
    if (it->mIndex == destination) {
      weight = it->mWeight;
      return true;
    }
    ++it;
  }
  weight = 0;
  return false;
}


void Graph::printGraph(std::ostream &os) const {
  for (size_t i = 0; i < mNumNodes; ++i) {
    const auto &current = mHead[i];
    os << "This node is " << current.cbegin()->mIndex;
    if (current.size() < 2)
      continue;
    os << ", linked to ";
    auto it_current = current.cbegin();
    std::advance(it_current, 1);
    while (it_current != current.cend()) {
      os << "(" << it_current->mIndex << " weight " << it_current->mWeight
         << ") ";
      std::advance(it_current, 1);
    }
    os << "\n";
  }
}


void Graph::summary(std::ostream &os) const {
  std::ios_base::fmtflags f(os.flags());
  os << "Summary of the graph: " << '\n';
  os << "Number of nodes: " << mHead.size() << '\n';
  os << "Number of edges: " << totalEdges() << '\n';
  os << std::boolalpha << "Is directed? " << mIsDirected << '\n';
  os.flags(f);
}

size_t Graph::totalEdges() const {
  size_t count = 0;
  for (size_t i = 0; i < mHead.size(); ++i) {
    for (size_t j = 1; j < mHead[i].size(); ++j) {
      ++count;
    }
  }
  return count;
}

void Graph::DFS(size_t start, std::function<void(const Node &)> func) const {
  std::vector<bool> visited(mNumNodes, false);
  DFSHelper(start, visited, func);
}

void Graph::DFS(size_t start, size_t end, std::function<void(const Node &)> func) const {
  std::vector<bool> visited(mNumNodes, false);
  DFSHelper(start, end, visited, func);
}

bool Graph::setEdgeHelper(size_t source, size_t destination, double weight) {
  if (source >= mNumNodes || destination >= mNumNodes || source == destination)
    return false;
  // access the source vertex
  auto &headList = mHead[source];
  // iterate over the list
  auto it = headList.begin();
  while (it != headList.end()) {
    if (it->mIndex == destination) {
      // if the edge is already in the list
      // just set the weight
      it->mWeight = weight;
      break;
    }
    ++it;
  }
  if (it == headList.end()) {
    // the destination vertex is not in the list
    // create a new edge
    headList.push_back(Node{destination, weight});
  }
  return true;
}

void Graph::DFSHelper(size_t i, std::vector<bool> &visited,
                      std::function<void(const Node &)> func) const {
  // mark this vertex as visited
  visited[i] = true;
  auto it = mHead[i].cbegin();
  func(*it);
  while (it != mHead[i].end()) {
    const size_t to_visit = it->mIndex;
    if (!visited[to_visit]) {
      // if this node is never visited, visit it
      DFSHelper(to_visit, visited, func);
    }
    ++it;
  }
}

void Graph::DFSHelper(size_t i, size_t j, std::vector<bool> &visited,
                      std::function<void(const Node &)> func) const {
  visited[i] = true;
  auto it = mHead[i].cbegin();
  func(*it);
  if (i == j) return;
  while (it != mHead[i].end()) {
    const size_t to_visit = it->mIndex;
    if (!visited[to_visit]) {
      // if this node is never visited, visit it
      DFSHelper(to_visit, j, visited, func);
    }
    ++it;
  }
}

Graph::FindPathResult Graph::Dijkstra(size_t start, size_t end,
                                      Graph::FindPathMode mode) const {
  switch (mode) {
  case Graph::FindPathMode::SumOfEdges: {
    const double dist_inf = std::numeric_limits<double>::max();
    return Dijkstra<double>(
        start, end, 0, dist_inf,
        [](const double &x, const double &y) { return x + y; });
    break;
  }
  default: {
    return FindPathResult();
    break;
  }
  }
}

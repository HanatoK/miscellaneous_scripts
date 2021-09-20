#include "Graph.h"
#include "Helper.h"
#include "Pathway.h"

#include <cmath>
#include <utility>
#include <map>
#include <fmt/format.h>

void read_pathway(std::istream& input_stream, std::vector<Image>& images) {
  std::vector<std::string> tmp_fields;
  size_t image_index = 0;
  if (!images.empty()) {
    image_index = images.back().mImageIndex + 1;
  }
  for (std::string line; std::getline(input_stream, line); ) {
    splitString(line, " ", tmp_fields);
    Image img;
    img.mImageIndex = image_index;
    for (const auto& elem: tmp_fields) {
      img.mPosition.push_back(std::stod(elem));
    }
    images.push_back(img);
    image_index = image_index  + 1;
    tmp_fields.clear();
  }
}

double image_distance(const Image& img1, const Image& img2) {
  double sum = 0;
  for (size_t i = 0; i < img1.mPosition.size(); ++i) {
    const double diff = img1.mPosition[i] - img2.mPosition[i];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

double average_distance(const std::vector<Image>& images) {
  if (images.size() < 2) return 0;
  double dist = 0;
  for (size_t i = 1; i < images.size(); ++i) {
    dist += image_distance(images[i], images[i-1]);
  }
  dist /= images.size() - 1;
  return dist;
}

// remove self loops
void remove_self_loops(std::vector<Image>& images, const double& distance_threshold) {
  bool has_self_loops = false;
  if (images.size() < 2) return;
  for (size_t i = 1; i < images.size(); ++i) {
    const double dist = image_distance(images[i], images[i-1]);
    if (dist < distance_threshold) {
      has_self_loops = true;
      images.erase(images.begin() + i);
      // break;
    }
  }
  if (has_self_loops) {
    remove_self_loops(images, distance_threshold);
  } else {
    return;
  }
  // return has_self_loops;
}

void print_pathway(const std::vector<Image>& images, std::ostream& os) {
  for (size_t i = 0; i < images.size(); ++i) {
    os << fmt::format("{:3d}", i);
    for (size_t j = 0; j < images[i].mPosition.size(); ++j) {
      os << fmt::format("{:15.7f}", images[i].mPosition[j]);
    }
    if (i > 0) {
      const double dist = image_distance(images[i], images[i-1]);
      os << fmt::format("{:15.7f}", dist);
    }
    os << "\n";
  }
}

std::tuple<Graph, std::vector<size_t>, size_t>
pathway_to_graph(const std::vector<Image>& images, const double& distance_threshold) {
  std::vector<size_t> nodes;
  std::map<size_t, size_t> image_to_node;
  image_to_node[0] = 0;
  std::vector<std::pair<size_t, size_t>> image_connections;
  // append the first image
  nodes.push_back(0);
  size_t previous_image_index = 0;
  size_t last_image_index = images.back().mImageIndex;
  for (size_t i = 1; i < images.size(); ++i) {
    bool add_image = true;
    size_t current_image_index = images[i].mImageIndex;
    for (size_t j = 0; j < nodes.size(); ++j) {
      const auto& exist_image = images[nodes[j]];
      const auto dist = image_distance(images[i], exist_image);
      if (dist < distance_threshold) {
        add_image = false;
        if (current_image_index == last_image_index) {
          last_image_index = nodes[j];
        }
        current_image_index = nodes[j];
        break;
      }
    }
    if (add_image) {
      nodes.push_back(current_image_index);
#ifdef DEBUG
      fmt::print("current_image_index = {} ; node_index = {}\n", current_image_index, nodes.size() - 1);
#endif
      image_to_node[current_image_index] = nodes.size() - 1;
    }
#ifdef DEBUG
    fmt::print("Edge from image {:02d} to {:02d}\n", previous_image_index, current_image_index);
#endif
    image_connections.push_back(std::make_pair(previous_image_index, current_image_index));
    previous_image_index = current_image_index;
  }
  Graph graph(nodes.size(), false);
  std::vector<Graph::Edge> edges;
  // image connections to edges
  for (size_t i = 0; i < nodes.size(); ++i) {
    const size_t source_image_index = nodes[i];
    // find source image index in connections
    for (size_t j = 0; j < image_connections.size(); ++j) {
      const size_t left = image_connections[j].first;
      const size_t right = image_connections[j].second;
      if (left == right) continue;
      if (source_image_index == left) {
#ifdef DEBUG
        fmt::print("i = {} ; right = {}\n", i, right);
#endif
        edges.push_back({i, image_to_node.at(right), 1.0});
      }
    }
  }
  graph.setEdges(edges);
#ifdef DEBUG
  graph.printGraph(std::cout);
  graph.summary(std::cout);
  fmt::print("Node to image:\n");
  for (size_t j = 0; j < nodes.size(); ++j) {
    const size_t image_index = nodes[j];
    fmt::print("Node {:02d} is image {:02d}, position = ", j, image_index);
    for (size_t k = 0; k < images[image_index].mPosition.size(); ++k) {
      fmt::print("{:15.7f}", (images[image_index].mPosition)[k]);
    }
    fmt::print("\n");
  }
#endif
  fmt::print("last_image_index = {:2d} \n", last_image_index);
  // fix the last index back to the origin endpoint
  if (last_image_index != images.back().mImageIndex) {
    size_t last_node_index = image_to_node.at(last_image_index);
    nodes[last_node_index] = images.back().mImageIndex;
    image_to_node[images.back().mImageIndex] = last_node_index;
  }
  return std::make_tuple(graph, nodes, image_to_node.at(last_image_index));
}

std::vector<Image> remove_loops_graph(const std::vector<Image>& images, const double& distance_threshold_factor) {
  const double avg_dist = average_distance(images);
  const auto g = pathway_to_graph(images, avg_dist * distance_threshold_factor);
  auto& graph = std::get<0>(g);
  const auto& node_image_map = std::get<1>(g);
  const auto& last_image_index = std::get<2>(g);
  const auto& path_find_result = graph.Dijkstra(0, last_image_index, Graph::FindPathMode::SumOfEdges);
  const auto& path_nodes = path_find_result.mPathNodes;
  std::vector<Image> pathway_new;
  for (size_t i = 0; i < path_nodes.size(); ++i) {
    const size_t image_index = node_image_map.at(path_nodes[i]);
    // fmt::print("Node {:2d}, origin image {:2d}\n", path_nodes[i], image_index);
    pathway_new.push_back(images[image_index]);
  }
  return pathway_new;
}

std::vector<Image> remove_loops_benoit(const std::vector<Image>& images) {
  std::vector<Image> results = images;
  bool has_loop = false;
  do {
    has_loop = false;
    for (size_t i = 1; i < results.size() - 1; ++i) {
      const auto neighbor_dist = image_distance(results[i-1], results[i]);
      size_t j = i + 1;
      for (; j < results.size(); ++j) {
        const auto dist = image_distance(results[i-1], results[j]);
        if (dist < neighbor_dist) {
          has_loop = true;
          break;
        }
      }
      if (has_loop) {
        results.erase(results.begin() + i, results.begin() + j);
        break;
      }
    }
  } while (has_loop);
  return results;
}
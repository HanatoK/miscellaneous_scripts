#include "Graph.h"
#include "Helper.h"
#include "Pathway.h"

#include <iostream>
#include <fmt/format.h>
#include <fmt/ostream.h>

int main(int argc, char* argv[]) {
  if (argc < 3) {
    fmt::print("Usage: {} pathfile_input pathfile_output\n", argv[0]);
    return 1;
  }
  std::vector<Image> pathway_origin;
  std::ifstream ifs(argv[1]);
  if (!ifs.is_open()) {
    std::cerr << "Cannot open input file " << argv[1] << "!\n";
    return 1;
  }
  read_pathway(ifs, pathway_origin);
  print_pathway(pathway_origin, std::cout);
  const double avg_dist = average_distance(pathway_origin);
  fmt::print("# Average distance = {:15.7f}\n", avg_dist);
  // remove_self_loops(pathway_origin, avg_dist * 0.5);
  // print_pathway(pathway_origin, std::cout);
  const auto g = pathway_to_graph(pathway_origin, avg_dist * 0.7);
  const auto graph = std::get<0>(g);
  const auto node_image_map = std::get<1>(g);
  const auto last_image_index = std::get<2>(g);
  const auto path_find_result = graph.Dijkstra(0, last_image_index, Graph::FindPathMode::SumOfEdges);
  const auto& path_nodes = path_find_result.mPathNodes;
  // construct the final pathway
  std::vector<Image> pathway_new;
  for (size_t i = 0; i < path_nodes.size(); ++i) {
    const size_t image_index = node_image_map.at(path_nodes[i]);
    fmt::print("Node {:2d}, origin image {:2d}\n", path_nodes[i], image_index);
    pathway_new.push_back(pathway_origin[image_index]);
  }
  print_pathway(pathway_new, std::cout);
  std::ofstream ofs(argv[2]);
  if (!ofs.is_open()) {
    std::cerr << "Cannot open output file " << argv[2] << "!\n";
    return 1;
  }
  for (size_t i = 0; i < pathway_new.size(); ++i) {
    const auto& position = pathway_new[i].mPosition;
    for (size_t j = 0; j < position.size(); ++j) {
      fmt::print(ofs, "{:15.7f}", position[j]);
    }
    fmt::print(ofs, "\n");
  }
  return 0;
}

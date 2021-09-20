#include "Graph.h"
#include "Helper.h"
#include "Pathway.h"

#include <iostream>
#include <fmt/format.h>
#include <fmt/ostream.h>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    fmt::print("Usage: {} pathfile_input pathfile_output_graph pathfile_output_simple\n", argv[0]);
    fmt::print("pathfile_output_graph and pathfile_output_simple are optional\n");
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
  // construct the final pathway
  const std::vector<Image> pathway_new = remove_loops_graph(pathway_origin, 0.7);
  print_pathway(pathway_new, std::cout);
  if (argc == 3) {
    std::ofstream ofs(argv[2]);
    if (ofs.is_open()) {
      for (size_t i = 0; i < pathway_new.size(); ++i) {
        const auto& position = pathway_new[i].mPosition;
        for (size_t j = 0; j < position.size(); ++j) {
          fmt::print(ofs, "{:15.7f}", position[j]);
        }
        fmt::print(ofs, "\n");
      }
    } else {
      std::cerr << "Cannot open output file " << argv[2] << "!\n";
    }
  }
  const std::vector<Image> pathway_benoit = remove_loops_benoit(pathway_origin);
  print_pathway(pathway_benoit, std::cout);
  if (argc == 4) {
    std::ofstream ofs(argv[3]);
    if (ofs.is_open()) {
      for (size_t i = 0; i < pathway_benoit.size(); ++i) {
        const auto& position = pathway_benoit[i].mPosition;
        for (size_t j = 0; j < position.size(); ++j) {
          fmt::print(ofs, "{:15.7f}", position[j]);
        }
        fmt::print(ofs, "\n");
      }
    } else {
      std::cerr << "Cannot open output file " << argv[3] << "!\n";
    }
  }
  return 0;
}

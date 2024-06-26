#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <execution>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#include <sys/types.h>

class cvc {
public:
  cvc(int x): data(x) {
    name = "cvc" + std::to_string(data);
  };
  void calc_value() {
    pid_t tid = gettid();
    usleep(1000);
    std::cout << "Calculate CVC " << name << " on thread id " << tid << std::endl;
  }
  std::vector<std::weak_ptr<cvc>> parents;
  std::vector<std::shared_ptr<cvc>> children;
  int data;
  std::string name;
};

std::vector<std::shared_ptr<cvc>> create_graph() {
  std::vector<std::shared_ptr<cvc>> results;
  auto root_node = std::make_shared<cvc>(cvc(100));
  root_node->children.push_back(std::make_shared<cvc>(cvc(101)));
  root_node->children.push_back(std::make_shared<cvc>(cvc(102)));
  for (auto it = root_node->children.begin(); it != root_node->children.end(); ++it) {
    (*it)->parents.push_back(root_node);
  }
  auto& tmp1 = root_node->children[0]->children;
  tmp1.push_back(std::make_shared<cvc>(cvc(103)));
  tmp1.push_back(std::make_shared<cvc>(cvc(104)));
  for (auto it = tmp1.begin(); it != tmp1.end(); ++it) {
    (*it)->parents.push_back(root_node->children[0]);
  }
  tmp1[1]->parents.push_back(root_node);
  root_node->children.push_back(tmp1[1]);
  auto& tmp2 = root_node->children[1]->children;
  tmp2.push_back(std::make_shared<cvc>(cvc(105)));
  tmp2.push_back(std::make_shared<cvc>(cvc(106)));
  tmp2.push_back(std::make_shared<cvc>(cvc(107)));
  tmp2.push_back(std::make_shared<cvc>(cvc(108)));
  for (auto it = tmp2.begin(); it != tmp2.end(); ++it) {
    (*it)->parents.push_back(root_node->children[1]);
  }
  tmp2[0]->parents.push_back(tmp1[1]);
  tmp1[1]->children.push_back(tmp2[0]);
  results.push_back(root_node);
  return results;
}

void node_max_depth(
  const std::vector<std::shared_ptr<cvc>>& nodes,
  std::unordered_map<std::shared_ptr<cvc>, int>& node_map,
  int current_level = 0) {
  for (auto it = nodes.begin(); it != nodes.end(); ++it) {
    auto map_it = node_map.find(*it);
    if (map_it == node_map.end()) {
      node_map.insert({*it, current_level});
    } else {
      if (map_it->second < current_level) {
        map_it->second = current_level;
      }
    }
    node_max_depth((*it)->children, node_map, current_level+1);
  }
}

void write_link(
  const std::vector<std::shared_ptr<cvc>>& nodes,
  const std::unordered_map<std::shared_ptr<cvc>, size_t>& node_to_index,
  std::ofstream& ofs) {
  for (auto it = nodes.begin(); it != nodes.end(); ++it) {
    const size_t parent_index = node_to_index.at(*it);
    ofs << "  node" << parent_index << " -> { ";
    auto children = (*it)->children;
    for (auto it_c = children.begin(); it_c != children.end(); ++it_c) {
      const size_t child_index = node_to_index.at(*it_c);
      ofs << "node" <<child_index << " ";
    }
    ofs << "} [label = \"reuse\" ]\n";
    write_link(children, node_to_index, ofs);
  }
}

void write_dot(
  const std::vector<std::shared_ptr<cvc>>& nodes,
  const std::unordered_map<std::shared_ptr<cvc>, int>& node_map,
  const std::string& filename) {
  std::ofstream ofs(filename.c_str());
  ofs << "digraph \" CVC graph dependency\" \n{\n";
  ofs << "  rankdir=\"LR\"\n{\n node [shape=box, style=\"rounded\"]";
  size_t node_index = 0;
  std::unordered_map<std::shared_ptr<cvc>, size_t> node_to_index;
  for (auto it = node_map.begin(); it != node_map.end(); ++it) {
    ofs << "    node" << node_index << " [label = \"" << it->first->name
        << "\\n" << "DEPTH: " << it->second << "\"]\n";
    node_to_index.insert({it->first, node_index});
    ++node_index;
  }
  ofs << "}\n";
  write_link(nodes, node_to_index, ofs);
  ofs << "}\n";
}

auto process_cvcs(const std::vector<std::shared_ptr<cvc>>& nodes) {
  std::unordered_map<std::shared_ptr<cvc>, int> node_map;
  // Get the max depth of each node
  node_max_depth(nodes, node_map, 0);
  int max_depth = 0;
  // Find the max depth
  for (auto it = node_map.begin(); it != node_map.end(); ++it) {
    if (it->second > max_depth) max_depth = it->second;
  }
  // Sort the nodes according to node depths
  std::vector<std::vector<std::shared_ptr<cvc>>> sorted_nodes(max_depth+1);
  for (auto it = node_map.begin(); it != node_map.end(); ++it) {
    sorted_nodes[it->second].push_back(it->first);
  }
  // Sort the nodes by name within the same depth
  for (auto it = sorted_nodes.begin(); it != sorted_nodes.end(); ++it) {
    std::sort(it->begin(), it->end(),
        [](const auto& a, const auto& b){
                return a->name < b->name;});
  }
  // Run nodes in each depth for parallel
  int current_depth = max_depth;
  for (auto it = sorted_nodes.rbegin(); it != sorted_nodes.rend(); ++it) {
    std::cout << "Depth " << current_depth-- << ":\n";
    // Mimic smp_colvars_loop
    std::for_each(
      std::execution::par_unseq, it->begin(), it->end(),
      [](auto& node){node->calc_value();});
    std::cout << std::endl;
  }
  return node_map;
}

int main() {
  auto nodes = create_graph();
  auto node_map = process_cvcs(nodes);
  // dot -Tpng test.dot -Gdpi=300 -o test.png
  write_dot(nodes, node_map, "test.dot");
  return 0;
}

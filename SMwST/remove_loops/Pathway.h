#ifndef PATHWAY_H
#define PATHWAY_H

#include <vector>
#include <iostream>
#include <tuple>

#include "Graph.h"

struct Image {
  size_t mImageIndex;
  std::vector<double> mPosition;
};

void read_pathway(std::istream& input_stream, std::vector<Image>& images);
double image_distance(const Image& img1, const Image& img2);
double average_distance(const std::vector<Image>& images);
void remove_self_loops(std::vector<Image>& images, const double& distance_threshold);
void print_pathway(const std::vector<Image>& images, std::ostream& os);

/**
 * @param images the images of the pathway from A to B
 * @param distance_threshold a threshold for checking if two images belong the same node
 * @return a graph constructed from the pathway
 * @return a vector indexed by the nodes contains the index of images
 * @return the node index of the last image
 */
std::tuple<Graph, std::vector<size_t>, size_t>
pathway_to_graph(const std::vector<Image>& images, const double& distance_threshold);

std::vector<Image> remove_loops_graph(const std::vector<Image>& images, const double& distance_threshold_factor);
std::vector<Image> remove_loops_simple(const std::vector<Image>& images);

#endif // PATHWAY_H

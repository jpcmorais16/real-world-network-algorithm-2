#pragma once

#include <vector>
#include <unordered_set>
#include <tuple>

std::vector<std::vector<int>> convertToFastGraph(std::vector<std::unordered_set<int>>& graph);

std::tuple<double, double, int> calculatePowerLaw(const std::vector<std::vector<int>>& fastGraph);
double calculateClusteringCoefficient(const std::vector<std::vector<int>>& fastGraph);
double calculateShortestPathLength(const std::vector<std::vector<int>>& fastGraph);


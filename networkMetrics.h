#pragma once

#include <vector>
#include <unordered_set>
#include <tuple>
#include <chrono>

struct NetworkMetrics {
    double avgPathLength;
    double clusteringCoefficient;
    double powerLawCoefficient;
    double pearsonR;
    int maxDegree;
    std::chrono::milliseconds powerLawDuration;
    std::chrono::milliseconds avgPathLengthDuration;
    std::chrono::milliseconds clusteringCoefficientDuration;
    size_t numThreads; 
};

std::vector<std::vector<int>> convertToFastGraph(std::vector<std::unordered_set<int>>& graph);

std::tuple<double, double, int> calculatePowerLaw(const std::vector<std::vector<int>>& fastGraph);
double calculateClusteringCoefficient(const std::vector<std::vector<int>>& fastGraph);
double calculateShortestPathLength(const std::vector<std::vector<int>>& fastGraph);


NetworkMetrics calculateAllMetrics(std::vector<std::unordered_set<int>>& graph);


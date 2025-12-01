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
    double powerLawCumulativeCoefficient;
    double pearsonRCumulative;
    int maxDegree;
    size_t numEdges;
    std::chrono::milliseconds powerLawDuration;
    std::chrono::milliseconds avgPathLengthDuration;
    std::chrono::milliseconds clusteringCoefficientDuration;
    size_t numThreads; 
};

// Network generation function (defined in main.cpp)
std::vector<std::unordered_set<int>> generate_network(int N, double m, double p, double fp);

// Network saving functions
void save_network_to_file(const std::vector<std::unordered_set<int>>& graph, const std::string& filename);
std::vector<std::unordered_set<int>> generate_network_with_save(int N, double m, double p, double fp, const std::string& filename = "");

std::vector<std::vector<int>> convertToFastGraph(std::vector<std::unordered_set<int>>& graph);

std::tuple<double, double, int> calculatePowerLaw(const std::vector<std::vector<int>>& fastGraph);
std::tuple<double, double, int> calculatePowerLawCumulative(const std::vector<std::vector<int>>& fastGraph);
double calculateClusteringCoefficient(const std::vector<std::vector<int>>& fastGraph);
double calculateShortestPathLength(const std::vector<std::vector<int>>& fastGraph);

// New functions for different approaches to calculate approximate average path length
double calculateApproximateAveragePathLengthDegreeMethod(const std::vector<std::vector<int>>& fastGraph);
double calculateApproximateAveragePathLengthSamplingMethod(const std::vector<std::vector<int>>& fastGraph, size_t sampleSize = 10000);
double calculateApproximateAveragePathLengthNodeSamplingMethod(const std::vector<std::vector<int>>& fastGraph, size_t sampleSize = 10000);

NetworkMetrics calculateAllMetrics(std::vector<std::unordered_set<int>>& graph);

std::tuple<double, double, double, double> compareApproximateAndExactAveragePathLength(const std::vector<std::vector<int>>& fastGraph);

// New function to compare all methods (exact vs three approximate methods)
struct MethodComparisonResult {
    double exactResult;
    double exactTime;
    double degreeMethodResult;
    double degreeMethodTime;
    double pairSamplingResult;
    double pairSamplingTime;
    double nodeSamplingResult;
    double nodeSamplingTime;
    double degreeMethodError;
    double pairSamplingError;
    double nodeSamplingError;
    size_t numNodes;
    size_t numEdges;
};

MethodComparisonResult compareAllMethods(int N, double p, double m, double friendship_probability, size_t pairSampleSize = 10000, size_t nodeSampleSize = 10000);


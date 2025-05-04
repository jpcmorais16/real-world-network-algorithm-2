#include <vector>
#include <unordered_set>
#include <queue>
#include <numeric>
#include <thread>
#include <atomic>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <tuple>
#include <chrono>

#include "networkMetrics.h"

using namespace std;

std::tuple<double, double, int> calculatePowerLaw(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    
    size_t maxDegree = 0;
    
    for (const auto& neighbors : fastGraph) {
        size_t degree = neighbors.size();
        maxDegree = std::max(maxDegree, degree);
    }
    
    vector<int> degreeFrequency(maxDegree + 1, 0);
    
    for (const auto& neighbors : fastGraph) {
        size_t degree = neighbors.size();
        degreeFrequency[degree]++;
    }
    
    vector<double> logDegree;
    vector<double> logProbability;
    
    for (size_t k = 1; k < degreeFrequency.size(); k++) {
        if (degreeFrequency[k] > 0) {
            double probability = static_cast<double>(degreeFrequency[k]) / n;
            logDegree.push_back(log(static_cast<double>(k)));
            logProbability.push_back(log(probability));
        }
    }
    
    if (logDegree.size() < 2) {
        return {0.0, 0.0, static_cast<int>(maxDegree)};
    }
    
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
    int m = logDegree.size();
    
    for (int i = 0; i < m; i++) {
        sumX += logDegree[i];
        sumY += logProbability[i];
        sumXY += logDegree[i] * logProbability[i];
        sumX2 += logDegree[i] * logDegree[i];
        sumY2 += logProbability[i] * logProbability[i];
    }
    
    double slope = (m * sumXY - sumX * sumY) / (m * sumX2 - sumX * sumX);
    
    double r = (m * sumXY - sumX * sumY) / 
               sqrt((m * sumX2 - sumX * sumX) * (m * sumY2 - sumY * sumY));
    
    return {-slope, r, static_cast<int>(maxDegree)};
}

std::vector<std::vector<int>> convertToFastGraph(std::vector<std::unordered_set<int>>& graph) {
    const size_t n = graph.size();
    std::vector<std::vector<int>> fastGraph(n);
    
    for (size_t i = 0; i < n; ++i) {
        const size_t degree = graph[i].size();
        fastGraph[i].reserve(degree);
        fastGraph[i].assign(graph[i].begin(), graph[i].end());
    }
    
    return fastGraph;
}

double calculateClusteringCoefficient(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    double totalSum = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const size_t degree = fastGraph[i].size();
        if (degree < 2) continue;

        unordered_set<int> neighborSet(fastGraph[i].begin(), fastGraph[i].end());
        
        int connections = 0;
        for (size_t j = 0; j < fastGraph[i].size(); ++j) {
            int neighbor1 = fastGraph[i][j];
            for (size_t k = j + 1; k < fastGraph[i].size(); ++k) {
                int neighbor2 = fastGraph[i][k];
                
                for (int adj : fastGraph[neighbor1]) {
                    if (adj == neighbor2) {
                        connections++;
                        break;
                    }
                }
            }
        }
        
        int maxPossibleConnections = degree * (degree - 1) / 2;
        
        double localCoefficient = connections / static_cast<double>(maxPossibleConnections);
        totalSum += localCoefficient;
    }

    return totalSum / n;
}

double calculateShortestPathLength(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    double totalSum = 0.0;

    const size_t maxThreads = std::thread::hardware_concurrency();
    const size_t numThreads = std::min(maxThreads ? maxThreads : 2, n);
    
    std::cout << "Using " << numThreads << " threads for path length calculation" << std::endl;
    
    std::vector<double> threadResults(numThreads, 0.0);
    std::vector<std::thread> threads;
    
    std::atomic<size_t> nextNode(0);
    
    auto processNodes = [&](size_t threadId) {
        double threadSum = 0.0;
        
        std::vector<char> visited(n, 0);
        std::vector<int> distances(n);
        std::queue<int> q;
        
        while (true) {
            size_t source = nextNode.fetch_add(1);
            if (source >= n) break;
            
            std::fill(visited.begin(), visited.end(), 0);
            std::fill(distances.begin(), distances.end(), 0);
            
            while (!q.empty()) q.pop();
            
            visited[source] = 1;
            q.push(source);
            
            while (!q.empty()) {
                int current = q.front();
                q.pop();
                
                const std::vector<int>& neighbors = fastGraph[current];
                const int currentDist = distances[current];
                
                for (int neighbor : neighbors) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = 1;
                        distances[neighbor] = currentDist + 1;
                        q.push(neighbor);
                    }
                }
            }
            
            int distanceSum = 0;
            int validPaths = 0;
            for (size_t i = 0; i < n; ++i) {
                if (i != source && visited[i]) {
                    distanceSum += distances[i];
                    validPaths++;
                }
            }
            
            if (validPaths > 0) {
                threadSum += static_cast<double>(distanceSum) / validPaths;
            }
        }
        
        threadResults[threadId] = threadSum;
    };
    
    for (size_t i = 0; i < numThreads; ++i) {
        threads.emplace_back(processNodes, i);
    }
    
    for (auto& thread : threads) {
        thread.join();
    }
    
    totalSum = std::accumulate(threadResults.begin(), threadResults.end(), 0.0);
    
    return totalSum / n;
}

NetworkMetrics calculateAllMetrics(std::vector<std::unordered_set<int>>& graph) {
    NetworkMetrics metrics;
    
    auto fastGraph = convertToFastGraph(graph);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    auto power_law_result = calculatePowerLaw(fastGraph);
    metrics.powerLawCoefficient = std::get<0>(power_law_result);
    metrics.pearsonR = std::get<1>(power_law_result);
    metrics.maxDegree = std::get<2>(power_law_result);
    auto end_time = std::chrono::high_resolution_clock::now();
    metrics.powerLawDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    const size_t maxThreads = std::thread::hardware_concurrency();
    const size_t n = fastGraph.size();
    metrics.numThreads = std::min(maxThreads ? maxThreads : 2, n);
    
    start_time = std::chrono::high_resolution_clock::now();
    metrics.avgPathLength = calculateShortestPathLength(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    metrics.avgPathLengthDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    metrics.clusteringCoefficient = calculateClusteringCoefficient(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    metrics.clusteringCoefficientDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    return metrics;
}

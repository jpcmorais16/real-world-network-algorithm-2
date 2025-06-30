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
#include <climits>

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
    // const size_t maxThreads = 5;
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

std::priority_queue<std::pair<int, int>> rankNodesByDegree(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    
    std::priority_queue<std::pair<int, int>> pq;
    
    for (size_t i = 0; i < n; ++i) {
        int degree = fastGraph[i].size();
        pq.push({degree, static_cast<int>(i)});
    }
    
    return pq;
}

std::vector<int> selectLandmarkNodes(const vector<vector<int>>& fastGraph, const size_t numLandmarks) {
    auto rankedNodes = rankNodesByDegree(fastGraph);
    std::vector<int> landmarkNodes;
    landmarkNodes.reserve(numLandmarks);
    
    std::unordered_set<int> prohibitedNodes;
    
    while (!rankedNodes.empty() && landmarkNodes.size() < numLandmarks) {
        auto top = rankedNodes.top();
        int nodeIndex = top.second;
        rankedNodes.pop();
        
        if (prohibitedNodes.count(nodeIndex) > 0) {
            continue;
        }
        
        landmarkNodes.push_back(nodeIndex);
        
        for (int neighbor : fastGraph[nodeIndex]) {
            prohibitedNodes.insert(neighbor);
        }
    }
    
    return landmarkNodes;
}

double calculateAproximateAveragePathLength(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    auto landmarks = selectLandmarkNodes(fastGraph, 100);
    cout << "Landmarks: " << landmarks.size() << endl;
    
    // Create result matrix: result[i] contains distances from node i to all landmarks
    std::vector<std::vector<int>> result(n, std::vector<int>(landmarks.size()));
    
    // For each landmark, run BFS to calculate distances to all nodes
    for (size_t landmarkIdx = 0; landmarkIdx < landmarks.size(); ++landmarkIdx) {
        int landmark = landmarks[landmarkIdx];
        
        // BFS setup
        std::vector<char> visited(n, 0);
        std::vector<int> distances(n);
        std::queue<int> q;
        
        // Start BFS from landmark
        visited[landmark] = 1;
        distances[landmark] = 0;
        q.push(landmark);
        
        while (!q.empty()) {
            int current = q.front();
            q.pop();
            
            // For each neighbor of current node
            for (int neighbor : fastGraph[current]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = 1;
                    distances[neighbor] = distances[current] + 1;
                    q.push(neighbor);
                }
            }
        }
        
        // Store distances in result matrix
        for (size_t i = 0; i < n; ++i) {
            result[i][landmarkIdx] = distances[i];
        }
    }
    
    cout << "BFS done" << endl;
    
    // Calculate approximate average path length using multiple threads
    const size_t maxThreads = std::thread::hardware_concurrency();
    const size_t numThreads = std::min(maxThreads ? maxThreads : 2, n);
    
    std::cout << "Using " << numThreads << " threads for path length calculation" << std::endl;
    
    std::vector<double> threadMinSumResults(numThreads, 0.0);
    std::vector<double> threadMaxDiffResults(numThreads, 0.0);
    std::vector<std::thread> threads;
    std::atomic<size_t> nextNode(0);
    
    auto processNodes = [&](size_t threadId) {
        double threadMinSum = 0.0;
        double threadMaxDiff = 0.0;
        
        while (true) {
            size_t s = nextNode.fetch_add(1);
            if (s >= n) break;
            
            for (size_t t = s + 1; t < n; ++t) {
                int minPathLength = INT_MAX;
                int maxAbsDiff = 0;
                
                // Find minimum s->landmark->t path length and maximum absolute difference
                for (size_t i = 0; i < landmarks.size(); ++i) {
                    int pathLength = result[s][i] + result[t][i];
                    minPathLength = std::min(minPathLength, pathLength);
                    
                    int absDiff = std::abs(result[s][i] - result[t][i]);
                    maxAbsDiff = std::max(maxAbsDiff, absDiff);
                }
                
                threadMinSum += minPathLength;
                threadMaxDiff += maxAbsDiff;
            }
        }
        
        threadMinSumResults[threadId] = threadMinSum;
        threadMaxDiffResults[threadId] = threadMaxDiff;
    };
    
    for (size_t i = 0; i < numThreads; ++i) {
        threads.emplace_back(processNodes, i);
    }
    
    for (auto& thread : threads) {
        thread.join();
    }
    
    double totalMinSum = std::accumulate(threadMinSumResults.begin(), threadMinSumResults.end(), 0.0);
    double totalMaxDiff = std::accumulate(threadMaxDiffResults.begin(), threadMaxDiffResults.end(), 0.0);
    
    double avgMinSum = totalMinSum / (n * (n - 1) / 2);
    double avgMaxDiff = totalMaxDiff / (n * (n - 1) / 2);
    
    cout << "Average minimum sum path length: " << avgMinSum << endl;
    cout << "Average maximum absolute difference: " << avgMaxDiff << endl;
    
    return avgMinSum;
}

NetworkMetrics calculateAllMetrics(std::vector<std::unordered_set<int>>& graph) {
    NetworkMetrics metrics;
    
    auto fastGraph = convertToFastGraph(graph);
    
    // Count number of edges
    size_t totalDegree = 0;
    for (const auto& neighbors : fastGraph) {
        totalDegree += neighbors.size();
    }
    metrics.numEdges = totalDegree / 2;  // Each edge is counted twice
    
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
    // metrics.avgPathLength = calculateShortestPathLength(fastGraph);
    metrics.avgPathLength = calculateAproximateAveragePathLength(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    metrics.avgPathLengthDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    // metrics.clusteringCoefficient = calculateClusteringCoefficient(fastGraph);
    metrics.clusteringCoefficient = 0.0;
    end_time = std::chrono::high_resolution_clock::now();
    metrics.clusteringCoefficientDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    return metrics;
}

pair<double, double> compareApproximateAndExactAveragePathLength(const vector<vector<int>>& fastGraph) {
    
    double exactAveragePathLength = calculateShortestPathLength(fastGraph);
    double approximateAveragePathLength = calculateAproximateAveragePathLength(fastGraph);
    
    return {exactAveragePathLength, approximateAveragePathLength};
}

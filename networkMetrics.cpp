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
#include <iomanip>
#include <random>
#include <algorithm>

#include "networkMetrics.h"

using namespace std;

std::tuple<double, double, int> calculatePowerLawCumulative(const vector<vector<int>>& fastGraph) {
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
    
    // Calculate complementary cumulative distribution P(X >= k)
    vector<double> complementaryCumulativeProbability(maxDegree + 1, 0.0);
    double runningSum = 0.0;
    
    // Start from maxDegree and work backwards
    for (int k = maxDegree; k >= 1; k--) {
        runningSum += static_cast<double>(degreeFrequency[k]) / n;
        complementaryCumulativeProbability[k] = runningSum;
    }

    
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
    int m = 0;
    
    const double MIN_PROBABILITY = 1e-10; // Minimum probability threshold
    
    for (size_t k = 1; k < degreeFrequency.size(); k++) {
        if (complementaryCumulativeProbability[k] > MIN_PROBABILITY) {
            double logK = log(static_cast<double>(k));
            double logProb = log(complementaryCumulativeProbability[k]);
            sumX += logK;
            sumY += logProb;
            sumXY += logK * logProb;
            sumX2 += logK * logK;
            sumY2 += logProb * logProb;
            m++;
        } else {
            // If complementaryCumulativeProbability[k] is too small, we can break early
            // since all subsequent values will also be very small (due to the cumulative nature)
            break;
        }
        
        // Additional safety check: if we've processed too many iterations, break
        if (k > 7500) {
            break;
        }
    }
    
    if (m < 2) {
        return {0.0, 0.0, static_cast<int>(maxDegree)};
    }
    
    double slope = (m * sumXY - sumX * sumY) / (m * sumX2 - sumX * sumX);
    
    double r = (m * sumXY - sumX * sumY) / 
               sqrt((m * sumX2 - sumX * sumX) * (m * sumY2 - sumY * sumY));
    
    // For complementary cumulative distribution: slope = -α, so α = 1-slope
    double alpha = 1-slope;

    return {alpha, r, static_cast<int>(maxDegree)};
}

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
    // Use the sampling method by default
    return calculateApproximateAveragePathLengthSamplingMethod(fastGraph, 10000);
}

double calculateApproximateAveragePathLengthDegreeMethod(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    auto landmarks = selectLandmarkNodes(fastGraph, 200);
    cout << "Landmarks: " << landmarks.size() << endl;
    vector<double> distance_sums(landmarks.size(), 0.0);
    
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
            distance_sums[landmarkIdx] += distances[i];
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
    
    std::atomic<size_t> processedPairs(0);
    const size_t totalPairs = n * (n - 1) / 2;
    const size_t progressInterval = std::max(static_cast<size_t>(1), totalPairs / 100); // Print every 1%
    
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
                
                // Progress tracking
                size_t currentPairs = processedPairs.fetch_add(1) + 1;
                if (currentPairs % progressInterval == 0) {
                    double progress = (static_cast<double>(currentPairs) / totalPairs) * 100.0;
                    std::cout << "Progress: " << std::fixed << std::setprecision(1) 
                              << progress << "% (" << currentPairs << "/" << totalPairs 
                              << " pairs processed)" << std::endl;
                }
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

double calculateApproximateAveragePathLengthSamplingMethod(const vector<vector<int>>& fastGraph, size_t sampleSize) {
    const size_t n = fastGraph.size();
    auto landmarks = selectLandmarkNodes(fastGraph, 200);
    cout << "Landmarks: " << landmarks.size() << endl;
    vector<double> distance_sums(landmarks.size(), 0.0);
    
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
            distance_sums[landmarkIdx] += distances[i];
        }
    }
    
    cout << "BFS done" << endl;
    
    // Sample pairs to calculate average distance
    cout << "Using sampling approach with " << sampleSize << " pairs" << endl;
    
    // Generate random pairs for sampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> nodeDist(0, n - 1);
    
    double totalSampledDistance = 0.0;
    size_t validPairs = 0;
    
    for (size_t sample = 0; sample < sampleSize; ++sample) {
        size_t s = nodeDist(gen);
        size_t t = nodeDist(gen);
        
        // Ensure s != t
        while (s == t) {
            t = nodeDist(gen);
        }
        
        // Find minimum path length through landmarks
        int minPathLength = INT_MAX;
        for (size_t i = 0; i < landmarks.size(); ++i) {
            int pathLength = result[s][i] + result[t][i];
            minPathLength = std::min(minPathLength, pathLength);
        }
        
        // Only count if we found a valid path
        if (minPathLength != INT_MAX) {
            totalSampledDistance += minPathLength;
            validPairs++;
        }
    }
    
    double sampledAveragePathLength = 0.0;
    if (validPairs > 0) {
        sampledAveragePathLength = totalSampledDistance / validPairs;
        cout << "Sampled average path length: " << sampledAveragePathLength << endl;
        cout << "Valid pairs in sample: " << validPairs << "/" << sampleSize << endl;
    }
    
    return sampledAveragePathLength;
}

double calculateApproximateAveragePathLengthNodeSamplingMethod(const vector<vector<int>>& fastGraph, size_t sampleSize) {
    const size_t n = fastGraph.size();
    
    // Limit sample size to the number of nodes
    sampleSize = std::min(sampleSize, n);
    
    cout << "Using node sampling approach with " << sampleSize << " source nodes" << endl;
    
    // Generate random source nodes for sampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> nodeDist(0, n - 1);
    
    std::vector<size_t> sourceNodes;
    sourceNodes.reserve(sampleSize);
    
    // Sample source nodes (avoid duplicates)
    std::unordered_set<size_t> sampledSet;
    while (sourceNodes.size() < sampleSize) {
        size_t node = nodeDist(gen);
        if (sampledSet.insert(node).second) {
            sourceNodes.push_back(node);
        }
    }
    
    double totalAveragePathLength = 0.0;
    size_t validSources = 0;
    
    // For each sampled source node, run BFS and calculate average path length
    for (size_t sourceIdx = 0; sourceIdx < sourceNodes.size(); ++sourceIdx) {
        size_t source = sourceNodes[sourceIdx];
        
        // BFS setup
        std::vector<char> visited(n, 0);
        std::vector<int> distances(n);
        std::queue<int> q;
        
        // Start BFS from source
        visited[source] = 1;
        distances[source] = 0;
        q.push(source);
        
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
        
        // Calculate average path length from this source to all other nodes
        int distanceSum = 0;
        int validPaths = 0;
        
        for (size_t i = 0; i < n; ++i) {
            if (i != source && visited[i]) {
                distanceSum += distances[i];
                validPaths++;
            }
        }
        
        // Only count if we found valid paths
        if (validPaths > 0) {
            double avgPathFromSource = static_cast<double>(distanceSum) / validPaths;
            totalAveragePathLength += avgPathFromSource;
            validSources++;
        }
    }
    
    double finalAveragePathLength = 0.0;
    if (validSources > 0) {
        finalAveragePathLength = totalAveragePathLength / validSources;
        cout << "Node sampling average path length: " << finalAveragePathLength << endl;
        cout << "Valid source nodes: " << validSources << "/" << sampleSize << endl;
    }
    
    return finalAveragePathLength;
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
    
    auto power_law_cumulative_result = calculatePowerLawCumulative(fastGraph);
    metrics.powerLawCumulativeCoefficient = std::get<0>(power_law_cumulative_result);
    metrics.pearsonRCumulative = std::get<1>(power_law_cumulative_result);
    auto end_time = std::chrono::high_resolution_clock::now();
    metrics.powerLawDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    const size_t maxThreads = std::thread::hardware_concurrency();
    const size_t n = fastGraph.size();
    metrics.numThreads = std::min(maxThreads ? maxThreads : 2, n);
    
    start_time = std::chrono::high_resolution_clock::now();
    metrics.avgPathLength = 0.0;
    metrics.avgPathLength = calculateShortestPathLength(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    metrics.avgPathLengthDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    metrics.clusteringCoefficient = calculateClusteringCoefficient(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    metrics.clusteringCoefficientDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    return metrics;
}

std::tuple<double, double, double, double> compareApproximateAndExactAveragePathLength(const vector<vector<int>>& fastGraph) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    double exactAveragePathLength = calculateShortestPathLength(fastGraph);
    auto end_time = std::chrono::high_resolution_clock::now();
    double exactTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    
    start_time = std::chrono::high_resolution_clock::now();
    double approximateAveragePathLength = calculateAproximateAveragePathLength(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    double approximateTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    
    return {exactAveragePathLength, approximateAveragePathLength, exactTime, approximateTime};
}

MethodComparisonResult compareAllMethods(int N, double p, double m, double friendship_probability, size_t pairSampleSize, size_t nodeSampleSize) {
    MethodComparisonResult result;
    
    cout << "\n=== Generating Network ===" << endl;
    cout << "Parameters: N=" << N << ", p=" << p << ", m=" << m << ", fp=" << friendship_probability << endl;
    
    // Generate the network
    auto graph = generate_network(N, m, p, friendship_probability);
    auto fastGraph = convertToFastGraph(graph);
    
    const size_t n = fastGraph.size();
    
    // Count number of edges
    size_t totalDegree = 0;
    for (const auto& neighbors : fastGraph) {
        totalDegree += neighbors.size();
    }
    result.numNodes = n;
    result.numEdges = totalDegree / 2;
    
    cout << "Network generated: " << n << " nodes, " << result.numEdges << " edges" << endl;
    cout << "\n=== Comparing All Methods ===" << endl;
    
    // 1. Exact Method
    cout << "\n--- Exact Method ---" << endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    result.exactResult = calculateShortestPathLength(fastGraph);
    auto end_time = std::chrono::high_resolution_clock::now();
    result.exactTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    cout << "Exact result: " << result.exactResult << " (took " << result.exactTime << "s)" << endl;
    
    // 2. Degree Method (All Pairs)
    cout << "\n--- Degree Method (All Pairs) ---" << endl;
    start_time = std::chrono::high_resolution_clock::now();
    result.degreeMethodResult = calculateApproximateAveragePathLengthDegreeMethod(fastGraph);
    end_time = std::chrono::high_resolution_clock::now();
    result.degreeMethodTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    result.degreeMethodError = std::abs(result.exactResult - result.degreeMethodResult) / result.exactResult;
    cout << "Degree method result: " << result.degreeMethodResult << " (took " << result.degreeMethodTime << "s)" << endl;
    cout << "Relative error: " << (result.degreeMethodError * 100) << "%" << endl;
    
    // 3. Pair Sampling Method
    cout << "\n--- Pair Sampling Method ---" << endl;
    start_time = std::chrono::high_resolution_clock::now();
    result.pairSamplingResult = calculateApproximateAveragePathLengthSamplingMethod(fastGraph, pairSampleSize);
    end_time = std::chrono::high_resolution_clock::now();
    result.pairSamplingTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    result.pairSamplingError = std::abs(result.exactResult - result.pairSamplingResult) / result.exactResult;
    cout << "Pair sampling result: " << result.pairSamplingResult << " (took " << result.pairSamplingTime << "s)" << endl;
    cout << "Relative error: " << (result.pairSamplingError * 100) << "%" << endl;
    
    // 4. Node Sampling Method
    cout << "\n--- Node Sampling Method ---" << endl;
    start_time = std::chrono::high_resolution_clock::now();
    result.nodeSamplingResult = calculateApproximateAveragePathLengthNodeSamplingMethod(fastGraph, nodeSampleSize);
    end_time = std::chrono::high_resolution_clock::now();
    result.nodeSamplingTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    result.nodeSamplingError = std::abs(result.exactResult - result.nodeSamplingResult) / result.exactResult;
    cout << "Node sampling result: " << result.nodeSamplingResult << " (took " << result.nodeSamplingTime << "s)" << endl;
    cout << "Relative error: " << (result.nodeSamplingError * 100) << "%" << endl;
    
    // Summary
    cout << "\n=== Summary ===" << endl;
    cout << "Method\t\t\tResult\t\tTime(s)\t\tError(%)" << endl;
    cout << "Exact\t\t\t" << std::fixed << std::setprecision(4) << result.exactResult << "\t\t" << result.exactTime << "\t\t0.00" << endl;
    cout << "Degree Method\t\t" << result.degreeMethodResult << "\t\t" << result.degreeMethodTime << "\t\t" << (result.degreeMethodError * 100) << endl;
    cout << "Pair Sampling\t\t" << result.pairSamplingResult << "\t\t" << result.pairSamplingTime << "\t\t" << (result.pairSamplingError * 100) << endl;
    cout << "Node Sampling\t\t" << result.nodeSamplingResult << "\t\t" << result.nodeSamplingTime << "\t\t" << (result.nodeSamplingError * 100) << endl;
    
    // Speedup analysis
    cout << "\n=== Speedup Analysis ===" << endl;
    cout << "Degree method speedup: " << (result.exactTime / result.degreeMethodTime) << "x" << endl;
    cout << "Pair sampling speedup: " << (result.exactTime / result.pairSamplingTime) << "x" << endl;
    cout << "Node sampling speedup: " << (result.exactTime / result.nodeSamplingTime) << "x" << endl;
    
    // Write results to file
    std::ofstream file("method_comparison_results_2.txt", std::ios::app);
    if (file.is_open()) {
        file << std::fixed << std::setprecision(6);
        file << "\n=== Method Comparison Results ===" << endl;
        file << "Timestamp: " << std::chrono::system_clock::now().time_since_epoch().count() << endl;
        file << "Network Parameters: N=" << N << ", p=" << p << ", m=" << m << ", fp=" << friendship_probability << endl;
        file << "Network: " << result.numNodes << " nodes, " << result.numEdges << " edges" << endl;
        file << "Sample Sizes: pair=" << pairSampleSize << ", node=" << nodeSampleSize << endl;
        file << endl;
        
        file << "Method\t\t\tResult\t\tTime(s)\t\tError(%)" << endl;
        file << "Exact\t\t\t" << result.exactResult << "\t\t" << result.exactTime << "\t\t0.00" << endl;
        file << "Degree Method\t\t" << result.degreeMethodResult << "\t\t" << result.degreeMethodTime << "\t\t" << (result.degreeMethodError * 100) << endl;
        file << "Pair Sampling\t\t" << result.pairSamplingResult << "\t\t" << result.pairSamplingTime << "\t\t" << (result.pairSamplingError * 100) << endl;
        file << "Node Sampling\t\t" << result.nodeSamplingResult << "\t\t" << result.nodeSamplingTime << "\t\t" << (result.nodeSamplingError * 100) << endl;
        file << endl;
        
        file << "Speedup Analysis:" << endl;
        file << "Degree method speedup: " << (result.exactTime / result.degreeMethodTime) << "x" << endl;
        file << "Pair sampling speedup: " << (result.exactTime / result.pairSamplingTime) << "x" << endl;
        file << "Node sampling speedup: " << (result.exactTime / result.nodeSamplingTime) << "x" << endl;
        file << "==========================================" << endl;
        file.close();
        cout << "\nResults written to method_comparison_results.txt" << endl;
    } else {
        cout << "\nWarning: Could not open file for writing results" << endl;
    }
    
    return result;
}

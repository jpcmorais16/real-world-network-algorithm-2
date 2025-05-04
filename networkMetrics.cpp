#include <vector>
#include <unordered_set>
#include <queue>
#include <numeric>
#include <thread>
#include <atomic>
#include <iostream>
#include <string>
#include <fstream>

#include "networkMetrics.h"

using namespace std;


double calculateCPUShortestPath(const vector<vector<int>>& fastGraph) {
    const size_t n = fastGraph.size();
    double totalSum = 0.0;

    const size_t maxThreads = std::thread::hardware_concurrency();
    const size_t numThreads = std::min(maxThreads ? maxThreads : 2, n);
    
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

double calculateAverageShortestPathLength(std::vector<std::unordered_set<int>>& graph) {
    const size_t n = graph.size();
    
    std::vector<std::vector<int>> fastGraph(n);
    for (size_t i = 0; i < n; ++i) {
        const size_t degree = graph[i].size();
        fastGraph[i].reserve(degree);
        fastGraph[i].assign(graph[i].begin(), graph[i].end());
    }
    
    return calculateCPUShortestPath(fastGraph);
}

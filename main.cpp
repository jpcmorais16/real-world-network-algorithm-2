#include <iostream>
#include <vector>
#include <unordered_set>
#include <random>
#include <numeric>
#include <chrono>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include "networkMetrics.h"

using namespace std;

std::random_device rd;
std::mt19937 gen(rd());

bool eventOccurs(double p) {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen) < p;
}

int get_random_integer(int min, int max) {
    std::uniform_int_distribution<> distr(min, max);
    return distr(gen);
}

std::unordered_map<double, std::vector<double>> p_cache;

double P(double p, int x) {
    if (p == 0.0) {
        return 0.1;
    }
    double denominator = 1.0 - pow(1.0 - p, 10);
    return (p * pow(1.0 - p, x - 1)) / denominator;
}

int get_number_of_steps(double p) {
    if (p_cache.find(p) == p_cache.end()) {
        std::vector<double> weights;
        for (int x = 1; x <= 10; ++x) {
            weights.push_back(P(p, x));
        }
        p_cache[p] = weights;
    }

    std::discrete_distribution<> dist(p_cache[p].begin(), p_cache[p].end());
    return dist(gen) + 1;
}

int random_walk(const std::vector<std::unordered_set<int>>& graph, int start, int number_of_steps) {
    int current_index = start;

    for (int i = 0; i < number_of_steps; i++) {
        const auto& current = graph[current_index];
        if (current.empty()) {
            break;
        }
        
        std::vector<int> neighbors(current.begin(), current.end());
        int set_index = get_random_integer(0, neighbors.size() - 1);
        current_index = neighbors[set_index];
    }

    return current_index;
}

std::vector<std::unordered_set<int>> generate_network(int N, int m, double p, double fp) {
    std::vector<std::unordered_set<int>> graph(10);

    for (int i = 0; i < 10; i++) {
        graph[i].insert((i + 1) % 10);
        if (i > 0) {
            graph[i].insert((i - 1));
        }
    }
    graph[0].insert(9);

    for (int i = 0; i < N; i++) {
        int current = get_random_integer(0, i + 9);
        std::unordered_set<int> marked;

        for (int j = 0; j < m; j++) {
            int number_of_steps = get_number_of_steps(p);
            current = random_walk(graph, current, number_of_steps);
            marked.insert(current);
        }

        graph.push_back(std::unordered_set<int>());
        
        for (int val : marked) {
            graph[i+10].insert(val);
            graph[val].insert(i+10);
        }
        
        for (int val1 : marked) {
            for (int val2 : marked) {
                if (val1 != val2 && eventOccurs(fp)) {
                    graph[val1].insert(val2);
                    graph[val2].insert(val1);
                }
            }
        }
    }

    std::cout << "Network generated successfully\n";

    return graph;
}

void writeNetworkToFile(const std::vector<std::unordered_set<int>>& graph, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    for (size_t i = 0; i < graph.size(); i++) {
        for (const int neighbor : graph[i]) {
            if (i < static_cast<size_t>(neighbor)) {
                file << i << " " << neighbor << "\n";
            }
        }
    }
    
    file.close();
    std::cout << "Network written to file: " << filename << std::endl;
}

int main() {
    int n;
    std::cin >> n;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    auto graph = generate_network(n, 2, 0.2, 0.5);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Network generation time: " << duration.count() << " milliseconds" << std::endl;
    
    writeNetworkToFile(graph, "network.txt");
    
    start_time = std::chrono::high_resolution_clock::now();
    double avg_path_length = calculateAverageShortestPathLength(graph);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "Average shortest path length: " << avg_path_length << std::endl;
    std::cout << "Metrics calculation time: " << duration.count() << " milliseconds" << std::endl;
}

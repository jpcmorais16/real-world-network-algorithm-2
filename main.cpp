#include <iostream>
#include <vector>
#include <unordered_set>
#include <random>
#include <numeric>
#include <chrono>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <iomanip>
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

    int n_nodes = 10;

    for (int i = 0; i < N; i++) {
        int start = get_random_integer(0, n_nodes - 1);
        int current = start;

        std::vector<int> marked;
        marked.push_back(start);

        for (int j = 0; j < m-1; j++) {
            int number_of_steps = get_number_of_steps(p);
            current = random_walk(graph, current, number_of_steps);
            marked.push_back(current);
        }

        graph.push_back(std::unordered_set<int>());
        
        for (int val : marked) {
            graph[n_nodes].insert(val);
            graph[val].insert(n_nodes);
        }
        
        for (size_t i = 0; i < marked.size(); i++) {
            for (size_t j = i + 1; j < marked.size(); j++) {
                if (marked[i] != marked[j] && eventOccurs(fp)) {
                    graph[marked[i]].insert(marked[j]);
                    graph[marked[j]].insert(marked[i]);
                }
            }
        }

        n_nodes++;
    }

    std::cout << "Network generated successfully\n";

    return graph;
}

struct SimulationResult {
    double avgClust;
    double transitivity;
    double avgPathLength;
    double pearsonR;
    double powerLawCoefficient;
    double maxDegree;
    double constructionTime;
    double pathCalcTime;
    double totalTime;
    size_t numThreads;
};

SimulationResult run_simulation(int N, double p, int m, double friendship_probability) {
    SimulationResult result;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    auto graph = generate_network(N, m, p, friendship_probability);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.constructionTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    
    auto metrics = calculateAllMetrics(graph);
    
    result.avgClust = metrics.clusteringCoefficient;
    result.transitivity = metrics.clusteringCoefficient;
    result.avgPathLength = metrics.avgPathLength;
    result.pearsonR = metrics.pearsonR;
    result.powerLawCoefficient = metrics.powerLawCoefficient;
    result.maxDegree = metrics.maxDegree;
    result.pathCalcTime = metrics.avgPathLengthDuration.count() / 1000.0;
    result.totalTime = result.constructionTime + result.pathCalcTime;
    result.numThreads = metrics.numThreads;
    
    std::cout << "Simulation completed - Clustering: " << result.avgClust 
              << ", Path length: " << result.avgPathLength 
              << ", Construction time: " << result.constructionTime << "s"
              << ", Path calc time: " << result.pathCalcTime << "s"
              << ", Total time: " << result.totalTime << "s"
              << ", Threads: " << metrics.numThreads << "\n";
              
    return result;
}

double calculate_std_dev(const std::vector<double>& values, double mean) {
    double variance = 0.0;
    for (double value : values) {
        double diff = value - mean;
        variance += diff * diff;
    }
    return std::sqrt(variance / values.size());
}

void create_point_sequential(const std::string& filepath, int N, double p, int m, double friendship_probability, int n_simulations) {
    std::cout << "Running " << n_simulations << " sequential simulations...\n";
    
    std::vector<SimulationResult> results;
    
    std::ofstream partial_file("partial.txt", std::ios_base::app);
    if (!partial_file.is_open()) {
        std::cerr << "Warning: Could not open partial.txt for writing individual results" << std::endl;
    } else {
        partial_file.seekp(0, std::ios::end);
        if (partial_file.tellp() == 0) {
            partial_file << "Simulation\tN\tp\tm\tfriendship_probability\tavg_clustering\ttransitivity\tavg_path\tr\tslope\tmax_degree\tconstruction_time\tpath_calc_time\ttotal_time\tthreads\n";
        }
        partial_file << std::fixed << std::setprecision(4);
    }
    
    for (int i = 0; i < n_simulations; i++) {
        std::cout << "Simulation " << (i+1) << "/" << n_simulations << "...\n";
        SimulationResult result = run_simulation(N, p, m, friendship_probability);
        results.push_back(result);
        
        if (partial_file.is_open()) {
            partial_file << (i+1) << "\t" 
                         << N << "\t" 
                         << p << "\t" 
                         << m << "\t" 
                         << friendship_probability << "\t"
                         << result.avgClust << "\t"
                         << result.transitivity << "\t"
                         << result.avgPathLength << "\t"
                         << result.pearsonR << "\t"
                         << result.powerLawCoefficient << "\t"
                         << result.maxDegree << "\t"
                         << result.constructionTime << "\t"
                         << result.pathCalcTime << "\t"
                         << result.totalTime << "\t"
                         << result.numThreads << "\n";
        }
    }
    
    if (partial_file.is_open()) {
        partial_file.close();
        std::cout << "Individual simulation results written to partial.txt" << std::endl;
    }
    
    // Calculate averages
    double avg_clust = 0.0, avg_trans = 0.0, avg_path = 0.0, avg_r = 0.0;
    double avg_slope = 0.0, avg_max_degree = 0.0, avg_construction_time = 0.0;
    double avg_shortest_path_calculation_time = 0.0;
    double avg_total_time = 0.0;
    
    std::vector<double> clust_values, trans_values, path_values, r_values;
    std::vector<double> slope_values, max_degree_values, construction_time_values, path_calc_time_values;
    std::vector<double> total_time_values;
    
    for (const auto& r : results) {
        avg_clust += r.avgClust;
        avg_trans += r.transitivity;
        avg_path += r.avgPathLength;
        avg_r += r.pearsonR;
        avg_slope += r.powerLawCoefficient;
        avg_max_degree += r.maxDegree;
        avg_construction_time += r.constructionTime;
        avg_shortest_path_calculation_time += r.pathCalcTime;
        avg_total_time += r.totalTime;
        
        clust_values.push_back(r.avgClust);
        trans_values.push_back(r.transitivity);
        path_values.push_back(r.avgPathLength);
        r_values.push_back(r.pearsonR);
        slope_values.push_back(r.powerLawCoefficient);
        max_degree_values.push_back(r.maxDegree);
        construction_time_values.push_back(r.constructionTime);
        path_calc_time_values.push_back(r.pathCalcTime);
        total_time_values.push_back(r.totalTime);
    }
    
    avg_clust /= n_simulations;
    avg_trans /= n_simulations;
    avg_path /= n_simulations;
    avg_r /= n_simulations;
    avg_slope /= n_simulations;
    avg_max_degree /= n_simulations;
    avg_construction_time /= n_simulations;
    avg_shortest_path_calculation_time /= n_simulations;
    avg_total_time /= n_simulations;
    
    double std_clust = calculate_std_dev(clust_values, avg_clust);
    double std_trans = calculate_std_dev(trans_values, avg_trans);
    double std_path = calculate_std_dev(path_values, avg_path);
    double std_r = calculate_std_dev(r_values, avg_r);
    double std_slope = calculate_std_dev(slope_values, avg_slope);
    double std_max_degree = calculate_std_dev(max_degree_values, avg_max_degree);
    double std_construction_time = calculate_std_dev(construction_time_values, avg_construction_time);
    double std_shortest_path_calculation_time = calculate_std_dev(path_calc_time_values, avg_shortest_path_calculation_time);
    double std_total_time = calculate_std_dev(total_time_values, avg_total_time);
    
    std::ofstream file(filepath, std::ios_base::app);
    if (file.is_open()) {
        file << std::fixed << std::setprecision(4);
        file << N << "\t" << p << "\t" << m << "\t" << friendship_probability << "\t"
             << avg_clust << " // " << std_clust << "\t"
             << avg_trans << " // " << std_trans << "\t"
             << avg_path << " // " << std_path << "\t"
             << avg_r << " // " << std_r << "\t"
             << avg_slope << " // " << std_slope << "\t"
             << avg_max_degree << " // " << std_max_degree << "\t"
             << avg_construction_time << " // " << std_construction_time << "\t"
             << avg_shortest_path_calculation_time << " // " << std_shortest_path_calculation_time << "\t"
             << avg_total_time << " // " << std_total_time << "\n";
        file.close();
        std::cout << "Results written to " << filepath << std::endl;
    } else {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
    }
}

struct PathComparisonResult {
    int N;
    double p;
    int m;
    double friendship_probability;
    double exactAveragePathLength;
    double approximateAveragePathLength;
    double error;
    size_t numEdges;
};

void create_comparison_sequential(const std::string& filepath, int N, double p, int m, double friendship_probability, int n_simulations) {
    std::cout << "Running " << n_simulations << " sequential simulations...\n";
    
    std::vector<PathComparisonResult> results;
    results.reserve(n_simulations);
    
    for (int i = 0; i < n_simulations; ++i) {
        std::cout << "Simulation " << (i + 1) << "/" << n_simulations << std::endl;
        
        // Create network
        auto network = generate_network(N, m, p, friendship_probability);
        
        // Convert to fast graph format
        auto fastGraph = convertToFastGraph(network);
        
        // Calculate both exact and approximate average path lengths
        auto pathLengths = compareApproximateAndExactAveragePathLength(fastGraph);
        double exactLength = pathLengths.first;
        double approximateLength = pathLengths.second;
        
        // Get number of edges
        size_t totalDegree = 0;
        for (const auto& neighbors : fastGraph) {
            totalDegree += neighbors.size();
        }
        size_t numEdges = totalDegree / 2;
        
        // Create result
        PathComparisonResult result;
        result.N = N;
        result.p = p;
        result.m = m;
        result.friendship_probability = friendship_probability;
        result.exactAveragePathLength = exactLength;
        result.approximateAveragePathLength = approximateLength;
        result.error = std::abs(exactLength - approximateLength) / exactLength; // Relative error
        result.numEdges = numEdges;
        
        results.push_back(result);
    }
    
    // Calculate statistics
    double avgExactLength = 0.0;
    double avgApproximateLength = 0.0;
    double avgError = 0.0;
    double maxError = 0.0;
    double minError = std::numeric_limits<double>::max();
    size_t avgNumEdges = 0;
    
    for (const auto& result : results) {
        avgExactLength += result.exactAveragePathLength;
        avgApproximateLength += result.approximateAveragePathLength;
        avgError += result.error;
        maxError = std::max(maxError, result.error);
        minError = std::min(minError, result.error);
        avgNumEdges += result.numEdges;
    }
    
    avgExactLength /= n_simulations;
    avgApproximateLength /= n_simulations;
    avgError /= n_simulations;
    avgNumEdges /= n_simulations;
    
    // Print results
    std::cout << "\nResults for N=" << N << ", p=" << p << ", m=" << m 
              << ", friendship_probability=" << friendship_probability << ":\n";
    std::cout << "Average Number of Edges: " << avgNumEdges << "\n";
    std::cout << "Average Exact Path Length: " << avgExactLength << "\n";
    std::cout << "Average Approximate Path Length: " << avgApproximateLength << "\n";
    std::cout << "Average Relative Error: " << (avgError * 100) << "%\n";
    std::cout << "Min Relative Error: " << (minError * 100) << "%\n";
    std::cout << "Max Relative Error: " << (maxError * 100) << "%\n\n";
    
    // Save results to file
    std::ofstream file(filepath, std::ios::app);
    if (file.is_open()) {
        file << std::fixed << std::setprecision(8);
        file << N << ", " << p << ", " << m << ", " << friendship_probability << "\t\t\t"
             << avgExactLength << ", " << avgApproximateLength << "\t\t\t"
             << avgError << ", " << minError << ", " << maxError << ", " << avgNumEdges << "\n\n";
        file.close();
    } else {
        std::cerr << "Error opening file: " << filepath << std::endl;
    }
}

int main() {
    // std::string filepath;
    // int N, m, n_simulations;
    // double p, friendship_probability;
    
    // std::cout << "Enter output filepath: ";
    // std::cin >> filepath;
    
    // std::cout << "Enter number of nodes (N): ";
    // std::cin >> N;
    
    // std::cout << "Enter step length probability (p): ";
    // std::cin >> p;
    
    // std::cout << "Enter number of marked nodes (m): ";
    // std::cin >> m;
    
    // std::cout << "Enter friendship probability: ";
    // std::cin >> friendship_probability;
    
    // std::cout << "Enter number of simulations: ";
    // std::cin >> n_simulations;

    // for (int N = 1000; N <= 100000; N=N+10000){
    //     for (int m = 2; m <= 10; m=m+8){
    //         for (double fp = 0.1; fp <= 0.9; fp=fp+0.8){
    //             create_comparison_sequential("approximation_results.txt", N, 0.5, m, fp, 5);
    //         }
    //     }
    // }   
    
    //create_comparison_sequential("approximation_results.txt", 200000, 0.5, 10, 0.9, 1);

    // for (int m = 2; m <= 10; m=m+1){
    //     for (double fp = 0.1; fp <= 0.9; fp=fp+0.1){
    //         for (double p = 0.1; p <= 1.0; p=p+0.9){
    //             create_comparison_sequential("edges_vs_error.txt", 10000, p, m, fp, 5);
    //         }
    //     }
    // }


    // for (double fp = 0.25; fp <= 0.75; fp=fp+0.25){
    //     for (int m = 2; m <= 10; m=m+1){
    //         create_comparison_sequential("edges_vs_error_feriado.txt", 50000, 0.5, m, fp, 5);
    //     }
    // }

    for(int i = 0; i < 2; i++){
        create_comparison_sequential("200000.txt", 200000, 0.5, 10, 0.9, 1);
    }

    // create_comparison_sequential("approximation_results.txt", 61000, 0.5, 10, 0.9, 5);
    // create_comparison_sequential("approximation_results.txt", 71000, 0.5, 2, 0.1, 5);
    
    // for (int N = 81000; N <= 101000; N=N+10000){
    //     create_comparison_sequential("approximation_results.txt", N, 0.5, 2, 0.1, 5);
    //     create_comparison_sequential("approximation_results.txt", N, 0.5, 10, 0.9, 5);
    // }

    // create_point_sequential(filepath, N, p, m, friendship_probability, n_simulations);
    
    return 0;
}

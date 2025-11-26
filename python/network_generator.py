"""
Core network generation functionality.
"""

import random
import math
import time
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
import powerlaw
import igraph

# Global random number generator
random.seed(time.time())

def event_occurs(p):
    """Returns True with probability p, False with probability 1-p"""
    return random.random() < p

def get_random_integer(min_val, max_val):
    """Returns a random integer between min_val and max_val (inclusive)"""
    return random.randint(min_val, max_val)

# Cache for P function calculations
p_cache = {}

def P(p, x):
    """Calculates the probability distribution for step lengths"""
    if p == 0.0:
        return 0.1
    denominator = 1.0 - pow(1.0 - p, 10)
    return (p * pow(1.0 - p, x - 1)) / denominator

def get_number_of_steps(p):
    """Returns a random number of steps based on probability distribution P"""
    if p not in p_cache:
        weights = []
        for x in range(1, 11):  # 1 to 10
            weights.append(P(p, x))
        p_cache[p] = weights
    
    # Simple weighted random selection without numpy
    total_weight = sum(p_cache[p])
    r = random.uniform(0, total_weight)
    cumulative_weight = 0
    for i, weight in enumerate(p_cache[p]):
        cumulative_weight += weight
        if r <= cumulative_weight:
            return i + 1
    return 10  # fallback

def random_walk(graph, start, number_of_steps):
    """Performs a random walk on the graph for the specified number of steps"""
    current_index = start
    
    for _ in range(number_of_steps):
        current_neighbors = graph[current_index]
        if not current_neighbors:
            break
        
        # Convert set to list for indexing
        neighbors_list = list(current_neighbors)
        set_index = get_random_integer(0, len(neighbors_list) - 1)
        current_index = neighbors_list[set_index]
    
    return current_index

def average_shortest_path_length(graph):
    n = len(graph)
    if n <= 1:
        return 0.0
    
    total_path_length = 0
    path_count = 0
    
    # For each node, find shortest paths to all other nodes
    for start in range(n):
        # BFS to find shortest distances from start to all other nodes
        distances = [-1] * n  # -1 means unreachable
        distances[start] = 0
        queue = [start]
        
        while queue:
            current = queue.pop(0)
            for neighbor in graph[current]:
                if distances[neighbor] == -1:  # Not visited yet
                    distances[neighbor] = distances[current] + 1
                    queue.append(neighbor)
        
        # Add all reachable distances to total
        for i in range(n):
            if i != start and distances[i] != -1:
                total_path_length += distances[i]
                path_count += 1
    
    # Return average path length, or 0 if no paths found (disconnected graph)
    return total_path_length / path_count if path_count > 0 else 0.0

def generate_network(N, m, p, fp=0.3) -> list[set]:
    # Input validation
    if N < 0:
        raise ValueError("N must be non-negative")
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must be between 0.0 and 1.0")
    if not 0.0 <= fp <= 1.0:
        raise ValueError("fp must be between 0.0 and 1.0")
    
    # Initialize graph with 10 nodes
    graph = [set() for _ in range(10)]
    
    # Create initial ring structure
    for i in range(10):
        graph[i].add((i + 1) % 10)
        if i > 0:
            graph[i].add(i - 1)
    graph[0].add(9)
    
    # Create random initial connections between the first 10 nodes
    for i in range(10):
        # Each node connects to 2-4 random other nodes
        num_connections = get_random_integer(2, 4)
        for _ in range(num_connections):
            target = get_random_integer(0, 9)
            if target != i:
                graph[i].add(target)
                graph[target].add(i)
    
    n_nodes = 10
    
    
    # Add N new nodes
    for _ in range(N):
        start = get_random_integer(0, n_nodes - 1)
        current = start
        
        marked = [start]
        
        # Perform random walks to find m-1 additional marked nodes
        for _ in range(int(m) - 1):
            number_of_steps = get_number_of_steps(p)
            current = random_walk(graph, current, number_of_steps)
            marked.append(current)

        if event_occurs(m - int(m)):
            number_of_steps = get_number_of_steps(p)
            current = random_walk(graph, current, number_of_steps)
            marked.append(current)
        
        # Add new node to graph
        graph.append(set())
        
        # Connect new node to all marked nodes
        for val in marked:
            graph[n_nodes].add(val)
            graph[val].add(n_nodes)
        
        # Create friendships between marked nodes with probability fp
        for i in range(len(marked)):
            for j in range(i + 1, len(marked)):
                if marked[i] != marked[j] and event_occurs(fp):
                    graph[marked[i]].add(marked[j])
                    graph[marked[j]].add(marked[i])
        
        n_nodes += 1
    
    return graph 

def calculate_alpha(graph):
    degree_distribution = [len(graph[i]) for i in range(len(graph))]
    count, bins = np.histogram(degree_distribution, bins=range(max(degree_distribution) + 1))
    pdf = count / sum(count)
    ccdf = 1 - np.cumsum(pdf)
    degrees = np.arange(len(pdf))
    valid = (degrees >= 1) & (ccdf > 0)
    slope, intercept = np.polyfit(np.log(degrees[valid]), np.log(ccdf[valid]), 1)
    alpha_estimate = 1.0 - slope
    return alpha_estimate

def graph_to_igraph(graph):
    edges = []
    for i, neighbors in enumerate(graph):
        for neighbor in neighbors:
            # Only add edge once (undirected graph)
            if i < neighbor:
                edges.append((i, neighbor))
    
    g = igraph.Graph(edges=edges, directed=False)
    return g

def get_last_completed_combination(output_file="output2.txt"):
    """Reads the output file and returns the last completed parameter combination.
    Returns (N, m, p, fp) tuple or None if file doesn't exist or is empty."""
    try:
        with open(output_file, "r") as f:
            lines = f.readlines()
            if not lines:
                return None
            
            # Get the last non-empty line
            last_line = None
            for line in reversed(lines):
                line = line.strip()
                if line:
                    last_line = line
                    break
            
            if not last_line:
                return None
            
            # Parse the line: "N, m, p, fp, avg_alpha, avg_path_length, avg_clustering"
            parts = [x.strip() for x in last_line.split(",")]
            if len(parts) < 4:
                return None
            
            N = int(parts[0])
            m = float(parts[1]) if '.' in parts[1] else int(parts[1])
            p = float(parts[2])
            fp = float(parts[3])
            
            return (N, m, p, fp)
    except FileNotFoundError:
        return None
    except (ValueError, IndexError):
        return None

def is_combination_before_or_equal(current, target):
    """Returns True if current combination is before or equal to target in iteration order."""
    if target is None:
        return False
    
    curr_N, curr_m, curr_p, curr_fp = current
    targ_N, targ_m, targ_p, targ_fp = target
    
    if curr_N < targ_N:
        return True
    if curr_N > targ_N:
        return False
    
    if curr_m < targ_m:
        return True
    if curr_m > targ_m:
        return False
    
    if curr_p < targ_p:
        return True
    if curr_p > targ_p:
        return False
    
    if curr_fp <= targ_fp:
        return True
    return False

last_completed = get_last_completed_combination("output2.txt")
if last_completed:
    print(f"Resuming from last completed combination: N={last_completed[0]}, m={last_completed[1]}, p={last_completed[2]}, fp={last_completed[3]}")
    print("Skipping already completed combinations...")
else:
    print("Starting from the beginning (no previous output found)")

for N in [500, 2500, 10000, 50000, 200000]:
    for m in [2, 4, 5, 6, 8, 10]:
        for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
            for fp in [0.1, 0.3, 0.5, 0.7, 0.9]:
                current_combination = (N, m, p, fp)
                if last_completed and is_combination_before_or_equal(current_combination, last_completed):
                    continue
                sum_alpha = 0
                sum_path_length = 0
                sum_clustering = 0

                n = 30
                
                for i in range(n):
                    brazil_time = datetime.now(ZoneInfo("America/Sao_Paulo"))
                    print(brazil_time.strftime("%Y-%m-%d %H:%M:%S %Z"), i, "/", n)

                    graph = generate_network(N, m, p, fp)

                    alpha_estimate = calculate_alpha(graph)
                    sum_alpha += alpha_estimate

                    g = graph_to_igraph(graph)
                    
                    if N < 50000 or i % 6 == 0:
                        path_length = g.average_path_length()
                        sum_path_length += path_length

                    clustering = g.transitivity_avglocal_undirected()
                    sum_clustering += clustering
                    
                avg_alpha = sum_alpha / n
                avg_path_length = sum_path_length / (n if N != 200000 else 5)
                avg_clustering = sum_clustering / n
                result_line = f"{N}, {m}, {p}, {fp}, {avg_alpha:.4f}, {avg_path_length:.4f}, {avg_clustering:.4f}\n"
                print(result_line.strip())
                
                with open("output2.txt", "a") as f:
                    f.write(result_line)


import networkx as nx
import igraph as ig
import time
import numpy as np

def calculate_avg_path_length_networkx(G):
    """Calculate average shortest path length using NetworkX."""
    # If the graph is not connected, calculate average path length for each component
    if not nx.is_connected(G):
        components = list(nx.connected_components(G))
        avg_paths = []
        nodes_count = []
        
        for component in components:
            if len(component) > 1:  # Only consider components with at least 2 nodes
                subgraph = G.subgraph(component)
                avg_paths.append(nx.average_shortest_path_length(subgraph))
                nodes_count.append(len(component))
                
        if not avg_paths:
            return 0
        
        # Weighted average based on number of nodes in each component
        return sum(p * n for p, n in zip(avg_paths, nodes_count)) / sum(nodes_count)
    else:
        return nx.average_shortest_path_length(G)

def calculate_avg_path_length_igraph(G):
    """Calculate average shortest path length using igraph."""
    # If the graph is not connected, calculate average path length for each component
    return G.average_path_length()

def main():
    # Read the network file
    print("Reading network file...")
    
    # NetworkX
    start_time = time.time()
    nx_graph = nx.read_edgelist('network.txt', nodetype=int)
    nx_load_time = time.time() - start_time
    print(f"NetworkX load time: {nx_load_time:.4f} seconds")
    print(f"NetworkX graph has {nx_graph.number_of_nodes()} nodes and {nx_graph.number_of_edges()} edges")
    
    # igraph
    start_time = time.time()
    ig_graph = ig.Graph.Read_Ncol('network.txt', names=True, directed=False)
    ig_load_time = time.time() - start_time
    print(f"igraph load time: {ig_load_time:.4f} seconds")
    print(f"igraph graph has {ig_graph.vcount()} nodes and {ig_graph.ecount()} edges")
    
    # Calculate average shortest path length with NetworkX
    # print("\nCalculating average shortest path length...")
    # start_time = time.time()
    # nx_avg_path = calculate_avg_path_length_networkx(nx_graph)
    # nx_calc_time = time.time() - start_time
    # print(f"NetworkX average shortest path length: {nx_avg_path:.6f}")
    # print(f"NetworkX calculation time: {nx_calc_time:.4f} seconds")
    
    # Calculate average shortest path length with igraph
    start_time = time.time()
    ig_avg_path = calculate_avg_path_length_igraph(ig_graph)
    ig_calc_time = time.time() - start_time
    print(f"igraph average shortest path length: {ig_avg_path:.6f}")
    print(f"igraph calculation time: {ig_calc_time:.4f} seconds")
    
    # Print summary
    print("\nSummary:")
    # print(f"NetworkX total time: {nx_load_time + nx_calc_time:.4f} seconds")
    print(f"igraph total time: {ig_load_time + ig_calc_time:.4f} seconds")
    # print(f"Speed ratio (NetworkX/igraph): {(nx_calc_time/ig_calc_time):.2f}x")

if __name__ == "__main__":
    main() 
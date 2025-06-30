import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

def parse_comparison_data(filename):
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    print(f"Read {len(lines)} lines from file")
    
    # Initialize list to store data
    data = []
    
    # Parse each line
    for line in lines:
        line = line.strip()
        if line:  # Skip empty lines
            # Parse the data line based on tab structure
            parts = line.split('\t\t\t')
            if len(parts) >= 3:  # We need three parts: params, metrics, and errors
                # First part contains parameters: N, p, m, fp
                params = parts[0].split(', ')
                # Second part contains exact and approx lengths
                metrics = parts[1].split(', ')
                # Third part contains errors: avg_error, min_error, max_error, num_edges
                errors = parts[2].split(', ')
                
                if len(params) >= 4 and len(errors) >= 4:
                    try:
                        N = int(params[0])
                        p = float(params[1])
                        m = int(params[2])
                        fp = float(params[3])
                        avg_error = float(errors[0])  # First value is average error
                        num_edges = int(errors[3])    # Fourth value is number of edges
                        
                        data.append({
                            'N': N,
                            'p': p,
                            'm': m,
                            'fp': fp,
                            'avg_error': avg_error,
                            'num_edges': num_edges
                        })
                        print(f"Successfully parsed data point: N={N}, m={m}, fp={fp}, edges={num_edges}, error={avg_error}")
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing line: {e}")
    
    print(f"Collected {len(data)} data points")
    df = pd.DataFrame(data)
    print(f"DataFrame columns: {df.columns}")
    print(f"DataFrame shape: {df.shape}")
    return df

def plot_m_vs_error(df):
    if df.empty:
        print("Error: DataFrame is empty!")
        return
        
    print("DataFrame contents:")
    print(df)
    
    # Create single plot
    plt.figure(figsize=(10, 6))
    
    # Filter data by friendship probability
    df_025 = df[df['fp'] == 0.25]
    df_075 = df[df['fp'] == 0.75]
    
    print(f"Data points with fp=0.25: {len(df_025)}")
    print(f"Data points with fp=0.75: {len(df_075)}")
    
    # Plot for fp=0.25
    if not df_025.empty:
        # Plot points
        plt.scatter(df_025['m'], df_025['avg_error'], c='blue', alpha=0.6, label='fp = 0.25')
        
        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(df_025['m'], df_025['avg_error'])
        line = slope * df_025['m'] + intercept
        #plt.plot(df_025['m'], line, 'b--', alpha=0.7, label=f'fp=0.25: y = {slope:.2e}x + {intercept:.2e} (R² = {r_value**2:.3f})')
    
    # Plot for fp=0.75
    if not df_075.empty:
        # Plot points
        plt.scatter(df_075['m'], df_075['avg_error'], c='red', alpha=0.6, label='fp = 0.75')
        
        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(df_075['m'], df_075['avg_error'])
        line = slope * df_075['m'] + intercept
        #plt.plot(df_075['m'], line, 'r--', alpha=0.7, label=f'fp=0.75: y = {slope:.2e}x + {intercept:.2e} (R² = {r_value**2:.3f})')
    
    # Customize the plot
    plt.xlabel('Number of Marked Nodes (m)')
    plt.ylabel('Average Relative Error')
    plt.title('Average Relative Error vs Number of Marked Nodes')
    plt.grid(True)
    plt.legend()
    
    # Save the plot
    plt.savefig('m_vs_error_feriado.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_log_edges_vs_error(df):
    if df.empty:
        print("Error: DataFrame is empty!")
        return
        
    print("DataFrame contents:")
    print(df)
    
    # Create single plot
    plt.figure(figsize=(10, 6))
    
    # Calculate log10 of edges for all data
    log_edges = np.log10(df['num_edges'])
    
    # Plot all points together
    plt.scatter(log_edges, df['avg_error'], c='blue', alpha=0.6, label='Data points')
    
    # Perform linear regression on all data
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_edges, df['avg_error'])
    line = slope * log_edges + intercept
    plt.plot(log_edges, line, 'r--', alpha=0.7, label=f'Linear regression: y = {slope:.2e}x + {intercept:.2e} (R² = {r_value**2:.3f})')
    
    # Customize the plot
    plt.xlabel('log₁₀(Number of Edges)')
    plt.ylabel('Average Relative Error')
    plt.title('Average Relative Error vs log₁₀(Number of Edges)')
    plt.grid(True)
    plt.legend()
    
    # Save the plot
    plt.savefig('log_edges_vs_error_feriado.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Parse the data
    df = parse_comparison_data('edges_vs_error_feriado.txt')
    
    # Create the plots
    plot_m_vs_error(df)
    plot_log_edges_vs_error(df)

if __name__ == "__main__":
    main() 
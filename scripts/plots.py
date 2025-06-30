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
                # Third part contains errors: avg_error, min_error, max_error
                errors = parts[2].split(', ')
                
                if len(params) >= 4 and len(errors) >= 3:
                    try:
                        N = int(params[0])
                        p = float(params[1])
                        m = int(params[2])
                        fp = float(params[3])
                        avg_error = float(errors[0])  # First value is average error
                        
                        data.append({
                            'N': N,
                            'p': p,
                            'm': m,
                            'fp': fp,
                            'avg_error': avg_error
                        })
                        print(f"Successfully parsed data point: N={N}, m={m}, fp={fp}, error={avg_error}")
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing line: {e}")
    
    print(f"Collected {len(data)} data points")
    df = pd.DataFrame(data)
    print(f"DataFrame columns: {df.columns}")
    print(f"DataFrame shape: {df.shape}")
    return df

def plot_n_vs_error(df):
    if df.empty:
        print("Error: DataFrame is empty!")
        return
        
    print("DataFrame contents:")
    print(df)
    
    # Filter data for m=2 with fp=0.1 and m=10 with fp=0.9, and N >= 11000
    df_m2 = df[(df['m'] == 2) & (df['fp'] == 0.1) & (df['N'] >= 11000)]
    df_m10 = df[(df['m'] == 10) & (df['fp'] == 0.9) & (df['N'] >= 11000)]
    
    print(f"Found {len(df_m2)} points for m=2 (fp=0.1) and {len(df_m10)} points for m=10 (fp=0.9)")
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Plot for m=2
    plt.plot(df_m2['N'], df_m2['avg_error'], 'bo-', label='m=2, fp=0.1')
    
    # Plot for m=10
    plt.plot(df_m10['N'], df_m10['avg_error'], 'ro-', label='m=10, fp=0.9')
    
    # Perform linear regression for m=2
    slope_m2, intercept_m2, r_value_m2, p_value_m2, std_err_m2 = stats.linregress(df_m2['N'], df_m2['avg_error'])
    line_m2 = slope_m2 * df_m2['N'] + intercept_m2
    plt.plot(df_m2['N'], line_m2, 'b--', alpha=0.5)
    equation_m2 = f'y = {slope_m2:.2e}x + {intercept_m2:.2e}\nR² = {r_value_m2**2:.3f}'
    
    # Perform linear regression for m=10
    slope_m10, intercept_m10, r_value_m10, p_value_m10, std_err_m10 = stats.linregress(df_m10['N'], df_m10['avg_error'])
    line_m10 = slope_m10 * df_m10['N'] + intercept_m10
    plt.plot(df_m10['N'], line_m10, 'r--', alpha=0.5)
    equation_m10 = f'y = {slope_m10:.2e}x + {intercept_m10:.2e}\nR² = {r_value_m10**2:.3f}'
    
    # Add equations to the plot
    plt.text(0.02, 0.45, f'm=2, fp=0.1:\n{equation_m2}', 
             transform=plt.gca().transAxes, verticalalignment='top', color='blue')
    plt.text(0.02, 0.65, f'm=10, fp=0.9:\n{equation_m10}', 
             transform=plt.gca().transAxes, verticalalignment='top', color='red')
    
    # Customize the plot
    plt.xlabel('Number of Nodes (N)')
    plt.ylabel('Average Relative Error')
    plt.title('Average Relative Error vs Number of Nodes\nfor m=2 (fp=0.1) and m=10 (fp=0.9)\nN ≥ 11000')
    plt.grid(True)
    plt.legend()
    
    # Save the plot
    plt.savefig('n_vs_error_feriado.png')
    plt.close()

def main():
    # Parse the data
    df = parse_comparison_data('edges_vs_error_feriado.txt')
    
    # Create the plot
    plot_n_vs_error(df)

if __name__ == "__main__":
    main()

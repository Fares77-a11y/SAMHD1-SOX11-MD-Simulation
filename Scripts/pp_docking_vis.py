import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.kaleido.scope.default_format = "svg"
pio.kaleido.scope.default_width = 1200
pio.kaleido.scope.default_height = 600
pio.kaleido.scope.default_scale = 2

pio.templates.default = 'ggplot2'

def cleaning(data, name):
    '''Cleaning MDS data.
    Inputs:
        data: a multi line string
        name: name of the file
    Outputs:
        a clean list of numbers only
    '''
    lines = data.strip().split('\n')
    final_data = []
    
    # Skip header lines that contain text
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # Skip lines with text headers
        if any(word in line.lower() for word in ['identifier', 'time', 'frame', '#', 'step']):
            continue
        
        # Handle space-separated data
        if ' ' in line:
            values = line.split()
            # Filter out non-numeric values
            numeric_values = []
            for val in values:
                try:
                    float(val)
                    numeric_values.append(val)
                except ValueError:
                    continue
            
            if not numeric_values:
                continue
                
            # Apply original logic for different file types
            if 'bonds' in name or 'RMSF' in name or 'comcom' in name or 'contacts' in name or 'interface' in name:
                # Take every second value starting from index 1
                for idx, val in enumerate(numeric_values):            
                    if idx % 2 == 1 and idx < len(numeric_values):
                        final_data.append(val)
            elif 'rog' in name.lower():
                # Skip first two columns if available
                if len(numeric_values) >= 3:
                    final_data.extend(numeric_values[2:])
                elif len(numeric_values) >= 2:
                    final_data.extend(numeric_values[1:])
                else:
                    final_data.extend(numeric_values)
            else:
                # For other files, take the last column
                final_data.append(numeric_values[-1])
        else:
            # Handle single-column data
            try:
                float(line)
                final_data.append(line)
            except ValueError:
                continue
    
    return final_data

def update_layout(fig, title):
    fig.update_layout(title={
                             'text': f'{title}',
                             'y':0.95,
                             'x':0.5,
                             'xanchor': 'center',
                             'yanchor': 'top'
                             }, 
                      legend=dict(
                                  y=-0.4,
                                  xanchor="center",
                                  x=0.5,
                                  orientation='h'
                                    ),
                      font_family="Times New Roman",
                      font_color="black",
                      font_size = 24,
                      title_font_family="Times New Roman",
                      title_font_color="black")     
    return fig

def plot(y, title, x=None):
    fig = go.Figure()
    
    # Updated color scheme - consistent across all figures
    chain_colors = {
        'A': '#1f77b4',      # Blue
        'B': '#ff7f0e',      # Orange  
        'C': '#2ca02c',      # Green
        'D': '#d62728',      # Red
        'E': '#9467bd',      # Purple
        'combined': '#8c564b', # Brown
        'complex': '#e377c2'   # Pink (changed from black for better visibility)
    }
    
    # Method-based line styles - consistent across all figures
    method_styles = {
        'Centroid': 'solid',
        'Random': 'dash'
    }
    
    opacity = 0.7 if 'Contact' in title or 'Hydrogen' in title else 0.8
    line_width = 3 if 'Contact' in title or 'Hydrogen' in title else 2
        
    if len(y) >= 1:
        if x is None:
            x = list(np.arange(len(y[0][1]))/10)  # Convert to ns
            
        for idx, n_data in enumerate(y):
            label, data = n_data
            words = label.split()
            
            # Initialize default values
            color = '#17becf'  # Default cyan color
            dash = 'solid'
            
            # Parse the label to determine color and line style
            if len(words) >= 3:
                # Format: "SAMHD1 A Centroid" or "SOX11 E Random" or "Complex complex Centroid"
                protein_or_complex = words[0]
                chain_or_identifier = words[1]
                method = words[2]
                
                if protein_or_complex in ['SAMHD1', 'SOX11']:
                    color = chain_colors.get(chain_or_identifier, color)
                elif protein_or_complex == 'Complex':
                    color = chain_colors.get('complex', color)
                
                dash = method_styles.get(method, 'solid')
                
            elif len(words) >= 2:
                # Format: "SAMHD1 A" or "Complex combined" or method names like "Centroid"
                first_word = words[0]
                second_word = words[1]
                
                if first_word in ['SAMHD1', 'SOX11']:
                    color = chain_colors.get(second_word, color)
                elif first_word == 'Complex':
                    color = chain_colors.get('complex', color)
                elif first_word in method_styles:
                    # Just method name
                    dash = method_styles.get(first_word, 'solid')
                    color = '#17becf'  # Default color for method-only labels
                elif second_word in method_styles:
                    # Something like "combined Centroid"
                    if first_word == 'combined':
                        color = chain_colors.get('combined', color)
                    dash = method_styles.get(second_word, 'solid')
                
            elif len(words) == 1:
                # Single word - could be method or identifier
                word = words[0]
                if word in method_styles:
                    dash = method_styles.get(word, 'solid')
                elif word in chain_colors:
                    color = chain_colors.get(word, color)
            
            fig.add_trace(go.Scatter(
                x=x,
                y=data,
                name=label,
                mode='lines',
                opacity=opacity,
                line=dict(
                    color=color,
                    width=line_width,
                    dash=dash
                )
            ))
            
    fig = update_layout(fig, title)
    return fig

def modify_limits(fig, name, data_series):
    if name == 'comcom':
        try:
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 1
            fig.update_yaxes(range=[0, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title="Å")
    
    elif name == 'hbonds':
        try:
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 2
            fig.update_yaxes(range=[0, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title="Number of H-bonds")
    
    elif name.startswith('contacts'):
        try:
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 5
            fig.update_yaxes(range=[0, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title="Number of Contacts")
    
    elif name == 'rog':
        try:
            lower_limit = min([np.array(s[1]).min() for s in data_series]) - 0.5
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 0.5
            fig.update_yaxes(range=[lower_limit, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title="Å")
    
    elif name == 'sasa' or name == 'interface_sasa':
        try:
            lower_limit = min([np.array(s[1]).min() for s in data_series]) - 1000
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 1000
            fig.update_yaxes(range=[lower_limit, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title="Å<sup>2</sup>")
    
    elif name.startswith('rmsd'):
        try:
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 1
            fig.update_yaxes(range=[0, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title="Å")
    
    elif name == 'rmsf':
        try:
            upper_limit = max([np.array(s[1]).max() for s in data_series]) + 1
            fig.update_yaxes(range=[0, upper_limit])
        except:
            pass
        fig.update_layout(xaxis_title="Amino Acid Number", yaxis_title="Å")
    
    return fig

def process_and_plot_docking(config):
    protein1 = config['protein1']
    protein2 = config['protein2']
    
    if not os.path.exists(config['output_dir']):
        os.makedirs(config['output_dir'])
    
    samhd1_chains = ['A', 'B', 'C', 'D']
    sox11_chains = ['E']
    
    # Interface metrics
    interface_metrics = {
        'comcom': f'comcom_{protein1}_{protein2}',
        'hbonds': f'hbonds_{protein1}_{protein2}',
        'interface_sasa': f'interface_SASA_{protein1}_{protein2}',
        'contacts_3.5A': f'contacts_{protein1}_{protein2}_3.5A',
        'contacts_5.0A': f'contacts_{protein1}_{protein2}_5.0A',
        'contacts_8.0A': f'contacts_{protein1}_{protein2}_8.0A'
    }
    
    # Plot interface metrics (if they exist)
    for metric, file_base in interface_metrics.items():
        data_series = []
        for data_dir, method in [(config['data_dir_centroid'], 'Centroid'), 
                               (config['data_dir_random'], 'Random')]:
            file_path = os.path.join(data_dir, f"{file_base}.dat")
            if os.path.exists(file_path):
                with open(file_path) as f:
                    data = f.read()
                cleaned_data = list(map(float, cleaning(data, file_path)))
                if cleaned_data:  # Only add if we have data
                    data_series.append((method, cleaned_data))
        
        if data_series:
            title = {
                'comcom': f"Distance between {protein1} and {protein2} Centers of Mass",
                'hbonds': f"Hydrogen Bonds between {protein1} and {protein2}",
                'interface_sasa': f"Interface SASA between {protein1} and {protein2}",
                'contacts_3.5A': f"Contacts between {protein1} and {protein2} (3.5Å)",
                'contacts_5.0A': f"Contacts between {protein1} and {protein2} (5.0Å)",
                'contacts_8.0A': f"Contacts between {protein1} and {protein2} (8.0Å)"
            }.get(metric, metric.replace('_', ' ').title())
            
            fig = plot(data_series, title)
            fig = modify_limits(fig, metric, data_series)
            output_path = os.path.join(config['output_dir'], f"{metric}")
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.svg", engine="kaleido")
    
    # Plot SAMHD1 chains - Updated to handle both 'rog' and 'RoG' naming
    for metric in ['rmsd', 'sasa', 'rog']:
        data_series = []
        for chain in samhd1_chains:
            for data_dir, method in [(config['data_dir_centroid'], 'Centroid'), 
                                   (config['data_dir_random'], 'Random')]:
                # Try both naming conventions
                file_bases = [f"{metric.upper()}_SAMHD1_{chain}"]
                if metric == 'rog':
                    file_bases.append(f"RoG_SAMHD1_{chain}")
                
                for file_base in file_bases:
                    file_path = os.path.join(data_dir, f"{file_base}.dat")
                    if os.path.exists(file_path):
                        with open(file_path) as f:
                            data = f.read()
                        cleaned_data = list(map(float, cleaning(data, file_path)))
                        if cleaned_data:  # Only add if we have data
                            data_series.append((f"SAMHD1 {chain} {method}", cleaned_data))
                        break  # Found the file, no need to try other naming
        
        if data_series:
            title = f"SAMHD1 Chains {metric.upper()}"
            fig = plot(data_series, title)
            fig = modify_limits(fig, metric, data_series)
            output_path = os.path.join(config['output_dir'], f"SAMHD1_chains_{metric}")
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.svg", engine="kaleido")
    
    # Plot SOX11 chain - Updated to handle both 'rog' and 'RoG' naming
    for metric in ['rmsd', 'sasa', 'rog']:
        data_series = []
        for chain in sox11_chains:
            for data_dir, method in [(config['data_dir_centroid'], 'Centroid'), 
                                   (config['data_dir_random'], 'Random')]:
                # Try both naming conventions
                file_bases = [f"{metric.upper()}_SOX11_{chain}"]
                if metric == 'rog':
                    file_bases.append(f"RoG_SOX11_{chain}")
                
                for file_base in file_bases:
                    file_path = os.path.join(data_dir, f"{file_base}.dat")
                    if os.path.exists(file_path):
                        with open(file_path) as f:
                            data = f.read()
                        cleaned_data = list(map(float, cleaning(data, file_path)))
                        if cleaned_data:  # Only add if we have data
                            data_series.append((f"SOX11 {chain} {method}", cleaned_data))
                        break  # Found the file, no need to try other naming
        
        if data_series:
            title = f"SOX11 {metric.upper()}"
            fig = plot(data_series, title)
            fig = modify_limits(fig, metric, data_series)
            output_path = os.path.join(config['output_dir'], f"SOX11_{metric}")
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.svg", engine="kaleido")
    
    # Plot combined SAMHD1 - Updated to handle both 'rog' and 'RoG' naming
    for metric in ['rmsd', 'sasa', 'rog']:
        data_series = []
        for data_dir, method in [(config['data_dir_centroid'], 'Centroid'), 
                               (config['data_dir_random'], 'Random')]:
            # Try both naming conventions
            file_bases = [f"{metric.upper()}_SAMHD1_combined"]
            if metric == 'rog':
                file_bases.append(f"RoG_SAMHD1_combined")
            
            for file_base in file_bases:
                file_path = os.path.join(data_dir, f"{file_base}.dat")
                if os.path.exists(file_path):
                    with open(file_path) as f:
                        data = f.read()
                    cleaned_data = list(map(float, cleaning(data, file_path)))
                    if cleaned_data:  # Only add if we have data
                        data_series.append((f"SAMHD1 combined {method}", cleaned_data))
                    break  # Found the file, no need to try other naming
        
        if data_series:
            title = f"SAMHD1 Combined {metric.upper()}"
            fig = plot(data_series, title)
            fig = modify_limits(fig, metric, data_series)
            output_path = os.path.join(config['output_dir'], f"SAMHD1_combined_{metric}")
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.svg", engine="kaleido")
    
    # Plot complex - Updated to handle both 'rog' and 'RoG' naming
    for metric in ['rmsd', 'sasa', 'rog']:
        data_series = []
        for data_dir, method in [(config['data_dir_centroid'], 'Centroid'), 
                               (config['data_dir_random'], 'Random')]:
            # Try both naming conventions
            file_bases = [f"{metric.upper()}_Complex"]
            if metric == 'rog':
                file_bases.append(f"RoG_Complex")
            
            for file_base in file_bases:
                file_path = os.path.join(data_dir, f"{file_base}.dat")
                if os.path.exists(file_path):
                    with open(file_path) as f:
                        data = f.read()
                    cleaned_data = list(map(float, cleaning(data, file_path)))
                    if cleaned_data:  # Only add if we have data
                        data_series.append((f"Complex complex {method}", cleaned_data))
                    break  # Found the file, no need to try other naming
        
        if data_series:
            title = f"Complex {metric.upper()}"
            fig = plot(data_series, title)
            fig = modify_limits(fig, metric, data_series)
            output_path = os.path.join(config['output_dir'], f"Complex_{metric}")
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.svg", engine="kaleido")
    
    # Plot RMSF for individual chains
    for chain in samhd1_chains + sox11_chains:
        data_series = []
        protein = 'SAMHD1' if chain in samhd1_chains else 'SOX11'
        for data_dir, method in [(config['data_dir_centroid'], 'Centroid'), 
                               (config['data_dir_random'], 'Random')]:
            file_base = f"RMSF_{protein}_{chain}"
            file_path = os.path.join(data_dir, f"{file_base}.dat")
            if os.path.exists(file_path):
                with open(file_path) as f:
                    data = f.read()
                cleaned_data = list(map(float, cleaning(data, file_path)))
                if cleaned_data:  # Only add if we have data
                    data_series.append((f"{protein} {chain} {method}", cleaned_data))
        
        if data_series:
            title = f"{protein} {chain} RMSF"
            x = list(map(str, config[f'aa_num_{protein}']))
            fig = plot(data_series, title, x=x)
            fig = modify_limits(fig, 'rmsf', data_series)
            output_path = os.path.join(config['output_dir'], f"{protein}_{chain}_rmsf")
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.svg", engine="kaleido")

if __name__ == "__main__":
    config = {
        'data_dir_centroid': '/Users/omara.soliman/Downloads/Master_work/SAMHD_SOX_Centroid',  # Update with your centroid data path
        'data_dir_random': '/Users/omara.soliman/Downloads/Master_work/SAMHD_SOX_Random',      # Update with your random data path
        'output_dir': '/Users/omara.soliman/Downloads/Master_work/output',          # Update with your output path
        'protein1': 'SAMHD1',
        'protein2': 'SOX11',
        'aa_num_SAMHD1': list(range(113, 600)),        # Adjust based on your SAMHD1 sequence
        'aa_num_SOX11': list(range(49, 196)),          # Adjust based on your SOX11 sequence
    }
    
    process_and_plot_docking(config)
#!/usr/bin/env python3
"""
Create PHR visualization filtered pafs. Methods section has full information on filtering. 
"""

import os
import sys
import glob
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np

def summarize_paf_file(paf_file):
    """
    Summarize alignments from PAF file
    Returns total bases and average identity
    """
    total_bases = 0
    identity_sum = 0
    alignment_count = 0
    
    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            # Get alignment length
            ref_start = int(fields[7])
            ref_end = int(fields[8])
            alignment_length = ref_end - ref_start
            
            # Get percent identity
            percent_identity = None
            for field in fields[12:]:
                if field.startswith('id:f:'):
                    percent_identity = float(field.split(':')[2]) * 100
                    break
                elif field.startswith('pi:f:'):
                    percent_identity = float(field.split(':')[2])
                    break
            
            if percent_identity is None:
                continue
            
            total_bases += alignment_length
            identity_sum += percent_identity * alignment_length
            alignment_count += 1
    
    avg_identity = identity_sum / total_bases if total_bases > 0 else 0
    
    return {
        'total_bases': total_bases,
        'avg_identity': avg_identity,
        'alignment_count': alignment_count
    }

def parse_chromosome_pair(filename):
    """Extract ref and query chromosome from filename"""
    basename = os.path.basename(filename).replace('.paf', '')
    parts = basename.split('_vs_')
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None

def create_network_plot(input_dir, output_file, min_bases=0):
    """
    Create network plot showing PHRs
    """
    # Get all PAF files and summarize
    paf_files = glob.glob(f"{input_dir}/*.paf")
    
    if len(paf_files) == 0:
        print(f"ERROR: No PAF files found in {input_dir}")
        sys.exit(1)
    
    print(f"Analyzing {len(paf_files)}")
    
    results = []
    for paf_file in paf_files:
        ref_chr, query_chr = parse_chromosome_pair(paf_file)
        if ref_chr is None or query_chr is None:
            continue
        
        summary = summarize_paf_file(paf_file)
        
        if summary['total_bases'] > min_bases:
            results.append({
                'ref_chr': ref_chr,
                'query_chr': query_chr,
                'total_bases': summary['total_bases'],
                'avg_identity': summary['avg_identity'],
                'alignment_count': summary['alignment_count']
            })
    
    if len(results) == 0:
        print("No PHRs found above minimum threshold.")
        sys.exit(1)
    
    df = pd.DataFrame(results)
    
    print(f"Found {len(df)} chromosome pairs with PHRs")
    print()
    
    # Create graph
    G = nx.Graph()
    
    # Add all chromosomes as nodes
    all_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    for chr_name in all_chromosomes:
        G.add_node(chr_name)
    
    # Add edges with attributes
    for _, row in df.iterrows():
        ref = row['ref_chr']
        query = row['query_chr']
        total_bases = row['total_bases']
        avg_identity = row['avg_identity']
        
        # Calculate edge properties
        line_width = (total_bases / 2_500_000) * 3 # Scale for visibility
        
        # Line color: normalize identity from 99-100% to 0-1
        if avg_identity >= 99:
            color_value = (avg_identity - 99.0) / 1.0
            color_value = max(0, min(1, color_value))
        else:
            color_value = 0
        
        G.add_edge(ref, query, 
                   weight=line_width,
                   color_value=color_value,
                   bases=total_bases,
                   identity=avg_identity)
    
    # Create circular layout
    def chr_sort_key(chr_name):
        """Sort chromosomes numerically"""
        if chr_name == 'chrX':
            return (23, 0)
        elif chr_name == 'chrY':
            return (24, 0)
        else:
            num = int(chr_name.replace('chr', ''))
            return (num, 0)
    
    sorted_nodes = sorted(G.nodes(), key=chr_sort_key)
    
    # Create circular layout
    pos = {}
    n = len(sorted_nodes)
    radius = 1.0
    for i, node in enumerate(sorted_nodes):
        angle = 2 * np.pi * i / n - np.pi / 2
        pos[node] = (radius * np.cos(angle), radius * np.sin(angle))
    
    # Prepare edge colors and widths
    edges = G.edges()
    edge_colors = [G[u][v]['color_value'] for u, v in edges] if len(edges) > 0 else []
    edge_widths = [G[u][v]['weight'] for u, v in edges] if len(edges) > 0 else []
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Draw edges
    if len(edges) > 0:
        cmap = cm.get_cmap('YlOrRd')
        nx.draw_networkx_edges(G, pos, 
                               edge_color=edge_colors, 
                               width=edge_widths,
                               edge_cmap=cmap,
                               edge_vmin=0,
                               edge_vmax=1,
                               alpha=0.6,
                               ax=ax)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, 
                          node_color='lightblue',
                          node_size=500,
                          ax=ax)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, 
                           font_size=9,
                           font_weight='bold',
                           ax=ax)
    
    # Add colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=99, vmax=100))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Average % Identity', rotation=270, labelpad=20)
    
    # Add title
    ax.set_title('PHRs\n'
                 f'Edge width = Shared sequence (scaled)\n'
                 f'Edge color = Avg % identity (99-100%)',
                 fontsize=12, pad=15)
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"PHR network plot saved to {output_file}")
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"  Nodes (chromosomes): {G.number_of_nodes()}")
    print(f"  Edges (chr pairs with PHRs): {G.number_of_edges()}")
    
    if G.number_of_edges() > 0:
        total_bases_all = sum([G[u][v]['bases'] for u, v in edges])
        avg_identity_all = np.mean([G[u][v]['identity'] for u, v in edges])
        print(f"  Total PHR sequence: {total_bases_all:,} bp ({total_bases_all/1e6:.1f} Mb)")
        print(f"  Average identity: {avg_identity_all:.2f}%")
        
        # Top 10 pairs
        print(f"\nTop 10 chromosome pairs by PHR content:")
        pairs = [(u, v, G[u][v]['bases'], G[u][v]['identity']) for u, v in edges]
        pairs.sort(key=lambda x: x[2], reverse=True)
        for i, (u, v, bases, identity) in enumerate(pairs[:10], 1):
            print(f"  {i}. {u} - {v}: {bases:,} bp ({bases/1e6:.2f} Mb), {identity:.2f}% identity")

def main():
    input_dir = "phrs"
    output_plot = "network.png"
    
    # Check if input directory exists
    if not os.path.exists(input_dir):
        print(f"Error: Directory {input_dir} not found!")
        sys.exit(1)
    
    # Create network plot
    create_network_plot(input_dir, output_plot, min_bases=10000)
    
    print("\nDone!")

if __name__ == "__main__":
    main()

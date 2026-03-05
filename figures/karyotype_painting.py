#!/usr/bin/env python3
"""
Generate karyotype painting of GCA_049354715.1 chromosomes colored by CHM13 chromosome alignments.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

def parse_paf(paf_file):
    """Parse PAF file and extract alignment information"""
    columns = ['query_name', 'query_len', 'query_start', 'query_end', 'strand',
               'target_name', 'target_len', 'target_start', 'target_end',
               'num_matches', 'aln_block_len', 'mapq']
    
    df = pd.read_csv(paf_file, sep='\t', usecols=range(12), names=columns, header=None)
    return df

def get_chromosome_order(chr_names):
    """Sort chromosomes"""
    def sort_key(chr_name):
        if chr_name.startswith('CM'):
            try:
                num = int(chr_name.replace('CM', '').split('.')[0])
                return num
            except:
                return float('inf')
        return float('inf')
    
    return sorted(chr_names, key=sort_key)

def get_hs1_color_scheme():
    """Define color scheme for human chromosomes"""
    colors = {
        'chr1':  '#332288',
        'chr2':  '#117733',
        'chr3':  '#44AA99',
        'chr4':  '#88CCEE',
        'chr5':  '#DDCC77',
        'chr6':  '#CC6677',
        'chr7':  '#AA4499',
        'chr8':  '#882255',
        'chr9':  '#E69F00',
        'chr10': '#56B4E9',
        'chr11': '#009E73',
        'chr12': '#F0E442',
        'chr13': '#0072B2',
        'chr14': '#D55E00',
        'chr15': '#CC79A7',
        'chr16': '#999933',
        'chr17': '#661100',
        'chr18': '#6699CC',
        'chr19': '#DDDDDD',
        'chr20': '#404040',
        'chr21': '#114477',
        'chr22': '#774411',
        'chrX':  '#EE3377',
        'chrY':  '#33BBEE',
        'chrM':  '#EE7733',
    }
    return colors

def create_karyotype_painting(paf_file, output_file, min_aln_len=10000, figsize=(20, 12), merge_gaps=False):
    """
    Create karyotype painting visualization
    """
    
    # Parse PAF file
    print("Parsing PAF file...")
    df = parse_paf(paf_file)
    
    # Filter by alignment length
    df = df[df['aln_block_len'] >= min_aln_len]
    print(f"Found {len(df)} alignments >= {min_aln_len} bp")
    
    # Get chromosome information
    query_chrs = get_chromosome_order(df['query_name'].unique())
    print(f"Found {len(query_chrs)} query chromosomes")
    
    # Get chromosome lengths
    chr_lengths = df.groupby('query_name')['query_len'].first().to_dict()
    
    # Get color scheme
    colors = get_hs1_color_scheme()
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot each chromosome
    chr_width = 0.8
    spacing = 0.2
    
    # Track maximum height for merged view
    max_height = 0
    
    for idx, chr_name in enumerate(query_chrs):
        chr_len = chr_lengths[chr_name]
        x_pos = idx
        
        # Get all alignments for this chromosome
        chr_alns = df[df['query_name'] == chr_name].sort_values('query_start')
        
        if merge_gaps:
            # Merge gaps: stack alignments vertically without whitespace
            current_y = 0
            
            for _, aln in chr_alns.iterrows():
                color = colors.get(aln['target_name'], '#cccccc')
                aln_len = aln['query_end'] - aln['query_start']
                
                # Draw alignment block at compressed position
                ax.add_patch(mpatches.Rectangle(
                    (x_pos - chr_width/2, current_y),
                    chr_width,
                    aln_len,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.8
                ))
                
                current_y += aln_len
            
            # Draw chromosome outline for merged view
            ax.add_patch(mpatches.Rectangle(
                (x_pos - chr_width/2, 0), chr_width, current_y,
                fill=False, edgecolor='black', linewidth=0.5
            ))
            
            max_height = max(max_height, current_y)
        else:
            # Original view: show actual genomic positions
            # Draw chromosome outline
            ax.add_patch(mpatches.Rectangle(
                (x_pos - chr_width/2, 0), chr_width, chr_len,
                fill=False, edgecolor='black', linewidth=0.5
            ))
            
            # Draw alignments
            for _, aln in chr_alns.iterrows():
                color = colors.get(aln['target_name'], '#cccccc')
                
                # Draw alignment block
                ax.add_patch(mpatches.Rectangle(
                    (x_pos - chr_width/2, aln['query_start']),
                    chr_width,
                    aln['query_end'] - aln['query_start'],
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.8
                ))
            
            max_height = max(max_height, chr_len)
    
    # Set axis properties
    ax.set_xlim(-0.5, len(query_chrs) - 0.5)
    ax.set_ylim(0, max_height * 1.02)
    
    # Set labels
    ax.set_xticks(range(len(query_chrs)))
    ax.set_xticklabels(query_chrs, rotation=90, ha='right', fontsize=8)
    ylabel = 'Cumulative aligned length (bp)' if merge_gaps else 'Position (bp)'
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel('GCA_049354715.1 Chromosomes', fontsize=12)
    title = 'Karyotype Painting: GCA_049354715.1 colored by hs1 alignments'
    if merge_gaps:
        title += ' (gaps removed)'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Format y-axis
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1e6)}Mb'))
    
    # Create legend
    target_chrs = sorted(df['target_name'].unique(), 
                        key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25').zfill(2)))
    legend_elements = [mpatches.Patch(facecolor=colors.get(chr, '#cccccc'), 
                                     edgecolor='black', 
                                     label=chr, alpha=0.8)
                      for chr in target_chrs if chr in colors]
    
    ax.legend(handles=legend_elements, 
             loc='center left', 
             bbox_to_anchor=(1.01, 0.5),
             ncol=1, 
             fontsize=9,
             title='hs1 chromosomes',
             title_fontsize=10,
             frameon=True)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    print(f"Saving figure to {output_file}...")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    print("Done!")
    
    # Print statistics
    print("\nAlignment statistics:")
    print(f"Total alignments: {len(df)}")
    print(f"Total aligned bases: {df['aln_block_len'].sum():,} bp")
    print(f"\nAlignments per hs1 chromosome:")
    for chr_name in target_chrs:
        count = len(df[df['target_name'] == chr_name])
        total_len = df[df['target_name'] == chr_name]['aln_block_len'].sum()
        print(f"  {chr_name}: {count} alignments ({total_len:,} bp)")

def main():
    parser = argparse.ArgumentParser(description='Generate karyotype painting from PAF alignments')
    parser.add_argument('paf_file', help='Input PAF file')
    parser.add_argument('-o', '--output', default='karyotype_painting.png',
                       help='Output figure file (default: karyotype_painting.png)')
    parser.add_argument('-m', '--min-length', type=int, default=10000,
                       help='Minimum alignment length to include (default: 10000)')
    parser.add_argument('-f', '--figsize', type=int, nargs=2, default=[20, 12],
                       help='Figure size in inches (width height) (default: 20 12)')
    parser.add_argument('--merge-gaps', action='store_true',
                       help='Remove whitespace by stacking alignments without gaps')
    
    args = parser.parse_args()
    
    create_karyotype_painting(
        args.paf_file,
        args.output,
        min_aln_len=args.min_length,
        figsize=tuple(args.figsize),
        merge_gaps=args.merge_gaps
    )

if __name__ == '__main__':
    main()

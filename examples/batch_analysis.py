"""
Batch Protein Structure Analysis Example

This script demonstrates how to analyze multiple protein structures
and compare their properties.
"""

from protein_structure_toolkit import (
    PDBParser,
    SecondaryStructureAnalyzer,
    StructureQualityChecker
)
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

print("""
╔═══════════════════════════════════════════════════════════════════╗
║     BATCH PROTEIN STRUCTURE ANALYSIS                              ║
╚═══════════════════════════════════════════════════════════════════╝
""")

# List of PDB files to analyze
pdb_files = [
    'data/example_protein.pdb',
    # Add more PDB files here
    # 'data/protein2.pdb',
    # 'data/protein3.pdb',
]

# Storage for results
results = []

# Analyze each structure
for pdb_file in pdb_files:
    print(f"\n{'='*70}")
    print(f"Analyzing: {pdb_file}")
    print('='*70)
    
    try:
        # Parse structure
        parser = PDBParser(pdb_file)
        summary = parser.get_structure_summary()
        
        # Secondary structure analysis
        ss_analyzer = SecondaryStructureAnalyzer(parser)
        ss_data = ss_analyzer.analyze_structure()
        ss_counts = ss_data['secondary_structure'].value_counts()
        
        # Quality metrics
        quality_checker = StructureQualityChecker(parser)
        b_metrics = quality_checker.analyze_temperature_factors()
        
        # Store results
        result = {
            'filename': pdb_file,
            'pdb_id': summary.get('pdb_id', 'N/A'),
            'total_atoms': summary['total_atoms'],
            'num_chains': summary['num_chains'],
            'num_ligands': summary['num_ligands'],
            'helix_count': ss_counts.get('H', 0),
            'sheet_count': ss_counts.get('E', 0),
            'coil_count': ss_counts.get('C', 0),
            'mean_b_factor': b_metrics['mean_b_factor'],
            'quality': 'High' if b_metrics['mean_b_factor'] < 30 else 
                      'Medium' if b_metrics['mean_b_factor'] < 50 else 'Low'
        }
        
        results.append(result)
        
        print(f"✅ Successfully analyzed {pdb_file}")
        
    except Exception as e:
        print(f"❌ Error analyzing {pdb_file}: {e}")
        continue

# Create results DataFrame
df_results = pd.DataFrame(results)

print("\n\n" + "="*70)
print("BATCH ANALYSIS SUMMARY")
print("="*70 + "\n")
print(df_results.to_string(index=False))

# Save results
df_results.to_csv('results/batch_analysis_summary.csv', index=False)
print("\n✅ Results saved to: results/batch_analysis_summary.csv")

# Create comparison plots if multiple structures
if len(results) > 1:
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Secondary structure comparison
    ax1 = axes[0, 0]
    ss_data = df_results[['helix_count', 'sheet_count', 'coil_count']]
    ss_data.plot(kind='bar', ax=ax1, color=['red', 'blue', 'gray'])
    ax1.set_xlabel('Structure', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Residue Count', fontsize=11, fontweight='bold')
    ax1.set_title('Secondary Structure Comparison', fontsize=12, fontweight='bold')
    ax1.legend(['α-Helix', 'β-Sheet', 'Coil/Loop'])
    ax1.set_xticklabels(df_results['pdb_id'], rotation=45)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: B-factor comparison
    ax2 = axes[0, 1]
    ax2.bar(range(len(df_results)), df_results['mean_b_factor'], 
           color='steelblue', edgecolor='black')
    ax2.axhline(y=30, color='green', linestyle='--', label='High Quality')
    ax2.axhline(y=50, color='orange', linestyle='--', label='Medium Quality')
    ax2.set_xlabel('Structure', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Mean B-factor (Å²)', fontsize=11, fontweight='bold')
    ax2.set_title('Structure Flexibility Comparison', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(df_results)))
    ax2.set_xticklabels(df_results['pdb_id'], rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Size comparison
    ax3 = axes[1, 0]
    ax3.scatter(df_results['total_atoms'], df_results['num_chains'],
               s=200, alpha=0.6, c=df_results['num_ligands'], 
               cmap='viridis', edgecolors='black', linewidth=2)
    for i, txt in enumerate(df_results['pdb_id']):
        ax3.annotate(txt, (df_results['total_atoms'].iloc[i], 
                          df_results['num_chains'].iloc[i]),
                    fontsize=10, fontweight='bold')
    ax3.set_xlabel('Total Atoms', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Number of Chains', fontsize=11, fontweight='bold')
    ax3.set_title('Structure Size Comparison', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Quality distribution
    ax4 = axes[1, 1]
    quality_counts = df_results['quality'].value_counts()
    colors_quality = {'High': 'green', 'Medium': 'orange', 'Low': 'red'}
    ax4.pie(quality_counts.values, labels=quality_counts.index,
           autopct='%1.1f%%', colors=[colors_quality[x] for x in quality_counts.index],
           startangle=90)
    ax4.set_title('Overall Quality Distribution', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('results/batch_comparison.png', dpi=300, bbox_inches='tight')
    print("✅ Comparison plots saved to: results/batch_comparison.png")
    plt.show()

print("\n" + "="*70)
print("✅ Batch Analysis Complete!")
print("="*70 + "\n")
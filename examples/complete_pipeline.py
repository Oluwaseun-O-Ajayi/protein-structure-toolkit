"""
Complete Protein Structure Analysis Pipeline

This script demonstrates the complete analysis pipeline that runs
all analyses at once.
"""

from protein_structure_toolkit import ProteinStructureAnalysisPipeline

print("""
╔═══════════════════════════════════════════════════════════════════╗
║     COMPLETE PROTEIN STRUCTURE ANALYSIS PIPELINE                  ║
╚═══════════════════════════════════════════════════════════════════╝

This will run all analyses:
1. Structure summary and parsing
2. Secondary structure analysis
3. Binding site identification
4. Structure quality assessment

All plots and data will be saved to the 'results/' directory.
""")

# Run complete analysis
pipeline = ProteinStructureAnalysisPipeline(
    pdb_file='data/example_protein.pdb',
    output_dir='results'
)

results = pipeline.run_complete_analysis()

# Access results
print("\n" + "="*70)
print("DETAILED RESULTS")
print("="*70)

if 'summary' in results:
    print("\n📊 Structure Summary:")
    for key, value in results['summary'].items():
        print(f"  {key}: {value}")

if 'secondary_structure' in results:
    print("\n🧬 Secondary Structure:")
    ss_counts = results['secondary_structure']['secondary_structure'].value_counts()
    total = ss_counts.sum()
    for ss_type, count in ss_counts.items():
        percentage = (count / total) * 100
        ss_name = {'H': 'α-Helix', 'E': 'β-Sheet', 'C': 'Coil/Loop'}[ss_type]
        print(f"  {ss_name}: {count} residues ({percentage:.1f}%)")

if 'binding_sites' in results:
    print("\n🎯 Binding Sites:")
    print(f"  Total residues: {len(results['binding_sites'])}")
    print(f"  Closest contact: {results['binding_sites']['distance'].min():.2f} Å")
    print(f"  Average distance: {results['binding_sites']['distance'].mean():.2f} Å")

if 'quality' in results:
    print("\n✅ Quality Metrics:")
    if 'b_factors' in results['quality']:
        b_metrics = results['quality']['b_factors']
        print(f"  Mean B-factor: {b_metrics['mean_b_factor']:.2f} Å²")
        if b_metrics['mean_b_factor'] < 30:
            print("  Assessment: Excellent quality (rigid structure)")
        elif b_metrics['mean_b_factor'] < 50:
            print("  Assessment: Good quality (moderate flexibility)")
        else:
            print("  Assessment: Lower quality (high flexibility)")

print("\n" + "="*70)
print("✅ PIPELINE COMPLETE!")
print("📁 All results saved to: results/")
print("="*70 + "\n")
"""
Basic Protein Structure Analysis Example

This script demonstrates basic usage of the Protein Structure Toolkit.
"""

from protein_structure_toolkit import (
    PDBParser,
    SecondaryStructureAnalyzer,
    BindingSiteAnalyzer,
    StructureQualityChecker
)

print("""
╔═══════════════════════════════════════════════════════════════════╗
║     BASIC PROTEIN STRUCTURE ANALYSIS                              ║
╚═══════════════════════════════════════════════════════════════════╝
""")

# Example 1: Parse PDB file and get summary
print("\n=== Example 1: PDB Parsing ===\n")

parser = PDBParser('data/example_protein.pdb')
summary = parser.get_structure_summary()

print("Structure Summary:")
for key, value in summary.items():
    print(f"  {key}: {value}")

# Calculate center of mass
com = parser.calculate_center_of_mass()
print(f"\nCenter of Mass: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})")


# Example 2: Secondary Structure Analysis
print("\n\n=== Example 2: Secondary Structure Analysis ===\n")

ss_analyzer = SecondaryStructureAnalyzer(parser)
ss_data = ss_analyzer.analyze_structure()

print(f"Total residues analyzed: {len(ss_data)}")
print("\nSecondary structure distribution:")
print(ss_data['secondary_structure'].value_counts())

# Generate plot
ss_analyzer.plot_secondary_structure(save_path='results/secondary_structure.png')


# Example 3: Binding Site Analysis
print("\n\n=== Example 3: Binding Site Analysis ===\n")

binding_analyzer = BindingSiteAnalyzer(parser)
binding_sites = binding_analyzer.identify_binding_sites(distance_threshold=5.0)

if not binding_sites.empty:
    print(f"Found {len(binding_sites)} residues in binding site")
    print("\nClosest contacts (< 4Å):")
    close_contacts = binding_sites[binding_sites['distance'] < 4.0]
    print(close_contacts[['residue_name', 'residue_number', 'distance']])
    
    # Generate plot
    binding_analyzer.plot_binding_site_analysis(save_path='results/binding_sites.png')
else:
    print("No ligands found in structure")


# Example 4: Structure Quality Check
print("\n\n=== Example 4: Structure Quality Assessment ===\n")

quality_checker = StructureQualityChecker(parser)
quality_checker.generate_quality_report()

# Generate plot
quality_checker.plot_quality_metrics(save_path='results/quality_metrics.png')


print("\n\n" + "="*70)
print("✅ Analysis Complete!")
print("📁 Results saved in 'results/' directory")
print("="*70 + "\n")
# Protein Structure Analysis Toolkit 🧬

A comprehensive Python toolkit for analyzing protein structures from PDB files. Designed for structural biology, drug discovery, and protein engineering research.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🎯 Overview

This toolkit provides robust analysis of protein structures with publication-quality visualizations, focusing on:

- **PDB File Parsing** - Extract structural information from Protein Data Bank files
- **Secondary Structure Analysis** - Identify α-helices, β-sheets, and loops
- **Binding Site Identification** - Locate and characterize ligand binding regions
- **Protein-Ligand Interactions** - Analyze contacts and distances
- **Structure Quality Assessment** - Evaluate B-factors, occupancy, and overall quality

## ✨ Features

- 📂 **Robust PDB parsing** - Extract atoms, residues, chains, and ligands
- 🔬 **Secondary structure prediction** - Classify residues into helix, sheet, and coil
- 🎯 **Binding site detection** - Identify residues near ligands
- 📊 **Comprehensive visualization** - Publication-quality plots
- ✅ **Quality metrics** - B-factor analysis, occupancy checks
- 🚀 **Batch processing** - Analyze multiple structures
- 📈 **Statistical analysis** - Distance distributions, residue composition

## 🚀 Quick Start

### Running Examples
```bash
# Navigate to project directory
cd protein-structure-toolkit

# Run basic analysis
python -m examples.basic_analysis

# Run complete pipeline
python -m examples.complete_pipeline

# Run batch analysis
python -m examples.batch_analysis
```

### Installation

```bash
# Clone the repository
git clone https://github.com/Oluwaseun-O-Ajayi/protein-structure-toolkit.git
cd protein-structure-toolkit

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
from protein_structure_toolkit import ProteinStructureAnalysisPipeline

# Run complete analysis on a PDB file
pipeline = ProteinStructureAnalysisPipeline('protein.pdb', output_dir='results')
results = pipeline.run_complete_analysis()
```

This single command will:
- Parse the PDB file
- Analyze secondary structure
- Identify binding sites
- Assess structure quality
- Generate all plots and reports

## 📋 Requirements

- Python 3.8+
- numpy >= 1.23.0
- pandas >= 1.5.0
- matplotlib >= 3.6.0
- seaborn >= 0.12.0

## 📊 Available Analyzers

### 1. PDBParser

Parse and extract information from PDB files.

```python
from protein_structure_toolkit import PDBParser

# Parse PDB file
parser = PDBParser('1abc.pdb')

# Get structure summary
summary = parser.get_structure_summary()
print(f"Total atoms: {summary['total_atoms']}")
print(f"Chains: {summary['chains']}")
print(f"Ligands: {summary['ligands']}")

# Get atoms as DataFrame for custom analysis
df = parser.get_atoms_dataframe()

# Calculate center of mass
com = parser.calculate_center_of_mass()
```

**Output:**
```
{
    'pdb_id': '1ABC',
    'total_atoms': 2456,
    'protein_atoms': 2312,
    'hetero_atoms': 144,
    'num_chains': 2,
    'chains': ['A', 'B'],
    'num_ligands': 1,
    'ligands': ['ATP:A:301']
}
```

### 2. SecondaryStructureAnalyzer

Analyze and visualize protein secondary structure.

```python
from protein_structure_toolkit import PDBParser, SecondaryStructureAnalyzer

parser = PDBParser('protein.pdb')
ss_analyzer = SecondaryStructureAnalyzer(parser)

# Analyze secondary structure
ss_data = ss_analyzer.analyze_structure()

# Visualize
ss_analyzer.plot_secondary_structure(save_path='secondary_structure.png')
```

**Features:**
- Classifies residues as α-helix (H), β-sheet (E), or coil (C)
- Shows distribution of secondary structure elements
- Plots secondary structure along the sequence
- Calculates percentages of each type

### 3. BindingSiteAnalyzer

Identify and characterize protein binding sites.

```python
from protein_structure_toolkit import PDBParser, BindingSiteAnalyzer

parser = PDBParser('protein_ligand.pdb')
binding_analyzer = BindingSiteAnalyzer(parser)

# Identify binding sites (residues within 5Å of ligand)
binding_sites = binding_analyzer.identify_binding_sites(distance_threshold=5.0)

# Visualize binding site composition
binding_analyzer.plot_binding_site_analysis(save_path='binding_sites.png')

# Export binding site residues
binding_sites.to_csv('binding_site_residues.csv', index=False)
```

**Output:**
- List of binding site residues
- Contact distances to ligand
- Residue composition in binding pocket
- Distance distribution histogram

### 4. StructureQualityChecker

Assess protein structure quality.

```python
from protein_structure_toolkit import PDBParser, StructureQualityChecker

parser = PDBParser('protein.pdb')
quality_checker = StructureQualityChecker(parser)

# Generate quality report
quality_checker.generate_quality_report()

# Analyze temperature factors
b_factor_metrics = quality_checker.analyze_temperature_factors()

# Check occupancy
occupancy_metrics = quality_checker.check_occupancy()

# Visualize quality metrics
quality_checker.plot_quality_metrics(save_path='quality_metrics.png')
```

**Metrics:**
- B-factor (temperature factor) analysis
- Occupancy values
- Structure flexibility assessment
- Quality classification (Excellent/Good/Poor)

## 📁 Project Structure

```
protein-structure-toolkit/
├── protein_structure_toolkit.py  # Main toolkit
├── README.md                      # This file
├── requirements.txt               # Dependencies
├── LICENSE                        # MIT License
├── .gitignore                    # Git ignore rules
├── examples/                      # Example scripts
│   ├── basic_analysis.py
│   ├── batch_analysis.py
│   └── custom_analysis.py
├── data/                          # Example PDB files
│   └── README.md
└── results/                       # Analysis outputs
    └── README.md
```

## 🔬 Example Workflows

### Workflow 1: Quick Structure Assessment

```python
from protein_structure_toolkit import ProteinStructureAnalysisPipeline

# Analyze downloaded PDB file
pipeline = ProteinStructureAnalysisPipeline('1abc.pdb')
results = pipeline.run_complete_analysis()
```

### Workflow 2: Binding Site Analysis for Drug Discovery

```python
from protein_structure_toolkit import PDBParser, BindingSiteAnalyzer

# Parse protein-ligand complex
parser = PDBParser('protein_drug_complex.pdb')

# Identify binding site
binding = BindingSiteAnalyzer(parser)
sites = binding.identify_binding_sites(distance_threshold=4.5)

# Analyze key residues
print("Key binding site residues:")
print(sites[sites['distance'] < 3.5])  # Close contacts

# Export for further analysis
sites.to_csv('drug_binding_site.csv')
```

### Workflow 3: Comparing Multiple Structures

```python
from protein_structure_toolkit import PDBParser, StructureQualityChecker
import pandas as pd

structures = ['wild_type.pdb', 'mutant1.pdb', 'mutant2.pdb']
quality_data = []

for pdb_file in structures:
    parser = PDBParser(pdb_file)
    checker = StructureQualityChecker(parser)
    metrics = checker.analyze_temperature_factors()
    
    quality_data.append({
        'structure': pdb_file,
        'mean_b_factor': metrics['mean_b_factor'],
        'flexibility': 'High' if metrics['mean_b_factor'] > 50 else 'Low'
    })

# Compare structures
comparison = pd.DataFrame(quality_data)
print(comparison)
```

### Workflow 4: Secondary Structure Comparison

```python
from protein_structure_toolkit import PDBParser, SecondaryStructureAnalyzer

proteins = {
    'Wild Type': 'wt.pdb',
    'Mutant A': 'mut_a.pdb',
    'Mutant B': 'mut_b.pdb'
}

for name, pdb_file in proteins.items():
    parser = PDBParser(pdb_file)
    ss = SecondaryStructureAnalyzer(parser)
    ss_data = ss.analyze_structure()
    
    print(f"\n{name}:")
    print(ss_data['secondary_structure'].value_counts())
```

## 🎯 Real-World Applications

### Drug Discovery
- **Target validation**: Analyze protein structure quality before screening
- **Binding site analysis**: Identify druggable pockets
- **Lead optimization**: Compare binding modes of different compounds

### Protein Engineering
- **Mutation effects**: Compare wild-type and mutant structures
- **Stability analysis**: Assess B-factors and flexibility
- **Structure validation**: Quality control for designed proteins

### Structural Biology Research
- **Publication figures**: Generate high-quality structure visualizations
- **Comparative analysis**: Batch process multiple structures
- **Data mining**: Extract structural features for machine learning

### Biotechnology
- **Enzyme engineering**: Analyze active site residues
- **Antibody design**: Study CDR regions and binding interfaces
- **Quality control**: Validate protein production batches

## 🛠️ Customization

### Adjust Distance Thresholds

```python
# Strict binding site definition
binding.identify_binding_sites(distance_threshold=3.5)

# Broad binding site definition
binding.identify_binding_sites(distance_threshold=6.0)
```

### Custom Plotting

```python
import matplotlib.pyplot as plt

# Customize figure size and resolution
ss_analyzer = SecondaryStructureAnalyzer(parser)
ss_analyzer.analyze_structure()

# High-resolution for publication
plt.rcParams['figure.dpi'] = 600
ss_analyzer.plot_secondary_structure(save_path='high_res_ss.png')
```

### Filter Specific Chains

```python
df = parser.get_atoms_dataframe()

# Analyze only chain A
chain_a = df[df['chain_id'] == 'A']
```

## 📖 PDB File Format

This toolkit expects standard PDB format files. You can download PDB files from:
- **RCSB PDB**: https://www.rcsb.org/
- **PDBe**: https://www.ebi.ac.uk/pdbe/
- **PDBj**: https://pdbj.org/

Example PDB file structure:
```
HEADER    HYDROLASE                               01-JAN-00   1ABC
ATOM      1  N   MET A   1      10.123  15.456  20.789  1.00 25.34      N
ATOM      2  CA  MET A   1      11.234  16.567  21.890  1.00 26.45      C
...
HETATM 2001  C1  ATP A 301      45.678  50.123  35.456  1.00 30.12      C
```

## ⚠️ Limitations & Notes

- **Simplified Secondary Structure**: Uses temperature factors for classification. For production use, consider BioPython's DSSP implementation for accurate phi-psi angle-based assignment.
- **No Explicit Hydrogen Bonds**: Binding site identification is distance-based. Advanced analysis might require hydrogen bond calculation.
- **Memory Usage**: Very large structures (>100,000 atoms) may require optimization.

## 🔄 Integration with Other Tools

This toolkit complements:
- **PyMOL**: For 3D visualization
- **BioPython**: For advanced structural analysis
- **MDAnalysis**: For MD trajectory analysis
- **VMD**: For molecular visualization
- **Chimera/ChimeraX**: For interactive structure exploration

## 🤝 Contributing

Contributions welcome! Areas for enhancement:
- DSSP-based secondary structure
- Hydrogen bond detection
- Salt bridge analysis
- Hydrophobic interaction finder
- Ramachandran plot generation
- Surface area calculations

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Developed for structural biology and drug discovery research
- Inspired by BioPython and MDAnalysis
- Thanks to the RCSB PDB for structural data

## 📧 Contact

**Oluwaseun O. Ajayi**  
PhD Researcher, Chemistry  
University of Georgia

- **GitHub**: [@Oluwaseun-O-Ajayi](https://github.com/Oluwaseun-O-Ajayi)
- **Academic Email**: oluwaseun.ajayi@uga.edu
- **Personal Email**: seunolanikeajayi@gmail.com

## 📖 Citation

If you use this toolkit in your research:

```bibtex
@software{protein_structure_toolkit,
  author = {Oluwaseun O. Ajayi},
  title = {Protein Structure Analysis Toolkit},
  year = {2024},
  url = {https://github.com/Oluwaseun-O-Ajayi/protein-structure-toolkit}
}
```

---

**Made with ❤️ for structural biology and computational chemistry research**
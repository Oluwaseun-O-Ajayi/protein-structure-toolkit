# Results Directory

Analysis outputs will be saved here.

## Generated Files

When you run the analysis scripts, the following files will be created:

### Plots
- `secondary_structure.png` - Secondary structure distribution and sequence plot
- `binding_sites.png` - Binding site composition and distance distribution
- `quality_metrics.png` - B-factor distribution and quality assessment
- `batch_comparison.png` - Comparison plots for multiple structures

### Data Files
- `batch_analysis_summary.csv` - Summary table for batch analysis

## Note

Generated output files (PNG, PDF, CSV) are not tracked by git to keep the repository size small. See `.gitignore`.

Run the example scripts to generate these files:
```bash
python examples/basic_analysis.py
python examples/batch_analysis.py
```
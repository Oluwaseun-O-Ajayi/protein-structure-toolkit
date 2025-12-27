# Data Directory

Place your PDB files here for analysis.

## Included Example

- `example_protein.pdb` - A simple example protein with ATP ligand for testing the toolkit

## Getting PDB Files

You can download PDB files from:

- **RCSB PDB**: https://www.rcsb.org/
- **PDB Europe**: https://www.ebi.ac.uk/pdbe/
- **PDB Japan**: https://pdbj.org/

## Example Downloads

```bash
# Download a specific PDB file (e.g., 1ABC)
wget https://files.rcsb.org/download/1ABC.pdb

# Or using curl
curl -O https://files.rcsb.org/download/1ABC.pdb
```

## File Format

PDB files should be in standard PDB format with ATOM and HETATM records.

Typical structure:
```
HEADER    ...
TITLE     ...
ATOM      1  N   MET A   1      ...
ATOM      2  CA  MET A   1      ...
...
HETATM  ...  (ligands)
END
```

## Note

Large PDB files (>10MB) are not tracked by git. See `.gitignore`.
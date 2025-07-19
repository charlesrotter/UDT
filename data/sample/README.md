# Sample SPARC Galaxy Data

This directory contains a sample galaxy rotation curve data file for testing and demonstration purposes.

## NGC3198_rotmod.dat

NGC 3198 is a well-studied spiral galaxy with high-quality rotation curve data, making it ideal for:
- Testing UDT analysis scripts
- Running the example Jupyter notebook
- Quick validation checks

### File Format
The file contains columns:
1. Radius (kpc)
2. Velocity (km/s)
3. Velocity error (km/s)
4. Additional model parameters

### Usage Example
```python
from udt.utils.data_loader import load_sparc_galaxy
data = load_sparc_galaxy('NGC3198', 'data/sample')
```

### Full Dataset
The complete SPARC database (~175 galaxies) is available at:
- Repository: `data/sparc_database/`
- Original source: http://astroweb.cwru.edu/SPARC/

This sample file is provided for quick testing without downloading the full dataset.
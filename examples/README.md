# UDT Examples

Interactive examples demonstrating UDT validation and analysis workflows.

## Notebooks

### `01_basic_udt_validation.ipynb`
Introduction to UDT validation framework:
- Loading real observational data
- Running basic galaxy and supernova analyses
- Understanding UDT predictions vs. ΛCDM

### `02_artifact_correction_demo.ipynb`
Demonstration of artifact correction methodology:
- ΛCDM contamination detection
- Bias-free correction application
- Validation of correction neutrality

### `03_data_integrity_workflow.ipynb`
Complete data integrity verification:
- SHA-256 manifest checking
- Synthetic data detection
- Quality assurance protocols

### `04_comprehensive_validation.ipynb`
Full UDT validation across all scales:
- Galaxy rotation curves (SPARC)
- Supernova cosmology (Pantheon+)
- CMB power spectrum analysis
- Multi-scale consistency checks

## Running the Examples

### Prerequisites
```bash
# Install UDT package
pip install -e .[dev]

# Start Jupyter
jupyter notebook examples/
```

### Data Requirements
Examples use sample data included in the repository. For full analyses:
```bash
# Download complete datasets (see docs/data.md)
udt-validate integrity --data-dir data/ --verify-manifest
```

## Example Outputs

All examples include:
- Clear explanations of UDT theory
- Real data visualization
- Statistical analysis results
- Comparison with ΛCDM predictions
- Discussion of limitations and uncertainties

## Scientific Standards

These examples demonstrate:
- ✅ **Real data only**: No synthetic data generation
- ✅ **Artifact correction**: Mandatory ΛCDM bias removal
- ✅ **Bias testing**: Validation of methodological neutrality
- ✅ **Uncertainty quantification**: Proper error propagation
- ✅ **Reproducibility**: Exact steps for replication

## Troubleshooting

### Common Issues
- **Import errors**: Check package installation
- **Data not found**: Verify data integrity manifest
- **Plotting issues**: Update matplotlib/jupyter

### Getting Help
1. Check notebook error messages
2. Verify data integrity: `udt-validate integrity --verify-manifest`
3. Test package imports: `python -c "from src.udt.models.core import UDTCosmology"`
4. Open GitHub issue with error details
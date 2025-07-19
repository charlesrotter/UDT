# Universal Distance Dilation Theory - Automated Results Generation
# Implements reviewer requirement: "make audit-all # reproduces all figures in manuscript"

.PHONY: all audit-all clean test figures tables data manuscript pdf
.DEFAULT_GOAL := audit-all

# Project configuration
PYTHON := python
RESULTS_DIR := results
FIG_DIR := fig
DOCS_DIR := docs
MANUSCRIPT := Universal_Distance_Dilation_Theory_Manuscript.md
OUTPUT_PDF := manuscript.pdf

# Data dependencies
SPARC_DATA := data/sparc_database/SPARC_Lelli2016c.mrt
LIGO_DATA := data/quantum_physics/GW150914_H1_32s.hdf5
MUON_DATA := data/quantum_physics/muon_g2_fermilab_data.json
PANTHEON_DATA := data/Pantheon_SH0ES.dat

# Analysis scripts
SPARC_SCRIPT := scripts/analyze_sparc_galaxies.py
LIGO_SCRIPT := quantum_validation/udt_ligo_final_analysis.py
MUON_SCRIPT := quantum_validation/pure_geometric_muon_g2_test.py
SUPERNOVA_SCRIPT := scripts/analyze_supernovae_raw.py
CMB_SCRIPT := scripts/analyze_cmb_planck.py
VALIDATION_SUITE := scripts/run_udt_validation_suite.py

# Output files
SPARC_RESULTS := $(RESULTS_DIR)/sparc_analysis/sparc_udt_results.csv
LIGO_RESULTS := $(RESULTS_DIR)/udt_ligo_final_assessment.json
MUON_RESULTS := $(RESULTS_DIR)/pure_geometric_muon_g2_results.json
SUPERNOVA_RESULTS := $(RESULTS_DIR)/supernova_raw_analysis/supernova_raw_results.csv

# Figures for manuscript
FIG1 := $(FIG_DIR)/fig1_sparc_rotation.png
FIG2 := $(FIG_DIR)/fig2_gw_residuals.png
FIG3 := $(FIG_DIR)/fig3_muon_g2_comparison.png
FIG4 := $(FIG_DIR)/fig4_multi_scale_validation.png

# Tables for manuscript
TABLES_DIR := tables
SPARC_TABLE := $(TABLES_DIR)/sparc_validation_statistics.csv
LIGO_TABLE := $(TABLES_DIR)/gravitational_wave_statistics.csv
MUON_TABLE := $(TABLES_DIR)/muon_g2_statistics.csv
VALIDATION_TABLE := $(TABLES_DIR)/multi_scale_validation_summary.csv

# Main target: reproduce all results
audit-all: figures tables data validation-report
	@echo "=== UDT Complete Audit Complete ==="
	@echo "All figures, tables, and data regenerated successfully"
	@echo "Results available in: $(RESULTS_DIR)/"
	@echo "Figures available in: $(FIG_DIR)/"
	@echo "Tables available in: $(TABLES_DIR)/"
	@echo ""
	@echo "To view manuscript: open $(MANUSCRIPT)"
	@echo "To generate PDF: make pdf"

# Generate all figures
figures: $(FIG1) $(FIG2) $(FIG3) $(FIG4)
	@echo "All manuscript figures generated"

# Generate all tables
tables: $(TABLES_DIR) $(SPARC_TABLE) $(LIGO_TABLE) $(MUON_TABLE) $(VALIDATION_TABLE)
	@echo "All manuscript tables generated"

# Generate all data
data: $(SPARC_RESULTS) $(LIGO_RESULTS) $(MUON_RESULTS) $(SUPERNOVA_RESULTS)
	@echo "All analysis data generated"

# Create directories
$(RESULTS_DIR):
	mkdir -p $(RESULTS_DIR)

$(FIG_DIR):
	mkdir -p $(FIG_DIR)

$(TABLES_DIR):
	mkdir -p $(TABLES_DIR)

# SPARC galaxy analysis
$(SPARC_RESULTS): $(SPARC_DATA) $(SPARC_SCRIPT) | $(RESULTS_DIR)
	$(PYTHON) $(SPARC_SCRIPT) --output-dir $(RESULTS_DIR)/sparc_analysis

$(FIG1): $(SPARC_RESULTS) | $(FIG_DIR)
	$(PYTHON) scripts/generate_sparc_figure.py --input $(SPARC_RESULTS) --output $(FIG1)

$(SPARC_TABLE): $(SPARC_RESULTS) | $(TABLES_DIR)
	$(PYTHON) scripts/extract_sparc_statistics.py --input $(SPARC_RESULTS) --output $(SPARC_TABLE)

# LIGO gravitational wave analysis
$(LIGO_RESULTS): $(LIGO_DATA) $(LIGO_SCRIPT) | $(RESULTS_DIR)
	$(PYTHON) $(LIGO_SCRIPT) --output-dir $(RESULTS_DIR)

$(FIG2): $(LIGO_RESULTS) | $(FIG_DIR)
	$(PYTHON) scripts/generate_ligo_figure.py --input $(LIGO_RESULTS) --output $(FIG2)

$(LIGO_TABLE): $(LIGO_RESULTS) | $(TABLES_DIR)
	$(PYTHON) scripts/extract_ligo_statistics.py --input $(LIGO_RESULTS) --output $(LIGO_TABLE)

# Muon g-2 analysis
$(MUON_RESULTS): $(MUON_DATA) $(MUON_SCRIPT) | $(RESULTS_DIR)
	$(PYTHON) $(MUON_SCRIPT) --output-dir $(RESULTS_DIR)

$(FIG3): $(MUON_RESULTS) | $(FIG_DIR)
	$(PYTHON) scripts/generate_muon_figure.py --input $(MUON_RESULTS) --output $(FIG3)

$(MUON_TABLE): $(MUON_RESULTS) | $(TABLES_DIR)
	$(PYTHON) scripts/extract_muon_statistics.py --input $(MUON_RESULTS) --output $(MUON_TABLE)

# Supernova analysis
$(SUPERNOVA_RESULTS): $(PANTHEON_DATA) $(SUPERNOVA_SCRIPT) | $(RESULTS_DIR)
	$(PYTHON) $(SUPERNOVA_SCRIPT) --output-dir $(RESULTS_DIR)/supernova_raw_analysis

# Multi-scale validation figure
$(FIG4): $(SPARC_RESULTS) $(LIGO_RESULTS) $(MUON_RESULTS) | $(FIG_DIR)
	$(PYTHON) scripts/generate_multiscale_figure.py --sparc $(SPARC_RESULTS) --ligo $(LIGO_RESULTS) --muon $(MUON_RESULTS) --output $(FIG4)

# Validation summary table
$(VALIDATION_TABLE): $(SPARC_RESULTS) $(LIGO_RESULTS) $(MUON_RESULTS) | $(TABLES_DIR)
	$(PYTHON) scripts/generate_validation_summary.py --sparc $(SPARC_RESULTS) --ligo $(LIGO_RESULTS) --muon $(MUON_RESULTS) --output $(VALIDATION_TABLE)

# Run comprehensive validation suite
validation-report: | $(RESULTS_DIR)
	$(PYTHON) $(VALIDATION_SUITE) --output-dir $(RESULTS_DIR)/validation_suite
	@echo "Comprehensive validation report generated"

# Generate PDF manuscript (requires pandoc)
pdf: $(MANUSCRIPT) figures tables
	@if command -v pandoc >/dev/null 2>&1; then \
		pandoc $(MANUSCRIPT) -o $(OUTPUT_PDF) --citeproc --number-sections --toc; \
		echo "PDF generated: $(OUTPUT_PDF)"; \
	else \
		echo "Warning: pandoc not found. Install pandoc to generate PDF"; \
		echo "Manuscript available as Markdown: $(MANUSCRIPT)"; \
	fi

# Quick validation (subset of data)
quick-audit: | $(RESULTS_DIR) $(FIG_DIR) $(TABLES_DIR)
	$(PYTHON) $(VALIDATION_SUITE) --max-galaxies 25 --output-dir $(RESULTS_DIR)/quick_validation
	$(PYTHON) scripts/generate_quick_figures.py --output-dir $(FIG_DIR)
	@echo "Quick audit complete (25 galaxies)"

# Test suite
test:
	$(PYTHON) -m pytest tests/ -v

# Data integrity check
check-data:
	@echo "Checking data integrity..."
	@if [ -f "$(SPARC_DATA)" ]; then echo "✓ SPARC data found"; else echo "✗ SPARC data missing"; fi
	@if [ -f "$(LIGO_DATA)" ]; then echo "✓ LIGO data found"; else echo "✗ LIGO data missing"; fi
	@if [ -f "$(MUON_DATA)" ]; then echo "✓ Muon g-2 data found"; else echo "✗ Muon g-2 data missing"; fi
	@if [ -f "$(PANTHEON_DATA)" ]; then echo "✓ Pantheon data found"; else echo "✗ Pantheon data missing"; fi

# Clean generated files
clean:
	rm -rf $(RESULTS_DIR)/sparc_analysis
	rm -rf $(RESULTS_DIR)/validation_suite
	rm -rf $(RESULTS_DIR)/quick_validation
	rm -f $(FIG_DIR)/*.png
	rm -f $(TABLES_DIR)/*.csv
	rm -f $(OUTPUT_PDF)
	@echo "Cleaned generated files"

# Development helpers
dev-setup:
	pip install -r requirements.txt
	@echo "Development environment setup complete"

# Repository integration check
repo-check:
	@echo "=== Repository Integration Check ==="
	@echo "Manuscript: $(shell [ -f "$(MANUSCRIPT)" ] && echo "✓ Found" || echo "✗ Missing")"
	@echo "Analysis scripts: $(shell [ -f "$(SPARC_SCRIPT)" ] && echo "✓ Found" || echo "✗ Missing")"
	@echo "Data files: $(shell [ -d "data" ] && echo "✓ Found" || echo "✗ Missing")"
	@echo "Results directory: $(shell [ -d "$(RESULTS_DIR)" ] && echo "✓ Found" || echo "✗ Missing")"
	@echo "Figure directory: $(shell [ -d "$(FIG_DIR)" ] && echo "✓ Found" || echo "✗ Missing")"

# Help target
help:
	@echo "UDT Automated Results Generation System"
	@echo ""
	@echo "Main targets:"
	@echo "  audit-all        - Reproduce all figures and tables in manuscript"
	@echo "  quick-audit      - Quick validation with subset of data"
	@echo "  figures          - Generate all manuscript figures"
	@echo "  tables           - Generate all manuscript tables"
	@echo "  data             - Generate all analysis data"
	@echo "  pdf              - Generate PDF manuscript (requires pandoc)"
	@echo "  test             - Run test suite"
	@echo ""
	@echo "Utility targets:"
	@echo "  clean            - Remove generated files"
	@echo "  check-data       - Verify data file availability"
	@echo "  repo-check       - Check repository structure"
	@echo "  dev-setup        - Setup development environment"
	@echo "  help             - Show this help"
	@echo ""
	@echo "Output locations:"
	@echo "  Results: $(RESULTS_DIR)/"
	@echo "  Figures: $(FIG_DIR)/"
	@echo "  Tables: $(TABLES_DIR)/"
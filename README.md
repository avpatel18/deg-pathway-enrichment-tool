
# DEG Pathway Enrichment Tool

A comprehensive Python package for differential gene expression (DEG) pathway enrichment analysis. This tool provides an easy-to-use interface for analyzing gene expression data and identifying enriched biological pathways.

## Features

- 🧬 **Comprehensive Pathway Analysis**: Supports multiple databases (KEGG, GO Biological Process, Reactome, MSigDB Hallmark)
- 📊 **Rich Visualizations**: Static and interactive plots with customizable DPI
- 🔍 **Gene Family Analysis**: Specialized analysis for keratin, claudin, and other gene families
- 📈 **Multiple Plot Types**: Bar plots, dot plots, volcano plots, and comprehensive summaries
- 📝 **Detailed Reports**: Markdown reports with analysis summaries
- 🎯 **Flexible Thresholds**: Customizable log fold change and p-value cutoffs
- 💻 **Command-line Interface**: Easy-to-use CLI for batch processing

## Installation

### From PyPI (Recommended)

```bash
pip install deg-pathway-enrichment-tool
```

### From Source

```bash
git clone https://github.com/yourusername/deg-pathway-enrichment-tool.git
cd deg-pathway-enrichment-tool
pip install -e .
```

## Quick Start

### Command Line Usage

```bash
# Basic analysis
deg-pathway-analysis your_deg_data.csv

# Custom output directory and thresholds
deg-pathway-analysis your_deg_data.csv -o results/ --logfc-threshold 2.0 --pval-threshold 0.001

# High-resolution figures
deg-pathway-analysis your_deg_data.csv --dpi 1200
```

### Python API Usage

```python
from deg_pathway_enrichment_tool import DEGPathwayAnalyzer

# Initialize analyzer
analyzer = DEGPathwayAnalyzer(
    deg_file="your_deg_data.csv",
    output_dir="./results",
    logfc_threshold=1.5,
    pval_threshold=0.01,
    dpi=600
)

# Run complete analysis
analyzer.run_complete_analysis()
```

## Input Data Format

Your CSV file must contain the following columns:

| Column Name | Description | Example |
|-------------|-------------|---------|
| `names` | Gene names/symbols | GAPDH, TP53, MYC |
| `logfoldchanges` | Log fold change values | 2.5, -1.8, 3.2 |
| `pvals_adj` | Adjusted p-values | 0.001, 0.05, 1e-10 |

### Example CSV format:
```csv
names,logfoldchanges,pvals_adj
KRT7,8.72,1e-300
CLDN10,8.52,3.35e-197
TP53,-2.1,0.001
GAPDH,1.2,0.05
```

## Output Files

The tool generates comprehensive results including:

### Visualizations
- `pathway_barplot.png` - Bar plot of top enriched pathways
- `pathway_dotplot.png` - Dot plot showing pathway significance vs effect size
- `interactive_pathway_barplot.html` - Interactive pathway visualization
- `keratin_expression.png` - Keratin gene family analysis
- `claudin_expression.png` - Claudin gene family analysis
- `comprehensive_summary.png` - Multi-panel summary figure

### Data Files
- `pathway_enrichment_results.csv` - Complete pathway enrichment results
- `keratin_genes.csv` - Keratin gene analysis results
- `claudin_genes.csv` - Claudin gene analysis results
- `analysis_report.md` - Comprehensive analysis report

## Command Line Options

```bash
deg-pathway-analysis --help
```

| Option | Description | Default |
|--------|-------------|---------|
| `input_file` | Path to CSV file containing DEG results | Required |
| `-o, --output-dir` | Output directory for results | `./deg_analysis_results` |
| `--logfc-threshold` | Log fold change threshold for significance | `1.5` |
| `--pval-threshold` | Adjusted p-value threshold for significance | `0.01` |
| `--dpi` | DPI for saved figures | `600` |
| `--databases` | Pathway databases to use | KEGG, GO-BP, Reactome, MSigDB |

## Supported Pathway Databases

- **KEGG_2021_Human**: KEGG pathway database
- **GO_Biological_Process_2021**: Gene Ontology Biological Process
- **Reactome_2022**: Reactome pathway database
- **MSigDB_Hallmark_2020**: MSigDB Hallmark gene sets

## Advanced Usage

### Custom Database Selection

```bash
deg-pathway-analysis input.csv --databases KEGG_2021_Human GO_Biological_Process_2021
```

### Python API - Step by Step

```python
from deg_pathway_enrichment_tool import DEGPathwayAnalyzer

# Initialize
analyzer = DEGPathwayAnalyzer("data.csv", output_dir="results")

# Run individual steps
pathway_results = analyzer.run_pathway_enrichment()
family_results = analyzer.analyze_gene_families()
analyzer.create_pathway_plots()
analyzer.create_comprehensive_summary()
analyzer.generate_report()
```

## Requirements

- Python ≥ 3.8
- pandas ≥ 2.0.0
- numpy ≥ 1.24.0
- matplotlib ≥ 3.7.0
- seaborn ≥ 0.12.0
- plotly ≥ 5.15.0
- gseapy ≥ 1.0.4
- scipy ≥ 1.10.0

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
DEG Pathway Enrichment Tool (2025). 
Available at: https://github.com/yourusername/deg-pathway-enrichment-tool
```

## Support

- 📖 Documentation: [GitHub Wiki](https://github.com/yourusername/deg-pathway-enrichment-tool/wiki)
- 🐛 Bug Reports: [GitHub Issues](https://github.com/yourusername/deg-pathway-enrichment-tool/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/yourusername/deg-pathway-enrichment-tool/discussions)

## Changelog

### v1.0.0 (2025-07-30)
- Initial release
- Comprehensive pathway enrichment analysis
- Multiple visualization options
- Gene family analysis
- Command-line interface
- Python API

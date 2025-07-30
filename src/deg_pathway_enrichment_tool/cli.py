
#!/usr/bin/env python3
"""
Command-line interface for DEG pathway enrichment analysis.
"""

import argparse
import sys
from pathlib import Path
from .core import DEGPathwayAnalyzer


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description="DEG Pathway Enrichment Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  deg-pathway-analysis input.csv

  # Custom output directory and thresholds
  deg-pathway-analysis input.csv -o results/ --logfc-threshold 2.0 --pval-threshold 0.001

  # High-resolution figures
  deg-pathway-analysis input.csv --dpi 1200

Required CSV columns:
  - names: Gene names/symbols
  - logfoldchanges: Log fold change values
  - pvals_adj: Adjusted p-values

For more information, visit: https://github.com/yourusername/deg-pathway-enrichment-tool
        """
    )
    
    parser.add_argument(
        "input_file",
        help="Path to CSV file containing DEG results"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        default="./deg_analysis_results",
        help="Output directory for results (default: ./deg_analysis_results)"
    )
    
    parser.add_argument(
        "--logfc-threshold",
        type=float,
        default=1.5,
        help="Log fold change threshold for significance (default: 1.5)"
    )
    
    parser.add_argument(
        "--pval-threshold",
        type=float,
        default=0.01,
        help="Adjusted p-value threshold for significance (default: 0.01)"
    )
    
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="DPI for saved figures (default: 600)"
    )
    
    parser.add_argument(
        "--databases",
        nargs="+",
        default=None,
        help="Pathway databases to use (default: KEGG, GO-BP, Reactome, MSigDB)"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="deg-pathway-enrichment-tool 1.0.0"
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    
    try:
        # Initialize analyzer
        analyzer = DEGPathwayAnalyzer(
            deg_file=args.input_file,
            output_dir=args.output_dir,
            logfc_threshold=args.logfc_threshold,
            pval_threshold=args.pval_threshold,
            dpi=args.dpi
        )
        
        # Run analysis
        analyzer.run_complete_analysis()
        
        print(f"\n✅ Analysis completed successfully!")
        print(f"📁 Results saved to: {Path(args.output_dir).absolute()}")
        
    except Exception as e:
        print(f"❌ Error during analysis: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()

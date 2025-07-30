
#!/usr/bin/env python3
"""
Core functionality for DEG pathway enrichment analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# import plotly.express as px
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
import gseapy as gp
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


class DEGPathwayAnalyzer:
    """Main class for DEG pathway enrichment analysis."""
    
    def __init__(self, deg_file, output_dir="./deg_analysis_results", 
                 logfc_threshold=1.5, pval_threshold=0.01, dpi=600):
        """
        Initialize the DEG pathway analyzer.
        
        Parameters:
        -----------
        deg_file : str
            Path to CSV file containing DEG results
        output_dir : str
            Directory to save output files
        logfc_threshold : float
            Log fold change threshold for significance
        pval_threshold : float
            Adjusted p-value threshold for significance
        dpi : int
            DPI for saved figures
        """
        self.deg_file = deg_file
        self.output_dir = Path(output_dir)
        self.logfc_threshold = logfc_threshold
        self.pval_threshold = pval_threshold
        self.dpi = dpi
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load and process data
        self.deg_df = None
        self.significant_genes = None
        self.pathway_results = {}
        
        self._load_data()
        
    def _load_data(self):
        """Load and validate DEG data."""
        try:
            self.deg_df = pd.read_csv(self.deg_file)
            print(f"Loaded {len(self.deg_df)} genes from {self.deg_file}")
            
            # Validate required columns
            required_cols = ['names', 'logfoldchanges', 'pvals_adj']
            missing_cols = [col for col in required_cols if col not in self.deg_df.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
            
            # Filter significant genes
            self.significant_genes = self.deg_df[
                (abs(self.deg_df['logfoldchanges']) > self.logfc_threshold) & 
                (self.deg_df['pvals_adj'] < self.pval_threshold)
            ].copy()
            
            print(f"Found {len(self.significant_genes)} significant genes")
            print(f"Upregulated: {len(self.significant_genes[self.significant_genes['logfoldchanges'] > 0])}")
            print(f"Downregulated: {len(self.significant_genes[self.significant_genes['logfoldchanges'] < 0])}")
            
        except Exception as e:
            raise ValueError(f"Error loading DEG data: {str(e)}")
    
    def run_pathway_enrichment(self, databases=None):
        """
        Run pathway enrichment analysis.
        
        Parameters:
        -----------
        databases : list
            List of databases to use for enrichment
        """
        if databases is None:
            databases = [
                'KEGG_2021_Human',
                'GO_Biological_Process_2021',
                'Reactome_2022',
                'MSigDB_Hallmark_2020'
            ]
        
        gene_list = self.significant_genes['names'].tolist()
        print(f"Running pathway enrichment with {len(gene_list)} genes...")
        
        all_results = []
        
        for db in databases:
            try:
                print(f"Analyzing {db}...")
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=db,
                    organism='Human',
                    cutoff=0.05,
                    no_plot=True
                )
                
                if not enr.results.empty:
                    results = enr.results.copy()
                    results['Database'] = db
                    all_results.append(results)
                    print(f"  Found {len(results)} significant pathways")
                else:
                    print(f"  No significant pathways found")
                    
            except Exception as e:
                print(f"  Error with {db}: {str(e)}")
                continue
        
        if all_results:
            combined_results = pd.concat(all_results, ignore_index=True)
            combined_results = combined_results.sort_values('Adjusted P-value')
            
            # Save results
            output_file = self.output_dir / 'pathway_enrichment_results.csv'
            combined_results.to_csv(output_file, index=False)
            print(f"Saved pathway results to {output_file}")
            
            self.pathway_results['all'] = combined_results
            return combined_results
        else:
            print("No pathway enrichment results found")
            return pd.DataFrame()
    
    def analyze_gene_families(self):
        """Analyze specific gene families (keratin, claudin, etc.)."""
        gene_families = {
            'Keratin': 'KRT',
            'Claudin': 'CLDN',
            'Tight Junction': 'TJP|OCLN|CLDN',
            'Cadherin': 'CDH'
        }
        
        family_results = {}
        
        for family_name, pattern in gene_families.items():
            family_genes = self.significant_genes[
                self.significant_genes['names'].str.contains(pattern, case=False, na=False)
            ].copy()
            
            if not family_genes.empty:
                family_genes = family_genes.sort_values('logfoldchanges', ascending=False)
                family_results[family_name] = family_genes
                
                # Save family results
                output_file = self.output_dir / f'{family_name.lower()}_genes.csv'
                family_genes.to_csv(output_file, index=False)
                print(f"Found {len(family_genes)} {family_name} genes")
        
        return family_results
    
    def create_pathway_plots(self):
        """Create pathway enrichment plots."""
        if 'all' not in self.pathway_results or self.pathway_results['all'].empty:
            print("No pathway results available for plotting")
            return
        
        results = self.pathway_results['all']
        
        # 1. Bar plot of top pathways
        self._create_pathway_barplot(results)
        
        # 2. Dot plot
        self._create_pathway_dotplot(results)
        
        # 3. Interactive plots (disabled for now)
        # self._create_interactive_plots(results)
    
    def _create_pathway_barplot(self, results, top_n=20):
        """Create bar plot of top pathways."""
        top_results = results.head(top_n)
        top_results['-log10_adjp'] = -np.log10(top_results['Adjusted P-value'])
        
        plt.figure(figsize=(12, 8))
        bars = plt.barh(range(len(top_results)), top_results['-log10_adjp'],
                       color=plt.cm.viridis(np.linspace(0, 1, len(top_results))))
        
        plt.yticks(range(len(top_results)), 
                  [f"{row['Database'].split('_')[0]}: {row['Term'][:50]}..." 
                   if len(row['Term']) > 50 else f"{row['Database'].split('_')[0]}: {row['Term']}"
                   for _, row in top_results.iterrows()])
        
        plt.xlabel('-log10(Adjusted P-value)')
        plt.title(f'Top {top_n} Enriched Pathways')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        
        output_file = self.output_dir / 'pathway_barplot.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        print(f"Saved pathway bar plot to {output_file}")
    
    def _create_pathway_dotplot(self, results, top_n=20):
        """Create dot plot of pathways."""
        top_results = results.head(top_n)
        top_results['Gene_Count'] = top_results['Genes'].str.split(';').str.len()
        top_results['-log10_adjp'] = -np.log10(top_results['Adjusted P-value'])
        
        plt.figure(figsize=(10, 8))
        scatter = plt.scatter(top_results['Combined Score'], 
                            range(len(top_results)),
                            s=top_results['Gene_Count'] * 10,
                            c=top_results['-log10_adjp'],
                            cmap='viridis', alpha=0.7)
        
        plt.yticks(range(len(top_results)),
                  [f"{row['Term'][:40]}..." if len(row['Term']) > 40 else row['Term']
                   for _, row in top_results.iterrows()])
        
        plt.xlabel('Combined Score')
        plt.title(f'Top {top_n} Pathways - Dot Plot')
        plt.colorbar(scatter, label='-log10(Adjusted P-value)')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        output_file = self.output_dir / 'pathway_dotplot.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        print(f"Saved pathway dot plot to {output_file}")
    
    def _create_interactive_plots(self, results):
        """Create interactive plotly visualizations."""
        print("Interactive plots disabled in this version")
    
    def create_gene_family_plots(self, family_results):
        """Create plots for gene families."""
        for family_name, genes in family_results.items():
            if genes.empty:
                continue
                
            plt.figure(figsize=(10, 6))
            genes_sorted = genes.sort_values('logfoldchanges', ascending=True)
            
            bars = plt.barh(genes_sorted['names'], genes_sorted['logfoldchanges'],
                           color='darkred' if family_name == 'Keratin' else 'darkblue',
                           alpha=0.8)
            
            # Add value labels
            for bar in bars:
                width = bar.get_width()
                plt.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                        f'{width:.1f}', ha='left', va='center')
            
            plt.xlabel('Log Fold Change')
            plt.title(f'{family_name} Gene Expression')
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            
            output_file = self.output_dir / f'{family_name.lower()}_expression.png'
            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            plt.close()
            print(f"Saved {family_name} plot to {output_file}")
    
    def create_comprehensive_summary(self):
        """Create comprehensive analysis summary."""
        # Create summary figure with multiple subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # 1. Top pathways
        if 'all' in self.pathway_results and not self.pathway_results['all'].empty:
            top_10 = self.pathway_results['all'].head(10)
            top_10['-log10_adjp'] = -np.log10(top_10['Adjusted P-value'])
            
            axes[0,0].barh(range(len(top_10)), top_10['-log10_adjp'],
                          color=plt.cm.viridis(np.linspace(0, 1, len(top_10))))
            axes[0,0].set_yticks(range(len(top_10)))
            axes[0,0].set_yticklabels([f"{row['Database'].split('_')[0]}: {row['Term'][:25]}..."
                                      for _, row in top_10.iterrows()], fontsize=8)
            axes[0,0].set_xlabel('-log10(Adjusted P-value)')
            axes[0,0].set_title('Top 10 Enriched Pathways')
            axes[0,0].grid(axis='x', alpha=0.3)
        
        # 2. Gene expression distribution
        axes[0,1].hist(self.significant_genes['logfoldchanges'], bins=30, 
                      alpha=0.7, color='green', edgecolor='black')
        axes[0,1].axvline(self.significant_genes['logfoldchanges'].mean(), 
                         color='red', linestyle='--',
                         label=f'Mean: {self.significant_genes["logfoldchanges"].mean():.2f}')
        axes[0,1].set_xlabel('Log Fold Change')
        axes[0,1].set_ylabel('Number of Genes')
        axes[0,1].set_title('Distribution of Log Fold Changes')
        axes[0,1].legend()
        axes[0,1].grid(alpha=0.3)
        
        # 3. Top genes volcano-like plot
        top_genes = self.significant_genes.head(15)
        top_genes['-log10_adjp'] = -np.log10(top_genes['pvals_adj'].replace(0, 1e-300))
        
        scatter = axes[0,2].scatter(top_genes['logfoldchanges'], top_genes['-log10_adjp'],
                                  s=60, alpha=0.7, c=range(len(top_genes)), cmap='viridis')
        axes[0,2].set_xlabel('Log Fold Change')
        axes[0,2].set_ylabel('-log10(Adjusted P-value)')
        axes[0,2].set_title('Top 15 Most Significant Genes')
        axes[0,2].grid(alpha=0.3)
        
        # 4. Keratin genes (if available)
        keratin_genes = self.significant_genes[
            self.significant_genes['names'].str.contains('KRT', case=False, na=False)
        ]
        if not keratin_genes.empty:
            keratin_sorted = keratin_genes.sort_values('logfoldchanges', ascending=True)
            axes[1,0].barh(keratin_sorted['names'], keratin_sorted['logfoldchanges'],
                          color='darkred', alpha=0.8)
            axes[1,0].set_xlabel('Log Fold Change')
            axes[1,0].set_title('Keratin Gene Expression')
            axes[1,0].grid(axis='x', alpha=0.3)
        
        # 5. Claudin genes (if available)
        claudin_genes = self.significant_genes[
            self.significant_genes['names'].str.contains('CLDN', case=False, na=False)
        ]
        if not claudin_genes.empty:
            claudin_sorted = claudin_genes.sort_values('logfoldchanges', ascending=True)
            axes[1,1].barh(claudin_sorted['names'], claudin_sorted['logfoldchanges'],
                          color='darkblue', alpha=0.8)
            axes[1,1].set_xlabel('Log Fold Change')
            axes[1,1].set_title('Claudin Gene Expression')
            axes[1,1].grid(axis='x', alpha=0.3)
        
        # 6. Summary statistics
        stats_text = f"""Analysis Summary:
        
Total genes: {len(self.deg_df)}
Significant genes: {len(self.significant_genes)}
Upregulated: {len(self.significant_genes[self.significant_genes['logfoldchanges'] > 0])}
Downregulated: {len(self.significant_genes[self.significant_genes['logfoldchanges'] < 0])}

Thresholds:
|LogFC| > {self.logfc_threshold}
Adj. p-value < {self.pval_threshold}

Top gene: {self.significant_genes.iloc[0]['names']}
LogFC: {self.significant_genes.iloc[0]['logfoldchanges']:.2f}
"""
        axes[1,2].text(0.1, 0.9, stats_text, transform=axes[1,2].transAxes,
                      fontsize=10, verticalalignment='top',
                      bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        axes[1,2].set_xlim(0, 1)
        axes[1,2].set_ylim(0, 1)
        axes[1,2].axis('off')
        axes[1,2].set_title('Analysis Summary')
        
        plt.suptitle('Comprehensive DEG Pathway Enrichment Analysis', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.output_dir / 'comprehensive_summary.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        print(f"Saved comprehensive summary to {output_file}")
    
    def generate_report(self):
        """Generate markdown report."""
        report_lines = [
            "# DEG Pathway Enrichment Analysis Report",
            "",
            "## Analysis Parameters",
            f"- Input file: {self.deg_file}",
            f"- Log fold change threshold: {self.logfc_threshold}",
            f"- Adjusted p-value threshold: {self.pval_threshold}",
            f"- Output directory: {self.output_dir}",
            "",
            "## Summary Statistics",
            f"- Total genes analyzed: {len(self.deg_df)}",
            f"- Significant genes: {len(self.significant_genes)}",
            f"- Upregulated genes: {len(self.significant_genes[self.significant_genes['logfoldchanges'] > 0])}",
            f"- Downregulated genes: {len(self.significant_genes[self.significant_genes['logfoldchanges'] < 0])}",
            "",
        ]
        
        # Add pathway results if available
        if 'all' in self.pathway_results and not self.pathway_results['all'].empty:
            results = self.pathway_results['all']
            report_lines.extend([
                "## Pathway Enrichment Results",
                f"- Total significant pathways found: {len(results)}",
                "",
                "### Top 10 Pathways",
                ""
            ])
            
            for i, (_, row) in enumerate(results.head(10).iterrows(), 1):
                report_lines.append(f"{i}. **{row['Database']}**: {row['Term']}")
                report_lines.append(f"   - Adjusted P-value: {row['Adjusted P-value']:.2e}")
                report_lines.append(f"   - Combined Score: {row['Combined Score']:.1f}")
                report_lines.append("")
        
        # Add top genes
        report_lines.extend([
            "## Top 10 Most Significant Genes",
            ""
        ])
        
        for i, (_, row) in enumerate(self.significant_genes.head(10).iterrows(), 1):
            report_lines.append(f"{i}. **{row['names']}**")
            report_lines.append(f"   - Log Fold Change: {row['logfoldchanges']:.2f}")
            report_lines.append(f"   - Adjusted P-value: {row['pvals_adj']:.2e}")
            report_lines.append("")
        
        # Save report
        output_file = self.output_dir / 'analysis_report.md'
        with open(output_file, 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"Saved analysis report to {output_file}")
    
    def run_complete_analysis(self):
        """Run complete pathway enrichment analysis."""
        print("=== Starting Complete DEG Pathway Enrichment Analysis ===")
        
        # 1. Run pathway enrichment
        pathway_results = self.run_pathway_enrichment()
        
        # 2. Analyze gene families
        family_results = self.analyze_gene_families()
        
        # 3. Create plots
        if not pathway_results.empty:
            self.create_pathway_plots()
        
        if family_results:
            self.create_gene_family_plots(family_results)
        
        # 4. Create comprehensive summary
        self.create_comprehensive_summary()
        
        # 5. Generate report
        self.generate_report()
        
        print(f"\n=== Analysis Complete ===")
        print(f"All results saved to: {self.output_dir}")
        print("\nGenerated files:")
        for file in sorted(self.output_dir.glob('*')):
            print(f"  - {file.name}")

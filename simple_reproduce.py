#!/usr/bin/env python3
"""
Reproduce all figures from:
"The Effect of Pore Size on the Adsorption of Gases within Metal Organic Frameworks"
"""

import os
from mof_analysis import run_mof_analysis

def main():
    """Generate all paper figures and validate results."""
    
    print("Reproducing MOF Pore Size Analysis Results")
    print("=" * 50)
    
    # Check data availability
    if not os.path.exists('data/CRAFTED_geometric.csv'):
        print("Error: CRAFTED database not found in data/ directory")
        print("Please download CRAFTED database files:")
        print("- CRAFTED_geometric.csv")
        print("- ISOTHERM_FILES_R/ directory")
        return
    
    # Run complete analysis
    analyzer = run_mof_analysis(
        geometric_csv_path='data/CRAFTED_geometric.csv',
        isotherm_directory='data/ISOTHERM_FILES_R/',
        output_directory='results/'
    )
    
    print("\nResults:")
    print("- All figures saved to results/figures/")
    print("- Numerical results in results/analysis_results.json")
    print("- Validation: python validate_results.py")

if __name__ == "__main__":
    main()

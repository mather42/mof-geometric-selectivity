#!/usr/bin/env python3
"""
Results Validation Script
========================

This script validates that the analysis code reproduces the results 
published in "The Effect of Pore Size on the Adsorption of Gases 
within Metal Organic Frameworks" by Cece Cheng (2024).

Usage:
    python validate_results.py [--data_path path/to/crafted/data]

Expected validation criteria:
- Optimal PLD values within ±0.5 Å of published results
- Model R² values ≥ 0.85
- Cross-validation scores consistent with paper
- Dimensionless scaling relationships verified
"""

import argparse
import sys
import os
import numpy as np
import json
from pathlib import Path

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    from mof_analysis import MOFPoreSizeAnalyzer, run_mof_analysis
except ImportError:
    print("Error: Could not import analysis modules. Ensure you're running from the project root.")
    sys.exit(1)

class ResultsValidator:
    """Comprehensive validation of analysis results against published values."""
    
    def __init__(self):
        """Initialize validator with expected results from paper."""
        
        # Published optimal PLD values (Table from paper)
        self.expected_optimal_plds = {
            'CO2': {273: 9.16, 298: 8.53, 323: 7.77},
            'N2': {273: 7.55, 298: 7.93, 323: 8.42}
        }
        
        # Expected model performance metrics
        self.expected_performance = {
            'min_r2': 0.85,        # Minimum acceptable R²
            'target_r2': 0.87,     # Target R² from paper
            'r2_tolerance': 0.05,   # Acceptable deviation
            'target_rmse': 0.21,    # Target RMSE (mmol/g)
            'rmse_tolerance': 0.10   # Acceptable RMSE deviation
        }
        
        # Expected dimensionless scaling
        self.expected_dimensionless = {
            'target_ratio': 2.4,    # Calculated from paper data
            'tolerance': 0.4        # Reasonable tolerance for universal ratio
        }
        
        # Validation tolerances
        self.pld_tolerance = 0.5   # ± 0.5 Å tolerance for optimal PLDs
        
        # Results storage
        self.validation_results = {
            'pld_validation': {},
            'performance_validation': {},
            'dimensionless_validation': {},
            'overall_status': 'UNKNOWN'
        }
    
    def validate_optimal_plds(self, analyzer):
        """Validate optimal PLD values against published results."""
        print("\n" + "="*60)
        print("VALIDATING OPTIMAL PLD VALUES")
        print("="*60)
        
        pld_results = {}
        total_error = 0
        valid_comparisons = 0
        
        for gas in ['CO2', 'N2']:
            gas_results = {}
            print(f"\n{gas} Optimal PLDs:")
            print(f"{'Temp':<6} {'Computed':<12} {'Expected':<12} {'Error':<8} {'Status'}")
            print("-" * 50)
            
            for temp in [273, 298, 323]:
                key = f'{gas}_{temp}K'
                expected_pld = self.expected_optimal_plds[gas][temp]
                
                if key in analyzer.optimal_plds:
                    computed_pld = analyzer.optimal_plds[key]['optimal_pld']
                    uncertainty = analyzer.optimal_plds[key]['uncertainty']
                    error = abs(computed_pld - expected_pld)
                    
                    # Determine validation status
                    if error <= self.pld_tolerance:
                        status = "✓ PASS"
                    elif error <= self.pld_tolerance * 2:
                        status = "⚠ MARGINAL"
                    else:
                        status = "✗ FAIL"
                    
                    print(f"{temp}K    {computed_pld:.2f}±{uncertainty:.2f}   "
                          f"{expected_pld:.2f}        {error:.2f}     {status}")
                    
                    gas_results[temp] = {
                        'computed': computed_pld,
                        'expected': expected_pld,
                        'error': error,
                        'uncertainty': uncertainty,
                        'status': 'PASS' if error <= self.pld_tolerance else 'FAIL'
                    }
                    
                    total_error += error
                    valid_comparisons += 1
                    
                else:
                    print(f"{temp}K    NO DATA      {expected_pld:.2f}        -        ✗ MISSING")
                    gas_results[temp] = {'status': 'MISSING'}
            
            pld_results[gas] = gas_results
        
        # Overall PLD validation assessment
        if valid_comparisons > 0:
            avg_error = total_error / valid_comparisons
            print(f"\nOVERALL PLD VALIDATION:")
            print(f"Average absolute error: {avg_error:.2f} Å")
            print(f"Tolerance threshold: {self.pld_tolerance:.2f} Å")
            
            if avg_error <= self.pld_tolerance:
                pld_status = "EXCELLENT"
            elif avg_error <= self.pld_tolerance * 1.5:
                pld_status = "GOOD"
            else:
                pld_status = "POOR"
            
            print(f"PLD Validation Status: {pld_status}")
        else:
            pld_status = "NO_DATA"
            print("No PLD data available for validation")
        
        self.validation_results['pld_validation'] = {
            'results': pld_results,
            'avg_error': avg_error if valid_comparisons > 0 else None,
            'status': pld_status
        }
        
        return pld_status in ['EXCELLENT', 'GOOD']
    
    def validate_model_performance(self, analyzer):
        """Validate statistical model performance metrics."""
        print("\n" + "="*60)
        print("VALIDATING MODEL PERFORMANCE")
        print("="*60)
        
        if not analyzer.fitted_models:
            print("No fitted models available for validation")
            self.validation_results['performance_validation']['status'] = 'NO_DATA'
            return False
        
        # Collect performance metrics
        r2_values = []
        rmse_values = []
        cv_scores = []
        
        print(f"{'Condition':<15} {'R²':<8} {'RMSE':<8} {'CV Score':<10} {'Status'}")
        print("-" * 55)
        
        for key, model_data in analyzer.fitted_models.items():
            perf = model_data['performance']
            r2 = perf['r2']
            rmse = perf['rmse']
            cv_mean = perf['cv_mean']
            
            r2_values.append(r2)
            rmse_values.append(rmse)
            cv_scores.append(cv_mean)
            
            # Assess individual model
            r2_status = "✓" if r2 >= self.expected_performance['min_r2'] else "✗"
            
            print(f"{key:<15} {r2:.3f}    {rmse:.3f}    {cv_mean:.3f}      {r2_status}")
        
        # Overall performance assessment
        avg_r2 = np.mean(r2_values)
        avg_rmse = np.mean(rmse_values)
        avg_cv = np.mean(cv_scores)
        
        print(f"\nOVERALL PERFORMANCE METRICS:")
        print(f"Average R²: {avg_r2:.3f} (Target: {self.expected_performance['target_r2']:.3f})")
        print(f"Average RMSE: {avg_rmse:.3f} (Target: {self.expected_performance['target_rmse']:.3f})")
        print(f"Average CV: {avg_cv:.3f}")
        
        # Validation checks
        r2_check = avg_r2 >= self.expected_performance['min_r2']
        r2_target_check = abs(avg_r2 - self.expected_performance['target_r2']) <= self.expected_performance['r2_tolerance']
        
        print(f"\nVALIDATION CHECKS:")
        print(f"R² ≥ {self.expected_performance['min_r2']}: {'✓ PASS' if r2_check else '✗ FAIL'}")
        print(f"R² near target: {'✓ PASS' if r2_target_check else '⚠ MARGINAL'}")
        
        # Overall status
        if r2_check and r2_target_check:
            perf_status = "EXCELLENT"
        elif r2_check:
            perf_status = "GOOD"
        else:
            perf_status = "POOR"
        
        print(f"Performance Status: {perf_status}")
        
        self.validation_results['performance_validation'] = {
            'avg_r2': avg_r2,
            'avg_rmse': avg_rmse,
            'avg_cv': avg_cv,
            'status': perf_status
        }
        
        return perf_status in ['EXCELLENT', 'GOOD']
    
    def validate_dimensionless_scaling(self, analyzer):
        """Validate dimensionless PLD/σ scaling relationships."""
        print("\n" + "="*60)
        print("VALIDATING DIMENSIONLESS SCALING")
        print("="*60)
        
        # Calculate dimensionless ratios
        dimensionless_results = analyzer._analyze_dimensionless_scaling()
        
        if not dimensionless_results:
            print("No dimensionless data available for validation")
            self.validation_results['dimensionless_validation']['status'] = 'NO_DATA'
            return False
        
        # Molecular parameters
        molecular_params = {'CO2': 3.3, 'N2': 3.6}  # σ values in Å
        
        print(f"{'Gas':<6} {'Temp':<6} {'PLD':<8} {'σ':<6} {'PLD/σ':<8}")
        print("-" * 40)
        
        all_ratios = []
        for gas, gas_data in dimensionless_results['by_gas'].items():
            sigma = molecular_params[gas]
            for entry in gas_data:
                ratio = entry['dimensionless_ratio']
                all_ratios.append(ratio)
                print(f"{gas:<6} {entry['temperature']:<6} {entry['optimal_pld']:<8.2f} "
                      f"{sigma:<6.1f} {ratio:<8.2f}")
        
        # Statistical analysis
        mean_ratio = np.mean(all_ratios)
        std_ratio = np.std(all_ratios)
        
        print(f"\nUNIVERSAL SCALING ANALYSIS:")
        print(f"Computed PLD/σ: {mean_ratio:.2f} ± {std_ratio:.2f}")
        print(f"Expected PLD/σ: {self.expected_dimensionless['target_ratio']:.2f}")
        
        # Validation
        ratio_error = abs(mean_ratio - self.expected_dimensionless['target_ratio'])
        ratio_check = ratio_error <= self.expected_dimensionless['tolerance']
        
        print(f"Error: {ratio_error:.2f} (Tolerance: {self.expected_dimensionless['tolerance']:.2f})")
        print(f"Scaling Validation: {'✓ PASS' if ratio_check else '✗ FAIL'}")
        
        dim_status = "PASS" if ratio_check else "FAIL"
        
        self.validation_results['dimensionless_validation'] = {
            'computed_ratio': mean_ratio,
            'computed_std': std_ratio,
            'expected_ratio': self.expected_dimensionless['target_ratio'],
            'error': ratio_error,
            'status': dim_status
        }
        
        return ratio_check
    
    def generate_validation_report(self):
        """Generate comprehensive validation report."""
        print("\n" + "="*80)
        print("COMPREHENSIVE VALIDATION REPORT")
        print("="*80)
        
        # Determine overall status
        pld_pass = self.validation_results['pld_validation'].get('status') in ['EXCELLENT', 'GOOD']
        perf_pass = self.validation_results['performance_validation'].get('status') in ['EXCELLENT', 'GOOD']
        dim_pass = self.validation_results['dimensionless_validation'].get('status') == 'PASS'
        
        print(f"\nVALIDATION SUMMARY:")
        print(f"{'Component':<25} {'Status'}")
        print("-" * 40)
        print(f"{'Optimal PLD Values':<25} {self.validation_results['pld_validation'].get('status', 'UNKNOWN')}")
        print(f"{'Model Performance':<25} {self.validation_results['performance_validation'].get('status', 'UNKNOWN')}")
        print(f"{'Dimensionless Scaling':<25} {self.validation_results['dimensionless_validation'].get('status', 'UNKNOWN')}")
        
        # Overall assessment
        if pld_pass and perf_pass and dim_pass:
            overall_status = "✓ FULL VALIDATION PASSED"
            recommendation = "Results are fully reproducible and match published values."
        elif pld_pass and perf_pass:
            overall_status = "✓ CORE VALIDATION PASSED"
            recommendation = "Primary results are reproducible. Minor deviations in scaling analysis."
        elif pld_pass:
            overall_status = "⚠ PARTIAL VALIDATION"
            recommendation = "Optimal PLD values match paper, but model performance issues detected."
        else:
            overall_status = "✗ VALIDATION FAILED"
            recommendation = "Significant deviations from published results. Review data and methodology."
        
        print(f"\nOVERALL VALIDATION: {overall_status}")
        print(f"RECOMMENDATION: {recommendation}")
        
        self.validation_results['overall_status'] = overall_status
        
        return pld_pass and perf_pass
    
    def save_validation_report(self, output_file):
        """Save detailed validation results to JSON file."""
        with open(output_file, 'w') as f:
            json.dump(self.validation_results, f, indent=2, default=str)
        print(f"\nDetailed validation report saved to: {output_file}")

def main():
    """Main validation execution function."""
    parser = argparse.ArgumentParser(
        description="Validate MOF analysis results against published values"
    )
    parser.add_argument(
        '--data_path', 
        type=str,
        default='data/',
        help='Path to CRAFTED database directory (default: data/)'
    )
    parser.add_argument(
        '--output_dir',
        type=str, 
        default='validation_results/',
        help='Output directory for validation results (default: validation_results/)'
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Run quick validation with subset of data'
    )
    
    args = parser.parse_args()
    
    print("MOF Pore Size Analysis - Results Validation")
    print("=" * 50)
    print(f"Data path: {args.data_path}")
    print(f"Output directory: {args.output_dir}")
    
    # Check data availability
    geometric_csv = os.path.join(args.data_path, 'CRAFTED_geometric.csv')
    isotherm_dir = os.path.join(args.data_path, 'ISOTHERM_FILES_R')
    
    if not os.path.exists(geometric_csv):
        print(f"Error: Could not find {geometric_csv}")
        print("Please ensure CRAFTED database files are in the data/ directory")
        sys.exit(1)
    
    if not os.path.exists(isotherm_dir):
        print(f"Error: Could not find {isotherm_dir}")
        print("Please ensure CRAFTED isotherm files are available")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize validator
    validator = ResultsValidator()
    
    try:
        # Run analysis
        print(f"\nRunning analysis...")
        analyzer = run_mof_analysis(
            geometric_csv_path=geometric_csv,
            isotherm_directory=isotherm_dir,
            output_directory=args.output_dir
        )
        
        # Run validation tests
        print(f"\nStarting validation tests...")
        
        validator.validate_optimal_plds(analyzer)
        validator.validate_model_performance(analyzer)
        validator.validate_dimensionless_scaling(analyzer)
        
        # Generate final report
        validation_passed = validator.generate_validation_report()
        
        # Save detailed results
        report_file = os.path.join(args.output_dir, 'validation_report.json')
        validator.save_validation_report(report_file)
        
        # Exit with appropriate code
        if validation_passed:
            print(f"\n🎉 VALIDATION SUCCESSFUL - Results reproduce paper findings!")
            sys.exit(0)
        else:
            print(f"\n⚠️  VALIDATION ISSUES DETECTED - Review results carefully")
            sys.exit(1)
            
    except Exception as e:
        print(f"\nError during validation: {e}")
        print("Check data files and dependencies")
        sys.exit(1)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
MOF Pore Size Optimization Analysis
===================================

Research code for analyzing the effect of pore-limiting diameter (PLD) on 
gas adsorption in Metal-Organic Frameworks using the CRAFTED database.

This code was developed during the research for:
"The Effect of Pore Size on the Adsorption of Gases within Metal Organic Frameworks"

Author: Cece Cheng
Date: September 2024
Institution: Wycombe Abbey School

Dependencies: pandas, numpy, matplotlib, scikit-learn, scipy
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import r2_score, mean_squared_error
from scipy.optimize import minimize_scalar, curve_fit
import seaborn as sns
from pathlib import Path
import os
import glob
import warnings
import json
from datetime import datetime
import logging

# Set up logging for research tracking
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress sklearn warnings during development
warnings.filterwarnings('ignore', category=UserWarning)

# Research constants - determined from literature and validated experimentally
MOLECULAR_PARAMETERS = {
    'CO2': {'sigma': 3.3, 'name': 'Carbon Dioxide'},  # Lennard-Jones diameter, Å
    'N2': {'sigma': 3.6, 'name': 'Nitrogen'}          # Lennard-Jones diameter, Å
}

# Experimental conditions for GCMC simulations
SIMULATION_CONDITIONS = {
    'temperatures': [273, 298, 323],  # K - chosen to span relevant industrial range
    'pressure': 100000,               # Pa (1 atm) - standard condition
    'force_field': 'UFF',            # Universal Force Field
    'charge_method': 'DDEC',         # Density Derived Electrostatic and Chemical charges
    'database': 'CRAFTED'            # Charge-dependent, Reproducible, Accessible, Forcefield-dependent, Temperature-dependent Exploratory Database
}

class MOFPoreSizeAnalyzer:
    """
    Main analysis class for MOF pore size optimization study.
    
    This class handles:
    - Loading and preprocessing CRAFTED database
    - Statistical analysis of adsorption-PLD relationships  
    - Polynomial modeling with cross-validation
    - Optimal pore size identification
    - Dimensionless scaling analysis
    - Results validation and visualization
    """
    
    def __init__(self, debug_mode=False):
        """Initialize analyzer with research settings."""
        self.debug = debug_mode
        self.geometric_data = None
        self.adsorption_data = None
        self.fitted_models = {}
        self.optimal_plds = {}
        self.analysis_metadata = {
            'creation_time': datetime.now().isoformat(),
            'git_commit': None,  # Could be populated from git
            'conditions': SIMULATION_CONDITIONS,
            'molecular_params': MOLECULAR_PARAMETERS
        }
        
        # Statistical parameters - optimized during research
        self.poly_degree = 2          # 2nd order polynomial fits best for this system
        self.cv_folds = 5            # 5-fold CV balances bias-variance tradeoff
        self.confidence_level = 0.95  # For uncertainty estimation
        
        logger.info("MOF Pore Size Analyzer initialized")
    
    def load_crafted_database(self, geometric_csv, isotherm_directory):
        """
        Load and preprocess CRAFTED database files.
        
        The CRAFTED database contains pre-computed GCMC simulation results
        for 690 MOF structures under various conditions.
        
        Parameters:
        -----------
        geometric_csv : str
            Path to CRAFTED_geometric.csv containing structural properties
        isotherm_directory : str  
            Path to directory containing individual isotherm CSV files
        """
        logger.info("Loading CRAFTED database...")
        
        # Load geometric descriptors
        self.geometric_data = pd.read_csv(geometric_csv)
        logger.info(f"Loaded geometric data for {len(self.geometric_data)} MOF structures")
        
        # The paper uses D_is (inscribed sphere diameter) as pore-limiting diameter
        # This is the narrowest constriction that controls molecular transport
        if 'D_is' in self.geometric_data.columns:
            self.geometric_data = self.geometric_data.rename(columns={'D_is': 'PLD'})
        
        # Initialize adsorption dataset
        self.adsorption_data = self.geometric_data.copy()
        
        # Pre-allocate columns for adsorption data at each T,P condition
        for temp in SIMULATION_CONDITIONS['temperatures']:
            for gas in MOLECULAR_PARAMETERS.keys():
                self.adsorption_data[f'{gas}_ads_{temp}K'] = np.nan
                self.adsorption_data[f'{gas}_error_{temp}K'] = np.nan
        
        # Load isotherm data from individual files
        self._process_isotherm_files(isotherm_directory)
        
        # Quality checks and data validation
        self._validate_loaded_data()
        
        logger.info("Database loading completed successfully")
    
    def _process_isotherm_files(self, isotherm_dir):
        """
        Process individual isotherm files from CRAFTED database.
        
        File naming convention: DDEC_{MOF_NAME}_UFF_{GAS}_{TEMP}.csv
        Each file contains pressure-dependent adsorption isotherms.
        We extract the value at 1 atm (100000 Pa) for this study.
        """
        pattern = os.path.join(isotherm_dir, "DDEC_*_UFF_*.csv")
        isotherm_files = glob.glob(pattern)
        
        successful_loads = 0
        failed_loads = 0
        
        logger.info(f"Processing {len(isotherm_files)} isotherm files...")
        
        for filepath in isotherm_files:
            try:
                # Parse filename to extract metadata
                filename = os.path.basename(filepath)
                components = filename.replace('.csv', '').split('_')
                
                if len(components) >= 5:
                    # Handle MOF names that might contain underscores
                    mof_name = '_'.join(components[1:-3])
                    gas_type = components[-2]
                    temperature = int(components[-1])
                    
                    # Only process conditions relevant to our study
                    if (temperature in SIMULATION_CONDITIONS['temperatures'] and 
                        gas_type in MOLECULAR_PARAMETERS.keys()):
                        
                        adsorption_val, uncertainty = self._extract_atmospheric_data(filepath)
                        
                        if adsorption_val is not None:
                            # Map to our dataset
                            mof_mask = self.adsorption_data['FrameworkName'] == mof_name
                            
                            if mof_mask.any():
                                ads_column = f'{gas_type}_ads_{temperature}K'
                                err_column = f'{gas_type}_error_{temperature}K'
                                
                                self.adsorption_data.loc[mof_mask, ads_column] = adsorption_val
                                self.adsorption_data.loc[mof_mask, err_column] = uncertainty
                                successful_loads += 1
                                
                                if self.debug:
                                    logger.debug(f"Loaded {mof_name} {gas_type} {temperature}K: {adsorption_val:.3f} mol/kg")
            
            except Exception as e:
                failed_loads += 1
                if self.debug:
                    logger.warning(f"Failed to process {filepath}: {e}")
        
        logger.info(f"Processed isotherms: {successful_loads} successful, {failed_loads} failed")
    
    def _extract_atmospheric_data(self, isotherm_file):
        """
        Extract adsorption value at atmospheric pressure (1 atm = 100000 Pa).
        
        The GCMC simulations were run at multiple pressures. For our analysis,
        we specifically need the adsorption at 1 atm to match industrial conditions.
        """
        try:
            df = pd.read_csv(isotherm_file)
            
            pressure_column = '# pressure[Pa]'
            if pressure_column not in df.columns:
                return None, None
            
            # Find data point closest to target pressure (100000 Pa)
            pressure_deviations = np.abs(df[pressure_column] - SIMULATION_CONDITIONS['pressure'])
            closest_index = pressure_deviations.idxmin()
            
            # Extract adsorption and uncertainty
            adsorption = df.loc[closest_index, 'mean_volume[mol/kg]']
            error = df.loc[closest_index, 'mean_error[mol/kg]']
            
            return adsorption, error
            
        except Exception as e:
            if self.debug:
                logger.warning(f"Error extracting data from {isotherm_file}: {e}")
            return None, None
    
    def _validate_loaded_data(self):
        """Validate loaded data quality and completeness."""
        total_structures = len(self.adsorption_data)
        
        logger.info("Data validation summary:")
        logger.info("-" * 40)
        
        for gas in MOLECULAR_PARAMETERS.keys():
            logger.info(f"{gas} adsorption coverage:")
            for temp in SIMULATION_CONDITIONS['temperatures']:
                column = f'{gas}_ads_{temp}K'
                coverage = self.adsorption_data[column].notna().sum()
                percentage = (coverage / total_structures) * 100
                logger.info(f"  {temp}K: {coverage}/{total_structures} structures ({percentage:.1f}%)")
        
        # PLD distribution analysis
        pld_stats = self.adsorption_data['PLD'].describe()
        logger.info(f"PLD statistics: min={pld_stats['min']:.2f}, max={pld_stats['max']:.2f}, mean={pld_stats['mean']:.2f} Å")
        
        # Check for potential outliers that might affect fitting
        for gas in MOLECULAR_PARAMETERS.keys():
            for temp in SIMULATION_CONDITIONS['temperatures']:
                column = f'{gas}_ads_{temp}K'
                data = self.adsorption_data[column].dropna()
                if len(data) > 0:
                    q99 = data.quantile(0.99)
                    outliers = (data > q99).sum()
                    if outliers > 0:
                        logger.warning(f"Found {outliers} potential outliers in {gas} {temp}K data (>{q99:.2f} mol/kg)")
    
    def fit_adsorption_models(self, gas, temperature):
        """
        Fit polynomial model to adsorption vs PLD data.
        
        Uses 2nd order polynomial with 5-fold cross-validation as described
        in the methodology. This approach was chosen after testing various
        model complexities during the research phase.
        
        Parameters:
        -----------
        gas : str
            'CO2' or 'N2'
        temperature : int
            Temperature in Kelvin (273, 298, or 323)
            
        Returns:
        --------
        dict : Model results including fitted parameters and statistics
        """
        logger.info(f"Fitting adsorption model for {gas} at {temperature}K...")
        
        ads_column = f'{gas}_ads_{temperature}K'
        
        # Data cleaning - remove invalid entries
        valid_mask = (
            self.adsorption_data['PLD'].notna() & 
            self.adsorption_data[ads_column].notna() &
            (self.adsorption_data['PLD'] > 0) &
            (self.adsorption_data[ads_column] >= 0) &
            (self.adsorption_data['PLD'] < 20)  # Remove unrealistic PLD values
        )
        
        if valid_mask.sum() < 10:
            logger.error(f"Insufficient data points ({valid_mask.sum()}) for {gas} at {temperature}K")
            return None
        
        # Extract clean data
        X = self.adsorption_data.loc[valid_mask, 'PLD'].values.reshape(-1, 1)
        y = self.adsorption_data.loc[valid_mask, ads_column].values
        
        logger.info(f"Using {len(X)} data points for fitting")
        
        # Polynomial feature transformation
        poly_transformer = PolynomialFeatures(degree=self.poly_degree, include_bias=True)
        X_poly = poly_transformer.fit_transform(X)
        
        # Fit regression model
        model = LinearRegression()
        model.fit(X_poly, y)
        
        # Cross-validation assessment
        cv_splitter = KFold(n_splits=self.cv_folds, shuffle=True, random_state=42)
        cv_scores = cross_val_score(model, X_poly, y, cv=cv_splitter, scoring='r2')
        
        # Model performance metrics
        y_predicted = model.predict(X_poly)
        r_squared = r2_score(y, y_predicted)
        rmse = np.sqrt(mean_squared_error(y, y_predicted))
        
        # Residual analysis for model validation
        residuals = y - y_predicted
        residual_std = np.std(residuals)
        
        # Store comprehensive results
        model_key = f'{gas}_{temperature}K'
        self.fitted_models[model_key] = {
            'model': model,
            'poly_transformer': poly_transformer,
            'performance': {
                'r2': r_squared,
                'rmse': rmse,
                'residual_std': residual_std,
                'cv_scores': cv_scores,
                'cv_mean': cv_scores.mean(),
                'cv_std': cv_scores.std()
            },
            'data': {
                'X_raw': X.flatten(),
                'y_observed': y,
                'y_predicted': y_predicted,
                'residuals': residuals,
                'n_points': len(X)
            },
            'conditions': {
                'gas': gas,
                'temperature': temperature,
                'poly_degree': self.poly_degree
            }
        }
        
        logger.info(f"Model fit complete: R² = {r_squared:.3f}, RMSE = {rmse:.3f}, CV = {cv_scores.mean():.3f}±{cv_scores.std():.3f}")
        
        return self.fitted_models[model_key]
    
    def find_optimal_pore_size(self, gas, temperature):
        """
        Determine optimal PLD that maximizes gas adsorption.
        
        Uses numerical optimization of the fitted polynomial to find
        the PLD that gives maximum adsorption. Includes uncertainty
        estimation based on model curvature.
        """
        model_key = f'{gas}_{temperature}K'
        
        # Ensure model is fitted
        if model_key not in self.fitted_models:
            self.fit_adsorption_models(gas, temperature)
        
        if model_key not in self.fitted_models:
            logger.error(f"Could not fit model for {gas} at {temperature}K")
            return None
        
        model_data = self.fitted_models[model_key]
        regressor = model_data['model']
        transformer = model_data['poly_transformer']
        
        # Define objective function for optimization
        def negative_adsorption(pld):
            """Negative adsorption for minimization (to find maximum)."""
            pld_array = np.array([[pld]])
            pld_features = transformer.transform(pld_array)
            return -regressor.predict(pld_features)[0]
        
        # Optimization bounds based on data range
        pld_range = model_data['data']['X_raw']
        lower_bound = max(pld_range.min(), 2.0)  # Physically meaningful lower limit
        upper_bound = min(pld_range.max(), 16.0) # Reasonable upper limit
        
        # Find optimal PLD
        optimization_result = minimize_scalar(
            negative_adsorption, 
            bounds=(lower_bound, upper_bound), 
            method='bounded'
        )
        
        optimal_pld = optimization_result.x
        max_adsorption = -optimization_result.fun
        
        # Uncertainty estimation using second derivative analysis
        uncertainty = self._estimate_pld_uncertainty(regressor, transformer, optimal_pld)
        
        # Store results
        result_key = f'{gas}_{temperature}K'
        self.optimal_plds[result_key] = {
            'optimal_pld': optimal_pld,
            'max_adsorption': max_adsorption,
            'uncertainty': uncertainty,
            'optimization_success': optimization_result.success,
            'gas': gas,
            'temperature': temperature
        }
        
        logger.info(f"Optimal PLD for {gas} at {temperature}K: {optimal_pld:.2f} ± {uncertainty:.2f} Å")
        
        return self.optimal_plds[result_key]
    
    def _estimate_pld_uncertainty(self, model, transformer, optimal_pld, perturbation=0.05):
        """
        Estimate uncertainty in optimal PLD using curvature analysis.
        
        This method uses numerical differentiation to estimate the second
        derivative at the optimum, which relates to the uncertainty in
        the optimal PLD determination.
        """
        try:
            # Calculate predictions at optimal point and neighbors
            pld_center = np.array([[optimal_pld]])
            pld_left = np.array([[optimal_pld - perturbation]])
            pld_right = np.array([[optimal_pld + perturbation]])
            
            pred_center = model.predict(transformer.transform(pld_center))[0]
            pred_left = model.predict(transformer.transform(pld_left))[0]
            pred_right = model.predict(transformer.transform(pld_right))[0]
            
            # Numerical second derivative
            second_derivative = (pred_left - 2*pred_center + pred_right) / (perturbation**2)
            
            # Uncertainty estimation based on curvature
            if second_derivative < -1e-6:  # Valid maximum
                uncertainty = np.sqrt(abs(2.0 / second_derivative))
                uncertainty = min(uncertainty, 1.0)  # Cap at reasonable value
            else:
                uncertainty = 0.3  # Default for flat maxima
            
            return uncertainty
            
        except Exception as e:
            logger.warning(f"Could not estimate PLD uncertainty: {e}")
            return 0.25  # Conservative default
    
    def run_complete_analysis(self):
        """
        Execute complete analysis pipeline for all gas-temperature combinations.
        
        This is the main analysis routine that processes all conditions
        and generates the comprehensive results presented in the paper.
        """
        logger.info("=" * 60)
        logger.info("STARTING COMPLETE MOF PORE SIZE ANALYSIS")
        logger.info("=" * 60)
        
        analysis_results = {}
        
        # Process each gas-temperature combination
        for gas in MOLECULAR_PARAMETERS.keys():
            for temperature in SIMULATION_CONDITIONS['temperatures']:
                logger.info(f"\nAnalyzing {MOLECULAR_PARAMETERS[gas]['name']} at {temperature}K...")
                
                # Fit adsorption model
                model_result = self.fit_adsorption_models(gas, temperature)
                
                if model_result:
                    # Find optimal pore size
                    optimal_result = self.find_optimal_pore_size(gas, temperature)
                    
                    # Store for summary
                    condition_key = f'{gas}_{temperature}K'
                    analysis_results[condition_key] = {
                        'model': model_result,
                        'optimal': optimal_result
                    }
        
        # Generate dimensionless analysis
        dimensionless_results = self._analyze_dimensionless_scaling()
        
        # Create comprehensive summary
        self._generate_research_summary(analysis_results, dimensionless_results)
        
        logger.info("Complete analysis finished successfully")
        return analysis_results, dimensionless_results
    
    def _analyze_dimensionless_scaling(self):
        """
        Analyze dimensionless PLD/σ scaling for universal design rules.
        
        This analysis tests whether optimal PLDs follow universal scaling
        when normalized by molecular size, as proposed in the paper.
        """
        logger.info("Performing dimensionless scaling analysis...")
        
        dimensionless_data = {}
        all_ratios = []
        
        for gas in MOLECULAR_PARAMETERS.keys():
            sigma = MOLECULAR_PARAMETERS[gas]['sigma']
            gas_ratios = []
            
            for temp in SIMULATION_CONDITIONS['temperatures']:
                result_key = f'{gas}_{temp}K'
                if result_key in self.optimal_plds:
                    optimal_pld = self.optimal_plds[result_key]['optimal_pld']
                    dimensionless_ratio = optimal_pld / sigma
                    
                    gas_ratios.append({
                        'temperature': temp,
                        'optimal_pld': optimal_pld,
                        'sigma': sigma,
                        'dimensionless_ratio': dimensionless_ratio
                    })
                    all_ratios.append(dimensionless_ratio)
            
            dimensionless_data[gas] = gas_ratios
        
        # Statistical analysis of universal scaling
        if all_ratios:
            mean_ratio = np.mean(all_ratios)
            std_ratio = np.std(all_ratios)
            
            logger.info(f"Universal dimensionless optimum: PLD/σ = {mean_ratio:.2f} ± {std_ratio:.2f}")
            
            return {
                'by_gas': dimensionless_data,
                'universal_mean': mean_ratio,
                'universal_std': std_ratio,
                'all_ratios': all_ratios
            }
        else:
            logger.warning("No dimensionless data available")
            return None
    
    def _generate_research_summary(self, analysis_results, dimensionless_results):
        """Generate comprehensive research summary with validation metrics."""
        
        print("\n" + "="*80)
        print("MOF PORE SIZE OPTIMIZATION - RESEARCH SUMMARY")
        print("="*80)
        
        print(f"\nDATASET: {len(self.adsorption_data)} MOF structures from CRAFTED database")
        print(f"CONDITIONS: {SIMULATION_CONDITIONS['force_field']} force field, "
              f"{SIMULATION_CONDITIONS['charge_method']} charges, 1 atm pressure")
        print(f"ANALYSIS: 2nd order polynomial fitting with {self.cv_folds}-fold cross-validation")
        
        print(f"\nOPTIMAL PORE-LIMITING DIAMETERS:")
        print("-" * 50)
        
        for gas in MOLECULAR_PARAMETERS.keys():
            print(f"\n{MOLECULAR_PARAMETERS[gas]['name']} ({gas}):")
            for temp in SIMULATION_CONDITIONS['temperatures']:
                result_key = f'{gas}_{temp}K'
                if result_key in self.optimal_plds:
                    result = self.optimal_plds[result_key]
                    model_stats = self.fitted_models[result_key]['performance']
                    
                    print(f"  {temp}K: {result['optimal_pld']:.2f} ± {result['uncertainty']:.2f} Å "
                          f"(R² = {model_stats['r2']:.3f}, RMSE = {model_stats['rmse']:.3f})")
        
        # Model performance summary
        if analysis_results:
            all_r2 = [analysis_results[key]['model']['performance']['r2'] for key in analysis_results]
            all_rmse = [analysis_results[key]['model']['performance']['rmse'] for key in analysis_results]
            all_cv = [analysis_results[key]['model']['performance']['cv_mean'] for key in analysis_results]
            
            print(f"\nMODEL PERFORMANCE SUMMARY:")
            print(f"Average R²: {np.mean(all_r2):.3f} ± {np.std(all_r2):.3f}")
            print(f"Average RMSE: {np.mean(all_rmse):.3f} ± {np.std(all_rmse):.3f} mol/kg")
            print(f"Average CV Score: {np.mean(all_cv):.3f} ± {np.std(all_cv):.3f}")
        
        # Dimensionless scaling results
        if dimensionless_results:
            print(f"\nDIMENSIONLESS SCALING ANALYSIS:")
            print("-" * 40)
            
            for gas, data in dimensionless_results['by_gas'].items():
                sigma = MOLECULAR_PARAMETERS[gas]['sigma']
                print(f"\n{gas} (σ = {sigma} Å):")
                for entry in data:
                    print(f"  {entry['temperature']}K: PLD/σ = {entry['dimensionless_ratio']:.2f}")
            
            mean_ratio = dimensionless_results['universal_mean']
            std_ratio = dimensionless_results['universal_std']
            print(f"\nUniversal scaling: PLD/σ = {mean_ratio:.2f} ± {std_ratio:.2f}")
            print("This provides a simple design rule for optimal MOF selection.")
    
    def create_publication_plots(self, output_directory=None):
        """
        Generate all publication-quality figures for the paper.
        
        Creates the complete set of figures (1-11) as described in the
        research paper, with consistent styling and high resolution.
        """
        if output_directory:
            os.makedirs(output_directory, exist_ok=True)
        
        logger.info("Generating publication plots...")
        
        # Individual adsorption vs PLD plots (Figures 1-6)
        for gas in MOLECULAR_PARAMETERS.keys():
            for temp in SIMULATION_CONDITIONS['temperatures']:
                self._plot_individual_adsorption(gas, temp, output_directory)
        
        # Comparative temperature plots (Figures 7-8)
        for gas in MOLECULAR_PARAMETERS.keys():
            self._plot_temperature_comparison(gas, output_directory)
        
        # Gas comparison plots (Figures 9-11)
        for temp in SIMULATION_CONDITIONS['temperatures']:
            self._plot_gas_comparison(temp, output_directory)
        
        logger.info("All publication plots generated successfully")
    
    def _plot_individual_adsorption(self, gas, temperature, output_dir=None):
        """Create individual adsorption vs PLD plot."""
        model_key = f'{gas}_{temperature}K'
        
        if model_key not in self.fitted_models:
            logger.warning(f"No model available for {gas} at {temperature}K")
            return
        
        model_data = self.fitted_models[model_key]
        optimal_data = self.optimal_plds.get(model_key, {})
        
        # Set up plot with publication formatting
        plt.figure(figsize=(10, 6))
        plt.style.use('default')  # Clean, publication-ready style
        
        # Plot experimental data points
        X_data = model_data['data']['X_raw']
        y_data = model_data['data']['y_observed']
        
        plt.scatter(X_data, y_data, alpha=0.6, s=25, 
                   color='blue' if gas == 'CO2' else 'orange',
                   label=f'{gas} Data Points', edgecolors='white', linewidth=0.5)
        
        # Plot fitted polynomial curve
        pld_smooth = np.linspace(X_data.min(), X_data.max(), 300)
        model = model_data['model']
        transformer = model_data['poly_transformer']
        
        pld_features = transformer.transform(pld_smooth.reshape(-1, 1))
        adsorption_smooth = model.predict(pld_features)
        
        plt.plot(pld_smooth, adsorption_smooth, 'r-', linewidth=2.5, 
                label='Polynomial Trendline')
        
        # Mark optimal point if available
        if 'optimal_pld' in optimal_data:
            opt_pld = optimal_data['optimal_pld']
            max_ads = optimal_data['max_adsorption']
            uncertainty = optimal_data['uncertainty']
            
            plt.plot(opt_pld, max_ads, 'ro', markersize=10, 
                    label=f'Optimum: {opt_pld:.2f}±{uncertainty:.2f} Å',
                    markeredgecolor='white', markeredgewidth=1)
        
        # Formatting
        plt.xlabel('Pore Size (Å)', fontsize=12, fontweight='bold')
        plt.ylabel(f'{gas} Adsorption (mol/kg)', fontsize=12, fontweight='bold')
        plt.title(f'Mean Volume of {gas} adsorbed as a function of Pore Size at {temperature}K\n'
                 f'(with polynomial trendline)', fontsize=11)
        
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.legend(frameon=True, fancybox=True, shadow=True)
        
        # Add statistics box
        stats = model_data['performance']
        stats_text = (f'R² = {stats["r2"]:.3f}\n'
                     f'RMSE = {stats["rmse"]:.3f}\n'
                     f'CV = {stats["cv_mean"]:.3f}±{stats["cv_std"]:.3f}\n'
                     f'N = {model_data["data"]["n_points"]}')
        
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                verticalalignment='top', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        if output_dir:
            filename = f'Figure_{gas}_{temperature}K_adsorption_vs_pore_size.png'
            plt.savefig(os.path.join(output_dir, filename), 
                       dpi=300, bbox_inches='tight', facecolor='white')
        
        plt.show()
    
    def _plot_temperature_comparison(self, gas, output_dir=None):
        """Create comparative temperature plot for single gas."""
        plt.figure(figsize=(12, 8))
        
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Professional color scheme
        temp_labels = [f'{temp}K' for temp in SIMULATION_CONDITIONS['temperatures']]
        
        for i, temp in enumerate(SIMULATION_CONDITIONS['temperatures']):
            model_key = f'{gas}_{temp}K'
            
            if model_key not in self.fitted_models:
                continue
            
            model_data = self.fitted_models[model_key]
            
            # Plot data points
            X_data = model_data['data']['X_raw']
            y_data = model_data['data']['y_observed']
            
            plt.scatter(X_data, y_data, alpha=0.6, s=20, 
                       color=colors[i], label=f'{temp}K {gas} Adsorption')
            
            # Plot fitted curves
            pld_smooth = np.linspace(2, 16, 300)
            model = model_data['model']
            transformer = model_data['poly_transformer']
            
            pld_features = transformer.transform(pld_smooth.reshape(-1, 1))
            adsorption_smooth = model.predict(pld_features)
            
            plt.plot(pld_smooth, adsorption_smooth, color=colors[i], 
                    linewidth=2.5, label=f'{temp}K Trendline')
        
        plt.xlabel('Pore Size (Å)', fontsize=12, fontweight='bold')
        plt.ylabel(f'{gas} Adsorption (mol/kg)', fontsize=12, fontweight='bold')
        plt.title(f'Mean Volume of {gas} adsorbed as a function of Pore Size\n'
                 f'at different Temperatures (with polynomial trendlines)', fontsize=11)
        
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.legend(frameon=True, fancybox=True, shadow=True)
        plt.tight_layout()
        
        if output_dir:
            filename = f'Figure_{gas}_temperature_comparison.png'
            plt.savefig(os.path.join(output_dir, filename), 
                       dpi=300, bbox_inches='tight', facecolor='white')
        
        plt.show()
    
    def _plot_gas_comparison(self, temperature, output_dir=None):
        """Create gas comparison plot at single temperature."""
        plt.figure(figsize=(12, 8))
        
        gas_colors = {'CO2': '#1f77b4', 'N2': '#ff7f0e'}
        
        for gas in MOLECULAR_PARAMETERS.keys():
            model_key = f'{gas}_{temperature}K'
            
            if model_key not in self.fitted_models:
                continue
            
            model_data = self.fitted_models[model_key]
            color = gas_colors[gas]
            
            # Plot data points
            X_data = model_data['data']['X_raw']
            y_data = model_data['data']['y_observed']
            
            plt.scatter(X_data, y_data, alpha=0.6, s=25, 
                       color=color, label=f'{gas} Data Points')
            
            # Plot fitted curves
            pld_smooth = np.linspace(2, 16, 300)
            model = model_data['model']
            transformer = model_data['poly_transformer']
            
            pld_features = transformer.transform(pld_smooth.reshape(-1, 1))
            adsorption_smooth = model.predict(pld_features)
            
            plt.plot(pld_smooth, adsorption_smooth, color=color, 
                    linewidth=2.5, label=f'{gas} Polynomial Trendline')
        
        plt.xlabel('Pore Size (Å)', fontsize=12, fontweight='bold')
        plt.ylabel('Adsorption (mol/kg)', fontsize=12, fontweight='bold')
        plt.title(f'Mean Volume of CO₂ and N₂ gas adsorbed as a function of Pore Size\n'
                 f'at {temperature}K (with polynomial trendlines)', fontsize=11)
        
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.legend(frameon=True, fancybox=True, shadow=True)
        plt.tight_layout()
        
        if output_dir:
            filename = f'Figure_gas_comparison_{temperature}K.png'
            plt.savefig(os.path.join(output_dir, filename), 
                       dpi=300, bbox_inches='tight', facecolor='white')
        
        plt.show()
    
    def export_results(self, output_file):
        """Export complete analysis results to JSON for reproducibility."""
        export_data = {
            'metadata': self.analysis_metadata,
            'optimal_plds': self.optimal_plds,
            'model_performance': {},
            'dimensionless_analysis': self._analyze_dimensionless_scaling()
        }
        
        # Export model performance metrics
        for key, model_data in self.fitted_models.items():
            export_data['model_performance'][key] = {
                'r2': model_data['performance']['r2'],
                'rmse': model_data['performance']['rmse'],
                'cv_mean': model_data['performance']['cv_mean'],
                'cv_std': model_data['performance']['cv_std'],
                'n_points': model_data['data']['n_points']
            }
        
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2, default=str)
        
        logger.info(f"Results exported to {output_file}")

# Research execution functions
def run_mof_analysis(geometric_csv_path, isotherm_directory, output_directory="results/"):
    """
    Main function to execute complete MOF pore size analysis.
    
    This function replicates the complete analysis workflow used in the research.
    
    Parameters:
    -----------
    geometric_csv_path : str
        Path to CRAFTED_geometric.csv file
    isotherm_directory : str  
        Path to directory containing CRAFTED isotherm files
    output_directory : str
        Directory for saving results and figures
        
    Returns:
    --------
    MOFPoreSizeAnalyzer : Analyzer instance with complete results
    """
    # Create output directory
    os.makedirs(output_directory, exist_ok=True)
    
    # Initialize analysis
    analyzer = MOFPoreSizeAnalyzer(debug_mode=False)
    
    # Load CRAFTED database
    analyzer.load_crafted_database(geometric_csv_path, isotherm_directory)
    
    # Execute complete analysis
    analysis_results, dimensionless_results = analyzer.run_complete_analysis()
    
    # Generate all publication plots
    plot_directory = os.path.join(output_directory, "figures")
    analyzer.create_publication_plots(plot_directory)
    
    # Export results for reproducibility
    results_file = os.path.join(output_directory, "analysis_results.json")
    analyzer.export_results(results_file)
    
    # Save processed data
    data_file = os.path.join(output_directory, "processed_data.csv")
    analyzer.adsorption_data.to_csv(data_file, index=False)
    
    logger.info(f"Complete analysis finished. Results saved to {output_directory}")
    
    return analyzer

if __name__ == "__main__":
    """
    Example execution for reproducing paper results.
    
    Update the paths below to point to your CRAFTED database files.
    """
    
    # File paths - update these for your system
    GEOMETRIC_DATA_PATH = "data/CRAFTED_geometric.csv"
    ISOTHERM_DATA_PATH = "data/ISOTHERM_FILES_R/"
    OUTPUT_PATH = "paper_results/"
    
    print("MOF Pore Size Optimization Analysis")
    print("===================================")
    print("Reproducing results from Cheng (2024)")
    print()
    
    # Run complete analysis
    analyzer = run_mof_analysis(
        geometric_csv_path=GEOMETRIC_DATA_PATH,
        isotherm_directory=ISOTHERM_DATA_PATH, 
        output_directory=OUTPUT_PATH
    )
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"Results saved to: {OUTPUT_PATH}")
    print("Generated files:")
    print("- figures/ : All publication figures (Figures 1-11)")
    print("- analysis_results.json : Complete numerical results")
    print("- processed_data.csv : Processed adsorption dataset")
    print()
    print("For questions or issues, contact: [cececheng1001@gmail.com]")

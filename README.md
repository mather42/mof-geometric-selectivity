# MOF Pore Size Optimization Analysis

Analysis code for: *"The Effect of Pore Size on the Adsorption of Gases within Metal Organic Frameworks"* (Cheng, 2024)

## Quick Start

```bash
# Install dependencies
pip install numpy pandas scipy scikit-learn matplotlib seaborn

# Download CRAFTED database to data/ directory
# Run analysis
python mof_analysis.py
```

## Key Results

**Optimal Pore-Limiting Diameters:**
- CO₂: 9.16 Å (273K), 8.53 Å (298K), 7.77 Å (323K)
- N₂: 7.55 Å (273K), 7.93 Å (298K), 8.42 Å (323K)

**Model Performance:** R² = 0.87 ± 0.03, RMSE = 0.21 mmol g⁻¹

## Validation

```bash
python validate_results.py
```

Expected output: All optimal PLDs within ±0.5 Å of paper values.


## 📂 Repository Structure

<repo-name>/
├── environment.yml           # Conda environment specification
├── requirements.txt          # pip install specification
├── LICENSE                   # License file
├── README.md                 # This file
├── data/                     # Input data (not checked in)
│   ├── CRAFTED_geometric.csv
│   └── ISOTHERM_FILES_R/
├── src/                      # Core analysis code
├── mof_analysis-2.py          
├── validate_results.py      # Results verification  
└── reproduce_figures.py     # Simple figure reproduction          


## 📊 Data & Code Availability
	•	CRAFTED database: https://doi.org/10.5281/zenodo.8190237
	•	CoRE MOF 2019: https://doi.org/10.5281/zenodo.3370144
	•	Environment specification:
	•	Conda: environment.yml
	•	pip: requirements.txt

⸻

## 📜 License

This project is licensed under the MIT License.

⸻

## 📬 Contact

Cece Cheng
Wycombe Abbey School, High Wycombe, UK
ORCID: 0009-0001-2481-328X
Email: cececheng1001@gmail.com

⸻

Feel free to open an issue or pull request on GitHub if you encounter any problems or have suggestions!


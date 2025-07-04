# MOF Pore-Size Optimization Analysis

**“Geometric Predictors of CO₂/N₂ Selectivity in Metal–Organic Frameworks” by Cece Cheng**  
Repository containing all scripts, data pointers, and instructions to reproduce the analyses and figures from the paper.

---

## 🚀 Quick Start

 ### 1. Clone the repository  

   git clone https://github.com/mather42/mof-geometric-selectivity.git  
cd mof-geometric-selectivity

### 2. Set up your environment 

Conda (recommended):

conda env create -f environment.yml
conda activate mof-adsorption

Or pip:

pip install -r requirements.txt

### 3.	Run the main analysis

python src/mof_analysis.py \
  --geometric data/CRAFTED_geometric.csv \
  --isotherms data/ISOTHERM_FILES_R/ \
  --out results/

### 4.	Reproduce paper figures

python examples/reproduce_figures.py



⸻

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
│   ├── mof_analysis.py


⸻

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


# Gaussian Process Spatial Interpolation — Code for Paper

This project uses Gaussian Process regression (scikit-learn) with geospatial data
(geopandas/shapely) to interpolate and visualize river measurement data. The main
analysis lives in the notebook **`main_csvInput.ipynb`**.

## Project structure

```
Code for Paper/
├── main_csvInput.ipynb   # Main analysis notebook
├── data/                 # Input CSV data files
├── requirements.txt      # Python dependencies
├── .gitignore
└── README.md
```

## Setup (Windows / PowerShell)

Requires Python 3.11 (any 3.10–3.12 should work). From this folder:

```powershell
# 1. Create the virtual environment (only needed once).
#    Keep it OUTSIDE synced folders like OneDrive — e.g. C:\venvs\gp-paper
py -3.11 -m venv C:\venvs\gp-paper

# 2. Activate it
C:\venvs\gp-paper\Scripts\Activate.ps1

# 3. Install dependencies (use requirements-lock.txt for exact versions)
python -m pip install -r requirements.txt
```

> **Tip:** If this project folder lives in OneDrive, do not create the
> environment inside it — OneDrive will try to sync ~20,000 package files,
> which slows Python down and wastes sync bandwidth.

If PowerShell blocks activation, run:
`Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`

### Setup (macOS / Linux)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Running the notebook

- **In Cursor / VS Code:** open `main_csvInput.ipynb` and select the kernel named
  `Python (.venv) GP Paper` (or point the kernel picker at
  `C:\venvs\gp-paper\Scripts\python.exe`).
- **In Jupyter:** with the environment activated, run `jupyter lab` (or
  `jupyter notebook`) and open the notebook.

The environment is registered as a named Jupyter kernel during setup, so it should
appear in the kernel list automatically. If not, register it manually:

```powershell
python -m ipykernel install --user --name gp-code-for-paper --display-name "Python (.venv) GP Paper"
```

## Data files

Input data is stored in `data/`:

| File | Description |
|---|---|
| `25_ecohydroData.csv` | 2025 ecohydrology measurements |
| `26_ecohydroData.csv` | 2026 ecohydrology measurements |
| `25_depthSensorData.csv` | 2025 depth sensor time series |
| `26_depthSensorData.csv` | 2026 depth sensor time series |
| `25_ADCPData.csv` | 2025 ADCP (SonTek) velocity data |

## Key dependencies

- **scikit-learn** — Gaussian Process regression (RBF, RationalQuadratic, White kernels)
- **geopandas / shapely / pyproj** — geospatial data handling and river geometry
- **pandas / numpy** — data wrangling
- **matplotlib** — plotting
- **ipywidgets** — interactive sliders in the notebook

Generated outputs (`complete_data_export.csv`, `gp_spatial_heatmap_results.csv`)
are written to the project folder and are git-ignored.

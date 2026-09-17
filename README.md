# Gaussian Process Spatial Interpolation — Code for Paper

This project uses Gaussian Process regression (scikit-learn) with geospatial data
(geopandas/shapely) to interpolate and visualize velocity measurements along three
Huron River reaches: **Upper**, **Middle**, and **Lower**.

The analysis is split into three stages, run in order:

1. **`dataSorting.ipynb`** — merge field data, attach upstream stage, trim, and write per-section CSVs / GeoJSONs
2. **`models/gpModel_{Section}.ipynb`** — train a 3-input GP and pickle the fitted model
3. **`model visualizations/modelVis_{Section}.ipynb`** — spatial maps, interactive sliders, and threshold analysis

## Project structure

```
Code for Paper/
├── dataSorting.ipynb              # Stage 1: clean and georeference field data
├── data/                          # Raw CSVs (ecohydro, depth sensors, ADCP)
├── sorted_data/                   # Per-section training tables from dataSorting
├── mapping/                       # Centerlines, survey points, bank/channel geometry
├── models/                        # Stage 2: GP training notebooks
│   ├── gpModel_Upper.ipynb
│   ├── gpModel_Middle.ipynb
│   └── gpModel_Lower.ipynb
├── model outputs/                 # Pickled GP models written by the training notebooks
├── model visualizations/          # Stage 3: maps and threshold figures
│   ├── modelVis_Upper.ipynb
│   ├── modelVis_Middle.ipynb
│   └── modelVis_Lower.ipynb
├── CrossSiteVisual/               # Separate network-scale WSE / rainfall figures
├── Claude outputs/                # Helper scripts and check figures for geometry
├── requirements.txt
├── requirements-lock.txt
├── .gitignore
└── README.md
```

Notebooks in `models/` and `model visualizations/` sit one level below the project
root, so they load files with `../` (for example `../sorted_data/...` and
`../model outputs/...`). Jupyter’s working directory is the notebook’s own folder,
so run those notebooks from where they live. `dataSorting.ipynb` stays at the
project root and uses paths like `data/` and `mapping/` with no `../`.

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

## Running the notebooks

- **In Cursor / VS Code:** open the notebook and select the kernel named
  `Python (.venv) GP Paper` (or point the kernel picker at
  `C:\venvs\gp-paper\Scripts\python.exe`).
- **In Jupyter:** with the environment activated, run `jupyter lab` (or
  `jupyter notebook`) and open the notebook.

The environment is registered as a named Jupyter kernel during setup, so it should
appear in the kernel list automatically. If not, register it manually:

```powershell
python -m ipykernel install --user --name gp-code-for-paper --display-name "Python (.venv) GP Paper"
```

### Suggested order

1. Run **`dataSorting.ipynb`** from the project root. This writes
   `sorted_data/complete_data_{Section}.csv` and `mapping/data_{Section}.geojson`,
   and computes distance-to-centerline / downstream distance.
2. Run **`models/gpModel_{Section}.ipynb`** for each reach. This trains the GP
   and writes `model outputs/gp_model_{Section}.pkl`.
3. Run **`model visualizations/modelVis_{Section}.ipynb`** to load that pickle
   and produce spatial maps, interactive upstream-depth sliders, and
   probabilistic threshold panels.

## Data files

Raw inputs live in `data/`:

| File | Description |
|---|---|
| `25_ecohydroData.csv` | 2025 ecohydrology measurements |
| `26_ecohydroData.csv` | 2026 ecohydrology measurements |
| `25_depthSensorData_Middle.csv` | 2025 upstream depth, Middle |
| `25_depthSensorData_Lower.csv` | 2025 upstream depth, Lower |
| `26_depthSensorData_Upper.csv` | 2026 upstream depth, Upper |
| `26_depthSensorData_Middle.csv` | 2026 upstream depth, Middle |
| `26_depthSensorData_Lower.csv` | 2026 upstream depth, Lower |
| `25_ADCPData.csv` | 2025 ADCP (SonTek) velocity data |

`sorted_data/` holds the per-section tables used for training
(`complete_data_Upper.csv`, `complete_data_Middle.csv`, `complete_data_Lower.csv`).

`mapping/` holds river geometry (NAD 1983 StatePlane Michigan South, ESRI:102690):

| File | Description |
|---|---|
| `Huron_RiverSystem.geojson` | Full river-system centerline |
| `centerline_union_{Section}.geojson` | Section main-stem centerline |
| `data_{Section}.geojson` | Survey points for that section |
| `bank_lines_{Section}.geojson` | Bank / waterline geometry |
| `channel_extent_{Section}.geojson` | Streambed / channel polygon |

`model outputs/gp_model_{Section}.pkl` is the trained 3-input GP (downstream
distance, distance to centerline, upstream depth → velocity) plus scalers and
cross-validation metrics.

## Other folders

- **`CrossSiteVisual/`** — network-scale water-surface elevation and rainfall
  figures (ridgeline, overlay, publication small-multiples). These scripts talk
  to InfluxDB / Open-Meteo and are separate from the GP interpolation.
- **`Claude outputs/`** — standalone helper scripts and check plots used while
  building bank lines, streambed extent, and related figures.

## Key dependencies

- **scikit-learn** — Gaussian Process regression (RBF, RationalQuadratic, White kernels)
- **geopandas / shapely / pyproj** — geospatial data handling and river geometry
- **pandas / numpy** — data wrangling
- **matplotlib / ipywidgets** — plotting and interactive sliders
- **scipy** — spatial / statistical helpers in the visualization notebooks

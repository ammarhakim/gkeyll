# kin-rad-gkeyll
**automated kinetic radiation parameter fitting for gkeyll.**

> *note: designed by jaime caballero for a phd project at differ and pppl. the scripts in this repository were created with the help of ai.*

this repository houses the complete data ingestion and non-linear optimization pipeline for calculating plasma radiation parameters. it takes raw bremsstrahlung and line radiation data from openadas, fits it to the analytical equations defined by j. roeltgen, and compiles the results into a lightweight database ready for the gkeyll kinetic simulation framework.

this is the **development and validation** repository. it contains the raw data, the master tracking databases, the physical plotting outputs, and the manual tuning sandboxes. the final, compiled output of this pipeline is exported to a lighter repository (`prod-kin-rad-gkeyll`) which is integrated directly into gkeyll.

---

### filesystem architecture

the repository is strictly separated into ingestion, batch processing, and testing environments:

```text
kin-rad-gkeyll/
├── raw_data/                  # raw openadas .dat files (git-ignored)
├── formatted_data/            # roeltgen-formatted 4-column .txt files
├── fits_data/                 # the master output directory
│   ├── roeltgen_data/         # legacy gkeyll database for visual benchmarking
│   ├── fit-db_<ID>/           # output folder for a specific batch run
│   │   ├── master-fit_db_<ID>.txt   # csv database of all fits and metadata
│   │   ├── gkeyll-db_fit_<ID>.txt   # the compiled gkeyll-ready text file
│   │   └── plots/             # dual log-log/linear-log pngs of every fit
│   └── individual_fits/       # sandbox directory for manually tuned fits (TBD)
├── src/                       # core python modules
│   ├── data_parser.py         # reads adas formats
│   ├── error_analysis.py      # handles zone-based fractional error physics
│   ├── optimizer_core.py      # handles the scipy/matlab integration wrapper
│   └── format_adas.py         # translator for raw .dat to formatted .txt
├── test/                      
│   └── check_fit.py           # cli viewer for instant plot/stat retrieval (TBD)
├── download_data.py           # pipeline entry point 1: ingestion
├── fit_batch.py               # pipeline entry point 2: mass optimization
├── fit_single.py              # pipeline entry point 3: manual sandbox (TBD)
├── swap_fits.py               # tool to overwrite batch fits with manual ones (TBD)
└── README.md
```

### pipeline usage

running the pipeline is a straightforward two-step process: fetch the data, then fire off the batch optimizer.

#### 1. data ingestion (`download_data.py`)
automatically fetches raw openadas radiation data and translates it into our strict formatted `.txt` grids.

```bash
# download and format all supported elements
python download_data.py
# or just target a specific element
python download_data.py --element He
```

#### 2. batch optimization 
fit_batch.py it sweeps through the $N_e$ grids, fits the radiation curves to roeltgen's analytical models, and dynamically saves the best parameters and plots to your fits_data/ directory.

```bash
# fit specific elements and assign a unique dataset id
python fit_batch.py --id run_1 --elements He,Li --optimizer trust-constr
```

```bash
# fit everything, and optionally chop off complex low-temperature tails (e.g., below 1.5 eV)
python fit_batch.py --id run_1 --elements all --min-te 1.5
```

#### 3. manual tuning & verification (TBD)

if automated fits aren't always perfect, this repository includes a "human-in-the-loop" sandbox to safely fix stubborn curves. sorry for the ai lingo this will change.

- `fit_single.py`: manually dial in weights and $V_0$ bounds for a highly specific $N_e$ / $T_e$ slice.
- `swap_fits.py`: surgically replace a bad batch fit with your handcrafted fit inside the master database.
- `test/check_fit.py`: instantly pull up error stats and visual plots for any fit ID directly from the terminal.
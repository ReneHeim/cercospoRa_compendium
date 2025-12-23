# Spatially explicit negative prognosis of Cercospora leaf spot

This repository is the open research compendium for the article *Spatially explicit negative prognosis of Cercospora leaf spot epidemics by process-based integration of leaf area index from remote sensing* (“cercospoRa” study). It contains the data, code, and documentation required to reproduce the analyses, figures, and tables in the manuscript.  

> **Status:** Under review at *[Journal name]*.  
> **Corresponding author:** Rene H. J. Heim (University of Göttingen).

---

## Quick links

- 📄 **Manuscript and supplements:** `docs/` (PDF and supporting material once available)  
- 📊 **Data:** `data/` (raw and example data; large files to be downloaded separately)  
- 🧮 **Analysis code:** `code/` (R and Python scripts / Rmds, ordered by prefix to indicate run order)  
- 📦 **R package:** [`cercospoRa`](https://cran.r-project.org/package=cercospoRa) – implementation of the CLS negative prognosis framework  

---

## Repository structure

The compendium follows a simplified research-compendium layout inspired by published guidelines on packaging data-analytical work reproducibly.  

- `data/`  
  - Raw and/or example data used in the analyses.  
  - Large files are not tracked in Git; download them into this directory following the instructions in `data/README.md` (if present).  

- `docs/`  
  - Manuscript, supplementary information, and any rendered HTML or PDF outputs.  

- `code/`  
  - R and Python scripts and R Markdown files to reproduce the analyses and visualisations.  
  - Files are numbered to indicate the recommended execution order (e.g. `01_*.R`, `02_*.R`, …).

---

## Reproducing the analysis

1. **Clone the repository**

```
git clone https://github.com/ReneHeim/cercospoRa_compendium.git
cd cercospoRa_compendium
```

2. **Obtain data**

Download the required raw data from HERE and copy them into `data/` (see the `data/` subfolder or the manuscript Data availability statement for links and access details).  

3. **Project root and paths**

All scripts assume the project root is the working directory, best achieved by creating an R project in the `cercospoRa_compendium/` directory, and use the `here` package (in R) to resolve paths for portability and reproducibility. 

<img src="https://github.com/user-attachments/assets/1c0e4b0e-0d93-4ce2-9fb0-b2a67a601a05"
     alt="rproject" style="width: 50%; height: auto;">

4. **Run the workflow**

Execute the scripts / Rmd files in `code/` in numerical order (e.g. `01_*.R` → `02_*.R` → …).  
Each script is documented with its inputs, outputs, and expected runtime.

---

## Software requirements

The analysis pipeline uses both R and Python.

### R environment

The analyses were run with **R 4.5.0**.  

Install the required packages (including the `cercospoRa` package implementing the CLS negative prognosis model) with:

```
install.packages(
c("cercospoRa", # CLS negative prognosis model
"readr",
"sf", # Spatial data manipulation
"terra", # Raster data manipulation
"RANN", # Nearest neighbour estimation
"flextable", # Tables
"officer", # Exporting tables/figures to Word
"data.table", # Data wrangling
"ggplot2", # Plotting
"mgcv", # GAMs
"here" # Project-rooted file paths
)
)

```

> For exact package versions, see the session information recorded in the analysis scripts / Rmd files.

### Python environment

Python code was run with **Python 3.11.0** in a Conda environment.  

Create and activate the environment with:

```
conda create -n cercospora_env -c conda-forge python=3.11
numpy pandas matplotlib seaborn geopandas rasterio rioxarray shapely
conda activate cercospora_env
```

Key Python libraries:

- `numpy` – Numerical computing  
- `pandas` – Tabular data handling  
- `matplotlib`, `seaborn` – Plotting and statistical visualisation  
- `geopandas`, `shapely` – Vector geospatial analysis  
- `rasterio`, `rioxarray` – Raster I/O and geospatial raster handling  

---

## Licenses

- **Manuscript text and figures:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)  
- **Code in this repository:** [MIT License](https://opensource.org/licenses/MIT)  
- **Data:** [CC0](https://creativecommons.org/publicdomain/zero/1.0/) (attribution requested when reusing).  

Please also cite the associated article when using this compendium in your own work (see below).

---

## How to cite

If you use this compendium or the `cercospoRa` package, please cite:

> Heim RHJ, *et al.* *Spatially explicit negative prognosis of Cercospora leaf spot epidemics by process-based integration of leaf area index from remote sensing.* *[Journal]*, [year], [volume]:[pages]. DOI: [DOI].

You may also wish to cite general guidance on research compendia, for example:

> Marwick B, Boettiger C, Mullen L (2018). “Packaging data analytical work reproducibly using R (and friends).” *The American Statistician* 72(1):80–88.  

---

## Acknowledgements

This compendium was inspired by existing research-compendium templates and workflows, including:

- Research compendium guidelines and examples by Ben Marwick and collaborators  
- The `manuscriptPackage` template by J. Hollister  
- Templates by Carl Boettiger and other open-science practitioners  
- *Reproducibility in Plant Pathology* materials  
- *R for Data Science*  

---

## Contact

Rene H. J. Heim  
Institute of Geography – Cartography, GIS and Remote Sensing  
University of Göttingen, Germany  

Email: [rene.heim@uni-goettingen.de](mailto:rene.heim@uni-goettingen.de)

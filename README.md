# cercospoRa_compendium
This repository contains data and code required to reproduce the analysis of the 
article *Spatially explicit negative prognosis of Cercospora leaf spot epidemics by process-based integration of leaf area index from remote sensing*. 

## Overview of contents

There are three basic folders:

- `data/` - raw and/or example data. - *Download large files to this location*
- `docs/` - published manuscript and supplementary information.
- `code/` - code to reproduce the analysis and visualizations of the associated manuscript.

## How to use?

Download the raw [data](https://github.com/ReneHeim/cercospoRa_compendium/tree/main/data) and run the provided [code](https://github.com/ReneHeim/cercospoRa_compendium/tree/main/code) in the sequence the files are named. 

## Requirements  
The analysis pipeline was undertaken in R and Python.

### R requirements 
The following R packages and was run on R version 4.5.0

#### CRAN packages  

```r
install.packages("cercospoRa", # cercospoRa negative prognosis model
                 "readr", 
                 "sf",         # spatial data manipulation
                 "terra",      # spatial data manipulation
                 "RANN",       # Nearest neighbour estimations
                 "flextable",  # presenting data in tables
                 "officer",    # exporting tables and figures to Microsoft Word
                 "data.table", # data rectangling / wrangling
                 "ggplot2,     # plotting data 
                 "mgcv",       # statistical analysis with GAM
                 "here")       # provides project folder path for reproducability
```

### Python requirements 
The following Python libraries were run on Python version 3.11.0

#### Python libraries (via conda)

```p
conda create -n my_analysis_env -c conda-forge python=3.11 \
    numpy pandas matplotlib seaborn geopandas rasterio rioxarray shapely
```

```p
seaborn         # Statistical data visualization
matplotlib      # Plotting and visualization
numpy           # Numerical computing
geopandas       # Geospatial data analysis
rasterio        # Raster data I/O and processing
rioxarray       # Geospatial raster data with xarray
shapely         # Geometry manipulation
pandas          # Data manipulation and analysis
```
## Licenses

Manuscript: [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)

Code: [MIT](https://opensource.org/licenses/MIT)

Data: [CC-0](https://creativecommons.org/publicdomain/zero/1.0/) attribution 
requested in reuse

## Credits

This compendium was developed based on the ideas and examples in other resources 
which suggest a research compendium structured as an R package. 
However, this is not an R package as we focused mainly on four main Rmd files 
that are used to generate the html output of the compendium, but there there is 
clear separation of the data and the output files.

[Marwick et al. 2017](https://doi.org/10.7287/peerj.preprints.3192v1)

[jhollist manuscriptpackage](https://github.com/jhollist/manuscriptPackage)

[cboettig template](https://github.com/cboettig/template)

[benmarwick researchcompendium](https://github.com/benmarwick/researchcompendium)

[Reproducibility_in_Plant_Pathology](Reproducibility_in_Plant_Pathology)

[R for Data Science](http://r4ds.had.co.nz/)

## Contact

Rene HJ Heim, Institute of Geography | Cartography, GIS and Remote Sensing | University of Goettingen | Germany
Email: rene.heim@uni-goettingen.de



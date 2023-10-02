# List of required packages
required_packages <- c(
  "calibrate", "dplyr", "FactoMineR", "factoextra", "geojsonio", 
  "GGally", "ggbiplot", "ggpubr", "ggrepel", "glmnet", "grid", 
  "gridExtra", "Hmisc", "kableExtra", "lctools", "lubridate", 
  "mctest", "olsrr", "parallel", "RandomFields","RColorBrewer", "reshape",
  "sf", "sp", "spatialreg", "spdep", "StatMeasures", "tidyverse"
)

# Check if each package is installed, and if not, install it
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Load the required packages
library("calibrate")
library("dplyr")
library("FactoMineR")
library("factoextra")
library("geojsonio")
library("GGally")
library("ggbiplot")
library("ggpubr")
library("ggrepel")
library("glmnet")
library("grid")
library("gridExtra")
library("Hmisc")
library("kableExtra")
library("lctools")
library("lubridate")
library("mctest")
library("olsrr")
library("parallel")
library("RandomFields")
library("RColorBrewer")
library("reshape")
library("sf")
library("sp")
library("spatialreg")
library("spdep")
library("StatMeasures")
library("tidyverse")

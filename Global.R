###########################################################################

# Application of synthetic touch DNA deposits in forensic research

# Emma Hermannsdóttir(1), Virginie Galais(2), Alexander Gray(2), Chris Gannicliffe(3), Niamh Nic Daeid(2), Agnieszka Kuffel(2)*

# (1) Department of Immunology, Genetics, and Pathology (IGP), Uppsala University, Uppsala, Sweden
# (2) Leverhulme Research Centre for Forensic Science, Department of Science and Engineering, University of Dundee, Dundee, DD1 4HN, UK
# (3) Scottish Police Authority Forensic Services, Aberdeen Laboratory, Aberdeen, AB24 5EQ, UK
# * Correspondence: akuffel001@dundee.ac.uk 

# Keywords: XXX

# Website: https://github.com/LRCFS/
# Contact: lrc@dundee.ac.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

###########################################################################
# To clean the Global environment
rm(list=ls()) 

#############################################################
#####                library requirement                #####
#############################################################

library(plyr)        # Tools for organising and summarising datasets (older package, mostly replaced by dplyr)
library(dplyr)       # Grammar of data wrangling for organising, filtering, and summarising datasets
library(tidyverse)   # Collection of packages for data analysis and visualisation (includes ggplot2, dplyr, tidyr, readr, etc.)
library(ggplot2)     # Create static data visualisations using the grammar of graphics
library(extrafont)   # Import, load, and use system fonts in plots
library(RColorBrewer)# Predefined colour palettes for plots
library(ggpubr)      # Publication-ready plots and tools to arrange multiple ggplots
library(gridExtra)   # Arrange multiple grid-based objects (including ggplots) on a page
library(grid)        # Low-level grid graphics system for building complex layouts
library(devtools)    # Tools for creating, testing, and installing R packages
library(ggrepel)     # Avoid overlapping text labels in ggplot2
library(reshape2)    # Convert data between wide and long formats
library(plotly)      # Create interactive plots or make ggplots interactive
library(readxl)      # Import Excel (.xls and .xlsx) files
library(rstatix)     # Easy pipe-friendly statistical tests and result formatting
library(ggtext)      # Use HTML/CSS-like text styling in ggplot2
library(ggbreak)     # Break or zoom plot axes for improved visualisation
library(cowplot)     # Publication-quality plot themes and multi-plot layouts

#############################################################
#####                   Functions                       #####
#############################################################

# Custom averaging function
custom_average <- function(values) {
  values <- na.omit(values)
  if (length(values) == 1) return(values)
  if (all(values == 0)) return(0)
  non_zero <- values[values != 0]
  mean(non_zero)
}

source("Functions/SearchAndReplace.R")

#############################################################
#####                Folder & Files                     #####
#############################################################

# where the generated figures are saved, create folder if not existing
# dir.create(Results.dir, recursive = TRUE)# will create folder if not already there.
# Results.dir <- "Results/"

#############################################################
#####                       Codes                       #####
#############################################################
# This codes can be run subsequently
source("Code/1 - Data upload.R")
source("Code/2 - Data Analysis samples.R")
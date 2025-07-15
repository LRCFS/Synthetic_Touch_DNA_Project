###########################################################################

# ARTICLE TITLE

# AUTHORS NAMES
# INSTITUTIONS

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

library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(grid)
library(devtools)
library(ggrepel)
library(gridExtra)
library(reshape2)
library(plotly)
library(readxl) #Library to read .xls files
library(rstatix)
library(ggtext) #library to allow inline formatting using HTML-like syntax 

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
#!/usr/bin/env Rscript
#
# Copyright 2011 Alex Zvoleff
#
# This file is part of the AccraABM agent-based model.
# 
# AccraABM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# AccraABM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# AccraABM.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact Alex Zvoleff in the Department of Geography at San Diego State 
# University with any comments or questions. See the README.txt file for 
# contact information.

###############################################################################
# Plots the LULC data from a model run.
###############################################################################

require(animation)
require(raster)

PLOT_WIDTH = 8.33
PLOT_HEIGHT = 5.53
DPI = 300
PLOT_WIDTH_PIXELS = round(PLOT_WIDTH * DPI)
PLOT_HEIGHT_PIXELS = round(PLOT_HEIGHT * DPI)

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]
DATA_PATH <- "R:/Data/AccraABM/Runs/Default/20120210-175757_azvoleff-THINK"

timesteps <- read.csv(paste(DATA_PATH, "/time.csv", sep=""))
    
plot_LULC <- function(DATA_PATH, timestep) {
    #world_mask <- raster(paste(DATA_PATH, "/AccraABM_world_mask.tif", sep=""))
    land_cover <- raster(paste(DATA_PATH, "/lulc_time_", timestep, ".tif", sep=""))
    land_cover_matrix <- as.matrix(land_cover)
    total_area <- sum(!is.na(land_cover_matrix) & (land_cover_matrix != 0))
    veg_pct <- format(round((sum(land_cover_matrix==2) / total_area) * 100), width=2)
    nonveg_pct <- format(round((sum(land_cover_matrix==1) / total_area) * 100), width=2)
    image(land_cover, col=c("white", "chocolate4", "chartreuse3"), 
          main=paste("Timestep:", format(timestep, width=2)), axes=FALSE, 
          xlab="", ylab="", cex.main=1.8, cex.sub=1.8, sub=paste("Vegetation: ", veg_pct, "%, Non-vegetation: ", nonveg_pct, "%", sep=""))
}

ani.options(ffmpeg=shQuote('C:/Program Files (x86)/ffmpeg-git-6833fe4-win32-static/bin/ffmpeg.exe'))
animation_file <- shQuote(paste(DATA_PATH, "/lulc_video.wmv", sep=""))
saveVideo({for (timestep in timesteps$timestep) plot_LULC(DATA_PATH, timestep)}, 
    interval=0.35, video.name=animation_file, other.opts="-s 800x800 -b 600K")
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
# Plots the health data from a model run.
###############################################################################

require(automap)
require(animation)
require(rgdal)
#require(raster)
#require(gstat)

PLOT_WIDTH = 8.33
PLOT_HEIGHT = 5.53
DPI = 300
PLOT_WIDTH_PIXELS = round(PLOT_WIDTH * DPI)
PLOT_HEIGHT_PIXELS = round(PLOT_HEIGHT * DPI)

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]

timesteps <- read.csv(paste(DATA_PATH, "/time.csv", sep=""))
    
new_mar = par("mar")
new_mar[1] <- new_mar[1] + .4
plot_health<- function(DATA_PATH, timestep) {
    # Load the grid on which to Krige. This GeoTIFF also will be used to mask 
    # the final kriging results.
    kriglocations <- readGDAL(paste(DATA_PATH, "AccraABM_world_mask.tif", sep="/"))
    kriglocations$band1[kriglocations$band1==min(kriglocations$band1)] <- 0
    kriglocations$band1[kriglocations$band1==max(kriglocations$band1)] <- 1
    persons <- read.csv(paste(DATA_PATH, "/psns_time_", timestep, ".csv", sep=""))
    persons.spatial <- SpatialPointsDataFrame(cbind(persons$x_utm30, persons$y_utm30), persons,
            coords.nrs=c(3,4), proj4string=CRS(proj4string(kriglocations)))
    krige_results <- autoKrige(health~1, persons.spatial, kriglocations)
    krigged_imagefile <- paste(DATA_PATH, "/health_ordinary_krig_END.tif", sep="")
    par("mar"=new_mar)
    image(krige_results$krige_output, main=paste("Timestep:", format(timestep, width=2)), 
          axes=FALSE, xlab="", ylab="", cex.main=5, cex.sub=4)
}
 
ani.options(convert=shQuote('C:/Program Files/ImageMagick-6.7.1-Q16/convert.exe'))
ani.options(outdir=DATA_PATH, ani.width=1200, ani.height=1200)
animation_file <- "health_animation.gif"
saveGIF({for (timestep in timesteps$timestep) plot_health(DATA_PATH, timestep)}, 
    interval=0.35, movie.name=animation_file)

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

require(animation)
require(raster)
require(gstat)

PLOT_WIDTH = 8.33
PLOT_HEIGHT = 5.53
DPI = 300
PLOT_WIDTH_PIXELS = round(PLOT_WIDTH * DPI)
PLOT_HEIGHT_PIXELS = round(PLOT_HEIGHT * DPI)

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]
DATA_PATH <- "R:/Data/AccraABM/Runs/Default/20120211-114145_azvoleff-THINK"

timesteps <- read.csv(paste(DATA_PATH, "/time.csv", sep=""))
    
plot_health<- function(DATA_PATH, timestep) {
    # Load the grid on which to Krige. This GeoTIFF also will be used to mask 
    # the final kriging results.
    kriglocations <- readGDAL(paste(DATA_PATH, "AccraABM_world_mask.tif", sep="/"))
    kriglocations$band1[kriglocations$band1==min(kriglocations$band1)] <- 0
    kriglocations$band1[kriglocations$band1==max(kriglocations$band1)] <- 1

    persons <- read.csv(paste(DATA_PATH, "/psns_time_", timestep, ".csv", sep=""))
    persons.spatial <- SpatialPointsDataFrame(cbind(persons$x_utm30, persons$y_utm30), persons,
            coords.nrs=c(3,4), proj4string=CRS(proj4string(kriglocations)))

    # Use ordinary kriging
    v <- variogram(health~1, persons.spatial)
    #v.fit <- fit.variogram(v, vgm(1, "Exp", 6000, .05))
    v.fit <- fit.variogram(v, vgm(1, "Sph", 6000, .05))
    krigged.ord <- krige(health~1, persons.spatial, kriglocations, v.fit)
    krigged.ord.pred <- krigged.ord["var1.pred"]
    # Mask out areas outside Accra using the study area mask. Set areas outside 
    # study area to -999
    krigged.ord.pred$var1.pred <- krigged.ord.pred$var1.pred * kriglocations$band1
    krigged.ord.pred$var1.pred[krigged.ord.pred$var1.pred==0] <- -1
    proj4string(krigged.ord.pred) <- CRS(proj4string(kriglocations))
    krigged_imagefile <- paste(DATA_PATH, "/health_ordinary_krig_END.tif", sep="")
    image(krigged.ord.pred, main=paste("Timestep:", format(timestep, width=2)), 
          axes=FALSE, xlab="", ylab="", cex.main=1.8, cex.sub=1.8, 
          sub=paste("Vegetation: ", veg_pct, "%, Non-vegetation: ", nonveg_pct, 
                    "%", sep=""))
}

ani.options(ffmpeg=shQuote('C:/Program Files (x86)/ffmpeg-git-6833fe4-win32-static/bin/ffmpeg.exe'))
animation_file <- shQuote(paste(DATA_PATH, "/health_video.wmv", sep=""))
saveVideo({for (timestep in timesteps$timestep) plot_health(DATA_PATH, timestep)}, 
    interval=0.35, video.name=animation_file, other.opts="-s 800x800 -b 600K")

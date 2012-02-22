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
require(ggplot2)
require(raster)
#require(gstat)

DPI <- 300
WIDTH <- 9
HEIGHT <- 5.67
theme_update(theme_grey(base_size=12))
update_geom_defaults("point", aes(size=3))
update_geom_defaults("line", aes(size=1))

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]

timesteps <- read.csv(paste(DATA_PATH, "/time.csv", sep=""))
time_Robj <- as.Date(paste(timesteps$time_date, "15", sep=","), 
                     format="%m/%Y,%d")
timesteps <- cbind(timesteps, time_Robj)
    
#new_mar = par("mar")
#new_mar[1] <- new_mar[1] + .4
plot_kriged_heath <- function(data_path, timestep) {
    # Load the grid on which to Krige. This GeoTIFF also will be used to mask 
    # the final kriging results.
    world_mask <- readGDAL(paste(data_path, "AccraABM_world_mask.tif", sep="/"))
    world_mask$band1[world_mask$band1==min(world_mask$band1)] <- 0
    world_mask$band1[world_mask$band1==max(world_mask$band1)] <- 1
    persons <- read.csv(paste(data_path, "/psns_time_", timestep, ".csv", sep=""))
    persons <- SpatialPointsDataFrame(cbind(persons$x_utm30, persons$y_utm30), persons,
            coords.nrs=c(3,4), proj4string=CRS(proj4string(world_mask)))
    mean_health <- format(round(mean(persons$health), 2), width=4)
    krige_results <- autoKrige(health~1, persons, world_mask)
    krige_results <- krige_results$krige_output
    krige_results$var1.pred <- krige_results$var1.pred * world_mask$band1
    krige_results$var1.pred[krige_results$var1.pred == 0] <- NA
    par("mar"=new_mar)
    # Convert to a raster here because image seems to plot rasters better (the 
    # titles and subtitles weren't showing up when plotting a 
    # SpatialGridDataFrame)
    image(raster(krige_results), main=paste("Timestep:", format(timestep, width=2)), 
          axes=FALSE, xlab="", ylab="", cex.main=5, cex.sub=4,
          sub=paste("Mean self-reported health:", mean_health))
}

calculate_mean_heath <- function(data_path, timesteps) {
    # Load the grid on which to Krige. This GeoTIFF also will be used to mask 
    # the final kriging results.
    mean_healths <- c()
    for (timestep in timesteps) {
        persons <- read.csv(paste(data_path, "/psns_time_", timestep, ".csv", sep=""))
        mean_healths <- c(mean_healths, mean(persons$health))
    }
    return(mean_healths)
}
#ani.options(convert=shQuote('C:/Program Files/ImageMagick-6.7.1-Q16/convert.exe'))
#ani.options(outdir=DATA_PATH, ani.width=1200, ani.height=1200)
#animation_file <- "health_animation.gif"
#saveGIF({for (timestep in timesteps$timestep) plot_kriged_health(DATA_PATH, timestep)}, 
#    interval=0.35, movie.name=animation_file)

mean_healths <- calculate_mean_heath(DATA_PATH, timesteps$timestep) 
mean_healths <- data.frame(health=mean_healths, timevalues=timesteps$time_Robj)
qplot(time_Robj, health, data=mean_healths, xlab="Time", ylab="Mean Self-reported Health")
ggsave(filename=paste(DATA_PATH, "/mean_health_plot", ".png", sep=""), width=WIDTH, height=HEIGHT, dpi=DPI)

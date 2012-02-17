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
require(ggplot2)

DPI <- 300
WIDTH <- 9
HEIGHT <- 5.67
theme_update(theme_grey(base_size=84))
update_geom_defaults("point", aes(size=3))

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]

timesteps <- read.csv(paste(DATA_PATH, "/time.csv", sep=""))
    
new_mar = par("mar")
new_mar[1] <- new_mar[1] + .4
plot_LULC <- function(timestep, data_path) {
    #world_mask <- raster(paste(data_path, "/AccraABM_world_mask.tif", sep=""))
    land_cover <- raster(paste(data_path, "/lulc_time_", timestep, ".tif", sep=""))
    land_cover_matrix <- as.matrix(land_cover)
    total_area <- sum(!is.na(land_cover_matrix) & (land_cover_matrix != 0))
    veg_pct <- format(round((sum(land_cover_matrix==2) / total_area) * 100), width=2)
    nonveg_pct <- format(round((sum(land_cover_matrix==1) / total_area) * 100), width=2)
    # Set margins (add a little to the default bottom margin)
    par("mar"=new_mar)
    image(land_cover, col=c("white", "chocolate4", "chartreuse3"), 
          main=paste("Timestep:", format(timestep, width=2)), axes=FALSE, 
          xlab="", ylab="", cex.main=5, cex.sub=4,
          sub=paste("Vegetation: ", veg_pct, "%, Non-vegetation: ", nonveg_pct, 
                    "%", sep=""))
}

plot_NBH_veg_fraction <- function(timestep, data_path) {
    persons <- read.csv(paste(data_path, "/psns_time_", timestep, ".csv", sep=""))
    # Convert fraction to a percentage
    persons <- cbind(persons, veg_percent=persons$veg_fraction * 100)
    p <- qplot(veg_percent, geom="histogram", binwidth=2, data=persons, 
          xlab="Percentage of NBH Covered in Vegetation",
          ylab="Number of Neighborhoods", xlim=c(0,25))
    p <- p + opts(plot.margin=unit(rep(2, 4),"lines"))
    return(p)
}
 
ani.options(convert=shQuote('C:/Program Files/ImageMagick-6.7.1-Q16/convert.exe'))
ani.options(outdir=DATA_PATH, ani.width=1200, ani.height=1200)

animation_file <- "lulc_animation.gif"
saveGIF({for (timestep in timesteps$timestep) plot_LULC(timestep, DATA_PATH)}, 
    interval=0.35, movie.name=animation_file)

animation_file <- "veg_fraction_animation.gif"
plot_list <- lapply(timesteps$timestep, plot_NBH_veg_fraction, data_path=DATA_PATH)
ani.options(ani.width=DPI*WIDTH, ani.height=DPI*HEIGHT)
saveGIF(lapply(plot_list, print), interval=0.35, movie.name=animation_file)

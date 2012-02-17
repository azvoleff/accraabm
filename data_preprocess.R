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
# This file preprocesses the CVFS data in R and cleans it so that it can be 
# used in initialize.py.
###############################################################################

require(rgdal)
require(rgeos) # needed for gBuffer function
require(raster)
require(ggplot2)

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]
IMAGERY_PATH <- commandArgs(trailingOnly=TRUE)[2]
WHSA1_FILE <- commandArgs(trailingOnly=TRUE)[3]
WHSA2_FILE <- "D:/Shared_Documents/SDSU/Ghana/AccraABM/whsa2_spdf.Rdata"
ACCRA_EA_PATH <- commandArgs(trailingOnly=TRUE)[4]
buffer_distance <- commandArgs(trailingOnly=TRUE)[5]

###############################################################################
# Load the data
###############################################################################

# Buffer distance specifies the area to include in the buffered neighborhood 
# (to expand the neighborhood to avoid boundary effects with neighborhoods are 
# calculated)

theme_update(theme_grey(base_size=18))
update_geom_defaults("smooth", aes(size=1))
update_geom_defaults("line", aes(size=1))
DPI <- 300
WIDTH <- 9
HEIGHT <- 5.67

# First load the imagery
imagery <- brick(IMAGERY_PATH)
layer_names <- c("NDVI_2001", "NDVI_2010", "VIS")
layerNames(imagery) <- layer_names

# Now load the human survey data
load(WHSA2_FILE)
whsa2 <- whsa2_spdf

EAs <- readOGR(ACCRA_EA_PATH, "Accra_EAs_Updated")
# These are the EAs IDs for clusters 1, 3, and 9, in order
#EA_clusters <- list(c(605017, 605029, 605030, 605014, 605006, 605039),
#                    c(506001, 505048),
#                    c(502023, 502012, 502009))
EA_clusters <- list(c(505038, 505001, 505013, 505029, 505030),
                    c(502023, 502012, 502009, 502011),
                    c(605003, 605016, 605005, 605006))
###############################################################################
# Clean the data
###############################################################################

# Define a function to replace NAs with resampling:
replace_nas <- function(input_vector) {
    na_loc <- is.na(input_vector)
    input_vector[na_loc] <- sample(input_vector[!na_loc], sum(na_loc), replace=TRUE)
    return(input_vector)
}

data_columns <- grep("^(id|x_utm30|y_utm30|w203_own_health|ses|hweight08|ea|major_ethnic|w116_religion|education|hhid)$", names(whsa2))
whsa2 <- whsa2[, data_columns]

whsa2$w203_own_health <- replace_nas(whsa2$w203_own_health)

###############################################################################
# Now output the clusters
###############################################################################

# Cut out each cluster from the raster and calculate transition matrices
for (clustnum in 1:length(EA_clusters)) {
    # First write out the respondent data
    EA <- EAs[EAs$EA %in% EA_clusters[[clustnum]],]
    EA <- gUnaryUnion(EAs[EAs$EA %in% EA_clusters[[clustnum]],], id=clustnum)
    sample_pop_rows <- gIntersects(whsa2, EA, byid=TRUE)
    sample_pop <- whsa2[as.vector(sample_pop_rows),]
    write.csv(sample_pop, file=paste(DATA_PATH, "/cluster_", clustnum, 
                                     "_sample.csv", sep=""), 
              row.names=FALSE)
    #writeOGR(sample_pop, DATA_PATH, paste("cluster_", clustnum, 
    #                                      "_sample_shapefile", sep=""),
    #        "ESRI Shapefile", overwrite_layer=TRUE, verbose=TRUE)
    
    clipped_imagery <- crop(imagery, EA)
    #plot(clipped_imagery)
    cluster_mask <- rasterize(EA, clipped_imagery)
    # Save this mask to use it later in the ABM as a "study area mask"
    writeRaster(cluster_mask, filename=paste(DATA_PATH, "/cluster_", clustnum, 
                                             "_area_mask.tif", sep=""), 
                format="Gtiff", overwrite=TRUE)
    clipped_imagery <- cluster_mask * clipped_imagery
    layerNames(clipped_imagery) <- layer_names

    # Also prepare a buffered version of the neighborhood
    cluster_poly <- gBuffer(EA, width=buffer_distance)
    clipped_imagery_buffered <- crop(imagery, cluster_poly)
    cluster_mask_buffered <- rasterize(cluster_poly, clipped_imagery_buffered)
    # Save this mask to use it later in the ABM as a "study area mask"
    writeRaster(cluster_mask_buffered, filename=paste(DATA_PATH, "/cluster_", 
                                                      clustnum, 
                                                      "_area_mask_buffered.tif", 
                                                      sep=""),
                format="Gtiff", overwrite=TRUE)
    clipped_imagery_buffered <- cluster_mask_buffered * clipped_imagery_buffered
    layerNames(clipped_imagery_buffered ) <- layer_names
    
    for (layernum in 1:nlayers(clipped_imagery)) {
        # First write out the raster clipped to this cluster
        writeRaster(subset(clipped_imagery, layernum), filename=paste(DATA_PATH, "/cluster_", clustnum, "_", layerNames(clipped_imagery)[layernum],".tif", sep=""), format="Gtiff", overwrite=TRUE)
        # Now write out the raster clipped to the buffered areas surrounding 
        # this cluster
        writeRaster(subset(clipped_imagery_buffered, layernum), filename=paste(DATA_PATH, "/cluster_", clustnum, "_", layerNames(clipped_imagery_buffered)[layernum], "_buffered_", buffer_distance,"m.tif", sep=""), format="Gtiff", overwrite=TRUE)
    }

    # NDVI classes:
    #   0 = NA
    #   1 = NONVEG
    #   2 = VEG
    classes <- c("NONVEG", "VEG")
    # For NDVI layers, calculate and write a matrix with the transition 
    # probabilities. NOTE: Need to convert to PER YEAR transition probabilities
    t1 <- subset(clipped_imagery, grep("NDVI_2001", layerNames(clipped_imagery)))
    t2 <- subset(clipped_imagery, grep("NDVI_2010", layerNames(clipped_imagery)))
    # Eliminate the unknown values (clouds), coded as 0
    t1[t1==0] <- NA
    t2[t2==0] <- NA
    transition_matrix <- crosstab(t1, t2, long=TRUE)
    for (value in unique(transition_matrix$first)) {
        value_rows <- transition_matrix$first==value
        value_total <- sum(transition_matrix$Freq[value_rows])
        transition_matrix$Freq[value_rows] <- transition_matrix$Freq[value_rows] / value_total
    }
    write.csv(crosstab(t1, t2), file=paste(DATA_PATH, "/cluster_", clustnum, "_transition_matrix_for_ppt.csv", sep=""))
    # Need to sort the transition matrix so the upper and lower bounds can be 
    # set for the random assignment (based on these probabilities) In the 
    new_col <- rep(NA, nrow(transition_matrix))
    # Round the frequencies now to avoid rounding errors later
    transition_matrix$Freq <- round(transition_matrix$Freq, 5)
    transition_matrix <- cbind(transition_matrix, lower=new_col, upper=new_col)
    transition_matrix <- transition_matrix[order(transition_matrix$first, transition_matrix$Freq),]
    for (value in unique(transition_matrix$first)) {
        value_rows <- which(transition_matrix$first == value)
        first_row <- value_rows[1]
        transition_matrix[first_row,]$lower <- 0
        transition_matrix[first_row,]$upper <- transition_matrix[first_row,]$Freq
        for (n in 2:length(value_rows)) {
            this_row <- value_rows[n]
            prev_row <- value_rows[n - 1]
            transition_matrix[this_row,]$lower <- transition_matrix[prev_row,]$upper
            transition_matrix[this_row,]$upper <- transition_matrix[this_row,]$lower + transition_matrix[this_row,]$Freq
        }
        transition_matrix[value_rows[length(value_rows)],]$upper <- 1
    }
    # Now drop the unneeded "Freq" column.
    transition_matrix <- transition_matrix[!names(transition_matrix)=="Freq"]
    write.csv(transition_matrix, file=paste(DATA_PATH, "/cluster_", clustnum, "_transition_matrix.csv", sep=""), row.names=FALSE)

    # Make plot of changes in composition (to use in powerpoints)
    t1 <- getValues(t1)
    t2 <- getValues(t2)
    year <- rep(c(2001, 2010), 2)
    comp_classes <- rep(classes, each=2)
    percent <- c(sum(t1==1, na.rm=T) / sum(!is.na(t1), na.rm=T),
                 sum(t2==1, na.rm=T) / sum(!is.na(t2), na.rm=T),
                 sum(t1==2, na.rm=T) / sum(!is.na(t1), na.rm=T),
                 sum(t2==2, na.rm=T) / sum(!is.na(t2), na.rm=T))
    changes <- data.frame(year, percent, class=comp_classes)
    qplot(year, percent, geom="line", colour=class, data=changes)
    ggsave(paste(DATA_PATH, "/cluster_", clustnum, "_composition_plot.png", sep=""), dpi=DPI, width=WIDTH, height=HEIGHT)
    write.csv(changes, file=paste(DATA_PATH, "/cluster_", clustnum, "_composition_data.csv", sep=""), row.names=FALSE)

    png(file=paste(DATA_PATH, "/cluster_", clustnum, "_map.png", sep=""), width=6.5, height=6.5, units="in",  res=300)
    brks <- c(0, 1, 2, 3, 4)
    nbrks <- length(brks) - 1
    plot(clipped_imagery, col=rev(rainbow(nbrks)), lab.breaks=brks)
    dev.off()
}

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
# This file preprocesses the Accra data in R and cleans it so that it can be 
# used in initialize.py.
###############################################################################

library(rgdal)
library(rgeos) # needed for gBuffer function
library(raster)
library(ggplot2)
library(expm)
library(sp)

DATA_PATH <- commandArgs(trailingOnly=TRUE)[1]
IMAGERY_PATH <- commandArgs(trailingOnly=TRUE)[2]
WHSA_FILE <- commandArgs(trailingOnly=TRUE)[3]
FMV_PATH <- commandArgs(trailingOnly=TRUE)[4]
buffer_distance <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
WINDOWED_MARKOV <- as.logical(commandArgs(trailingOnly=TRUE)[6])
MARKOV_WINDOW_SIZE <- as.numeric(commandArgs(trailingOnly=TRUE)[7])

DATA_PATH <- "R:/Data/Ghana/AccraABM/Initialization"
IMAGERY_PATH <- "R:/Data/Ghana/Imagery/Layer_Stack/NDVI2002_NDVI2010_VIS.tif"
WHSA_FILE <- "D:/Shared_Documents/SDSU/Ghana/AccraABM/whsa1_whsa2linked.Rdata"
FMV_PATH <- "R:/Data/Ghana/GIS/Accra_DB_Export"
buffer_distance <- 700
WINDOWED_MARKOV <- TRUE
MARKOV_WINDOW_SIZE <- 5

# The number of years included in the calibration dataset, so that the 
# transition matrix can be adjusted to a period of one year.
MARKOV_CALIBRATION_INTERVAL <- 8 

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
load(WHSA_FILE)
#whsa <- whsa2_spdf_FMV_ego
whsa <- whsa1_whsa2linked

FMVs <- readOGR(FMV_PATH, "FMVNeighborhoods")
# The FMVs are:
#   31: Dansoman Estate
#   67: Nima
#   70: North Kaneshie
Chosen_FMV_IDs <- c(31, 67, 70)
###############################################################################
# Clean the data
###############################################################################

# Define a function to replace NAs with resampling:
replace_nas <- function(input_vector) {
    na_loc <- is.na(input_vector)
    input_vector[na_loc] <- sample(input_vector[!na_loc], sum(na_loc), replace=TRUE)
    return(input_vector)
}

data_columns <- grep("^(id|x_utm30|y_utm30|rand_sf36_pf|ea|major_ethnic|age|education|Charcoal|OwnToilet|YrsInHouse_Cat|hhid)$", names(whsa))
whsa <- whsa[, data_columns]

whsa$education <- replace_nas(whsa$education)
whsa$age <- replace_nas(whsa$age)
whsa$major_ethnic <- replace_nas(whsa$major_ethnic)
whsa$OwnToilet <- replace_nas(whsa$OwnToilet)
whsa$Charcoal <- replace_nas(whsa$Charcoal)
whsa$YrsInHouse_Cat <- replace_nas(whsa$YrsInHouse_Cat)
whsa$rand_sf36_pf <- replace_nas(whsa$rand_sf36_pf)

###############################################################################
# Now output the clusters
###############################################################################

# Cut out each cluster from the raster and calculate Markov transition matrices
for (FMV_num in 1:length(Chosen_FMV_IDs)) {
    # First write out the respondent data
    FMV <- FMVs[FMVs$FMV_Number == Chosen_FMV_IDs[FMV_num],]
    FMV_name <- abbreviate(FMV$FMV_Name, minlength=6)
    sample_pop_rows <- gIntersects(whsa, FMV, byid=TRUE)
    #FMV_spdf <- SpatialPolygonsDataFrame(FMV, data=data.frame(ID=row.names(FMV)))
    #writeOGR(FMV_spdf, DATA_PATH, paste("cluster_", FMV_name, 
    #                                      "_FMV_shapefile", sep=""),
    #        "ESRI Shapefile", overwrite_layer=TRUE, verbose=TRUE)
    sample_pop <- whsa[as.vector(sample_pop_rows),]
    write.csv(sample_pop, file=paste(DATA_PATH, "/cluster_", FMV_name, 
                                     "_sample.csv", sep=""), row.names=FALSE)
    # Also write out summary stats:
    sample_pop_df <- as.data.frame(sample_pop)
    var_columns <- grep("^(rand_sf36_pf|major_ethnic|age|education|YrsInHouse_Cat|Charcoal|OwnToilet)$", names(sample_pop_df))
    mins <- apply(sample_pop_df[var_columns], 2, function(x) min(as.numeric(x), na.rm=TRUE))
    maxs <- apply(sample_pop_df[var_columns], 2, function(x) max(as.numeric(x), na.rm=TRUE))
    means <- apply(sample_pop_df[var_columns], 2, function(x) mean(as.numeric(x), na.rm=TRUE))
    sds <- apply(sample_pop_df[var_columns], 2, function(x) sd(as.numeric(x), na.rm=TRUE))
    samp_size <- apply(sample_pop_df[var_columns], 2, function(x) sum(!is.na(x)))
    summary_stats <- data.frame(mins, maxs, means, sds, samp_size, 
                                row.names=names(sample_pop_df[var_columns]))
    write.csv(summary_stats, file=paste(DATA_PATH, "/cluster_", FMV_name, 
                                        "_summary_stats.csv", sep=""))
    #sample_pop_ogr <- sample_pop
    #sample_pop_ogr$major_ethnic <- factor(sample_pop_ogr$major_ethnic, ordered=FALSE)
    #sample_pop_ogr$education <- factor(sample_pop_ogr$education, ordered=FALSE)
    #sample_pop_ogr$srh <- factor(sample_pop_ogr$srh, ordered=FALSE)
    #writeOGR(sample_pop_ogr, DATA_PATH, paste("cluster_", FMV_name, 
    #                                      "_sample_shapefile", sep=""),
    #        "ESRI Shapefile", overwrite_layer=TRUE, verbose=TRUE)
    
    clipped_imagery <- crop(imagery, FMV)
    #plot(clipped_imagery)
    cluster_mask <- rasterize(FMV, clipped_imagery)
    # Save this mask to use it later in the ABM as a "study area mask"
    writeRaster(cluster_mask, filename=paste(DATA_PATH, "/cluster_", FMV_name, 
                                             "_area_mask.tif", sep=""), 
                format="Gtiff", overwrite=TRUE)
    clipped_imagery <- cluster_mask * clipped_imagery
    layerNames(clipped_imagery) <- layer_names
    # Also prepare a buffered version of the neighborhood
    cluster_poly <- gBuffer(FMV, width=buffer_distance)
    clipped_imagery_buffered <- crop(imagery, cluster_poly)
    cluster_mask_buffered <- rasterize(cluster_poly, clipped_imagery_buffered)
    # Save this mask to use it later in the ABM as a "study area mask"
    writeRaster(cluster_mask_buffered, filename=paste(DATA_PATH, "/cluster_", 
                                                      FMV_name, 
                                                      "_area_mask_buffered.tif", 
                                                      sep=""),
                format="Gtiff", overwrite=TRUE)
    clipped_imagery_buffered <- cluster_mask_buffered * clipped_imagery_buffered
    layerNames(clipped_imagery_buffered ) <- layer_names
    for (layernum in 1:nlayers(clipped_imagery)) {
        # First write out the raster clipped to this cluster
        writeRaster(subset(clipped_imagery, layernum), filename=paste(DATA_PATH, "/cluster_", FMV_name, "_", layerNames(clipped_imagery)[layernum],".tif", sep=""), format="Gtiff", overwrite=TRUE)
        # Now write out the raster clipped to the buffered areas surrounding 
        # this cluster
        writeRaster(subset(clipped_imagery_buffered, layernum), filename=paste(DATA_PATH, "/cluster_", FMV_name, "_", layerNames(clipped_imagery_buffered)[layernum], "_buffered_", buffer_distance,"m.tif", sep=""), format="Gtiff", overwrite=TRUE)
    }

    # NDVI classes:
    #   0 = NA
    #   1 = NONVEG
    #   2 = VEG
    class_names <- c("NONVEG", "VEG")
    class_codes <- c(1, 2)
    # For NDVI layers, calculate and write a matrix with the Markov transition 
    # probabilities.
    t1 <- subset(clipped_imagery, grep("NDVI_2001", layerNames(clipped_imagery)))
    t2 <- subset(clipped_imagery, grep("NDVI_2010", layerNames(clipped_imagery)))
    # Eliminate the unknown values (clouds), coded as 0
    t1[t1==0] <- NA
    t2[t2==0] <- NA
    if (WINDOWED_MARKOV) {
        # TODO: this currently only works if there are only two classes.
        t1_veg_prob <- focal(t1, w=MARKOV_WINDOW_SIZE, function(x) sum(x==2) / 
                             sum(!is.na(x)))
        t2_veg_prob <- focal(t2, w=MARKOV_WINDOW_SIZE, function(x) sum(x==2) / 
                             sum(!is.na(x)))
        prob_nonveg2veg <- (1 - t1_veg_prob) * t2_veg_prob
        prob_nonveg2veg <- mean(as.matrix(prob_nonveg2veg), na.rm=TRUE)
        prob_veg2nonveg <- t1_veg_prob * (1 - t2_veg_prob)
        prob_veg2nonveg <- mean(as.matrix(prob_veg2nonveg), na.rm=TRUE)
        # Setup the transition matrix:
        trans_matrix_cal <- matrix(c(1 - prob_nonveg2veg, prob_nonveg2veg,
                                      prob_veg2nonveg, 1 - prob_veg2nonveg),
                                    byrow=TRUE, nrow=2, ncol=2)
    } else {
        trans_matrix_cal <- crosstab(t1, t2)
        trans_matrix_cal <- trans_matrix_cal / rowSums(trans_matrix_cal)
    }

    # Now convert to per-year transition probabilities:
    trans_matrix <- expm((1/MARKOV_CALIBRATION_INTERVAL) * logm(trans_matrix_cal))
    rownames(trans_matrix) <- colnames(trans_matrix) <- class_codes
    # Also save a version of the matrix with textual class names, for ppt, etc.
    trans_matrix_ppt <- trans_matrix
    rownames(trans_matrix_ppt) <- colnames(trans_matrix_ppt) <- class_names
    write.csv(trans_matrix_ppt, file=paste(DATA_PATH, "/cluster_", FMV_name, 
                                            "_trans_matrix_for_ppt.csv", 
                                            sep=""))
    trans_matrix_long <- data.frame(first=rep(colnames(trans_matrix), 2))
    trans_matrix_long$second <- rep(colnames(trans_matrix), each=2)
    trans_matrix_long$Freq <- rep(NA, nrow(trans_matrix_long))
    for (n in 1:nrow(trans_matrix_long)) {
        class1 <- trans_matrix_long$first[n]
        class2 <- trans_matrix_long$second[n]
        trans_matrix_long$Freq[n] <- trans_matrix[class1, class2]
    }
    for (value in unique(trans_matrix_long$first)) {
        value_rows <- trans_matrix_long$first==value
        value_total <- sum(trans_matrix_long$Freq[value_rows])
        trans_matrix_long$Freq[value_rows] <- 
            trans_matrix_long$Freq[value_rows] / value_total
    }
    # Need to sort the transition matrix so the upper and lower bounds can be 
    # set for the random assignment (based on these probabilities) In the 
    new_col <- rep(NA, nrow(trans_matrix_long))
    # Round the frequencies now to avoid rounding errors later
    trans_matrix_long$Freq <- round(trans_matrix_long$Freq, 5)
    trans_matrix_long <- cbind(trans_matrix_long, lower=new_col, upper=new_col)
    trans_matrix_long <- trans_matrix_long[order(trans_matrix_long$first, trans_matrix_long$Freq),]
    for (value in unique(trans_matrix_long$first)) {
        value_rows <- which(trans_matrix_long$first == value)
        first_row <- value_rows[1]
        trans_matrix_long[first_row,]$lower <- 0
        trans_matrix_long[first_row,]$upper <- trans_matrix_long[first_row,]$Freq
        for (n in 2:length(value_rows)) {
            this_row <- value_rows[n]
            prev_row <- value_rows[n - 1]
            trans_matrix_long[this_row,]$lower <- trans_matrix_long[prev_row,]$upper
            trans_matrix_long[this_row,]$upper <- trans_matrix_long[this_row,]$lower + trans_matrix_long[this_row,]$Freq
        }
        trans_matrix_long[value_rows[length(value_rows)],]$upper <- 1
    }
    # Now drop the unneeded "Freq" column.
    trans_matrix_long <- trans_matrix_long[!names(trans_matrix_long)=="Freq"]
    write.csv(trans_matrix_long, file=paste(DATA_PATH, "/cluster_", FMV_name, "_trans_matrix.csv", sep=""), row.names=FALSE)

    # Make plot of changes in composition (to use in powerpoints)
    t1 <- getValues(t1)
    t2 <- getValues(t2)
    year <- rep(c(2002, 2010), 2)
    comp_classes <- rep(class_names, each=2)
    percent <- c(sum(t1==1, na.rm=T) / sum(!is.na(t1), na.rm=T),
                 sum(t2==1, na.rm=T) / sum(!is.na(t2), na.rm=T),
                 sum(t1==2, na.rm=T) / sum(!is.na(t1), na.rm=T),
                 sum(t2==2, na.rm=T) / sum(!is.na(t2), na.rm=T))
    changes <- data.frame(year, percent, class=comp_classes)
    qplot(year, percent, geom="line", colour=class, data=changes)
    ggsave(paste(DATA_PATH, "/cluster_", FMV_name, "_composition_plot.png", sep=""), dpi=DPI, width=WIDTH, height=HEIGHT)
    write.csv(changes, file=paste(DATA_PATH, "/cluster_", FMV_name, "_composition_data.csv", sep=""), row.names=FALSE)
    png(file=paste(DATA_PATH, "/cluster_", FMV_name, "_map.png", sep=""), width=6.5, height=6.5, units="in",  res=300)
    brks <- c(0, 1, 2, 3, 4)
    nbrks <- length(brks) - 1
    plot(clipped_imagery, col=rev(rainbow(nbrks)), lab.breaks=brks)
    dev.off()
}

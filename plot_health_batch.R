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

require(ggplot2)
require(reshape) # For 'melt' command

DPI <- 300
WIDTH <- 9
HEIGHT <- 5.67
theme_update(theme_bw(base_size=18))
update_geom_defaults("point", aes(size=2))
update_geom_defaults("line", aes(size=.75))

DATA_PATH<- commandArgs(trailingOnly=TRUE)[1]
# Give the root of the FMV name without the "_vernacular" or "_egocentric" 
# suffix
FMV_NAME <- commandArgs(trailingOnly=TRUE)[2]

DATA_PATH <- "R:/Data/Ghana/AccraABM/Runs/"
#FMV_NAME <- "Nima"
#FMV_NAME <- "DnsmnE"
FMV_NAME <- "NrthKn"

calculate_mean_heath <- function(data_path) {
    mean_healths <- c()
    timesteps <- read.csv(paste(data_path, "/time.csv", sep=""))
    for (timestep in timesteps$timestep) {
        persons <- read.csv(paste(data_path, "/psns_time_", timestep, ".csv", sep=""))
        mean_healths <- c(mean_healths, mean(persons$health))
    }
    time_Robj <- as.Date(paste(timesteps$time_date, "15", sep=","), 
                         format="%m/%Y,%d")
    mean_healths <- data.frame(time_Robj, mean_healths)
    return(mean_healths)
}

scenario_names <- c("_none", "_egocentric", "_vernacular")
scenario_paths <- paste(DATA_PATH, FMV_NAME, scenario_names, sep="")
scenario_num <- 1
for (scenario_path in scenario_paths) {
    directories <- list.files(scenario_path)
    # Only match the model results folders - don't match any other folders or files 
    # in the directory, as trying to read results from these other files/folders 
    # would lead to an error.
    directories <- directories[grep("[0-9]{8}-[0-9]{6}", directories)]
    if (length(directories)<1) stop(paste("can't run plot_LULC_batch with", length(directories), "model runs."))
    if (length(directories)<5) warning(paste("Only", length(directories), "model runs found."))
    n <- 1
    for (directory in directories) {
        full_directory_path <- paste(scenario_path, directory, sep="/") 
        mean_healths <- calculate_mean_heath(full_directory_path)
        runname <- paste("run", n, sep="")
        names(mean_healths)[2] <- paste("mean_health_", runname, sep="")
        if (n==1) {
            mean_healths_all_runs <- mean_healths
        } else  {
            mean_healths_all_runs <- merge(mean_healths_all_runs, mean_healths, by="time_Robj")
        }
        n <- n + 1
    }
    mean_healths_cols <- grep('mean_health', names(mean_healths_all_runs))
    mean_healths_stats <- data.frame(time_Robj=mean_healths_all_runs$time_Robj)
    mean_healths_stats$meanhealth_FMV <- apply(mean_healths_all_runs[mean_healths_cols], 1, mean, na.rm=T)
    mean_healths_stats$sd_FMV <- apply(mean_healths_all_runs[mean_healths_cols], 1, sd, na.rm=T)
    mean_healths_stats$conf_lower_FMV <- mean_healths_stats$meanhealth - mean_healths_stats$sd
    mean_healths_stats$conf_upper_FMV <- mean_healths_stats$meanhealth + mean_healths_stats$sd
    names(mean_healths_stats) <- gsub('FMV', scenario_num, names(mean_healths_stats))
    if (scenario_num==1) {
        scenario_stats <- mean_healths_stats
    } else {
        scenario_stats <- merge(scenario_stats, mean_healths_stats, by="time_Robj")
    }
    scenario_num <- scenario_num + 1
}

#qplot(time_Robj, health, data=mean_healths, xlab="Time", ylab="Mean Self-reported Health")
#ggsave(filename=paste(DATA_PATH, "/mean_health_plot", ".png", sep=""), width=WIDTH, height=HEIGHT, dpi=DPI)

mean_cols <- grep('^(time_Robj|meanhealth)', names(scenario_stats))
meanhealth <- scenario_stats[mean_cols]
meanhealth <- meanhealth[-1,]
meanhealth <- melt(meanhealth, id='time_Robj')

p <- qplot(time_Robj, value, geom="line", linetype=variable, data=meanhealth)
p + labs(x="Time", y="Physical Functioning Score") +
    scale_linetype_discrete(name="Model Type",
                            breaks=c("meanhealth_1", "meanhealth_2", "meanhealth_3"), 
                            labels=c("No NBH Effects", "Egocentric NBHs", "Vernacular NBHs")) +
    opts(legend.position=c(.85, .85)) +
    geom_vline(xintercept=as.numeric(meanhealth$time_Robj[6]), linetype=4)
ggsave(paste(DATA_PATH, "/", FMV_NAME, "_mean_health.png", sep=""), 
       width=WIDTH,
        height=HEIGHT, dpi=DPI)

#### This file serves to perform the tortuosity analysis with varying parameters on the three types of blood vessels
#### in order to analyze the metrics generated.

setwd("/Users/kaitlinylim/Documents/angicart/data/images/mouse_brain/tCar/dat files")

source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/dat_file_reader.R")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/single_vessel_plotter.R")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/vessel_poly_fit.R")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/vessel_spline_fit.R")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/frenet_vectors.R")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/tortuosity_metrics.R")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/curve_torse_check_plotter.R")

library("rgl")
library(ggplot2)
library(grid)

filename <- "PT 12 IH001_vd_1242x1242x2200_sm_th_0.17.dat"

## Use the dat_file_reader() function to read in a .dat file for analysis.
vessels_slice <- dat_file_reader(dat_filename = filename)

## Extracting an individual vessel to perform the parameter analysis.
vessel <- vessels_slice[which(vessels_slice$ID == 50), 1:3]

## Global Variables
cur_tor_met <- data.frame(TC=double(), AC=double(), TT=double(), AT=double(), MC=double(), MT=double(), TCC=double(), ACC=double(),
                          stringsAsFactors = FALSE)
dist_met <- c()
inflec_met <- c()
soam <- c()

for (ind in 3:12) {

  ## We will begin with the polynomial fitting.
  smth_vessel <- vessel_poly_fit(vessel = vessel, number_samples = 10000, poly_degree = ind, plot = FALSE)
  
  ## Now we save the metrics generated with smth_vessel into the respective lists.
  cur_tor_met <- rbind(cur_tor_met, curvature_torsion_calculator(tangent = smth_vessel[[1]], normal = smth_vessel[[2]], 
                                                                 binormal = smth_vessel[[3]], vessel_coords = smth_vessel[[4]])[4:11])
  
  temp <- distance_metric(vessel_coords = smth_vessel[[4]])
  names(temp) <- NULL
  dist_met <- c(dist_met, temp[2])
  
  inflec_met <- c(inflec_met, inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = TRUE, return_coords = FALSE))
  
  soam <- c(soam, sum_of_all_angles_metric(smth_vessel[[4]]))
}

## Here, we graph these results. We will test using the soam list.
## Note to self: find a way to alter the dimensions of the axes.
for (i in 1:8) {
  df_cur_tor_met <- data.frame(degrees = c(3:12),
                               vals = cur_tor_met[,i])
  plot <- ggplot(data=df_cur_tor_met, aes(x=degrees, y=vals, group=1)) +
    geom_line() + geom_point() +
    labs(title=paste0(colnames(cur_tor_met)[i], " vs. Polynomial fitting parameter"), x="Degree",y=colnames(cur_tor_met)[i]) +
    theme_light() +
    theme(plot.title = element_text(size=12, hjust=0.5)) +
    scale_x_reverse()
  ggsave(filename = paste0(filename, "_poly_", colnames(cur_tor_met)[i], ".pdf"), plot=plot, height=7, width=7)
}

df_dist_met <- data.frame(degrees = c(3:12),
                          vals = dist_met)
plot <- ggplot(data=df_dist_met, aes(x=degrees, y=vals, group=1)) +
  geom_line() + geom_point() +
  labs(title="DM vs. Polynomial fitting parameter",x="Degree",y="DM") +
  theme_light() +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  scale_x_reverse()
ggsave(filename = paste0(filename, "_poly_DM.pdf"), plot=plot, height = 7, width = 7)

df_inflec_met <- data.frame(degrees = c(3:12),
                            vals = inflec_met)
plot <- ggplot(data=df_inflec_met, aes(x=degrees, y=vals, group=1)) +
  geom_line() + geom_point() +
  labs(title="ICM vs. Polynomial fitting parameter", x="Degree", y="ICM") +
  theme_light() +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  scale_x_reverse()
ggsave(filename = paste0(filename, "_poly_ICM.pdf"), plot=plot, height = 7, width = 7)

df_soam <- data.frame(degrees = c(3:12),
                      vals = soam)
plot <- ggplot(data=df_soam, aes(x=degrees, y=vals, group=1)) +
  geom_line() + geom_point() +
  labs(title="SOAM vs. Polynomial fitting parameter",x="Degree", y = "SOAM") +
  theme_light() +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  scale_x_reverse()
ggsave(filename = paste0(filename, "_poly_SOAM.pdf"), plot=plot, height = 7, width = 7)


## Now, we will use the spline fitting.
# First we must refresh the global variables.
cur_tor_met <- data.frame(TC=double(), AC=double(), TT=double(), AT=double(), MC=double(), MT=double(), TCC=double(), ACC=double(),
                          stringsAsFactors = FALSE)
dist_met <- c()
inflec_met <- c()
soam <- c()

## We must first save the metrics when the smoothing parameter is NULL.
smth_vessel <- vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = NULL, plot = FALSE)

cur_tor_met <- rbind(cur_tor_met, curvature_torsion_calculator(tangent = smth_vessel[[1]], normal = smth_vessel[[2]], 
                                                               binormal = smth_vessel[[3]], vessel_coords = smth_vessel[[4]])[4:11])

temp <- distance_metric(vessel_coords = smth_vessel[[4]])
names(temp) <- NULL
dist_met <- c(dist_met, temp[2])

inflec_met <- c(inflec_met, inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = TRUE, return_coords = FALSE))

soam <- c(soam, sum_of_all_angles_metric(smth_vessel[[4]]))

## Now, we can loop through the smoothing parameters.
for (ind in seq(from=0, to=1, by=0.1)) {
  smth_vessel <- vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = ind, plot = FALSE)
  
  cur_tor_met <- rbind(cur_tor_met, curvature_torsion_calculator(tangent = smth_vessel[[1]], normal = smth_vessel[[2]], 
                                                                 binormal = smth_vessel[[3]], vessel_coords = smth_vessel[[4]])[4:11])
  
  temp <- distance_metric(vessel_coords = smth_vessel[[4]])
  names(temp) <- NULL
  dist_met <- c(dist_met, temp[2])
  
  inflec_met <- c(inflec_met, inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = TRUE, return_coords = FALSE))
  
  soam <- c(soam, sum_of_all_angles_metric(smth_vessel[[4]]))
}

## Here, we graph these results.
for (i in 1:8) {
  df_cur_tor_met <- data.frame(smoothing = c(-0.1, seq(from=0, to=1, by=0.1)),
                               vals = cur_tor_met[,i])
  plot <- ggplot(data=df_cur_tor_met, aes(x=smoothing, y=vals, group=1)) +
    geom_line() + geom_point() +
    labs(title=paste0(colnames(cur_tor_met)[i], " vs. Spline fitting parameter"), x="Smoothing",y=colnames(cur_tor_met)[i]) +
    theme_light() +
    theme(plot.title = element_text(size=10, hjust=0.5)) +
    scale_x_reverse()
  ggsave(filename = paste0(filename, "_spline_", colnames(cur_tor_met)[i], ".pdf"), plot=plot, height=7, width=7)
}

df_dist_met <- data.frame(smoothing = c(-0.1, seq(from=0, to=1, by=0.1)),
                          vals = dist_met)
plot <- ggplot(data=df_dist_met, aes(x=smoothing, y=vals, group=1)) +
  geom_line() + geom_point() +
  labs(title="DM vs. Spline fitting parameter",x="Smoothing",y="DM") +
  theme_light() +
  theme(plot.title = element_text(size=10, hjust=0.5)) +
  scale_x_reverse()
ggsave(filename = paste0(filename, "_spline_DM.pdf"), plot=plot, height = 7, width = 7)

df_inflec_met <- data.frame(smoothing = c(-0.1, seq(from=0, to=1, by=0.1)),
                            vals = inflec_met)
plot <- ggplot(data=df_inflec_met, aes(x=smoothing, y=vals, group=1)) +
  geom_line() + geom_point() +
  labs(title="ICM vs. Spline fitting parameter", x="Smoothing", y="ICM") +
  theme_light() +
  theme(plot.title = element_text(size=10, hjust=0.5)) +
  scale_x_reverse()
ggsave(filename = paste0(filename, "_spline_ICM.pdf"), plot=plot, height = 7, width = 7)

df_soam <- data.frame(smoothing = c(-0.1, seq(from=0, to=1, by=0.1)),
                      vals = soam)
plot <- ggplot(data=df_soam, aes(x=smoothing, y=vals, group=1)) +
  geom_line() + geom_point() +
  labs(title="SOAM vs. Spline fitting parameter",x="Smoothing", y = "SOAM") +
  theme_light() +
  theme(plot.title = element_text(size=10, hjust=0.5)) +
  scale_x_reverse()
ggsave(filename = paste0(filename, "_spline_SOAM.pdf"), plot=plot, height = 7, width = 7)

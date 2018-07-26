################# Tortuosity Metrics #################

# This script contains functions for calculating tortuosity metrics for easy sourcing.
# Metrics included are: (1) distance metric; (2) inflection count metric; (3) sum of all angles metric; (4) curvature torsion metric

library("pracma")
source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/frenet_vectors.R")

################## Distance Metric #################

## This code is the function for calculating the standard distance metric (DM) of tortuosity of a curve.
## This method calculated the path length along a curve, and divides that value by the euclidean distance between the endpoints of the curve.
## For documentation on this method, see Bullitt et al., "Measuring Tortuosity of the Intracerebral Vasculature From MRA Images", IEEE Transactions on medical Imaging, Vol. 22, No. 9, September 2003.
## Note that this metric is impervious to differentiating between curvature and tortuosity...but we all must start somewhere.


euclidean_distance <- function(pt1, pt2){
  ## This function calculates the euclidean distance between two points in 3D, labeled as pt1 and pt2.
  distance <- sqrt((pt2[1] - pt1[1])**2 + (pt2[2] - pt1[2])**2 + (pt2[3] - pt1[3])**2)
  return(distance)
}

distance_metric <- function(vessel_coords, units){
  pointwise_distance <- mat.or.vec(length(vessel_coords[,1]), 1)
  for(i in 1:(length(pointwise_distance)-1)){
    pointwise_distance[i] <- euclidean_distance(vessel_coords[i,], vessel_coords[i+1,])
  }
  pathlength <- sum(pointwise_distance)
  end_to_end_distance <- euclidean_distance(vessel_coords[1,], vessel_coords[length(pointwise_distance),])
  return(c(pathlength, pathlength/end_to_end_distance))
}

################## Inflection Count Metric #################

## This code is the function for calculating the inflection count metric (ICM) of tortuosity of a curve.
## This method calculates the distance metric (DM) and multiples it by the number of inflection points along the curve.
## Inflection points can be identified by identifying local maxima of the quantity delta N dot delta N, where N is the unit vector
## representing the Frenet normal axis, and delta N represnts changes in this vector associated with points pt1 and pt2 along the curve.
## For documentation on this method, see Bullitt et al., "Measuring Tortuosity of the Intracerebral Vasculature From MRA Images", IEEE Transactions on medical Imaging, Vol. 22, No. 9, September 2003.
## Note that this metric is intended to differentiate between long broad curves with few oscillations and curves with high oscillations.
## Realize that here oscillations is to be interpreted as in plane oscillations, not spirals (or straight helices).
## This is because spirals do not actually exhibit inflection points.

inflection_count_metric <- function(vessel_coords, return_count = TRUE, return_coords = FALSE){
  ## Initialize normal_array, delta_N array, and delta_N squared vector
  normal_array <- mat.or.vec(length(vessel_coords[,1]), 3)
  delta_N <- mat.or.vec(length(vessel_coords[,1]), 3)
  delta_N_squared <- mat.or.vec(length(vessel_coords[,1]), 1)
  
  for(i in 1:(length(normal_array[,1]) - 2)){
    normal_array[i+1,] <- normal_vector(vessel_coords[i,], vessel_coords[i+1,], vessel_coords[i+2,])
  }
  
  for(i in 1:(length(delta_N[,1]) - 1)){
    delta_N[i+1,] <- normal_array[i+1,] - normal_array[i,]
  }
  
  for(i in 1:length(delta_N[,1])){
    delta_N_squared[i] <- sum(delta_N[i,]*delta_N[i,])
  }
  
  inflection_count <- length(which(delta_N_squared > 2))
  if(return_count == TRUE){
    return(inflection_count)
  }
  if(return_coords == TRUE){
    return(delta_N_squared)
  }
}

################## Sum of All Angles #################

## This code is the function for calculating the sum of all angles metric (SOAM) of tortuosity of a curve.
## This method calculates the both the angular changes within plane (a proxy measure for local curvature, 
## represented as IP_k for in-plane angle at point P_k) and out of plane (a proxy measure for local torsion, represented as TP_k for torsion angle at point P_k). 
## These two quantities are combined via a square root of the sum of squares to generate the total angle 
## (a proxy for total "curvature", represented as CP_k for curvature at point P_k).  
## Finally, the curvature CP_k is "integreated" by discrete sum along the lenght of the curve, and the total quantity is normalized by path length.  
## For documentation on this method, see Bullitt et al., "Measuring Tortuosity of the Intracerebral Vasculature From MRA Images", IEEE Transactions on medical Imaging, Vol. 22, No. 9, September 2003.
## Note that this metric is intended to better identify tight coils void of inflection points.

## Should explore how this metric performs in identifying in plane arc versus out of plane tortuous vessel.
## That is, output value does not uniquely distinguish between two cases as the angles are combined.  

## Consider an extention of this metric that calculates integrated curvature and torsion seperately.

cpk_calc <- function(pt1, pt2, pt3, pt4){
  
  T1 <- pt2 - pt1
  T1_norm <- T1/sqrt(sum(T1*T1))
  T2 <- pt3 - pt2
  T2_norm <- T2/sqrt(sum(T2*T2))
  
  ipk_value <- acos(sum(T1_norm*T2_norm))
  
  T3 <- pt4 - pt3
  T1_cross_T2 <- cross(T1, T2)
  T1_cross_T2_norm <- T1_cross_T2/sqrt(sum(T1_cross_T2*T1_cross_T2))
  T2_cross_T3 <- cross(T2, T3)
  T2_cross_T3_norm <- T2_cross_T3/sqrt(sum(T2_cross_T3*T2_cross_T3))
  
  tpk_value <- acos(sum(T1_cross_T2_norm*T2_cross_T3_norm))
  
  cpk_value <- sqrt(ipk_value**2 + tpk_value**2)
  return(cpk_value)
}

sum_of_all_angles_metric <- function(vessel_coords){
  
  num_points <- length(vessel_coords[,1])
  
  cpk_array <- mat.or.vec(num_points,1)
  
  for(i in 2:(num_points - 2)){
    pt1 <- vessel_coords[i-1,]
    pt2 <- vessel_coords[i,]
    pt3 <- vessel_coords[i+1,]
    pt4 <- vessel_coords[i+2,]
    
    cpk_array[i] <- cpk_calc(pt1, pt2, pt3, pt4)
  }
  
  path_length <- distance_metric(vessel_coords)[1]
  
  soam <- sum(cpk_array[-c(1, num_points, num_points-1, num_points-2)])/path_length
  # soam <- sum(cpk_array[-c(1, num_points, num_points-1, num_points-2)])
  return(soam)
}

################## Curvature Torsion Calculator #################

curvature_torsion_calculator <- function(tangent, normal, binormal, vessel_coords){
  
  ## This function calculates curvature at a point(s) given a set Frenet-Serrat frame vectors.
  ## The approach used to calculate curvature and torsion is rather simple, and is taken from RHB pps. 340-343, however it seems sufficient for the time being.
  ## Currently, this is done using the defining relations between curvature and torsion and the Frenet-Serrat frame vectors.
  ## Specifically, from dt-hat/ds = kappa*n-hat, we can find curvature kappa from kappa = n-hat dot dt-hat/ds.
  ## Similarly, from db-hat/ds = -tau*n-hat, we can find torsion tau from tau = -n-hat dot db-hat/ds.
  ## This approach seems considerably subject to noise in the data, thus the next version to try would be calculting curvature and torsion directly from the 
  ## r(s) parameterization, where we discretize derivatives of r(s) to extract curvature and torsion (this is taken fomr O'Flynn et al., Annals of Biomedical Engineering, 2007).
  ## Currently, we calculate the Frenet-Serrat frame vectors from r(s), then discretize derivatives of the frame vectors.
  
  # Coerce frenet arrays to matrices.
  
  tangent <- as.matrix(tangent)
  normal <- as.matrix(normal)
  binormal <- as.matrix(binormal)
  
  # Initialize curvature and torsion values along vessel
  curve_vec <- mat.or.vec(nrow(tangent),1)
  torse_vec <- mat.or.vec(nrow(tangent),1)
  combined_curvature <- mat.or.vec(nrow(tangent),1)
  curve_check_vec <- mat.or.vec(nrow(tangent), 1)
  torse_check_vec <- mat.or.vec(nrow(tangent), 1)
  
  curve_vec[] <- NaN
  torse_vec[] <- NaN
  combined_curvature[] <- NaN
  curve_check_vec[] <- NaN
  torse_check_vec[] <- NaN
  
  # Arc length is needed for the entire length of the curve.
  
  arc_length <- mat.or.vec(nrow(tangent),1)
  arc_length[] <- NaN
  
  pointwise_distance <- mat.or.vec(length(vessel_coords[,1]), 1)
  for(i in 1:(length(pointwise_distance)-1)){
    pointwise_distance[i] <- euclidean_distance(vessel_coords[i,], vessel_coords[i+1,])
    arc_length[i] <- sum(pointwise_distance)
  }
  
  diff_tan <- mat.or.vec(nrow(tangent),3)
  for(i in 1:(nrow(tangent)-2)){
    diff_tan[i+1,] <- (tangent[i+2,] - tangent[i,])/(arc_length[i+2]-arc_length[i])
  }
  
  diff_bin <- mat.or.vec(nrow(tangent),3)
  for(i in 1:(nrow(tangent)-2)){
    diff_bin[i+1,] <- (binormal[i+2,] - binormal[i,])/(arc_length[i+2]-arc_length[i])
  }
  
  diff_nor <- mat.or.vec(nrow(tangent),3)
  for(i in 1:(nrow(tangent)-2)){
    diff_nor[i+1,] <- (normal[i+2,] - normal[i,])/(arc_length[i+2] - arc_length[i])
  }
  
  for(i in 1:nrow(tangent)){
    curve_vec[i] <- sum(diff_tan[i,]*normal[i,])
    torse_vec[i] <- -sum(diff_bin[i,]*normal[i,])
    combined_curvature[i] <- sqrt(curve_vec[i]**2 + torse_vec[i]**2)
    curve_check_vec[i] <- sum(diff_nor[i,]*tangent[i,]) + curve_vec[i]*sum(tangent[i,]*tangent[i,])
    torse_check_vec[i] <- sum(diff_nor[i,]*binormal[i,]) - torse_vec[i]*sum(binormal[i,]*binormal[i,])
  }
  
  TC <- sum(abs(curve_vec), na.rm = T)
  AC <- mean(abs(curve_vec), na.rm = T)
  TT <- sum(torse_vec, na.rm = T)
  AT <- mean(torse_vec, na.rm = T)
  MC <- max(curve_vec, na.rm = T)
  MT <- max(abs(torse_vec), na.rm = T)
  TCC <- sum(abs(combined_curvature), na.rm = T)
  ACC <- mean(abs(combined_curvature), na.rm = T)
  
  return(list(curvature = curve_vec, torsion = torse_vec, arclength = arc_length, TC = TC, AC = AC, TT = TT, AT = AT, MC = MC, MT = MT, TCC = TCC, ACC = ACC, curvature_check = curve_check_vec, torsion_check = torse_check_vec))
}

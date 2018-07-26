vessel_spline_fit <- function(vessel, number_samples = 10000, smoothing = NULL, plot = TRUE){
  # This function performs a spline smoothing fit of the provided voxels for a given vessel.
  # The vessel parameter needs to be an N X 3 matrix of the (x, y, z) coordinates of the vessel voxels.
  # The number_samples parameter determines how many new points will be formed for defining the output vessel.
  # Currently we just define it as a single value, but in the future this should be modified to be representative of samples 
  # per unit arclength to maintain a constant number of samples through curves.  
  # The smoothing parameter is a place holder for now that will be used to toggle the extent of smoothing.  
  # This can be more complicated than just a single number.  
  # See ?smooth.spline for more information on the use of the smoothing parameter, which in the documentation is referred to as "spar".
  library("pracma")
  
  source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/single_vessel_plotter.R")
  
  
  vessel_steps <- nrow(vessel)
  
  if(plot){
    open3d()
    single_vessel_plotter(vessel_coords = vessel)
  }
  
  # filtering will happen were.  A default smoothing parameter value is estimated, however we can control that by varying the value of spar and removing the use of cv.  
  # Note that cv is a leave-one-out cross-validation procedure for fitting.  
  # Annoyingly, the use of cv = TRUE yields the best fits to known curves (common functions) even with noise, but it appears to grossly overfit real vessel data.
  
  if(!is.null(smoothing)){
    fitx <- smooth.spline(x = c(1:vessel_steps), y = vessel[,1], spar = smoothing)
    fity <- smooth.spline(x = c(1:vessel_steps), y = vessel[,2], spar = smoothing)
    fitz <- smooth.spline(x = c(1:vessel_steps), y = vessel[,3], spar = smoothing)
  }else{
    fitx <- smooth.spline(x = c(1:vessel_steps), y = vessel[,1], cv = TRUE)
    fity <- smooth.spline(x = c(1:vessel_steps), y = vessel[,2], cv = TRUE)
    fitz <- smooth.spline(x = c(1:vessel_steps), y = vessel[,3], cv = TRUE)
  }
  
  # Here is where we interpolate additional points to yield (ideally) more contiunuous measurements of vessel coordinates.
  fit_length <- number_samples
  
  x <- predict(fitx,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
  y <- predict(fity,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
  z <- predict(fitz,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
  
  smth_vessel <- cbind(x, y, z)
  
  if(plot){
    single_vessel_plotter(vessel_coords = smth_vessel, new = TRUE)
  }
  
  # Here we initialize the Frenet-Serrat frame vectors.
  tangent_array <- mat.or.vec(fit_length, 3)
  normal_array <- mat.or.vec(fit_length, 3)
  binormal_array <- mat.or.vec(fit_length, 3)
  
  tangent_array[,] <- NaN
  normal_array[,] <- NaN
  binormal_array[,] <- NaN
  
  for(i in 1:(length(smth_vessel[,1]) - 2)){
    tangent_array[i+1,] <- tangent_vector(smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
    normal_array[i+1,] <- normal_vector(smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
    binormal_array[i+1,] <- binormal_vector(smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
  }
  
  return(list(tangent = tangent_array, normal = normal_array, binormal = binormal_array, vssl_coords = smth_vessel))
}

vessel_poly_fit <- function(vessel, number_samples = 10000, poly_degree = 5, plot = TRUE){
  ## Generic comments about this function...
  
  library("pracma")
  
  source("/Users/kaitlinylim/Documents/angicart/tortuosity_R-master/code_and_sample_data/single_vessel_plotter.R")
  
  vessel_steps <- nrow(vessel)
  
  if(plot){
    open3d()
    single_vessel_plotter(vessel_coords = vessel)
  }
  
  ## Here begins the polynomial fitting
  ## We need to add a column indexing the voxels.  This appears to be needed to the way that the predict() function identifies and tracks data.  
  fit_length <- number_samples
  
  vessel <- data.frame(vessel, "voxel" = 1:vessel_steps)
  
  fitx <- lm(formula = x ~ poly(voxel, degree = poly_degree, raw = TRUE), data = vessel)
  fity <- lm(formula = y ~ poly(voxel, degree = poly_degree, raw = TRUE), data = vessel)
  fitz <- lm(formula = z ~ poly(voxel, degree = poly_degree, raw = TRUE), data = vessel)
  
  smth_vessel <- data.frame(voxel = seq(from = 1, to = vessel_steps, length.out = 1000))

  smth_vessel$predx <- predict(fitx, smth_vessel)
  smth_vessel$predy <- predict(fity, smth_vessel)
  smth_vessel$predz <- predict(fitz, smth_vessel)
  
  smth_vessel <- smth_vessel[,-1]
  smth_vessel <- as.matrix(smth_vessel)

  
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

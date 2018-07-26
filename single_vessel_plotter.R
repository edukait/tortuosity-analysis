single_vessel_plotter <- function(vessel_coords, centered = TRUE, frenet = FALSE, scale = 1, frenet_vectors, col_metric = NULL, ...){
  ## This function will plot vessel coordinates, whether they are source voxels or smoothed and interpolated coordinates.  The parameter center is needed if one wants to also plot the Frenet vectors in order to keep the vectors each at unit length.  The frenet parameter toggles whether or not the Frenet vectors are calculated and subsequently plotted.  The scale parameter is a multiplying factor for rescaling the Frenet vectors in case they are too large or small for viewing.  The frenet_vectors parameter determines the vessel coordinates at which to plot the Frenet vectors.  This object needs to be a vector that will be used to index the Frenet vectors calculated.  The col_metric parameter is used for color coding the vessels based on the values provided in col_metric.  col_metric needs to be a numeric vector with length equal to the number of rows in vessel_coords.  Examples for color coding are curvature and torsion.
  library(rgl)
  
  # Centered = TRUE performs a subroutine for centering the graph.  This is helpful if plotting the Frenet vectors for keeping them all of unit-length.
  if(centered){
    mean_x <- mean(vessel_coords[,1], na.rm = T)
    mean_y <- mean(vessel_coords[,2], na.rm = T)
    mean_z <- mean(vessel_coords[,3], na.rm = T)
    
    mean_vec <- c(mean_x, mean_y, mean_z)
    
    for(i in 1:nrow(vessel_coords)){
      vessel_coords[i,] <- vessel_coords[i,] - mean_vec
    }
  }
  
  x_range <- range(vessel_coords[,1])[2] - range(vessel_coords[,1])[1]
  y_range <- range(vessel_coords[,2])[2] - range(vessel_coords[,2])[1]
  z_range <- range(vessel_coords[,3])[2] - range(vessel_coords[,3])[1]
  
  bbox_range <- round(max(x_range, y_range, z_range)/2) + 1
  
  ## Here is the call to specify the colors of the vessel.  If col_metric is provided, then the vessel will be color coded based on the values of col_metric (for example, curvature values or torsion values).  Otherise, the vessel will be made black.
  if(!is.null(col_metric)){
    library(RColorBrewer)
    colors <- colorRampPalette(brewer.pal(20, "Spectral"))
    colors_key <- as.numeric(cut(col_metric, 20))
    col <- colors(20)[colors_key]
  }else{
    col <- "black"
  }
  
  ## Here is the call to plot only the vessel without Frenet vectors.
  if(!frenet){
    plot3d(vessel_coords[,1], vessel_coords[,2], vessel_coords[,3], xlim = c(-bbox_range, bbox_range), ylim = c(-bbox_range, bbox_range), zlim = c(-bbox_range, bbox_range), col = col, ...)
  }
  
  ## Here is the call to calculate the Frenet vectors for plotting with the vessel.
  if(frenet){
    tangent_array <- mat.or.vec(nrow(vessel_coords), 3)
    normal_array <- mat.or.vec(nrow(vessel_coords), 3)
    binormal_array <- mat.or.vec(nrow(vessel_coords), 3)
    
    tangent_array[,] <- NaN
    normal_array[,] <- NaN
    binormal_array[,] <- NaN
    
    for(i in 1:(nrow(vessel_coords) - 2)){
      tangent_array[i+1,] <- tangent_vector(vessel_coords[i,], vessel_coords[i+1,], vessel_coords[i+2,])
      normal_array[i+1,] <- normal_vector(vessel_coords[i,], vessel_coords[i+1,], vessel_coords[i+2,])
      binormal_array[i+1,] <- binormal_vector(vessel_coords[i,], vessel_coords[i+1,], vessel_coords[i+2,])
    }
    
    plot3d(vessel_coords[,1], vessel_coords[,2], vessel_coords[,3], xlim = c(-bbox_range, bbox_range), ylim = c(-bbox_range, bbox_range), zlim = c(-bbox_range, bbox_range), col = col, ...)
    for(i in frenet_vectors){
      if(!is.nan(binormal_array[i,1])){
        arrow3d(p0 = vessel_coords[i,], p1 = c(vessel_coords[i,1] + scale*tangent_array[i,1], vessel_coords[i,2] + scale*tangent_array[i,2], vessel_coords[i,3] + scale*tangent_array[i,3]), type = "rotation", col = "blue")
        arrow3d(p0 = vessel_coords[i,], p1 = c(vessel_coords[i,1] + scale*normal_array[i,1], vessel_coords[i,2] + scale*normal_array[i,2], vessel_coords[i,3] + scale*normal_array[i,3]), type = "rotation", col = "red")
        arrow3d(p0 = vessel_coords[i,], p1 = c(vessel_coords[i,1] + scale*binormal_array[i,1], vessel_coords[i,2] + scale*binormal_array[i,2], vessel_coords[i,3] + scale*binormal_array[i,3]), type = "rotation", col = "green")
      }
    }
  }
}
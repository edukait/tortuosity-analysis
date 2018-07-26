### This code is for calculating Frenet frame vectors given an array of spatial coordinates in 3D.  Results can be used for caclulating various tortuosity metrics given that rely first on having Frenet frame vectors at each point in space.  Note that our implementation requires three points to define a single frenet frame.  Thus, if data consists of N points, N-2 vectors will be generated.  We need to incorporate five point stencils for use instead of three-points stencils for higher accuracy.

## To calculate Frenet frame vectors, we will need vector operations of cross products and dot products.  Fortunately these can be implemented with the package "pracma"

library("pracma")

tangent_vector <- function(pt1, pt2, pt3){
  ## This function calculates the unit tangent vector of pt2.  The idea is to calculate the velocity vectory, then normalize by it's length.  That is, T = V/|V|. This can be done both directly, and by using the velocity vectory function.  We will write both, but comment out that which uses the functions.
  v <- pt3 - pt1
  mag_v <- sqrt(sum(v*v))
  t <- v/mag_v
  return(t)
}

normal_vector <- function(pt1, pt2, pt3){
  ## This function calculates the unit normal vector of pt2.  To do so, we will need both the velocity and acceleration vectors.  The formula for the unit normal vector is N = V X (A X V) / |V X (A X V)|.  Note that Bullitt et al. is regrettably vague about the order of operations of V X A X V, and it must either be determined from the provided Figure, or by looking up the definition elsewhere.  As usual, we can calculate either directly, or call upon previously defined functions.  We will write both, but comment out that which uses the functions.
  v <- pt3 - pt1
  a <- pt3 - 2*pt2 + pt1
  n <- cross(v, cross(a, v))
  mag_n <- sqrt(sum(n*n))
  N <- n/mag_n
  return(N)
}

binormal_vector <- function(pt1, pt2, pt3){
  ## This function calculates the unit binormal vector of pt2.  To do so, we will need the unit tangent and unit normal vectors.  The formula for the unit binormal vector is B = T X N.  As the unit tangent and normal vectors should already be properly normalized, we need not normalize this vector after calculation of it's direction.
  t <- tangent_vector(pt1, pt2, pt3)
  n <- normal_vector(pt1, pt2, pt3)
  B <- cross(t, n)
  return(B)
}

frenet_frame_calc <- function(vessel_coords){
  ## This function calculates the Frenet frame vectors along the length of a vessel.  It is simply a loop function running the frevet vector functions along the appropriate voxel coordinates.
  
  if(length(vessel_coords) >= 4){
    vessel_coords <- vessel_coords[,-c(4:ncol(vessel_coords))]
  }
  
  if(class(vessel_coords) != "matrix"){
    vessel_coords <- as.matrix(vessel_coords)
  }
  
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
  
  frenet_array <- cbind(tangent_array, normal_array, binormal_array)
  frenet_array <- as.data.frame(frenet_array)
  names(frenet_array) <- c("tanx", "tany", "tanz", "norx", "nory", "norz", "binx", "biny", "binz")
  
  
  return(frenet_array)
}
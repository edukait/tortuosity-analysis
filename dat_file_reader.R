dat_file_reader <- function(dat_filename){
  
  ## This function is intended for use as a dat file reader for voxel backbone output.  It will automatically convert the voxels to physical units based on the voxel dimensions provided in the filename.  Furthermore, it will attach a column of vessel IDs to individual dat files.  The function returns a data frame of the backbone (x,y,z) coordinates with ID factors.
  
  vessel_coords <- read.table(file = dat_filename, blank.lines.skip = FALSE)
  vessel_coords <- rbind(c(rep(NA, 3)), vessel_coords)
  names(vessel_coords) <- c("x", "y", "z")
  
  ## Convert pixel values to micrometer distances.  Use values from file name.
  xfactor <- as.numeric(strsplit(grep("\\d+[x]\\d+[x]\\d+", strsplit(dat_filename, "_")[[1]], value = TRUE), "x")[[1]][1])/1000
  yfactor <- as.numeric(strsplit(grep("\\d+[x]\\d+[x]\\d+", strsplit(dat_filename, "_")[[1]], value = TRUE), "x")[[1]][2])/1000
  zfactor <- as.numeric(strsplit(grep("\\d+[x]\\d+[x]\\d+", strsplit(dat_filename, "_")[[1]], value = TRUE), "x")[[1]][3])/1000
  
  vessel_coords$x <- vessel_coords$x*xfactor
  vessel_coords$y <- vessel_coords$y*yfactor
  vessel_coords$z <- vessel_coords$z*zfactor
  
  ## Indeces of line breaks between individual vessels.
  vessel_breaks <- which(is.na(vessel_coords[,1]))
  
  ## Number of inividual vessels.
  length(which(is.na(vessel_coords[,1])))
  
  ## Need to build ID tag column for vessels based on location of vessel_breaks.
  ID <- mat.or.vec(length(vessel_coords[,1]), 1)
  ID_tag <- 1
  for(i in 1:(length(vessel_breaks)-1)){
    for(j in (vessel_breaks[i]+1):(vessel_breaks[i+1]-1)){
      ID[j] <- ID_tag
    }
    ID_tag <- ID_tag + 1
  }
  
  ## Add ID tag column
  vessel_coords <- data.frame("x" = vessel_coords$x, "y" = vessel_coords$y, "z" = vessel_coords$z, "ID" = ID)
  
  ## Remove vessel line breaks
  vessel_coords <- vessel_coords[-vessel_breaks,]
  
  return(vessel_coords)
}
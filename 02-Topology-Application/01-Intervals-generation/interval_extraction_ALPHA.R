
"""
Purpose of this script is to generate the Rips and Alpha Complex intervals for each C1 atom.

Input: Excel file
Output: CSV file
"""

# Load package
# Install TDA package
  # TDA package: https://cran.r-project.org/web/packages/TDA/index.html
  # file_name_with_path = 'C:/Users/brandon.yong/Downloads/bio3d_2.3-4.zip'
  # install.packages(file_name_with_path, repos = NULL, type = 'source')
library("TDA") # Generate interval

# Install xlsx package
  # xlsx package: https://cran.r-project.org/web/packages/xlsx/index.html
  # file_name_with_path = 'C:/Users/brandon.yong/Downloads/xlsx_0.6.1.zip'
  # install.packages(file_name_with_path, repos = NULL, type = 'source')
library("xlsx") # Read excel file


# STEP 2
# All relevant pdb files are in a single folder
main_folder <- "D:/PHML B factor estimation/02 Project"
subfolder_input <- "01 Download and extract/03 B factor normalization/Distance - %d"
subfolder_output <- "02 Topology Application/01 Intervals generation/Distance - %d"

# For each distance threshold
for (i in seq(10, 11, by = 5)) {

  # Create path
  path <- paste(main_folder, sprintf(subfolder_input, i), sep = "/") 
  files <- list.files(path, full.names = TRUE)

  # For each file, generate interval
  for (file in files) {

    # Read excel file
    data <- read.xlsx2(file, sheetIndex = 1)
    coordinates <- data[c("x","y","z")]
    b_factor <- data[data["marker"] == 1][7]
    atom_name <- strsplit(file, split = "/")[[1]][[7]]
    file_name <- paste("Alpha", sub(".xlsx", "", atom_name), b_factor, sep = "=")
    file_name <- paste(file_name, ".csv", sep = "")

    # print status
    print(atom_name)

    # Generate and save intervals
    diagram <- alphaComplexDiag(X = data.matrix(coordinates), library = c("GUDHI","Dionysus"), printProgress = TRUE) # Generate interval
    interval <- diagram[["diagram"]] # Extract the intervals
    diagram.scaled <- cbind(interval[, 1, drop = F], 2*sqrt(interval[, c(2,3), drop = F]) )

    # Save the intervals 
    write.csv(diagram.scaled, file = paste(main_folder, sprintf(subfolder_output,i), file_name, sep = "/" ), row.names = F)

    # Generate jpeg diagram
    title <- paste("Alpha", sub(".xlsx", ".jpg", atom_name), sep = "=")
    jpeg(file = paste(main_folder, sprintf(subfolder_output,i), title, sep = "/" ))
    title <- sub(".xlsx", "", atom_name)
    try(plot(x = interval, barcode = T))
    dev.off()

  }
}

'''
Purpose is to pull the clean extracted pdb files and extract the coordinates and the B factor. The data are then saved as csv for future computation.

Input: PDB files
Output: CSV
'''


# Load package
library("bio3d")

letter_list <- c("C1","all atoms","N","O","P") # Choices are C1, all atoms or one of CNOP
for (letter in letter_list){
	main_folder = 'D:/PHML B factor estimation/02 Project'
	input_folder = '01 Download and extract/02 Processed PDB files/C1-%s'
	input_folder = sprintf(input_folder, letter)
	path = paste(main_folder, input_folder, sep = "/")

	# For PDB format files
	pdbFiles <- list.files(path = path, pattern = "Clean_pdb", full.names = FALSE)
	filename <- gsub('.{4}$', "", pdbFiles)
	filename <- gsub("Clean_pdb_", "Required_", filename)
	files <- file.path(path, c(pdbFiles))
	pdb.list <- list()
	pdb.list <- lapply(files, read.pdb)


	# Data will be in a dataframe when extracted
	for (i in 1:length(pdbFiles)){
	  	temp <- pdb.list[[i]]$atom[, c("elety","x","y","z","b")]  
	  	path <- paste(main_folder, input_folder, filename[i], sep = "/")
	  	write.csv(temp, file = paste(path, ".csv", sep = ""), row.names = F)	}


	path = paste(main_folder, input_folder, sep = "/")
	# For CIF format files
	pdbFiles <- list.files(path = path, pattern = "Clean_cif", full.names = FALSE)
	filename <- gsub('.{4}$', "", pdbFiles)
	filename <- gsub("Clean_cif_", "Required_", filename)
	files <- file.path(path, c(pdbFiles))
	pdb.list <- list()
	pdb.list <- lapply(files, read.cif)

	# Data will be in a dataframe when extracted
	for (i in 1:length(pdbFiles)){
	  	temp <- pdb.list[[i]]$atom[, c("elety","x","y","z","b")]  
	  	path <- paste(main_folder, input_folder, filename[i], sep = "/")
	  	write.csv(temp, file = paste(path, ".csv", sep = ""), row.names = F)	}

}
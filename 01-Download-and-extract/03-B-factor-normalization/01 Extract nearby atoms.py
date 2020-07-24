"""
Purpose is to normalize the B factors for each chain. 
After that, for each C1 atom, neighbouring atoms that are <= 10 Angstrom distance away are identified and saved into separate files.

Input: CSV file
Output: CSV file

"""
import time
#time.sleep(6*60*60)

# Import packages
import numpy as np
from statsmodels import robust
import math
import os
import pandas as pd
from scipy.spatial import distance

np.seterr(invalid='ignore')

#letter_list = ["all atoms","C","N","O","P"]
letter_list = ["C"]
#letter = "all atoms"
euclidean_distance_dict = {	"C1":range(45, 56, 10),
							"C": range(35, 41, 5)}
"""euclidean_distance_dict = {"C1":range(10,36,5),
							"C":range(5,31,5), 
							"N":range(5,31,5),
							"O":range(5,31,5),
							"P":range(5,31,5),
							"all atoms":range(5,31,5)}
"""
#letter = "all atoms" # Choices are C1, all atoms or one of CNOP
#euclidean_distances = range(5, 16, 5)

"""
Atoms 		Range
all atoms 	5 -15
C1 			10 - 30
C 			5 - 20
N 			5 - 20
O 			5 - 20
P 			5 - 20
"""

for letter in letter_list:
	euclidean_distances = euclidean_distance_dict[letter]

	#main_folder = r'D:\PHML B factor estimation\02 Project'
	#subfolder_input = r'01 Download and extract\02 Processed PDB files\C1-%s' %letter
	#subfolder_output = r'01 Download and extract\03 B factor normalization'

	main_folder = "/home/brandon.yong/02 Project"
	subfolder_input = '01 Download and extract/02 Processed PDB files/C1-%s' %letter
	subfolder_output = '01 Download and extract/03 B factor normalization'

	# Create folder if does not exist
	folders = os.listdir(os.path.join(main_folder, subfolder_output))
	if "C1-%s" %letter not in folders:
		os.mkdir(os.path.join(main_folder, subfolder_output, "C1-%s" %letter))

	# Create subfolder if does not exist
	subfolder_output = os.path.join(subfolder_output, "C1-%s" %letter)
	folders = os.listdir(os.path.join(main_folder, subfolder_output))
	for euclidean_distance in euclidean_distances:
		if "Distance - %s" %euclidean_distance not in folders:
			os.mkdir(os.path.join(main_folder, subfolder_output, "Distance - %s" %euclidean_distance))

	# List of input csv files
	files = [f for f in os.listdir(os.path.join(main_folder, subfolder_input)) if f.endswith(".csv")]

	# Load the dataset
	dataset = []
	for file in files:
		input_file = os.path.join(main_folder, subfolder_input, file)
		dataset.append(pd.read_csv(input_file))

	# B factor calculation 
	for j, full_data in enumerate(dataset):

		# Replace b factors of non-C1 atoms with np.nan
		criteria = full_data["elety"].apply(lambda x: "C1" in x)
		data_c1 = full_data.loc[criteria,:]
		data_others = full_data.loc[~criteria,:]

		# Compute median and median absolute displacement
		median = np.nanmedian(data_c1["b"])
		MAD = robust.mad(data_c1["b"])*0.6745

		# If the MAD is NaN, it is calculated manually
		if math.isnan(MAD):
			temp = [abs(n - median) for n in data_c1["b"]]
			MAD = np.nanmedian(temp)

		# Calculate the distance from the median and determine which are outliers
		data_c1["dist"] = 0.6745*(data_c1["b"] - median)/MAD
		data_c1["outlier_truth"] = data_c1["dist"].apply(lambda x: 1 if x < 3.5 else np.nan) # 1 if not outlier

		# Normalize the B factors
		data_c1["b_norm"] = data_c1["b"]*data_c1["outlier_truth"]
		mean = np.nanmean(data_c1["b_norm"])
		std = np.nanstd(data_c1["b_norm"])
		data_c1["b_norm"] = (data_c1["b_norm"]-mean)/std

		# concatenate the two dataframes
		dataset[j] = pd.concat([data_c1, data_others], axis = 0, sort = False).reset_index(drop = True)



	# For each dataset, calculate the distance of each atom from each other 
	# Local compactness
	for j, data in enumerate(dataset):
		chain_name = files[j].strip(".csv")[9:]

		# Calculate euclidean distance for each pair of points
		values_1 = data[["x","y","z"]].values
		#values_2 = values_1.reshape(values_1.shape[0], 1, values_1.shape[1])
		#distances = np.sqrt(np.einsum('ijk, ijk->ij', values_1 - values_2, values_1 - values_2)) # shape = n x n
		distances = distance.cdist(values_1, values_1, metric = "euclidean")

		c1_index = data['b_norm'].index[~data['b_norm'].apply(np.isnan)]

		for euclidean_distance in euclidean_distances:
			print("Processing C1-%s \tchain: %s \tEuclidean distance: %s" %(letter, chain_name, euclidean_distance))
			distances_1 = (distances <= euclidean_distance)*1

			for num, index in enumerate(c1_index):

				# Filter out nearby atoms for each C1
				data["neighbouring_atoms"] = distances_1[index]
				data["marker"] = 0
				data.loc[index, "marker"] = 1 # Identify the C1 atom in question so that I'll know which b_norm to use
				temp = data[data["neighbouring_atoms"] == 1]

				# Save the filtered data as xlsx
				file_name = os.path.join(main_folder, subfolder_output, r"Distance - %s\%s - %s.xlsx" %(euclidean_distance, chain_name, num+1))
				writer = pd.ExcelWriter(file_name)
				temp.to_excel(writer, index = False)
				writer.save()

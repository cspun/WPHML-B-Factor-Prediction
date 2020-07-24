"""
Purpose of this script is to generate the Rips and Alpha Complex intervals for each C1 atom.

Input: Excel file
Output: CSV file
"""

# STEP 1
# Import packages
import time
#time.sleep(3*60*60)

import pandas as pd
from ripser import ripser
import os
import matplotlib.pyplot as plt

# STEP 2
# All relevant pdb files are in a single folder
main_folder = r'C:\Users\User\Desktop\B factor'
subfolder_input = r"01 Download and extract\03 B factor normalization\C1-%s\Distance - %d"
subfolder_output = r"02 Topology Application\01 Intervals generation\C1-%s\Scale - 2d\Distance - %d"

#main_folder = "/home/brandon.yong/02 Project"
#subfolder_input = "01 Download and extract/03 B factor normalization/C1-%s/Distance - %d"
#subfolder_output = "02 Topology Application/01 Intervals generation/C1-%s/Scale - 2d/Distance - %d"


#Initial parameters
max_dimension = 1 # Maximum betti number
#letter = "O"
save_image_truth = False
group = 3

letter_list = ["C"]
for letter in letter_list:
	# For each euclidean distance
	for euclidean_distance in [45]:
		# Adjust the max scale such that it is twice that of the distance threshold
		max_scale = 1.5*euclidean_distance # Max filtration distance (Angstrom)

		# Save image
	#		if euclidean_distance == 30:
	#			save_image_truth = True
	#		else:
	#			save_image_truth = False
		

		# Create path
		path = os.path.join(main_folder, subfolder_input %(letter, euclidean_distance))
		files = os.listdir(path)

		# C, E = 45
		if letter == "C" and euclidean_distance == 45 and group == 1:
			files = [f for f in files if f.startswith('4lnt_RA - 1888') or f.startswith('4lnt_RA - 1889') or 
						f.startswith('4lnt_RA - 189') or f.startswith('4lnt_RA - 19') ]

		if letter == "C" and euclidean_distance == 45 and group == 3:
			files = [f for f in files if f.startswith('4lnt_RA - 774') ]
			files = reversed(files)

		# For each file, generate interval
		for file in files:

			# Read excel file
			full_path = os.path.join(main_folder, subfolder_input %(letter, euclidean_distance), file)
			data = pd.read_excel(full_path)

			# Remove excess carbons for NOP
			criteria_bfactor = data["marker"] == 1
			criteria_carbon = data["elety"].str.contains("C")
			b_factor = data[criteria_bfactor]["b_norm"].values[0]

			# Relevant information
			if letter in list("NOP"):
				data = data[~(~criteria_bfactor & criteria_carbon)]
			coordinates = data[["x","y","z"]].values
			
			# Rips complex
			if b_factor != "":
				atom_name = file.strip(r".xlsx")
				file_name = "Rips=" + atom_name + "=" + str(b_factor) + ".csv"
				file_output_path = os.path.join(main_folder, subfolder_output %(letter, euclidean_distance), file_name)

				# Print status
				print("Atoms - %s \tDistance : %s \tChain : %s" %(letter, euclidean_distance, atom_name))

				# Generate diagram
				diagram = ripser(coordinates, maxdim = max_dimension, thresh = max_scale)["dgms"]

				dimensions, births, deaths = [], [], []
				for i in range(max_dimension+1):
					temp = diagram[i]
					dimensions.extend([i]*len(temp))
					births.extend(temp[:, 0])
					deaths.extend(temp[:, 1])

				intervals = pd.DataFrame({"dimension":dimensions,"Birth":births, "Death":deaths})
				intervals.to_csv(file_output_path, index = False)
            
			del diagram, dimensions, births, deaths
				# Generate jpeg file
	#			if save_image_truth:
	#				title = "Rips=" + atom_name + r".jpg"
	#				image_output_path = os.path.join(main_folder, subfolder_output %(letter, euclidean_distance), title)
	#				gd.plot_persistence_barcode(diagram, alpha = 1)
	#				plt.title("Persistence barcode - %s" %atom_name)
	#				plt.xlim(0)
	#				plt.savefig(image_output_path, bbox_inches = "tight")
	#				plt.close()

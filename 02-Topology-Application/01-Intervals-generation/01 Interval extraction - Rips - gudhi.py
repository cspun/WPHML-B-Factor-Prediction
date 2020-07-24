"""
Purpose of this script is to generate the Rips and Alpha Complex intervals for each C1 atom.

Input: Excel file
Output: CSV file
"""

# STEP 1
# Import packages
import pandas as pd
import gudhi as gd
import os
import matplotlib.pyplot as plt
import time

# STEP 2
# All relevant pdb files are in a single folder
main_folder = r"C:\Users\User\Desktop\B factor"
subfolder_input = r"01 Download and extract\03 B factor normalization\C1-%s\Distance - %d"
subfolder_output = r"02 Topology Application\01 Intervals generation\C1-%s\Scale - 2d\Distance - %d"

#main_folder = "/home/brandon.yong/02 Project"
#subfolder_input = "01 Download and extract/03 B factor normalization/C1-%s/Distance - %d"
#subfolder_output = "02 Topology Application/01 Intervals generation/C1-%s/Scale - 2d/Distance - %d"


#Initial parameters
max_dimension = 1 # Maximum betti number
letter = "O"
save_image_truth = False


# For each euclidean distance
for euclidean_distance in reversed(range(30, 31, 5)):
	# Adjust the max scale such that it is twice that of the distance threshold
	max_scale = 2*euclidean_distance # Max filtration distance (Angstrom)

	# Save image
#		if euclidean_distance == 30:
#			save_image_truth = True
#		else:
#			save_image_truth = False
	

	# Create path
	path = os.path.join(main_folder, subfolder_input %(letter, euclidean_distance))
	files = os.listdir(path)
	files = [f for f in files if f.startswith('4v9o_BA')]

    
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
			rips_complex = gd.RipsComplex(points = coordinates, max_edge_length = max_scale).create_simplex_tree(max_dimension = max_dimension + 1)
			diagram = rips_complex.persistence()
			del rips_complex
            
			time.sleep(5)
			dimensions = [f[0] for f in diagram]
			births, deaths = map(list, zip(*[f[1] for f in diagram]))
			intervals = pd.DataFrame({"dimension":dimensions,"Birth":births, "Death":deaths})
			intervals.sort_values("dimension", inplace = True)
			intervals.to_csv(file_output_path, index = False)

			# Generate jpeg file
#			if save_image_truth:
#				title = "Rips" + atom_name + r".jpg"
#				image_output_path = os.path.join(main_folder, subfolder_output %(letter, euclidean_distance), title)
#				gd.plot_persistence_barcode(diagram, alpha = 1)
#				plt.title("Persistence barcode - %s" %atom_name)
#				plt.xlim(0)
#				plt.savefig(image_output_path, bbox_inches = "tight")
#				plt.close()

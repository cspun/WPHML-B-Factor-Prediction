"""
Purpose of this is to extract the relevant ATOM lines from the raw pdb/cif files and then save them into a separate pdb files.

Each pdb file corresponds to a specific chain in each RNA (as per the dataset).
"""

import os

main_folder = r'D:\PHML B factor estimation\02 Project'
input_folder = r'01 Download and extract\01 Raw PDB files'
output_folder = r'01 Download and extract\02 Processed PDB files\C1-C1'
list_folder = r'00 Basic information'

# Extract the list of pdb files downloaded
input_path = os.path.join(main_folder, input_folder)
pdb_files_extension = {}
for file in os.listdir(input_path):
	name, ext = file.split(".")
	pdb_files_extension[name] = ext

# Extract the list of RNA proteins required (along with the residue type)
sources = ["Test_list", "Train_list"]
RNAproteins = []
for src in sources:
	list_path = os.path.join(main_folder, list_folder, src)
	file = open(list_path, "r").read().split("\n")[:-1]
	RNAproteins.extend(file)

# Split the above proteins and store RNA proteins and its residues in a dictionary
RNAproteins_dict = {}
for chain in RNAproteins:
	RNA, residue = chain.split("_")
	if RNA in RNAproteins_dict:
		RNAproteins_dict[RNA].append(residue)
		# Certain RNA may have multiple residue. Hence, it is attached rather than overwritten
	else:
		RNAproteins_dict[RNA] = [residue]

# Extract the required ATOM lines from the pdb files and save them into a separate pdb files
for chain in list(RNAproteins_dict.keys()):
	print("Extracting ... %s" %(chain))

	# Check if the pdb file is downloaded. Else, generate an empty pdb file with suffix '_missing'
	if chain in pdb_files_extension.keys():

		# This portion is for pdb files
		if pdb_files_extension[chain] == "pdb":
			for j in RNAproteins_dict[chain]:

				# Generate the input and output file name
				input_filename = chain + ".pdb"
				output_filename = "Extract_pdb_" + chain + "_" + j + ".pdb"

				# Input and output path
				input_path = os.path.join(main_folder, input_folder, input_filename)
				output_path = os.path.join(main_folder, output_folder, output_filename)

				oldfile = open(input_path, "r")
				newfile = open(output_path,"w")

				# Extract the residue lines
				for line in oldfile:
					if line:
						if line[0:4] == "ATOM" and line[21] == j and line[13:16] == "C1'":
							newfile.write(line)

				oldfile.close()
				newfile.close()

		# This portion is for cif files
		else:
			for j in RNAproteins_dict[chain]:

				# Generate the input and output files
				input_filename = chain + ".cif"
				output_filename = "Extract_cif_" + chain + "_" + j + ".pdb"

				# Input and output path
				input_path = os.path.join(main_folder, input_folder, input_filename)
				output_path = os.path.join(main_folder, output_folder, output_filename)

				oldfile = open(input_path, "r")
				newfile = open(output_path,"w")

				# Extract the residue lines
				for line in oldfile:
					if line:
						if line[0:4] == "ATOM" and j in line[90:95].split(" ") and line[18:21] == "C1'":
							newfile.write(line)

				oldfile.close()
				newfile.close()

	# This portion for missing pdb files
	else:
		output_filename = "Missing_" + chain + ".pdb"
		newfile = open(output_filename,"w")
		newfile.close()

print("Extraction complete ...")
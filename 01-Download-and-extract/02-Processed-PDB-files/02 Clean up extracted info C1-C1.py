"""
Purpose of this script is to clean up the extracted information from the raw pdb files.

Clean up activities are:
1. Select C1 atom at the same position with the higher occupancy factor.
	If the occupancy factors are 0.5 each, then the first C1 atom is selected.

2. some thing about the insertion code

Inputs: 	Extract_<pdb/cif>_<chain name>.pdb
Outputs: 	Clean_<pdb/cif>_<chain name>.pdb
"""



# STEP 1
# Import packages
import os
import re

letter = "C1 test"
main_folder = r'D:\PHML B factor estimation\02 Project'
input_folder = r'01 Download and extract\02 Processed PDB files\C1-%s' %letter

# STEP 2
# Extract list of extraction file
extraction_files = []
input_path = os.path.join(main_folder, input_folder)
extraction_files = [f for f in os.listdir(input_path) if f.startswith("Extract_")]

# STEP 3
#extraction_files = ["Extract_pdb_1dk1_B.pdb"]
#for old_filename in example:
num = 1
ignore_insertion_code = re.compile('[b-zB-Z]').search

for old_filename in extraction_files:
	print("Cleaning up ... %s - %s" %(old_filename, num))
	num += 1
	chosen_line = "X"*100
	prev_line = "Y"*100

	# Process for pdb file
	if old_filename.split("_")[1] == "pdb":

		# Open and read/write file
		new_filename = "Clean_" + old_filename[8:]
		oldfile = open(os.path.join(main_folder, input_folder, old_filename), "r")
		newfile = open(os.path.join(main_folder, input_folder, new_filename), "w")

		# For each line in oldfile, check for occupancy factor
		for line in oldfile:

			# If previous line and current line are different positions, save previous line
			if line[23:27] != prev_line[23:27]:
				if prev_line[0] != "Y":
					newfile.write(prev_line)
				prev_line = line

			# If previous line and current line are of same position (ignoring insertion code), check for occupancy factor
			else:
				print("Multi line found")
				occupancy_1 = float(prev_line[56:60])
				occupancy_2 = float(line[56:60])

				# Line with higher occupancy factor is chosen
				if occupancy_2 > occupancy_1:
					prev_line = line

		# Write the very last line that would not be saved
		newfile.write(prev_line)

		# Save files
		oldfile.close()
		newfile.close()

	# Process for cif file
	else:
		# Open and read/write file
		new_filename = "Clean_" + old_filename[8:]
		oldfile = open(os.path.join(main_folder, input_folder, old_filename), "r")
		newfile = open(os.path.join(main_folder, input_folder, new_filename), "w")

		# For each line in oldfile, check for occupancy factor
		for line in oldfile:

			# If previous line and current line are different positions, save previous line
			if line[36:41] != prev_line[36:41]:
				if prev_line[0] != "Y":
					newfile.write(prev_line)
				prev_line = line

			# If previous line and current line are of same position (ignoring insertion code), check for occupancy factor
			else:
				print("Multi line found")
				occupancy_1 = float(prev_line[69:73])
				occupancy_2 = float(line[69:73])

				# Line with higher occupancy factor is chosen
				if occupancy_2 > occupancy_1:
					prev_line = line

		# Write the very last line that would not be saved
		newfile.write(prev_line)

		# Save files
		oldfile.close()
		newfile.close()


print("Clean up process completed ....")
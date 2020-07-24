import pandas as pd
import os

main_folder = r'C:\Users\User\Desktop\B factor'
subfolder = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\C1-%s'

letter_list = ["N","O","P"]

for letter in letter_list:
	# List of files to read
	path = os.path.join(main_folder, subfolder %letter)
	files = [f for f in os.listdir(path) if f.endswith("csv")]

	number_files = len(files)
	number = 1

	for file in files:
		# Read file
		file_path = os.path.join(path, file)
		data_edit = pd.read_csv(file_path)

		number_features = int(file.split("_")[-1].split(".")[0])

		# Extract and overwrite the necessary data only
		data_keep = data_edit.iloc[:, 2*number_features:]
		data_keep.to_csv(file_path, index = False)

		# Print progress
		print("Complete %d/%d : %s" %(number, number_files, file))
		number += 1
import pandas as pd
import os

main_folder = r'C:\Users\User\Desktop\B factor'
#main_folder = r'C:\Users\brandon.yong\Desktop\02 Project'
subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s'
subfolder_output = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\C1-%s'

letter_input = "C-N-O-P"
#letter_output = ["C","N","O","P"]
letter_output = ["C"]
input_path = os.path.join(main_folder, subfolder_input %letter_input)
files = [f for f in os.listdir(input_path) if f.endswith("csv")]
files = [f for f in files if f.startswith("%s_35_" %letter_input) or f.startswith("%s_25_" %letter_input) or
			f.startswith("%s_30_" %letter_input)]

num = len(files)
number_file = 1
for file in files:
	file_path = os.path.join(input_path, file)
	data = pd.read_csv(file_path)

	number_features = int(file.split("_")[-1].split(".")[0])
	data_bfactor = data.iloc[:, -2:]

	print("Performing file %s/%s: %s" %(number_file, num, file))
	number_file += 1
	for i, letter in enumerate(letter_output):
		data_temp = data.iloc[:, i*2*number_features : (i+1)*2*number_features]
		data_temp = pd.concat([data_temp, data_bfactor], axis = 1, ignore_index = True)

		output_file = file.replace(letter_input, "C1-%s" %letter)
		output_path = os.path.join(main_folder, subfolder_output %letter, output_file)
		data_temp.to_csv(output_path, index = False)
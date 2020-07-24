"""
Purpose of this script is to count the number of points at each interval (boundary line) for each dimension.

Input: CSV
Output: CSV
"""
#import time
#time.sleep(16*20*60) # Secs

import numpy as np
import os
import pandas as pd

# Main and subfolder
#main_folder = r"D:\PHML B factor estimation\02 Project"
#main_folder = r"C:\Users\brandon.yong\Desktop\02 Project"
main_folder = r"C:\Users\User\Desktop\B factor"
subfolder_input = r"02 Topology Application\01 Intervals generation\C1-%s\Scale - 2d\Distance - %d"
subfolder_output = r"02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s"

#main_folder = "/home/brandon.yong/02 Project"
#subfolder_input = "02 Topology Application/01 Intervals generation/C1-%s/Scale - 2d/Distance - %d"
#subfolder_output = "02 Topology Application/02 Topological features representation/ML inputs - Counts of points at interval/%s"

# Basic information
bin_sizes = [0.1, 0.05, 0.15]
#atom_types_list = [["C","N"],["C","O"],["C","P"],["N","O"],["N","P"],["O","P"]]
atom_types_list = [["C","P"]]
#atom_types = ["O", "P"]
max_dimension = 1 # Maximum betti number
float_range = np.linspace(.5, .6, 2)

for atom_types in atom_types_list:

    for euclidean_distance in [40, 45]:

        original_files, short_names = [], []
        # List of input files
        for atom_type in atom_types:
            path = os.path.join(main_folder, subfolder_input %(atom_type, euclidean_distance))
            original_file = [f for f in os.listdir(path) if f.endswith(".csv") and f.startswith("Rips")]
            short_name = [f.split("=")[1].split(" - ")[0] + "="+ f.split("=")[2].strip(r".csv")[:13] + r".csv" for f in original_file] # Chain=b_factor
            original_files.append(original_file)
            short_names.append(short_name)
        
        common_atoms = [f for f in short_names[0] if f in short_names[1]]
        atom_type = "-".join(f for f in atom_types)

        # The list of filtration distance to apply
        #filtration_distance_1 = filtration_distances[euclidean_distance]
        filtration_distance_1 = float_range*euclidean_distance

        for filtration_distance in filtration_distance_1:

            for bin_size in bin_sizes:

                # Create bin num, bin threshold, blank rows
                bin_num = int(np.ceil(filtration_distance/bin_size))
                bins = [f*bin_size for f in range(1, bin_num + 1)]
                bins = np.asarray(bins).reshape((1, len(bins)))
                blank_rows = [0]*bin_num
                
                # List of csv files
                feature_output = []
                
                for common_atom in common_atoms:
                    print("TFR \tatom %s \tEuclidean %s \tFiltration %s \tatom: %s"%(atom_type, euclidean_distance, filtration_distance, common_atom))
                    
                    # Extract the normalized b factor and chain name
                    b_factor = common_atom.split("=")[-1].strip(".csv")
                    chain_name = common_atom.split("=")[0]
                    feature_row = []

                    # Extract the file that is to be loaded
                    indexes = [short_names[0].index(common_atom), short_names[1].index(common_atom)]

                    for i, index in enumerate(indexes):

                        # Read the CSV file
                        path = os.path.join(main_folder, subfolder_input %(atom_types[i], euclidean_distance))
                        data = pd.read_csv(os.path.join(path, original_files[i][index]))
                        data.loc[data["Death"] == np.inf, "Death"] = filtration_distance*2
                        
                        # For each dimension
                        for betti_num in range(0, max_dimension + 1):
                            temp = data[data["dimension"] == betti_num].values
                            
                            # Check the subdataframe dimension
                            if len(temp) != 0:
                                counts = ((temp[:, 1].reshape((len(temp), 1)) <= bins) & (temp[:, 2].reshape((len(temp), 1)) >= bins))*1
                                counts = counts.sum(axis = 0)
                                feature_row.extend(counts) 
                                    
                            else:
                                # If the subdataframe is blank
                                feature_row.extend(blank_rows) 
                            
                    # Append b_factor and chain name
                    feature_row.extend([b_factor, chain_name])
                    
                    # Append row to dataframe
                    feature_output.append(feature_row)
                    
                # Save dataframe to csv
                feature_output = pd.DataFrame(feature_output)
                # file name format: <atom>_<euclidean distance>_<filtration distance>_<bin_size>_<bin_num>.csv
                file_name = "%s_%s_%s_%s_%s.csv" %(atom_type, euclidean_distance, filtration_distance, bin_size, bin_num)
                path_2 = os.path.join(main_folder, subfolder_output %atom_type, file_name)
                feature_output.to_csv(path_2, index = False)
"""
Purpose of this script is to train the model with kfold cross validation.
The evaluation metric is MSE and Pearson correlation.

Input: CSV file, train_list
Output: sav file (saved model), CSV file (PCC and MSE)

"""
import time
#time.sleep(2*60*60) # secs

import pandas as pd
import pickle
import os
import datetime
from random import Random

from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

# Variables
num_fold = 5
seed = 50
#num_tree = 1000
max_depth = 3
subsample = 0.8
colsample = 0.8
save_model = False
#letter = "C1-N"
range_value = "all"

num_tree_list = [50, 100, 250, 500]
letter_list = ["C1-P"]

for letter in letter_list:

    #main_folder = r'C:\Users\User\Desktop\B factor'
    main_folder = r'D:\PHML B factor estimation\02 Project'
    subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s' %letter
    subfolder_list = r'00 Basic information'
    subfolder_output_model = r'03 ML models\03 XGBoost\Trained models - Counts'
    subfolder_output_results = r'03 ML models\03 XGBoost\Results - Counts'
    writer = pd.ExcelWriter(os.path.join(main_folder, subfolder_output_results, 
                            str(datetime.datetime.today().date()) + " XGBoost results - Counts 4 %s - %s.xlsx" %(range_value, letter)))

    # List of input files
    path_input = os.path.join(main_folder, subfolder_input)
    input_files = [f for f in os.listdir(path_input) if f.endswith(".csv")]
    input_files = [f for f in input_files if "_0.15_" in f or "_0.5_" in f or "_1.0_" in f or "_1.5_" in f]
    input_files = [f for f in input_files if f.startswith("%s_10_10" %letter) or f.startswith("%s_15_15" %letter) or f.startswith("%s_20_20" %letter)]
    
    # Initial lists
    atoms, bin_sizes, bin_nums, euclidean_distances, filtration_distances = [], [], [], [], []
    num_folds, seeds, num_trees, max_depths, subsamples, colsamples = [], [], [], [], [], []
    mse_train, mse_cv, pcc_train, pcc_cv = [], [], [], []

    for file in input_files:

        # Read train list
        path_list = os.path.join(main_folder, subfolder_list, r'Train_list')
        train_list = open(path_list, "r").read().split("\n")[:-1]

        # Read the input data
        path_files = path_input + "\\" + file
        data = pd.read_csv(path_files)

        # Extract the train set
        criteria = data.iloc[:, -1].str.strip(" ").isin(train_list)
        data_set = data[criteria]
        data_set.reset_index(inplace = True, drop = True)

        # Dataset details
        file = file.strip(".csv")
        atom, euclidean_distance, filtration_distance, bin_size, bin_num = file.split("_")
        euclidean_distance, filtration_distance = int(euclidean_distance), int(float(filtration_distance))
        bin_size, bin_num = float(bin_size), int(bin_num)

        # kFold cross validation
        limit = round(len(train_list)/num_fold)
        train_index, cv_index = [], []
        for i in range(num_fold):
            Random(seed).shuffle(train_list)
            train_list_split = train_list[:limit]
            train_index_temp = data_set.iloc[:, -1].str.strip(" ").isin(train_list_split)
            cv_index_temp = ~data_set.iloc[:, -1].str.strip(" ").isin(train_list_split)
            train_index.append(train_index_temp)
            cv_index.append(cv_index_temp)

        # Split the dataset into x and y data
        data_set = data_set.iloc[:, :-1]
        data_x, data_y = data_set.iloc[:,:-1], data_set.iloc[:, -1]

        for num_tree in num_tree_list:
            # Train model 
            for i in range(num_fold):
                
                print(datetime.datetime.now(), "\tXGBoost Regressor - Training %s \tkFold = %d" %(file, i+1))
                print("Number of trees: %s" %num_tree)
                
                # Extract train and cv set
                train_x, train_y =  data_x.loc[train_index[i]], data_y[train_index[i]]
                cv_x,  cv_y =       data_x.loc[cv_index[i]],    data_y[cv_index[i]]

                # Training
                clf = XGBRegressor(max_depth = max_depth, n_estimators = num_tree, random_state = seed,
                                    subsample = subsample, colsample_bytree = colsample)
                clf.fit(train_x, train_y)
                y_pred = clf.predict(train_x)
                mse_train.append(mean_squared_error(train_y, y_pred))
                pcc_train.append(pearsonr(y_pred, train_y)[0])

                # Predict values
                y_pred = clf.predict(cv_x)
                
                # Calculate MSE and PCC
                mse_cv.append(mean_squared_error(cv_y, y_pred))
                pcc_cv.append(pearsonr(y_pred, cv_y)[0])

                # Save the trained model
                if save_model:
                    model_filename = os.path.join(main_folder, subfolder_output_model, file + "_model_%s" %(i+1) + ".sav")
                    pickle.dump(clf, open(model_filename, "wb"))

                # Print the latest results
                print("MSE train: %.5f \tPCC train: %.3f \nMSE CV: %.5f \tPCC CV: %.3f\n" %(mse_train[-1], pcc_train[-1], mse_cv[-1], pcc_cv[-1]))

            # Create lists for dataframe
            print("Creating dataframe details")
            ## Dataset details
            atoms.extend([atom]*num_fold)
            bin_sizes.extend([bin_size]*num_fold)
            bin_nums.extend([bin_num]*num_fold)
            euclidean_distances.extend([euclidean_distance]*num_fold)
            filtration_distances.extend([filtration_distance]*num_fold)

            ## Model details
            num_folds.extend([num_fold]*num_fold)
            seeds.extend([seed]*num_fold)
            num_trees.extend([num_tree]*num_fold)
            max_depths.extend([max_depth]*num_fold)
            subsamples.extend([subsample]*num_fold)
            colsamples.extend([colsample]*num_fold)

    # Create dataframe
    output_results = pd.DataFrame({
        "Atoms":atoms,
        "Euclidean distance":euclidean_distances,
        "Filtration distance":filtration_distances,
        "Bin size":bin_sizes,
        "Bin num":bin_nums,
        "Seed":seeds,
        "No trees":num_trees,
        "Max depth":max_depths,
        "Subsample":subsamples,
        "Colsample":colsamples,
        "MSE_train":mse_train,
        "MSE_CV":mse_cv,
        "PCC_train": pcc_train,
        "PCC_CV":pcc_cv
        }) 

    # Save dataframe to excel
    print("Saving results")
    output_results.to_excel(writer, index = False)
    writer.save()

    print("Completed at", datetime.datetime.now())


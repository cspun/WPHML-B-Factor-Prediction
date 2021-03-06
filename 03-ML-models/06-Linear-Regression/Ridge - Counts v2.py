"""
Purpose of this script is to train the model with kfold cross validation.
The evaluation metric is MSE and Pearson correlation.

Input: CSV file, train_list
Output: sav file (saved model), CSV file (PCC and MSE)

"""
import pandas as pd
import pickle
import os
import datetime
from random import Random

from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

# Variables
num_fold = 5
seed = 50
#alpha = 0.5
#alpha_list = [0.01, 0.05, 0.1, 0.5, 1, 5] 
alpha_list = [10, 50, 100, 500] 
save_model = False
#letter = "C1-P"
range_value = "all"
letter_list = ["C-N"]
#letter_list = ["C-N-O","C-N-P","C-O-P","N-O-P"]
for letter in letter_list:

    main_folder = r'C:\Users\User\Desktop\B factor'
    #main_folder = r'C:\Users\brandon.yong\Desktop\02 Project'
    subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s' %letter
    subfolder_list = r'00 Basic information'
    subfolder_output_model = r'03 ML models\06 Linear Regression\Trained models - Counts'
    subfolder_output_results = r'03 ML models\06 Linear Regression\Results - Counts'
    writer = pd.ExcelWriter(os.path.join(main_folder, subfolder_output_results, 
                            str(datetime.datetime.today().date()) + " Ridge results - Counts %s 2 - %s.xlsx" %(range_value, letter)))

    # List of input files
    path_input = os.path.join(main_folder, subfolder_input)
    input_files = [f for f in os.listdir(path_input) if f.endswith(".csv")]
    #input_files = [f for f in input_files if "_0.15_" in f ]
    input_files = [f for f in input_files if f.startswith("%s_45_" %letter)]

    # Initial lists
    atoms, bin_sizes, bin_nums, euclidean_distances, filtration_distances = [], [], [], [], []
    num_folds, seeds, alphas = [], [], []
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

        for alpha in alpha_list:
            # Train model 
            for i in range(num_fold):
                
                print(datetime.datetime.now(), "\tLinear Regression - Training %s \tkFold = %d" %(file, i+1))
                print("Dataset size: ", data_x.shape)

                # Extract train and cv set
                train_x, train_y =  data_x.loc[train_index[i]], data_y[train_index[i]]
                cv_x,  cv_y =       data_x.loc[cv_index[i]],    data_y[cv_index[i]]
                
                # Training
                clf = Ridge(alpha = alpha, random_state = seed)
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
            alphas.extend([alpha]*num_fold)

    # Create dataframe
    output_results = pd.DataFrame({
        "Atoms":atoms,
        "Euclidean distance":euclidean_distances,
        "Filtration distance":filtration_distances,
        "Bin size":bin_sizes,
        "Bin num":bin_nums,
        
        "Seed":seeds,
        "Alpha":alphas,
        
        "MSE_train":mse_train,
        "MSE_CV":mse_cv,
        "PCC_train": pcc_train,
        "PCC_CV":pcc_cv
        }) 

    cols = ["Atoms","Euclidean distance","Filtration distance","Bin size","Bin num",
            "Seed","Alpha",
            "MSE_train","MSE_CV","PCC_train","PCC_CV"]
    output_results = output_results.loc[:, cols]
    
    # Save dataframe to excel
    print("Saving results")
    output_results.to_excel(writer, index = False)
    writer.save()

    print("Completed at", datetime.datetime.now())

    del output_results

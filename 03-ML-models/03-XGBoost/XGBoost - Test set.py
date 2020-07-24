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
import numpy as np

from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

# Variables
num_fold = 1
seed = 50
num_tree = 50
max_depth = 3
subsample = 0.8 #0.8
colsample = 0.8 #0.8

# List of input files
letter = "C-N-O-P"
input_files = ["C-N-O-P_45_31.5_0.45_70.csv"]

# Boolean truth
save_model = True
BDP_directory = True
remove_correlated_features = True

# Directories
main_folder = r'D:\PHML B factor estimation\02 Project'
subfolder_list = r'00 Basic information'

if BDP_directory:
    subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - BDP\%s' %letter
    subfolder_output_model = r'03 ML models\03 XGBoost\Trained models - BDP'
    subfolder_output_results = r'03 ML models\03 XGBoost\Results - BDP'
else:
    subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s' %letter
    subfolder_output_model = r'03 ML models\03 XGBoost\Trained models - Counts'
    subfolder_output_results = r'03 ML models\03 XGBoost\Results - Counts'

writer = pd.ExcelWriter(os.path.join(main_folder, subfolder_output_results, 
                        str(datetime.datetime.today().date()) + " XGBoost results - Test set 4 - %s - remove.xlsx" %letter))

# Read train list
path_list = os.path.join(main_folder, subfolder_list, r'Train_list')
train_list = open(path_list, "r").read().split("\n")[:-1]
path_list = os.path.join(main_folder, subfolder_list, r'Test_list')
test_list = open(path_list, "r").read().split("\n")[:-1]


# Initial lists
atoms, bin_sizes, bin_nums, euclidean_distances, filtration_distances = [], [], [], [], []
num_folds, seeds, num_trees, max_depths, subsamples, colsamples = [], [], [], [], [], []
mse_train, mse_test, pcc_train, pcc_test = [], [], [], []

for file in input_files:

    # Read the input data
    path_files = os.path.join(main_folder, subfolder_input, file)
    data = pd.read_csv(path_files)

    # Extract the train set
    criteria = data.iloc[:, -1].str.strip(" ").isin(train_list)
    train_set = data[criteria]
    train_set = train_set.iloc[:, :-1]
    train_set_x, train_set_y = train_set.iloc[:,:-1], train_set.iloc[:, -1]
    train_set_x.reset_index(inplace = True, drop = True)
    train_set_y.reset_index(inplace = True, drop = True)

    # Extract the test set
    criteria = data.iloc[:, -1].str.strip(" ").isin(test_list)
    test_set = data[criteria]
    test_set = test_set.iloc[:, :-1]
    test_set_x, test_set_y = test_set.iloc[:,:-1], test_set.iloc[:, -1]
    test_set_x.reset_index(inplace = True, drop = True)
    test_set_y.reset_index(inplace = True, drop = True)

    # Dataset details
    file = file.strip(".csv")
    atom, euclidean_distance, filtration_distance, bin_size, bin_num = file.split("_")
    euclidean_distance, filtration_distance = int(euclidean_distance), int(float(filtration_distance))
    bin_size, bin_num = float(bin_size), int(bin_num)
    
    # Train model 
    print(datetime.datetime.now(), "\tXGBoost - Training %s" %(file))

    if remove_correlated_features:
        # Identify highly correlated columns in train set
        ## Create correlation matrix
        corr_matrix = train_set_x.corr().abs()

        ## Select upper triangle of correlation matrix
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k = 1).astype(np.bool))
        
        ## Find index of feature columns with correlation greater than 0.95
        to_drop = [col for col in upper.columns if any(upper[col] > 0.95)]

        # Remove highly correlated columns in train and cv set based on train set above
        train_set_x.drop(train_set_x[to_drop], axis = 1, inplace = True)
        test_set_x.drop(test_set_x[to_drop], axis = 1, inplace = True)

    print("Training")
    # Training
    clf = XGBRegressor(max_depth = max_depth, n_estimators = num_tree, random_state = seed,
                        subsample = subsample, colsample_bytree = colsample)
    clf.fit(train_set_x, train_set_y)
    y_pred = clf.predict(train_set_x)
    mse_train.append(mean_squared_error(train_set_y, y_pred))
    pcc_train.append(pearsonr(y_pred, train_set_y)[0])

    print("Predicting")
    # Predict values
    y_pred = clf.predict(test_set_x)
    print("Calculating MSE and PCC")
    # Calculate MSE and PCC
    mse_test.append(mean_squared_error(test_set_y, y_pred))
    pcc_test.append(pearsonr(y_pred, test_set_y)[0])

    # Save the trained model
    if save_model:
        model_filename = os.path.join(main_folder, subfolder_output_model, file + "_model.sav")
        pickle.dump(clf, open(model_filename, "wb"))

    # Print the latest results
    print("MSE train: %.5f \tPCC train: %.3f \nMSE test: %.5f \tPCC test: %.3f\n" %(mse_train[-1], pcc_train[-1], mse_test[-1], pcc_test[-1]))

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



    print("Predicting for each chain in test set")
    mse_individual_test, pcc_individual_test = [], []
    for chain in test_list:
        criteria = data.iloc[:, -1].str.strip(" ") == chain
        individual_test_data = data[criteria]
        individual_test_data = individual_test_data.iloc[:, :-1]
        individual_test_data_x, individual_test_data_y = individual_test_data.iloc[:, :-1], individual_test_data.iloc[:, -1]

        if remove_correlated_features:
            individual_test_data_x.drop(individual_test_data_x[to_drop], axis = 1, inplace = True)

        individual_test_data_pred = clf.predict(individual_test_data_x)
        mse_individual_test.append(mean_squared_error(individual_test_data_y, individual_test_data_pred))
        pcc_individual_test.append(pearsonr(individual_test_data_pred, individual_test_data_y)[0])


    print("Predicting for each chain in train set")
    mse_individual_train, pcc_individual_train = [], []
    for chain in train_list:
        criteria = data.iloc[:, -1].str.strip(" ") == chain
        individual_train_data = data[criteria]
        individual_train_data = individual_train_data.iloc[:, :-1]
        individual_train_data_x, individual_train_data_y = individual_train_data.iloc[:, :-1], individual_train_data.iloc[:, -1]

        if remove_correlated_features:
            individual_train_data_x.drop(individual_train_data_x[to_drop], axis = 1, inplace = True)

        individual_train_data_pred = clf.predict(individual_train_data_x)
        mse_individual_train.append(mean_squared_error(individual_train_data_y, individual_train_data_pred))
        pcc_individual_train.append(pearsonr(individual_train_data_pred, individual_train_data_y)[0])




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
    "MSE_test":mse_test,
    "PCC_train": pcc_train,
    "PCC_test":pcc_test
    }) 

test_results = pd.DataFrame({
    "Chain":test_list,
    "MSE":mse_individual_test,
    "PCC":pcc_individual_test
    })

train_results = pd.DataFrame({
    "Chain":train_list,
    "MSE":mse_individual_train,
    "PCC":pcc_individual_train
    })

# Save dataframe to excel
print("Saving results")
output_results.to_excel(writer, sheet_name = "Results", index = False)
test_results.to_excel(writer, sheet_name = "Test results", index = False)
train_results.to_excel(writer, sheet_name = "Train results", index = False)
writer.save()

print("Completed at", datetime.datetime.now())


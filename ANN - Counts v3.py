"""
Purpose of this script is to train the model with kfold cross validation.
The evaluation metric is MSE and Pearson correlation.

Input: CSV file, train_list
Output: sav file (saved model), CSV file (PCC and MSE)

"""
import time
#time.sleep(15/60*60*60)

import pandas as pd
import pickle
import os
import datetime
from random import Random

from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.layers import LeakyReLU
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

# Variables
num_fold = 5
seed = 50
range_value = 0.15

## Hidden layers parameter
#hidden_layer = 3
#hidden_unit = "75%" #512 # No of units in each hidden units
hidden_activation = "sigmoid" # Activation function for each hidden layer
output_activation = LeakyReLU(alpha = 1.0)
output_activation_1 = "Leaky ReLU"
dropout_rate = 0.2

hidden_layer_list = [2,3,4]
epoch_list = [10, 15, 20]

## Learning rate parameters
optimizer = "Adam"
loss_function = "mean_squared_error"
#epoch = 10
batch_size = None # Chose a value such as the remainder of len(train)%batch_size == 0 else leave as None
save_model = False
letter_list = ["C1-P"]

for letter in letter_list:

    #main_folder = r'C:\Users\User\Desktop\B factor'
    main_folder = r'D:\PHML B factor estimation\02 Project'
    subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s' %letter
    subfolder_list = r'00 Basic information'
    subfolder_output_model = r'03 ML models\05 ANN\Trained models - Counts'
    subfolder_output_results = r'03 ML models\05 ANN\Results - Counts'

    writer = pd.ExcelWriter(os.path.join(main_folder, subfolder_output_results, 
                            str(datetime.datetime.today().date()) + " ANN results - Counts %s 3 - %s.xlsx" %(range_value, letter)))

    # List of input files
    path_input = os.path.join(main_folder, subfolder_input)
    input_files = [f for f in os.listdir(path_input) if f.endswith(".csv")]
    input_files = [f for f in input_files if "_0.5_" in f or "_1.0_" in f or "_1.5_" in f ]
    input_files = [f for f in input_files if f.startswith("%s_35_35" %letter) or f.startswith("%s_40_40" %letter) or f.startswith("%s_45_45" %letter)]


    # Initial lists
    atoms, bin_sizes, bin_nums, euclidean_distances, filtration_distances = [], [], [], [], []
    num_folds, seeds = [], []
    hidden_layers, hidden_units, hidden_activations, output_activations  = [], [], [], []
    optimizers, epochs, batch_sizes, dropout_rates = [], [], [], []
    mse_train, mse_cv, pcc_train, pcc_cv = [], [], [], []

    for file in input_files:

        _, euclidean_distance, filtration_distance, _, _ = file.split("_")
        euclidean_distance, filtration_distance = int(euclidean_distance), int(float(filtration_distance))

        if filtration_distance <= euclidean_distance:

            # Read train list
            path_list = os.path.join(main_folder, subfolder_list, r'Train_list')
            train_list = open(path_list, "r").read().split("\n")[:-1]

            # Read the input data
            path_files = os.path.join(path_input, file)
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
            
            for hidden_layer in hidden_layer_list:
                for epoch in epoch_list:
                        # Train model 
                        for i in range(num_fold):
                            
                            print(datetime.datetime.now(), "\tANN - Training %s \tkFold = %d" %(file, i+1))

                            # Extract train and cv set
                            train_x, train_y =  data_x.loc[train_index[i]], data_y[train_index[i]]
                            cv_x,  cv_y =       data_x.loc[cv_index[i]],    data_y[cv_index[i]]
                            
                            # Normalization
                            ## Normalization parameters
                            mean_train = train_x.mean(axis = 0)
                            std_train = train_x.std(axis = 0).replace(0,1) # Replace 0 with 1 to prevent nan values

                            ## Normalizing data
                            train_x = (train_x - mean_train)/std_train
                            cv_x = (cv_x - mean_train)/std_train

                            # No of hidden units == 0.75 of input nodes
                            input_nodes = train_x.shape[1]
                            hidden_unit = round(0.75*input_nodes)

                            print("Layers: %s \tUnits: %s \tBatch: %s \tEpoch: %s \tDropout: %s" %(hidden_layer, hidden_unit, batch_size, epoch, dropout_rate))

                            # Create ANN model
                            clf = Sequential()
                            clf.add(Dense(hidden_unit, input_dim = input_nodes, batch_size = batch_size))
                            clf.add(Activation(hidden_activation))
                            clf.add(Dropout(rate = dropout_rate))

                            ## Add in the 2nd hidden layers onwards
                            if hidden_layer >= 2:
                                for j in range(2, hidden_layer + 1):
                                    clf.add(Dense(hidden_unit))
                                    clf.add(Activation(hidden_activation))
                                    clf.add(Dropout(dropout_rate))
                            
                            ## Output layer
                            clf.add(Dense(1))
                            clf.add(Activation(output_activation))
                            
                            ## Compile the learning rate
                            clf.compile(optimizer = optimizer, loss = loss_function)

                            # Training
                            clf.fit(train_x, train_y, batch_size = batch_size, epochs = epoch, verbose = 0)
                            y_pred = clf.predict(train_x)
                            y_pred = y_pred.reshape(len(y_pred))
                            mse_train.append(mean_squared_error(train_y, y_pred))
                            pcc_train.append(pearsonr(y_pred, train_y)[0])

                            # Predict values
                            y_pred = clf.predict(cv_x)
                            y_pred = y_pred.reshape(len(y_pred))

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
                        hidden_layers.extend([hidden_layer]*num_fold)
                        hidden_units.extend([hidden_unit]*num_fold)
                        hidden_activations.extend([hidden_activation]*num_fold)
                        output_activations.extend([output_activation_1]*num_fold)
                        
                        optimizers.extend([optimizer]*num_fold)
                        epochs.extend([epoch]*num_fold)
                        batch_sizes.extend([batch_size]*num_fold)
                        dropout_rates.extend([dropout_rate]*num_fold)

    # Create dataframe
    output_results = pd.DataFrame({
        "Atoms":atoms,
        "Euclidean distance":euclidean_distances,
        "Filtration distance":filtration_distances,
        "Bin size":bin_sizes,
        
        "Bin num":bin_nums,
        "Seed":seeds,
        "Hidden layers":hidden_layers,
        "Hidden units":hidden_units,
        "Hidden activation":hidden_activations,
        "Output activation":output_activations,
        
        "Optimizer":optimizers,
        "Epoch":epochs,
        "Batch size":batch_sizes,
        "Drouput rate":dropout_rates,
        
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


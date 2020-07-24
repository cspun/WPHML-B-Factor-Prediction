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

from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.layers import LeakyReLU
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

# Variables
num_fold = 1
seed = 50
hidden_activation = "sigmoid" # Activation function for each hidden layer
output_activation = LeakyReLU(alpha = 1.0)
output_activation_1 = "Leaky ReLU"
dropout_rate = 0.2
hidden_layer = 2
epoch = 10
optimizer = "Adam"
loss_function = "mean_squared_error"
batch_size = None # Chose a value such as the remainder of len(train)%batch_size == 0 else leave as None
save_model = False
letter = "C-N-O-P"

#main_folder = r'C:\Users\User\Desktop\B factor'
main_folder = r'D:\PHML B factor estimation\02 Project'
subfolder_input = r'02 Topology Application\02 Topological features representation\ML inputs - Counts of points at interval\%s' %letter
subfolder_list = r'00 Basic information'
subfolder_output_model = r'03 ML models\05 ANN\Trained models - Counts'
subfolder_output_results = r'03 ML models\05 ANN\Results - Counts'
writer = pd.ExcelWriter(os.path.join(main_folder, subfolder_output_results, 
                        str(datetime.datetime.today().date()) + " ANN results - Test set 1 - %s.xlsx" %letter))

# Read train list
path_list = os.path.join(main_folder, subfolder_list, r'Train_list')
train_list = open(path_list, "r").read().split("\n")[:-1]
path_list = os.path.join(main_folder, subfolder_list, r'Test_list')
test_list = open(path_list, "r").read().split("\n")[:-1]

# List of input files
path_input = os.path.join(main_folder, subfolder_input)
input_files = ["C-N-O-P_45_22.5_0.5_45.csv"]

# Initial lists
atoms, bin_sizes, bin_nums, euclidean_distances, filtration_distances = [], [], [], [], []
num_folds, seeds = [], []
hidden_layers, hidden_units, hidden_activations, output_activations  = [], [], [], []
optimizers, epochs, batch_sizes, dropout_rates = [], [], [], []
mse_train, mse_test, pcc_train, pcc_test = [], [], [], []

for file in input_files:

    # Read the input data
    path_files = path_input + "\\" + file
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

    # Normalization
    ## Normalization parameters
    mean_train = train_set_x.mean(axis = 0)
    std_train = train_set_x.std(axis = 0).replace(0,1) # Replace 0 with 1 to prevent nan values

    ## Normalizing data
    train_set_x = (train_set_x - mean_train)/std_train
    test_set_x = (test_set_x - mean_train)/std_train

    # No of hidden units == 0.75 of input nodes
    input_nodes = train_set_x.shape[1]
    hidden_unit = round(0.75*input_nodes)

    # Train model 
    print(datetime.datetime.now(), "\tANN - Training %s" %(file))
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
    clf.fit(train_set_x, train_set_y, batch_size = batch_size, epochs = epoch, verbose = 0)

    y_pred = clf.predict(train_set_x)
    y_pred = y_pred.reshape(len(y_pred))
    mse_train.append(mean_squared_error(train_set_y, y_pred))
    pcc_train.append(pearsonr(y_pred, train_set_y)[0])

    # Predict values
    y_pred = clf.predict(test_set_x)
    y_pred = y_pred.reshape(len(y_pred))

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
    hidden_layers.extend([hidden_layer]*num_fold)
    hidden_units.extend([hidden_unit]*num_fold)
    hidden_activations.extend([hidden_activation]*num_fold)
    output_activations.extend([output_activation_1]*num_fold)
    
    optimizers.extend([optimizer]*num_fold)
    epochs.extend([epoch]*num_fold)
    batch_sizes.extend([batch_size]*num_fold)
    dropout_rates.extend([dropout_rate]*num_fold)



    print("Predicting for each chain in test set")
    mse_individual_test, pcc_individual_test = [], []
    for chain in test_list:
        criteria = data.iloc[:, -1].str.strip(" ") == chain
        individual_test_data = data[criteria]
        individual_test_data = individual_test_data.iloc[:, :-1]
        individual_test_data_x, individual_test_data_y = individual_test_data.iloc[:, :-1], individual_test_data.iloc[:, -1]

        individual_test_data_pred = clf.predict(individual_test_data_x)
        individual_test_data_pred = individual_test_data_pred.reshape(len(individual_test_data_pred))
        mse_individual_test.append(mean_squared_error(individual_test_data_y, individual_test_data_pred))
        pcc_individual_test.append(pearsonr(individual_test_data_pred, individual_test_data_y)[0])

    print("Predicting for each chain in train set")
    mse_individual_train, pcc_individual_train = [], []
    for chain in train_list:
        criteria = data.iloc[:, -1].str.strip(" ") == chain
        individual_train_data = data[criteria]
        individual_train_data = individual_train_data.iloc[:, :-1]
        individual_train_data_x, individual_train_data_y = individual_train_data.iloc[:, :-1], individual_train_data.iloc[:, -1]

        individual_train_data_pred = clf.predict(individual_train_data_x)
        individual_train_data_pred = individual_train_data_pred.reshape(len(individual_train_data_pred))
        mse_individual_train.append(mean_squared_error(individual_train_data_y, individual_train_data_pred))
        pcc_individual_train.append(pearsonr(individual_train_data_pred, individual_train_data_y)[0])


# Create dataframe
output_results = pd.DataFrame({
    "Atoms":atoms,
    "Euclidean distance":euclidean_distances,
    "Filtration distance":filtration_distances,
    "Bin size":bin_sizes,
    "Bin num":bin_nums,

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
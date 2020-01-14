import keras
import numpy as np
import pandas as pd
from keras import layers, optimizers, losses, metrics
from keras.models import load_model
from keras.models import Model
import sys
import csv
import os.path
from os import path

# Check arguments
if len(sys.argv) != 3:
	sys.exit("Not enough arguments")

if sys.argv[1] != "-i":
	sys.exit("Invalid type of argument")

if not path.exists(sys.argv[2]):
	sys.exit("None existing file")

# Load pretrained model
complete_model = load_model('./WindDenseNN.h5')

# Read data
data = pd.read_csv(sys.argv[2], usecols = [i+1 for i in range(128)], header=None)
timestamps = pd.read_csv(sys.argv[2], usecols = [0], header=None)

# Predict the results of the first layer
intermediate_model = Model(inputs=complete_model.input, outputs =complete_model.layers[0].output)
result = intermediate_model.predict(data)

# Concatenate timestamps with result matrix
csv_contents = np.hstack((timestamps,result))

# Write final results to csv
with open('new_representation.csv', 'w') as file:
	writer = csv.writer(file)
	writer.writerows(csv_contents)
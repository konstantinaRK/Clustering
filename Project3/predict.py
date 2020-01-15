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

# Check if real results exist
if not path.exists("actual.csv"):
	sys.exit("You must provide the file of the actual results first")

# Read actual results
actual_results = pd.read_csv("actual.csv", usecols = [i+1 for i in range(7)], header=None)

# Load pretrained model
model = load_model('./WindDenseNN.h5')

# Read data
data = pd.read_csv(sys.argv[2], usecols = [i+1 for i in range(128)], header=None)
timestamps = pd.read_csv(sys.argv[2], usecols = [0], header=None)

# Predict model
result = model.predict(data)

# for line in range(len(actual_results)):
# 	for column in range(7):
		
# Find the mean absolute error
difference = np.subtract(actual_results, result)
abs_diff = abs(difference)
m_e_a = abs_diff.mean().mean()

# Find the mean absolute percentage error
abs_diff_perc = np.divide(abs_diff, actual_results, out=np.zeros_like(abs_diff), where=actual_results!=0)
m_e_p_a = abs_diff_perc.mean().mean() * 100

# Find the mean square error
square_difference = np.power(difference, 2)
m_s_e = square_difference.mean().mean()

# Concatenate timestamps with result matrix
csv_contents = np.hstack((timestamps,result))

# Write final results to csv
with open('predicted.csv', 'w') as file:
	writer = csv.writer(file)
	writer.writerow(["MAE:" + str(m_e_a), "MAPE:" + str(m_e_p_a) + "%", "MSE:" + str(m_s_e)])
	writer.writerows(csv_contents)
'''
krabbe, sm, cyclo
wt cyclo
gauche, sm, cyclo
npa, sm, cyclo
npc, cyclo
'''
import os
import csv
import sys
from dicom import read_file
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind
from skimage.io import imread
import seaborn as sns
import pandas as pd


cell_lines = ["fabry", "krabbe", "gauche", "npa", "npc", "wt"]
experiments = ["sm", "cyclo"]

def importing_images(dir):
	"""
	Import images in each cell line directory and save the pixel array data as numpy array

	Returns a dictionary with keys = cell_lines/treatment and values = list of arrays
	"""
	img_dict = {}
	for root, dirs, files in os.walk(dir):
		for file in files:
			if ".tif" in file:
				cell_type = os.path.basename(root)

				if cell_type not in img_dict.keys():
					img_dict[cell_type] = []

				img_array = imread(root + "/" + file)

				img_dict[cell_type].append(img_array)

	return img_dict

def lipid_to_cell_ratio(array):
	flat = array.flatten()
	segment_cells = flat[flat > 400]
	segment_lipid = flat[flat > 1000]
	ratio = len(segment_lipid)/len(segment_cells)
	return ratio

def trimmed_histogram(trimmed_dict):
	histogram_dict = {}
	for key, val in trimmed_dict.items():
		histogram_dict[key] = np.histogram(val, bins=100, normed=False, density=False)[0]
	return histogram_dict

def mean_lipid(array):
	flat = array.flatten()
	segment_lipid = flat[flat > 1000]
	return segment_lipid.mean()

def mean_cell_line_pixel_intensities(img_dict, function = lipid_to_cell_ratio):

	mean_std_dict = {}
	for file, arrays in img_dict.items():

		if file not in mean_std_dict.keys():
			mean_std_dict[file] = []

		accum_array = []
		for array in arrays:
			accum_array.append(function(array))

		mean = np.round(np.mean(accum_array), decimals = 6)
		std = np.round(np.std(accum_array), decimals = 6)

		mean_std_dict[file] = (mean, std)

	return mean_std_dict

def mean_image_arrays(img_dict, function = lipid_to_cell_ratio):

	mean_image_array = {}
	for file, arrays in img_dict.items():

		if file not in mean_image_array.keys():
			mean_image_array[file] = []

		accum_array = []
		for array in arrays:
			accum_array.append(function(array))

		mean_image_array[file] = accum_array

	return mean_image_array

def pairwise_t_test(img_dict, function = lipid_to_cell_ratio):
	array_dict = {}

	for file, arrays in img_dict.items():

		if file not in array_dict.keys():
			array_dict[file] = []

		accum_array = []
		for array in arrays:
			accum_array.append(function(array))
			array_dict[file] = accum_array 

	pairwise_dict = {}
	for cell_line1, accum_array1 in array_dict.items(): ### PAIRWISE comparison between cell lines and their expperiments
		for cell_line2, accum_array2 in array_dict.items():
			for cell_line in cell_lines:
				for exper in experiments:
					if (cell_line in cell_line1) and (cell_line in cell_line2): # matching cell lines
						if (exper not in cell_line1) and (exper in cell_line2): # non-matching comparisons
							p_val = ttest_ind(accum_array1, accum_array2, equal_var=False)[1] # only the p-values, leave out the t-statistics
							pairwise_dict[cell_line1 + "_vs_" + cell_line2] = (
								# int(np.mean(accum_array1) * 100), int(np.mean(accum_array2) * 100)), np.round(p_val, decimals = 4))
								(np.mean(accum_array1), np.mean(accum_array2)), np.round(p_val, decimals = 4))

	return pairwise_dict

def merge_dictionaries(dict1, dict2):
	merged_dict = {}

	for key1, val1 in dict1.items():
		if key1 not in merged_dict.keys():
			merged_dict[key1] = []
		merged_dict[key1].extend(val1)

	for key2, val2 in dict2.items():
		if key2 not in merged_dict.keys():
			merged_dict[key2] = []
			merged_dict[key2].extend(val2)

		elif key2 in merged_dict.keys():
			merged_dict[key2].extend(val2)

	return merged_dict


if __name__ == '__main__':

	img_dict = importing_images("/Users/toddhollon/Desktop/")
	mean_val_dict = mean_image_arrays(img_dict, function = lipid_to_cell_ratio)
	pd.DataFrame.from_dict(data=mean_val_dict, orient='index').to_csv('mark_cells_mark4.csv', header=False)
	


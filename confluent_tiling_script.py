

import os
from dicom import read_file
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from skimage.io import imsave, imread
import sys

tile_size = 250
step_size = 100

def cell_background_ratio(img):
	flat = img.flatten()
	segment_cells = flat[flat > 400]
	try:
		ratio = round(len(segment_cells)/len(flat), ndigits=2)
	except ZeroDivisionError:
		ratio = 0
	return ratio

def confluent_tiles(image):

	num_starts = int(image.shape[0]/step_size) - 3 # this will cut off last 50 pixels
	confluent_tiles = []

	for i in range(num_starts):
		for j in range(num_starts):
			y_start = i * step_size
			y_stop = y_start + tile_size
			x_start = j * step_size
			x_stop = x_start + tile_size

			tile = image[y_start:y_stop, x_start:x_stop]

			if cell_background_ratio(tile) > 0.9: # select FOVs with greater than 90% cell confluence
				confluent_tiles.append((i, j, cell_background_ratio(tile)))
			else:
				continue

	return confluent_tiles

def image_tile_dictionary(img_dict):
	tile_dict = {}
	for cell_line, image_list in img_dict.items():
		for image in image_list:

			conf_tiles = confluent_tiles(image) # list of confluent tiles

			if cell_line not in tile_dict.keys():
					tile_dict[cell_line] = []

			for i, j, ratio in conf_tiles:
				xstart = j * step_size  # "column start"
				ystart = i * step_size # "row start"
				xstop = xstart + tile_size
				ystop = ystart + tile_size

				tile = image[ystart:ystop, xstart:xstop]
 
				tile_dict[cell_line].append(tile)

	return tile_dict

def save_tiles(img_dict, dir):
	tile_dict = {}
	for dcm_file, tiles in img_dict.items():
		for tile_number, tile in enumerate(tiles):
			os.chdir(dir)
			imsave(cell_type + "_" + str(tile_number) + ".tif", tile)

if __name__ == '__main__':

	img_dict = {}
	cell_type = "npc_cyclo"

	os.chdir("/Users/toddhollon/Desktop/mark_dcm_files/mark4" + "/" + cell_type)
	file_list = os.listdir()

	img_dict[cell_type] = []

	for file in file_list:
		if "dcm" in file:
			dcm_file = read_file(file)
			img_dict[cell_type].append(dcm_file.pixel_array)

	tile_dict = image_tile_dictionary(img_dict)

	save_tiles(tile_dict, "/Users/toddhollon/Desktop/mark_tiles/mark4" + "/" + cell_type)


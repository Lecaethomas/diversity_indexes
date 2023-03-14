# Author : T.Lecae - INRAE - US-ODR
# Date : 09/03/2023
#######
## Code executed using Gis_1 environment
## The environment may be reproduced using $ conda create --name <env> --file requirements.txt
## requirement.txt created using :: conda list -e > requirements.txt
import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from rasterio.mask import mask
from scipy.stats import entropy
from scipy import ndimage
import os
import json


# Ignore "division by zeros" numpy errors 
np.seterr(all='ignore')

# Define directories
### Get the path to the parent dir
parent_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
#### Data dir
# build the path to the code dir
code_dir = os.path.join(parent_dir, "code")
#### Data dir
# build the path to the data dir
data_dir = os.path.join(parent_dir, "data")
#### Output dir
# build the path to the output dir
output_dir = os.path.join(parent_dir, "output")
# Create the dir if it doesn't exist
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

params_f = os.path.join(code_dir, 'param.json')
f = open(params_f)
params = json.load(f)

raster_name = params["raster_file_name"]
vector_name = params["rpg_complete_name"]
output_name = params["name_for_output_parcels"]


def compute_indices(polygons, src, params, output_dir, output_name):
    def count_landcover_patches(masked):
        # implementation of count_landcover_patches() function goes here
        pass
    
    def calc_ldi(masked):
        # implementation of calc_ldi() function goes here
        pass
    
    def calc_pci(masked):
        # implementation of calc_pci() function goes here
        pass
    
    def calc_contag(masked):
        # implementation of calc_contag() function goes here
        pass
    
    def calc_lsi(masked):
        # implementation of calc_lsi() function goes here
        pass
    
    def calc_pri(masked):
        # implementation of calc_pri() function goes here
        pass
    
    def calc_iji(masked):
        # implementation of calc_iji() function goes here
        pass
    
    for index, polygon in polygons.iterrows():
        # Extract the geometry of the polygon
        geom = polygon.geometry
        # Mask the raster with the polygon geometry
        masked, transform = rasterio.mask.mask(src, [geom], crop=True)
        # Compute basic informations
        cell_values = np.unique(masked.flatten())
        # counts = np.bincount(masked.flatten())
        # print('bincount : ', counts)
        unique, counts = np.unique(masked, return_counts=True)
        total_cells = np.sum(counts)
        patch_numb = count_landcover_patches(masked)
        filter_arr = unique > 0
        class_numb = unique[filter_arr].size

       # Calculate "basics" diversity indexes
        shannon_diversity = -np.sum([(count/total_cells) * np.log2(count/total_cells) for count in counts if count > 0])
        simpson_diversity = 1 - np.sum([(count/total_cells)**2 for count in counts])
        diversity_index = np.sum(np.square(counts/total_cells))
        # Calculate Shannon evenness index
        if (np.log2(class_numb)>1):
            shannon_evenness = shannon_diversity / np.log2(class_numb)
        else : 
            shannon_evenness = '1 seule classe'
        # Calculate dominance index
        dominance_index = np.max(counts) / total_cells
        # Calculate Landscape Division Index (LDI)
        ldi = calc_ldi(masked)
        ## Patch based indexes
        # Patch Cohesion index
        pci = calc_pci(masked)
        # Calc contagion index
        contag = calc_contag(masked)
        # Landscape Shape Index
        lsi = calc_lsi(masked)
        # Patch Richness Index
        pri = calc_pri(masked)
        # Interspersion and Juxtaposition Index
        iji = calc_iji(masked)
        
        ## Store computed data 
        # Store basic information about the mask (cell number, classes number)
        polygons.loc[index, 'cells_n'] = total_cells
        polygons.loc[index, 'class_n'] = class_numb
        polygons.loc[index, 'patch_n'] = patch_numb
        # Store the diversity indices in the new columns
        polygons.loc[index, 'shannon_d'] = shannon_diversity
        polygons.loc[index, 'simpson_d'] = simpson_diversity
        polygons.loc[index, 'class_d'] = diversity_index
        polygons.loc[index, 'shannon_e'] = shannon_evenness
        # Store the dominance index and LDIin the new columns
        polygons.loc[index, 'dominance_i'] = dominance_index
        polygons.loc[index, 'ldi'] = ldi
        # Patch based indexes
        polygons.loc[index,'pci'] = pci
        polygons.loc[index, 'lsi'] = lsi
        polygons.loc[index, 'pri'] = pri
        polygons.loc[index, 'iji'] = iji
        polygons.loc[index, 'contag_i'] = contag
    
    if (bool(params["save_all_fields"]) == True):
        print('true')
        polygons.to_file(os.path.join(output_dir, output_name), driver='ESRI Shapefile', index=True)

    else :
        print('false')
        light_polygons = polygons[[ 'geometry','cells_n', 'patch_n','class_n', 'shannon_d', 'simpson_d', 'class_d', 'shannon_e', 'dominance_i', 'ldi', 'contag_i',  'pci', 'lsi', 'pri', 'iji']] 
            # Save the polygon data to a GeoJSON file
        light_polygons.to_file(os.path.join(output_dir, output_name), driver='ESRI Shapefile', index=True)

    
    # Enlight polygons by keeping only interesting columns 
# Load the raster data and vector data
with rasterio.open(os.path.join(data_dir, raster_name)) as src:
    # print('CRS : ', src.crs, 'EPSG : ', rasterio.crs.CRS.from_string(src.crs))
    # print("connection_meta : " , src.meta)
    raster = src.read(1)
    r_crs = src.crs
    r_epsg = r_crs.to_epsg()
    print('Your raster file is using crs: ', r_crs, 'This coordinate system will be used to reproject vectorial data if not.')
    
    
    polygons = gpd.read_file(os.path.join(data_dir, vector_name))
    if (polygons.crs.to_epsg()!= r_epsg) :
        print('Your polygons EPSG is : ', polygons.crs, '. It will be reprojected using raster\'s EPSG.' )
        polygons = polygons.to_crs(r_epsg)
        #do nothing
    else :
        print('Your polygons EPSG is : ', polygons.crs, '. It won\'t be reprojected.' )

    compute_indices(polygons, src, params, output_dir, output_name)

    
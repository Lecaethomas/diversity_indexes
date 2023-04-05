import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from rasterio.mask import mask
from scipy.stats import entropy
from scipy import ndimage
import os
import json
# Local Imports
from indexes_func import *

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
with open(params_f, 'r') as f:
    params = json.load(f)

    if params['save_all_fields'].lower() == 'true':
        save_all_fields = True
    else:
        save_all_fields = False
    raster_name = params["raster_file_name"]
    vector_name = params["rpg_complete_name"]
    output_name = params["name_for_output_parcels"]

def _compute_indexes(polygons, src, output_dir, output_name):
    for index, polygon in polygons.iterrows():
        # Extract the geometry of the polygon
        geom = polygon.geometry
        # Mask the raster with the polygon geometry
        masked, transform = rasterio.mask.mask(src, [geom], crop=True)
        # Compute basic informations
        cell_values = np.unique(masked.flatten())
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
        
        ##### EDGES #####
        # compute the edges's geometry
        edges = calc_edges(polygon) 
        # simple diversity index under edges
        edges_d_index = calc_edge_diversity_index(edges, src)
        # this function returns two variables edges mode (most common lc class) and the number of classes
        edges_mode, edges_class_numb = calc_most_common_landcover_class(edges, src)#[0]
        
        #edges_class_numb = calc_most_common_landcover_class(edges, src)#[1]
        
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
        # Store the edge-diversity index in the new column
        polygons.loc[index, 'e_div_i'] = edges_d_index
        polygons.loc[index, 'e_mode'] = edges_mode
        polygons.loc[index, 'e_class_n'] = edges_class_numb 
        # Store difference between edge's number of classes and polygon total number of classes
        polygons.loc[index, 'e_p_class_n']= class_numb - edges_class_numb
        
    
    
    # Depending on the parameters file we keep all the fields or we only keep the original geometry 
    if save_all_fields:   
        polygons.to_file(os.path.join(output_dir, output_name), driver='ESRI Shapefile', index=True)

    else :
        print('false')
        light_polygons = polygons[[ 'geometry','cells_n', 'patch_n','class_n', 'shannon_d', 'simpson_d', 'class_d', 'shannon_e', 'dominance_i', 'ldi', 'contag_i',  'pci', 'lsi', 'pri', 'iji', 'e_div_i', 'e_mode', 'e_class_n', 'e_p_class_n']] 
            # Save the polygon data to a GeoJSON file
        light_polygons.to_file(os.path.join(output_dir, output_name), driver='ESRI Shapefile', index=True)
    # Enlight polygons by keeping only interesting columns 
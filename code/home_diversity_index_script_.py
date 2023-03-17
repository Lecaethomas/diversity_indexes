# Author : T.Lecae - INRAE - US-ODR
# Date : 09/03/2023
#######
## Code executed using gis_1 environment
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
from multiprocessing import Pool
from indexes_func import *


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
with open(params_f, 'r') as f:
    params = json.load(f)

    if params['save_all_fields'].lower() == 'true':
        save_all_fields = True
    else:
        save_all_fields = False
    raster_name = params["raster_file_name"]
    vector_name = params["rpg_complete_name"]
    output_name = params["name_for_output_parcels"]





def compute_indices(polygon, src, output_dir, output_name):
        print('compute_START')
        # Extract the geometry of the polygon
        geom = polygon.geometry
        
        index = polygon.index
        print('1')
        # Mask the raster with the polygon geometry
        masked, transform = rasterio.mask.mask(src, [geom], crop=True)
        print('1.5')
        # Compute basic informations
        cell_values = np.unique(masked.flatten())
        # counts = np.bincount(masked.flatten())
        # print('bincount : ', counts)
        unique, counts = np.unique(masked, return_counts=True)
        total_cells = np.sum(counts)
        print('1.5')
        patch_numb = count_landcover_patches(masked)
        print('1.6')
        filter_arr = unique > 0
        class_numb = unique[filter_arr].size
        print('2')
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
        
        indices = {
        'cells_n': total_cells,
        'class_n': class_numb,
        'patch_n': patch_numb,
        'shannon_d': shannon_diversity,
        'simpson_d': simpson_diversity,
        'class_d': diversity_index,
        'shannon_e': shannon_evenness,
        'dominance_i': dominance_index,
        'ldi': ldi,
        'pci': pci,
        'lsi': lsi,
        'pri': pri,
        'iji': iji,
        'contag_i': contag,
    }
        return indices

    # ## Here I try to get a boolean from the json parameters file. I couldn't parse it correctly so 
    #     if save_all_fields:
    #         print('true')
    #         polygons.to_file(os.path.join(output_dir, output_name), driver='ESRI Shapefile', index=True)

    #     else :
    #         print('false')
    #         light_polygons = polygons[[ 'geometry','cells_n', 'patch_n','class_n', 'shannon_d', 'simpson_d', 'class_d', 'shannon_e', 'dominance_i', 'ldi', 'contag_i',  'pci', 'lsi', 'pri', 'iji']] 
    #             # Save the polygon data to a GeoJSON file
    #         light_polygons.to_file(os.path.join(output_dir, output_name), driver='ESRI Shapefile', index=True)
    #     # Enlight polygons by keeping only interesting columns 
    #     print('compute_END')
   
def compute_indices_parallel(poly):
    poly_id, polygon = poly
    print('poly : ', polygon.index)
    indices = compute_indices(polygon, src, output_dir, output_name)
    print('output : ', poly_id, indices)
    return poly_id, indices
results_dict = {}

         
# Load the raster data and vector data
with rasterio.open(os.path.join(data_dir, raster_name)) as src:
    #Get infos about projection (in case data don't use the same projection the code will reproject vectorial data so that it matches raster projection)
    raster = src.read(1)
    r_crs = src.crs
    r_epsg = r_crs.to_epsg()
    # print('Your raster file is using crs: ', r_epsg, 'This coordinate system will be used to reproject vectorial data if not.')
    
    
    polygons = gpd.read_file(os.path.join(data_dir, vector_name))
    if (polygons.crs.to_epsg()!= r_epsg) :
        # print('Your polygons EPSG is : ', polygons.crs, '. It will be reprojected using raster\'s EPSG.' )
        polygons = polygons.to_crs(r_epsg)
        #do nothing
    else :
        pass
        # print('Your polygons EPSG is : ', polygons.crs, '. It won\'t be reprojected.' )


    if __name__ == '__main__':
        try:
            p = Pool(processes=5)
            for poly_id, indices in p.imap_unordered(compute_indices_parallel, polygons.iterrows()):
                results_dict[poly_id] = indices
                # do something with poly_id and indices
        except Exception as e:
            print(e)
        finally:
            p.close()
            p.join()
        
        results = polygons.copy()
        for poly_id, indices in results_dict.items():
            print('Merging')
            for key, value in indices.items():
                results.loc[poly_id, key] = value
    # compute_indices(polygons, src, output_dir, output_name)
    #for index, polygon in polygons.iterrows()
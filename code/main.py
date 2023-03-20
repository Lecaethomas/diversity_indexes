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
# local import
from indexes_func import *
import compute_indexes 



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


def main():
    with open(params_f, 'r') as f:
        params = json.load(f)

        if params['save_all_fields'].lower() == 'true':
            save_all_fields = True
        else:
            save_all_fields = False
        raster_name = params["raster_file_name"]
        vector_name = params["rpg_complete_name"]
        output_name = params["name_for_output_parcels"]
    # Load the raster data and vector data
    with rasterio.open(os.path.join(data_dir, raster_name)) as src:
        #Get infos about projection (in case data don't use the same projection the code will reproject vectorial data so that it matches raster projection)
        raster = src.read(1)
        r_crs = src.crs
        r_epsg = r_crs.to_epsg()
        print('Your raster file is using crs: ', r_epsg, 'This coordinate system will be used to reproject vectorial data if not.')
        polygons = gpd.read_file(os.path.join(data_dir, vector_name))
        if (polygons.crs.to_epsg()!= r_epsg) :
            print('Your polygons EPSG is : ', polygons.crs, '. It will be reprojected using raster\'s EPSG.' )
            polygons = polygons.to_crs(r_epsg)
            #do nothing
        else :
            print('Your polygons EPSG is : ', polygons.crs, '. It won\'t be reprojected.' )
        
        compute_indexes._compute_indexes(polygons, src, output_dir, output_name)
main()

print( 'The process finished successfully.')
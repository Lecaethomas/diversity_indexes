# Author : T.Lecae - INRAE - US-ODR
# Date : 09/03/2023
#######
## Code executed using Gis_1 environment
## The environment may be reproduced using $ conda create --name <env> --file requirements.txt
## requirement.txt created using :: conda list -e > requirements.txt

#Performances : 
  # using i7 32gb RAM 
  # #over 135k parcels and OSO (10*10m cells) finished in 15min without parallelization
  
import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
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
# build the path to the code dir
code_dir = os.path.join(parent_dir, "code")
# build the path to the data dir
data_dir = os.path.join(parent_dir, "data")
# build the path to the output dir
output_dir = os.path.join(parent_dir, "output")
# Create the output dir if it doesn't exist
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


params_f = os.path.join(code_dir, 'param.json')


def main():
    """main function used to launch the processing. It's kept without parameters for now because it may be
    amended with parameters like an administrative divisions list so that it may be parallelized using multiprocessing
    """    
    with open(params_f, 'r') as f: #Get parameters from params.json
        params = json.load(f)
        raster_name = params["raster_file_name"]
        vector_name = params["rpg_complete_name"]
        output_name = params["name_for_output_parcels"]
    # Load the raster data and vector data
    with rasterio.open("C:\_DATA\OCS\OSO\OCS_2021.tif",mode="r+",
                       nodata= 0) as src: 
        #Get infos about projection (in case data don't use the same projection the code will reproject vectorial data so that it matches raster projection)
        r_crs = src.crs
        r_epsg = r_crs.to_epsg()
        print('UNIT : ',src.crs.linear_units)
        if (src.crs.linear_units.strip in ('metre', 'meter')): #check what's the raster's projection in order to inform user that he must reproject it using a meter-based crs
            print('Your raster input is not using a projected CRS (measured in meters). In order to process correctly reproject it to your local projected CRS.' )
        print('Your raster file is using crs: ', r_epsg, 'This coordinate system will be used to reproject vectorial data if not.')
        polygons = gpd.read_file(os.path.join(data_dir, vector_name))
        
        if (polygons.crs.to_epsg()!= r_epsg) : # reprojecting depending on raster's projection
            print('Your polygons EPSG is : ', polygons.crs, '. It will be reprojected using raster\'s EPSG.' )
            polygons = polygons.to_crs(r_epsg)
            
        else : #do nothing
            print('Your polygons EPSG is : ', polygons.crs, '. It won\'t be reprojected.' )
            
        #Once parameters are defined let's launch it 
        compute_indexes.compute_indexes_(polygons, src,
                                         output_dir, output_name)
        
main()

print( 'The process finished successfully.')
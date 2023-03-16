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


## Definition de quelques fonctions permettant de calculer des indices de diversité. Diversité brute de pixels mais aussi pour certaines prenant en compte les patches. 
## Elles seront appliquées itterativement aux "morceaux" d'OSO découpés par les polygones de parcelles du RPG

def count_landcover_patches(masked):
    """
     @brief Count the patches (defined as contiguous cells having the same digital number) in a landcover raster
     @param raster An array
     @return The number of patches in the raster
    """
    # Get the unique landcover classes in the raster.
    classes = np.unique(masked)
    # Initialize a variable to store the total number of patches.
    patch_count = 0
    # Loop over each landcover class.
    # This function will label the connected regions of the current class.
    for c in classes:
        # Create a binary mask for the current class.
        binary_raster = np.where(masked == c, 1, 0)
        # Label the connected regions in the binary mask.
        labeled_raster, num_labels = ndimage.label(binary_raster)
        # Count the number of patches with contiguous cells of the same class.
        for i in range(1, num_labels + 1):
            patch_mask = labeled_raster == i
            # If patch_mask is not empty it will increment patch count.
            if patch_mask.any():
                patch_count += 1
    return patch_count

#Landscape Division Index
def calc_ldi(masked):
    """
     @brief Calculate the LDI of a set of cells. This is a measure of how much the cell is in the data set but not the area of the cell that is used to calculate the distance between the cell and its neighboring cells
     @param masked A 2D NumPy array with shape ( nb_cells num_cells )
     @return An LDI value between 0 and 1 ( inclusive ). The value is calculated as the standard deviation divided by the mean
    """
    # Calculate the standard deviation of cell values
    std = np.std(masked.flatten())
    # Calculate the mean of cell values
    mean = np.mean(masked.flatten())
    # Calculate the LDI
    ldi = std / mean
    return ldi

## Patch Cohesion Index
def calc_pci(masked):
    """
     @brief Calculate the PCI for a set of images. This is based on the number of patches and perimeter
     @param masked Masked image to calculate the PCI for
     @return pci The PCI for the image as a function of the number of patches and perimeter ( float
    """
    # Calculate the Patch Cohesion Index (PCI)
    # Returns the sum of masked masked values.
    if np.sum(masked) == 0:
        return 0.0

    # Calculate the number of patches
    label_im, nb_labels = ndimage.label(masked)
    # Returns the number of labels.
    if nb_labels == 1:
        return 1.0

    # Calculate the total perimeter
    perimeter = np.sum(ndimage.binary_dilation(masked) & ~masked)

    # Calculate the perimeter for each patch
    patch_perimeter = np.zeros(nb_labels, dtype=np.float64)
    # Calculate patch perimeter for each label.
    for label in range(1, nb_labels + 1):
        mask_label = label_im == label
        patch_perimeter[label - 1] = np.sum(ndimage.binary_dilation(mask_label) & ~mask_label)

    # Calculate the average patch perimeter
    avg_patch_perimeter = np.sum(patch_perimeter) / nb_labels

    # Calculate the PCI
    pci = (perimeter - avg_patch_perimeter) / perimeter

    return pci

## Contagion Index 
def calc_contag(masked):
    """
     @brief Calculate the contagion index of a mask. This is based on the proportion of the patches in the mask divided by the total area
     @param masked The mask to calculate the contagion index for
     @return The contagion index of the mask as a function of the area of the image and the proportion
    """
    # Calculate the Contagion Index (CONTAG)
    # Returns the sum of masked masked values.
    if np.sum(masked) == 0:
        return 0.0

    # Calculate the number of patches
    label_im, nb_labels = ndimage.label(masked)
    # Returns the number of labels.
    if nb_labels == 1:
        return 0.0

    # Calculate the total area
    total_area = masked.size

    # Calculate the area for each patch
    patch_area = np.zeros(nb_labels, dtype=np.float64)
    # Calculate the patch area for each label
    for label in range(1, nb_labels + 1):
        mask_label = label_im == label
        patch_area[label - 1] = np.sum(mask_label)

    # Calculate the proportion of each patch
    patch_prop = patch_area / total_area

    # Calculate the CONTAG
    contag = -np.sum(patch_prop * np.log(patch_prop))

    return contag

## Landscape Shape Index
def calc_lsi(masked):
    """
     @brief Calculate LSI for a landcover raster. This is based on the proportion of pixels in each row and column that have landcover.
     @param lc The landcover raster to calculate lsi for.
     @return The LSI for the raster as a 1 - ( var ( row ) + var ( column ))
    """
    # convert landcover raster to binary (1 = landcover class present, 0 = landcover class absent)
    binary_lc = np.where(masked > 0, 1, 0)
    
    # calculate the proportion of pixels in each row and column that have landcover
    row_prop = np.mean(binary_lc, axis=1)
    col_prop = np.mean(binary_lc, axis=0)
    
    # calculate row and column variance
    row_var = np.var(row_prop)
    col_var = np.var(col_prop)
    
    # calculate LSI
    lsi = 1 - (row_var + col_var)
    
    return lsi

## Patch Richness Index
def calc_pri(masked):
    """
     @brief Calculate pri for landcover. This is based on the number of patches and total area for each class
     @param lc numpy array of landcover classes
     @return pri ( float ) : pri for landcover in range 0.. 1 where 0 is no patches and 1 is
    """
    # get unique landcover classes
    classes = np.unique(masked)
    
    # calculate the number of patches and total area for each class
    num_patches = []
    total_area = []
    # Add the number of patches and labels for each class.
    for c in classes:
        binary_lc = np.where(masked == c, 1, 0)
        num, labels = ndimage.label(binary_lc)
        num_patches.append(num)
        total_area.append(np.sum(binary_lc))
    
    # calculate PRI
    pri = np.sum([num_patches[i] * total_area[i] for i in range(len(classes))]) / np.sum(total_area)
    
    return pri

## Interspersion and Juxtaposition Index
def calc_iji(masked):
    """
     @brief Calculate the number of neighboring pixels with different landcover classes. This is used to calculate the i. i. d.
     @param lc The landcover raster of the same shape as the grid
     @return The number of neighboring pixels with different landcover classes in the grid as a 1x1 np
    """
    # convert landcover raster to binary (1 = landcover class present, 0 = landcover class absent)
    binary_lc = np.where(masked > 0, 1, 0)
    
    # calculate the number of adjacent pixels with different landcover classes
    iji = np.sum(np.abs(np.diff(binary_lc, axis=0))) + np.sum(np.abs(np.diff(binary_lc, axis=1)))
    
    return iji


def compute_indices(polygons, src, output_dir, output_name):
    
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
    ## Here I try to get a boolean from the json parameters file. I couldn't parse it correctly so 
        if save_all_fields:
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
    if __name__ == '__main__':
        try:
            p = Pool(processes=5)
            p.map_async(lambda x: compute_indices(x[1], src, output_dir, output_name), polygons.iterrows())
            print('The process finished successfully.')
            p.close()
        except Exception as e:
            print(e)
            p.terminate()
            print('except')
        finally:
            p.join()
            print('join')
    # compute_indices(polygons, src, output_dir, output_name)
    #for index, polygon in polygons.iterrows()
    

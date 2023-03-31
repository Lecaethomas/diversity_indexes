import rasterio
import numpy as np
from scipy import stats as st
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString
from shapely import ops
from rasterio.mask import mask
from scipy.stats import entropy
from scipy import ndimage
# from collections import Counter

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

def calc_edges_(polygon, src):
    # Extract the edges of the polygon as line segments
    edges = [LineString([pt1, pt2]) for pt1, pt2 in zip(polygon.geometry.exterior.coords[:-1], polygon.geometry.exterior.coords[1:])]
    edges = [LineString([pt1, pt2]) for pt1, pt2 in zip(polygon.geometry.exterior.coords[:-1], polygon.geometry.exterior.coords[1:])]

    # Create a MultiLineString object from the edges
    mls = MultiLineString(edges)

    # Merge the lines into a single linestring
    merged_line = ops.linemerge(mls)
    # Clip the raster to the line segments
    edge_rasters = []
    # for edge in edges:
    try:
            clipped_raster, _ = mask(src, [merged_line], crop=True)
            edge_rasters.append(clipped_raster)
    except ValueError as e:
            print(e)
        
    # Calculate the diversity index for each edge raster
    edge_diversity_indices = []
    pol_unique_values = []
    for edge_raster in edge_rasters:
        unique_values = np.unique(edge_raster)
        total_cells = np.size(edge_raster)
        diversity_index = len(unique_values) / total_cells
        edge_diversity_indices.append(diversity_index)
    
    # Calculate the average edge-diversity index
    if edge_diversity_indices:
        edge_diversity_index = np.mean(edge_diversity_indices)
        pol_unique_values = np.sum(unique_values)
    else:
        edge_diversity_index = np.nan
        pol_unique_values = np.nan
    
    return edge_diversity_index, pol_unique_values

def calc_edge_diversity_index(polygon, src):
    # Extract the edges of the polygon as a single LineString
    edges = LineString(polygon.geometry.exterior.coords)
    # Clip the raster to the LineString
    try:
        clipped_raster, _ = mask(src, [edges], crop=True)
        # Calculate the diversity index for the clipped raster
        unique_values = np.unique(clipped_raster)
        total_cells = np.size(clipped_raster)
        diversity_index = len(unique_values) / total_cells
        return diversity_index
    except ValueError as e:
        print(f"Error: {e}. Skipping polygon...")
        return np.nan

def calc_most_common_landcover_class(polygon, src):
    # Extract the edges of the polygon as a single LineString
    edges = LineString(polygon.geometry.exterior.coords)
    # Clip the raster to the LineString
    try:
        clipped_raster, _ = mask(src, [edges], crop=True)
        # Calculate the mode of the clipped raster. The numpy array includes the rasterio no-data (0) so we still need to remove it.
        vals,counts = np.unique(clipped_raster[clipped_raster != src.nodata], return_counts=True)
        # find the index of the maximum count
        mode_index = np.argmax(counts)
        filter_arr = vals > 0
        edge_class_numb = vals[filter_arr].size
        # the mode is the value at the mode index
        mode = vals[mode_index]
        return mode, edge_class_numb
    
    except ValueError as e:
        print(f"Error: {e}. Skipping polygon...")
        return np.nan
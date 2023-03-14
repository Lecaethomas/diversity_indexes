# Diversity Indexes
This repo store codes used to process raster and vectorial data aiming at computing/testing few diversity indices. 
## Synthesis
### What for?
These focus on showing diversity of classes of landcover in a raster layer. That diversity may speak out interms of raw number of classes, but also may be thought as diversity of spatial distribution (number of patches).
### For who?
It may be useful for landscape ecologists or any other kind of domains aiming at assessing the diversity of raster cells contained in polygon data.
### How?
Despite the code had been parameterized, it still needs installs and little skills in python/tech (using conda for install, editing .json file...)
### Appended fields :
**Between parenthesis : the field name used**
#### Basics infos (computed for each polygons individually): 
- total cells number (*cells_n*)
- total patch number (patch defined as a contiguous set of cells belonging to the same class) (*patchs_n*)
- class number (*class_n*)
#### Diversity indexes
For more infos about diversity indexes and their meaning look at Fragstat docs
- Class diversity (*class_d*)
- Shannon diversity (*shannon_d*)
- Simpson diversity (*simpson_d*)
- Jaccard diversity (*jaccard_d*)
- Shannon eveness (*shannon_e*)
- Dominance index (*dominance_i*)
- Contagion index (*contag_i*)
- Landscape Division Index (LDI) (*ldi*)
- Landscape Shape Index (LSI) (*lsi*)
- Patch Richness Index (PRI) (*pri*)
- Interspersion end Juxtaposition Index (IJI) (*iji*)
## Usage
This tool take two inputs : 
- one one-band raster layer made of integers
- one vectorial layer (format readen by geopandas/fiona see below) made of polygons 
### Installs : 
- DL this repo
- Using console get to the repo directory and type : 

`conda create --name <env> --file requirements.txt`
### Parameters :
The file param.json allows you to precise the name of your inputs wich nevertheless must be stored in the "data" directory.
You can precise the name of your output. If the directory "output" doesn't exists it will be created.

##### List of supported vector format:
`'DXF',
 'CSV',
 'OpenFileGDB',
 'ESRIJSON',
 'ESRI Shapefile',
 'FlatGeobuf',
 'GeoJSON',
 'GeoJSONSeq',
 'GPKG',
 'GML',
 'OGR_GMT',
 'GPX',
 'Idrisi',
 'MapInfo File',
 'DGN',
 'PCIDSK',
 'OGR_PDS',
 'S57',
 'SQLite',
 'TopoJSON'` 
# Diversity Indexes
This repo store codes used to process raster and vectorial data aiming at computing/testing few diversity indices. 
## Synthesis
### What for?
These focus on showing diversity of classes of landcover in a raster layer. That diversity may speak out interms of raw number of classes, but also may be thought as diversity of spatial distribution (number of patches).
### For who?
It may be useful for landscape ecologists or any other kind of domains aiming at assessing the diversity of raster cells contained in polygon data.
### How?
Despite the code had been parameterized, it still needs installs and little skills in coding (using conda for install, editing .json file...)
### List of indexes computed :

## In depths
This tool take two inputs : 
- one one-band raster layer made of cells with integers
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
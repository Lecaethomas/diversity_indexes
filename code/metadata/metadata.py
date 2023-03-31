## A tool used to read a shapefile, and create a .json that store 
# metadata (field name, field type, create description for each field,
# that should be filled manually).
# This .json file may then be parsed by metadata_to_pdf.pdf to create a 
# .pdf file containing these info 
# There's the possibility to write the bbox as geojson or to create an image 
# showing the layer with a mapbox static layer as background but it 
# is not functionnal
import geopandas as gpd
import pandas as pd
import json
from shapely.geometry import shape, box
from io import BytesIO
import requests
from PIL import Image


def shapefile_to_json(shapefile_path, show_layer=False, bbox_file=True):
    """
     @brief Converts a shapefile to GeoJSON. This is a utility function to convert a GDS - style shapefile to a GeoJSON - style feature collection that can be used in an OCR query
     @param shapefile_path Path to the shapefile to convert
     @param show_layer Whether or not to show the layer
     @param bbox_file Whether or not to output the bounding box as a separate file
     @return GeoJSON - compatible feature collection of the shapefile and the bounding box in EPSG : 4
    """
    print('shp_to_json start')
    gdf = gpd.read_file(shapefile_path)
    feature_collection = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": None,
                "extent" : None,
                "properties": {
                    field: {"type": gdf.dtypes[field].name} for field in gdf.columns
                }
            }
        ]
    }
    # Show the layer of the shapefile.
    if show_layer :
        print('show_layer')
        
        # Get the bounding box of the shapefile
        bbox1 = gdf.to_crs('EPSG:4326').total_bounds
        bbox = box(*bbox1)
        geojson_bbox = {
        "type": "Feature",
        "geometry": bbox.__geo_interface__,
        "properties": {}
        }
        bbox_2 = []
        # Add bbox to bbox1.
        for p in bbox1 :
            bbox_2.append(p)
        print(bbox_2)
        # Output the bounding box as a separate file if bbox_file is provided
        # Display the bounding box on a map or a map
        if bbox_file:
            print('bbox_file')
            with open('bbox_file_.geojson', "w") as f:
                f.write(json.dumps(bbox.__geo_interface__))
                
        # Otherwise, display the bounding box on a map
        else:
            print('map_snapshot')
            url = "https://api.mapbox.com/styles/v1/mapbox/streets-v11/static/geojson("
            url += str(geojson_bbox)
            url += ")/auto/500x300?access_token="
            url += "pk.eyJ1IjoidGxlY2FlIiwiYSI6ImNreHJqOGFwaTAzN3Ayd281dTBmb3VzYTYifQ.8fwOWabWWbcfcUxi1rIxAQ"
            print(url)
            response = requests.get(url)
            print("BBOX : ",json.dumps(bbox.__geo_interface__),"RESPONSE : ", response)
            image = Image.open(BytesIO(response.content))
            image.format = "PNG"
            image.show()
    else : 
        gdf.to_crs('EPSG:4326')
    return feature_collection


shapefile_path = "C:/_DEV/test/landscape_metrics/landscape_metrics_home/output/parcels_w_indexes.shp"
json_output = shapefile_to_json(shapefile_path,show_layer= True, bbox_file=False)

# Convert json_output to a Python dictionary before writing it to

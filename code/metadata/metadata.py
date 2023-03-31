import geopandas as gpd
import json
from shapely.geometry import shape, box
from io import BytesIO
import requests
from PIL import Image


def shapefile_to_json(shapefile_path, show_layer=True, bbox_file=None):
    gdf = gpd.read_file(shapefile_path)
    feature_collection = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": None,
                "properties": {
                    field: {"type": gdf.dtypes[field].name} for field in gdf.columns
                }
            }
        ]
    }
    if show_layer:
        # Get the bounding box of the shapefile
        bbox = gdf.total_bounds
        bbox = box(*bbox)
        # Output the bounding box as a separate file if bbox_file is provided
        if bbox_file:
            with open(bbox_file, "w") as f:
                f.write(json.dumps(bbox.__geo_interface__))
        # Otherwise, display the bounding box on a map
        else:
            url = "https://api.mapbox.com/styles/v1/mapbox/streets-v11/static/geojson("
            url += bbox
            url += ")/auto/500x300?access_token="
            url += "pk.eyJ1IjoidGxlY2FlIiwiYSI6ImNreHJqOGFwaTAzN3Ayd281dTBmb3VzYTYifQ.8fwOWabWWbcfcUxi1rIxAQ"
            response = requests.get(url)
            print("BBOX : ",json.dumps(bbox.__geo_interface__),"RESPONSE : ", response)
            image = Image.open(BytesIO(response.content))
            image.format = "PNG"
            image.show()
    return json.dumps(feature_collection)



shapefile_path = "C:/_DEV/test/landscape_metrics/landscape_metrics_home/output/parcels_w_indexes.shp"
json_output = shapefile_to_json(shapefile_path, show_layer=True)
print(json_output)
with open("metadata.json", "w") as f:
    json.dump(json_output, f)
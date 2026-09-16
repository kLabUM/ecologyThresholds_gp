import csv
import json

esri_json = {
    "geometryType": "esriGeometryPoint",
    "spatialReference": {"wkid": 4326},
    "fields": [
        {"name": "Site", "type": "esriFieldTypeString"},
        {"name": "Label", "type": "esriFieldTypeString"},
    ],
    "features": [],
}

with open("CrossSiteVisual/siteLatLong.csv", mode="r", encoding="utf-8-sig") as file:
    reader = csv.DictReader(file)
    for row in reader:
        feature = {
            "geometry": {
                "x": float(row["Long"].strip()),
                "y": float(row["Lat"].strip()),
            },
            "attributes": {
                "Site": row["Site"].strip(),
                "Site_Label": row["Label"].strip(),
            },
        }
        esri_json["features"].append(feature)

with open("CrossSiteVisual/sensor_sites_esri.json", "w", encoding="utf-8") as file:
    json.dump(esri_json, file, indent=2)
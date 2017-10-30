from spacenetutilities.labeltools import coreLabelTools
import json
import glob
import argparse
from datetime import datetime
import os

def modifyTimeField(geoJson, geoJsonNew, featureItemsToAdd=['ingest_tim', 'edit_date'], featureKeyListToRemove=[]):
    now = datetime.today()
    with open(geoJson) as json_data:
        d = json.load(json_data)


    featureList = d['features']
    newFeatureList = []
    for feature in featureList:
        tmpFeature = dict(feature)
        for featureKey in featureKeyListToRemove:
            if featureKey in tmpFeature['properties']:
                del tmpFeature['properties'][featureKey]
        for featureKey in featureItemsToAdd:
            if not (featureKey in tmpFeature['properties']):
                print('inserting missing field')
                print(now.isoformat())
                tmpFeature['properties'][featureKey] = now.isoformat()
            else:
                if not tmpFeature['properties'][featureKey]:
                    print('filling empty field')

                    tmpFeature['properties'][featureKey] = now.isoformat()

        newFeatureList.append(tmpFeature)

    d['features']=newFeatureList

    if os.path.exists(geoJsonNew):
        os.remove(geoJsonNew)
    with open(geoJsonNew, 'w') as json_data:
        json.dump(d, json_data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Inserts Missing Time Fields in geojson which causes errors with certain tools')
    parser.add_argument("geojson", help="geoJsonToModify")
    parser.add_argument("--geojsonIDModifierer", help="modifier to extension of geojson",
                        default='noid')


    args = parser.parse_args()
    geojson = args.geojson
    geojsonNew = geojson.replace('.geojson', args.geojsonIDModifierer+".geojson")
    modifyTimeField(args.geojson,
                                           args.geo)




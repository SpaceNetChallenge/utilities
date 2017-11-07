from spacenetutilities.labeltools import coreLabelTools
import json
import glob
import argparse
from datetime import datetime
import os


def modifyTimeField(geoJson, geoJsonNew, featureItemsToAdd=['ingest_tim', 'ingest_time', 'edit_date'], featureKeyListToRemove=[]):
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


def removeIdFieldFromJsonEntries(geoJson, geoJsonNew, featureKeyListToRemove=['Id', 'id'], featureItemsToAdd={}):
    with open(geoJson) as json_data:
        d = json.load(json_data)


    featureList = d['features']
    newFeatureList = []
    for feature in featureList:
        tmpFeature = dict(feature)
        for featureKey in featureKeyListToRemove:
            if featureKey in tmpFeature['properties']:
                del tmpFeature['properties'][featureKey]

        tmpFeature.update(featureItemsToAdd)
        newFeatureList.append(tmpFeature)

    d['features']=newFeatureList

    if os.path.exists(geoJsonNew):
        os.remove(geoJsonNew)
    with open(geoJsonNew, 'w') as json_data:
        json.dump(d, json_data)


def removeIdinGeoJSONFolder(folder, modifier='noid'):

    geoJsonList = glob.glob(os.path.join(folder, '*.geojson'))

    for geojsonName in geoJsonList:
        removeIdFieldFromJsonEntries(geojsonName, geojsonName.replace('.geojson', '{}.geojson'.format(modifier)))


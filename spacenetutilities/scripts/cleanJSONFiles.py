import json
import glob
import os




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

    with open(geoJsonNew, 'w') as json_data:
        json.dump(d, json_data)


def removeIdinGeoJSONFolder(folder, modifier='noid'):

    geoJsonList = glob.glob(os.path.join(folder, '*.geojson'))

    for geojsonName in geoJsonList:
        removeIdFieldFromJsonEntries(geojsonName, geojsonName.replace('.geojson', '{}.geojson'.format(modifier)))




if __name__ == '__main__':


    folderList = ['AOI_2_Vegas/geojson/buildings',
                  'AOI_3_Paris/geojson/buildings',
                  'AOI_4_Shanghai/geojson/buildings',
                  'AOI_5_Khartoum/geojson/buildings']

    for folder in folderList:
        removeIdinGeoJSONFolder(folder)














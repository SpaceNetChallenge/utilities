import json
import glob
import geopandas as gpd
import os
import pandas as pd




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


    fillOpacity=0.3
    AOI_List = [
        {'AOI_Name':'AOI_1_Rio', 'aoi_color':'#800000'},
        {'AOI_Name': 'AOI_2_Vegas', 'aoi_color': '#FF0000'},
        {'AOI_Name': 'AOI_3_Paris', 'aoi_color': '#FFA500'},
        {'AOI_Name': 'AOI_4_Shanghai', 'aoi_color': '#FFFF00'},
        {'AOI_Name': 'AOI_5_Khartoum', 'aoi_color': '#808000'}
        ]




    baseLocation='./datasets/'

    #spaceNetAWSLocation = os.path.join('spacenet-dataset', AOI_Name)


    panLocation = '/srcData/rasterData/PAN/'
    mulLocation = '/srcData/rasterData/MUL/'
    rgbLocation = '/srcData/rasterData/RGB-PanSharpen/'
    mulLocation = '/srcData/rasterData/MUL-PanSharpen/'


    buildingLocation = './processedData/geojson/buildings'


    removeId = False
    if removeId:
        for AOI_DICT in AOI_List:
            folder = os.path.join(baseLocation, AOI_DICT['AOI_Name'], buildingLocation)
            removeIdinGeoJSONFolder(folder)



    totalGDFList = []
    for AOI_DICT in AOI_List:
        AOI_Name = AOI_DICT['AOI_Name']
        srcGeoJson = os.path.join(baseLocation,AOI_Name,'{}_SrcTindex.geojson'.format(AOI_Name))
        srcGDF = gpd.read_file(srcGeoJson)
        srcGDF['AOI']=AOI_Name

        spaceNetAWSLocation = os.path.join('spacenet-dataset', AOI_Name, 'srcData')
        srcGDF['AWSLocation'] = [os.path.join(spaceNetAWSLocation, x) for x in srcGDF['location']]
        srcGDF['fill-opacity']=fillOpacity
        srcGDF['fill']=AOI_DICT['aoi_color']
        srcGDF.to_file(srcGeoJson.replace('.geojson', 'ex.geojson'), driver='GeoJSON')
        #TODO finish merge of geoDataFrame
        #TODO Implement visual geojson features for github rendering
        totalGDFList.append(srcGDF)

    totalGDF = gpd.GeoDataFrame(pd.concat(totalGDFList))
    totalGDF.crs = srcGDF.crs
    totalGeoJsonName = os.path.join(baseLocation, 'SpaceNetSummaryTindex.geojson')
    if os.path.exists(totalGeoJsonName):
        os.remove(totalGeoJsonName)


    totalGDF.to_file(os.path.join(baseLocation, 'SpaceNetSummaryTindex.geojson'), driver='GeoJSON')





from spaceNetUtilities import geoTools as gT
import geopandas as gpd
import glob
import pandas as pd
import os
import csv
EPSGLIST = {'Rio': 32723,
            'AOI_2_Vegas_Train': 32612,
            'AOI_3_Paris_Train': 32632,
            'AOI_4_Shanghai_Train':32608,
            'AOI_5_Khartoum_Train': 32635
             }






row = ['ImageId', 'Count'
       'Area_Mean', 'Area_std', 'Area_min', 'Area_25per', 'Area_50per', 'Area_75per', 'Area_max',
       'PartialDec_Mean', 'PartialDec_std', 'PartialDec_min', 'PartialDec_25per', 'PartialDec_50per', 'PartialDec_75per', 'PartialDec_max',
       'partialBuilding_Mean', 'partialBuilding_std', 'partialBuilding_min', 'partialBuilding_25per', 'partialBuilding_50per', 'partialBuilding_75per', 'partialBuilding_max',
       ]
def processGeoJson(geoJsonFileName, pixelSizeM = 0.3, pixelSizeDeg=0.000002700000000):
    df = gpd.read_file(geoJsonFileName)
    ImageId = os.path.basename(geoJsonFileName)
    ImageId = os.path.splitext(ImageId)[0]
    ImageId = ImageId.split('_', 1)[1]
    dataRow = [ImageId]

    if not df.size == 0:
        dfArea = (df.area * (pixelSizeM ** 2) / (pixelSizeDeg ** 2))
        dataRow.append(dfArea.sum())
        dataRow.extend(dfArea.describe().values)
        dataRow.extend(df['partialDec'].describe().values[1:])
        dataRow.extend(df['partialBuilding'].describe().values[1:])

    return dataRow

if __name__ == "__main__":


    summaryFileLocation = 'imageIDSummaryFile.csv'
    inputDirectory = ''
    directoryList = [['AOI_2_Vegas_Train', 0.000002700000000],
                 ['AOI_3_Paris_Train', 0.000002700000000],
                 ['AOI_4_Shanghai_Train', 0.000002700000000],
                 ['AOI_5_Khartoum_Train', 0.000002700000000]]

    row = ['ImageId',  'Area_Sum', 'Count',
                      'Area_Mean', 'Area_std', 'Area_min', 'Area_25per', 'Area_50per', 'Area_75per', 'Area_max',
           'PartialDec_Mean', 'PartialDec_std', 'PartialDec_min', 'PartialDec_25per', 'PartialDec_50per',
           'PartialDec_75per', 'PartialDec_max',
           'partialBuilding_Mean', 'partialBuilding_std', 'partialBuilding_min', 'partialBuilding_25per',
           'partialBuilding_50per', 'partialBuilding_75per', 'partialBuilding_max',
           ]

    dataRowList = []
    for directory in directoryList:
        fullDirectory = os.path.join(inputDirectory, directory[0], 'geojson', 'buildings', '*.geojson')
        geojsonlist = glob.glob(fullDirectory)

        for geoJsonFileName in geojsonlist:
            print(geoJsonFileName)
            dataRowList.append(processGeoJson(geoJsonFileName, pixelSizeM=0.3, pixelSizeDeg=0.000002700000000))


    with open(summaryFileLocation, 'wb') as csvfile:
        rowwriter = csv.writer(csvfile, delimiter=',')

        rowwriter.writerow(row)

        for dataRow in dataRowList:
            rowwriter.writerow(dataRow)







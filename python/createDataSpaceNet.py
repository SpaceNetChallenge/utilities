import os
import sys
import glob
from osgeo import gdal
import random
sys.path.extend([os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../../utilities_test/python/')])
from spaceNetUtilities import labelTools as lT


dataPath = os.path.dirname(os.path.realpath(__file__))

spaceNetDirectory = '/data/TopCoder_Challenge2/spacenetV2_Train'
listOfAOIs = [subdir for subdir in os.listdir(spaceNetDirectory) if os.path.isdir(os.path.join(spaceNetDirectory, subdir))]

srcImageryDirectory = 'RGB-PanSharpen' # 'PAN', 'MUL, 'MUL-PanSharpen'
geojsonDirectory = 'geojson/buildings/'
annotationsDirectory = 'annotations'


entryList = []

for aoiSubDir in listOfAOIs:
    fullPathSubDir = os.path.join(spaceNetDirectory, aoiSubDir)

    ## Create Annotations directory
    fullPathAnnotationsDirectory = os.path.join(fullPathSubDir, annotationsDirectory)
    if not os.path.exists(fullPathAnnotationsDirectory):
        os.makedirs(fullPathAnnotationsDirectory)

    fullPathImageDirectory = os.path.join(fullPathSubDir, srcImageryDirectory)
    fullPathGeoJsonDirectory = os.path.join(fullPathSubDir, geojsonDirectory)

    listofRaster = sorted(glob.glob(os.path.join(fullPathImageDirectory, '*.tif')))
    listofgeojson = sorted(glob.glob(os.path.join(fullPathGeoJsonDirectory, '*.geojson')))


    if len(listofRaster) != len(listofgeojson):
        print('Error lists do not match fix source errors')

        break

    else:

        for rasterImage, geoJson in zip(listofRaster, listofgeojson):

            annotationName = os.path.basename(rasterImage)
            annotationName = annotationName.replace('.tif', '.xml')
            annotationName = os.path.join(fullPathAnnotationsDirectory, annotationName)


            lT.geoJsonToPascalVOC(annotationName, geoJson, rasterImage,
                                  dataset= 'spacenetV2',
                                  folder_name='spacenetV2',
                                  annotationStyle='PASCAL VOC2012',
                                  segment=True,
                                  bufferSizePix=2.5)

            gdalimage = gdal.Open(rasterImage)
            ## todo make this more workable with 8bit and 16bit
            outputRaster = annotationName.replace('.xml', '.jpg')
            outputRaster = outputRaster.replace('_img', '_8bit_img')

            entry= {'rasterFileName': outputRaster,
                    'geoJsonFileName': geoJson,
                    'annotationName': annotationName,
                    'width': gdalimage.RasterXSize,
                    'height': gdalimage.RasterYSize,
                    'basename': os.path.splitext(os.path.basename(rasterImage))[0]
                    }

            entryList.append(entry)

random.shuffle(entryList)

splitPoint=int(round(len(entryList)*0.8))
trainvalList = entryList[0:splitPoint]
testList     = entryList[splitPoint+1:]

print('creating trainval.txt {} entries'.format(len(trainvalList)))
with open(os.path.join(dataPath, 'trainval.txt'), 'w') as f:
    for entry in trainvalList:
        f.write('{} {}\n'.format(entry['rasterFileName'], entry['annotationName']))

print('creating test.txt {} entries'.format(len(testList)))
with open(os.path.join(dataPath, 'test.txt'), 'w') as f, open(os.path.join(dataPath, 'test_name_size.txt'), 'w') as fname:
    for entry in testList:
        f.write('{} {}\n'.format(entry['rasterFileName'], entry['annotationName']))
        fname.write('{} {} {}\n'.format(entry['basename'], entry['width'], entry['height']))






















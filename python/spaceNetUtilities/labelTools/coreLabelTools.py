import numpy as np
import os
from spaceNetUtilities import geoTools as gT
import math
import pickle
import csv
import glob
from PIL import Image
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.etree import ElementTree
from xml.dom import minidom
import subprocess
import scipy.io
from scipy.sparse import csr_matrix
import json
import re
import shapely
import fiona
import geopandas as gpd
import rasterio
from scipy.ndimage import morphology
from rasterio import features
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.linestring import LineString
from shapely.geometry.multilinestring import MultiLineString
from shapely.geometry import shape, box
from shapely import affinity
from osgeo import gdal, osr, ogr, gdalnumeric

def evaluateLineStringPlane(geom, label='Airplane'):

    poly = Polygon(geom)

    noseToTail = LineString((geom.coords[0], geom.coords[2]))
    wingLength = LineString((geom.coords[1], geom.coords[3]))



    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(geom)

    noseToTail = shapely.ops.tranform(transform_WGS84_To_UTM, noseToTail)
    wingLength = shapely.ops.tranform(transform_WGS84_To_UTM, wingLength)

    Length = noseToTail.length
    Width = wingLength.length
    Aspect = Length/Width
    Direction = (math.atan2(geom.coords[2][0]-geom.coords[0][0], geom.coords[2][1]-geom.coords[0][1])*180/math.pi) % 360

    return [poly, Length, Width, Aspect, Direction]

def evaluateLineStringBoat(geom, label='Boat', aspectRatio=3):


    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(geom)
    geom = shapely.ops.tranform(transform_WGS84_To_UTM, geom)

    pt0 = geom.coords[0] # Stern
    pt1 = geom.coords[1] # Bow
    Length = math.sqrt((pt1[0]-pt0[0])**2 + (pt1[1]-pt0[1])**2)
    Direction = (math.atan2(pt1[0]-pt0[0], pt1[1]-pt0[1])*180/math.pi) % 360


    poly, areaM, angRad, lengthM = gT.createBoxFromLine(geom, aspectRatio,
                                                              transformRequired=True,
                                                              transform_WGS84_To_UTM=transform_WGS84_To_UTM,
                                                              transform_UTM_To_WGS84=transform_UTM_To_WGS84)

    Width = Length/aspectRatio
    Aspect = aspectRatio

    return [poly, Length, Width, Aspect, Direction]

def evaluateLineStringFeature(geom, labelType='Boat', labelName='', aspectRatio=3):

    if labelName=='':
        labelName = labelType
    if labelType=='Boat':

        return evaluateLineStringBoat(geom, label=labelName, aspectRatio=aspectRatio)

    elif labelType=='Airplane':

        return evaluateLineStringPlane(geom, label=labelName)

    else: # Default treat like boat
        return evaluateLineStringBoat(geom, label=labelName, aspectRatio=aspectRatio)


def convertLabelStringToPoly(shapeFileSrc, outGeoJSon, labelType='Airplane', aspectRatio=3, labelName=''):

        source_layer_gdf = gpd.read_file(shapeFileSrc)
        # Create the output Layer
        if os.path.exists(outGeoJSon):
            fiona.remove(outGeoJSon, 'GeoJSON').DeleteDataSource(outGeoJSon)

        source_layer_gdf["geometry"], \
        source_layer_gdf["Length_m"], \
        source_layer_gdf["Width_m"], \
        source_layer_gdf["Aspect(L/W)"], \
        source_layer_gdf["compassDeg"] = source_layer_gdf.apply(
            lambda x: evaluateLineStringFeature(x['geometry'],
                                                labelType=labelType,
                                                aspectRatio=aspectRatio,
                                                labelName=labelName
                                                ),
            axis=1)

        source_layer_gdf.to_file(outGeoJSon, driver='GeoJSON')


def createNPPixArrayDist(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):

    print('Deprecated, use createDistanceTransform')

    return createDistanceTransform(rasterSrc, vectorSrc, npDistFileName=npDistFileName, units=units)


def polygonize(imageArray, transformAffineObject, maskValue=0):

    mask = imageArray!=maskValue

    featureGenerator = features.shapes(imageArray,
                             transform=transformAffineObject,
                             mask=mask)

    return featureGenerator

def createGDFfromShapes(featureGenerator, fieldName='rasterVal'):

    geomList = []
    rasterValList = []

    for feat in featureGenerator:
        geomList.append(shape(feat[0]))
        rasterValList.append(feat[1])

    featureGDF = gpd.GeoDataFrame({'geometry': geomList, fieldName: rasterValList})

    return featureGDF


def createGeoJSONFromRaster(geoJsonFileName,
                            imageArray,
                            geom,
                            crs,
                            maskValue=0,
                            fieldName="rasterVal"):

    featureGenerator = polygonize(imageArray,
                                  geom,
                                  maskValue=maskValue,
                                  fieldName=fieldName)

    featureGDF = createGDFfromShapes(featureGenerator,
                                     fieldName=fieldName)
    featureGDF.crs = crs

    gT.exporttogeojson(geoJsonFileName, featureGDF)

    return featureGDF


def createRasterFromGeoJson(srcGeoJson,
                            srcRasterFileName,
                            outRasterFileName,
                            burnValue=255,
                            burnValueField=''):

    srcGDF = gpd.read_file(srcGeoJson)
    if burnValueField == '':
        featureList = ((geom, value) for geom, value in zip(srcGDF.geometry, burnValue))
    else:
        featureList = ((geom, value) for geom, value in zip(srcGDF.geometry, srcGDF[burnValueField]))

    with rasterio.open(srcRasterFileName) as rst:
        meta = rst.meta.copy()
        meta.update(count=1)
        meta.update(dtype='uint8')

        with rasterio.open(
                outRasterFileName, 'w',
                **meta) as dst:


            burned = features.rasterize(shapes=featureList,
                                        out=rst.shape,
                                        transform=rst.transform)

            dst.write(burned, indexes=1)

    return srcGDF

def createCSVSummaryFile(chipSummaryList, outputFileName, rasterChipDirectory='', replaceImageID='',
                         createProposalsFile=False,
                         pixPrecision=2):


    with open(outputFileName, 'wb') as csvfile:
        writerTotal = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        if createProposalsFile:
            writerTotal.writerow(['ImageId', 'BuildingId', 'PolygonWKT_Pix', 'Confidence'])
        else:
            writerTotal.writerow(['ImageId', 'BuildingId', 'PolygonWKT_Pix', 'PolygonWKT_Geo'])

        for chipSummary in chipSummaryList:
            chipName = chipSummary['chipName']
            print(chipName)
            geoVectorName = chipSummary['geoVectorName']
            #pixVectorName = chipSummary['pixVectorName']
            buildingList = gT.convert_wgs84geojson_to_pixgeojson(geoVectorName,
                                                                 os.path.join(rasterChipDirectory, chipName),
                                                                 pixPrecision=pixPrecision)

            if len(buildingList) > 0:
                for building in buildingList:
                    imageId = os.path.basename(building['ImageId']).replace(replaceImageID, "")
                    if createProposalsFile:
                        writerTotal.writerow([imageId, building['BuildingId'],
                                              building['polyPix'], 1])
                    else:
                        writerTotal.writerow([imageId, building['BuildingId'],
                                          building['polyPix'], building['polyGeo']])
            else:
                imageId = os.path.splitext(os.path.basename(chipName))[0].replace(replaceImageID, "")
                if createProposalsFile:
                    writerTotal.writerow([imageId, -1,
                                      'POLYGON EMPTY', 1])
                else:
                    writerTotal.writerow([imageId, -1,
                                          'POLYGON EMPTY', 'POLYGON EMPTY'])

    return 1

def createCSVSummaryFileFromJsonList(geoJsonList, outputFileName, chipnameList=[],
                                     input='Geo',
                                     replaceImageID=''):

    # Note will not Create Pixle inputs.  No current input for Raster
    if chipnameList:
        pass
    else:
        for geoJson in geoJsonList:
            chipnameList.append(os.path.basename(os.path.splitext(geoJson)[0]))



    with open(outputFileName, 'wb') as csvfile:
        writerTotal = csv.writer(csvfile, delimiter=',', lineterminator='\n')

        writerTotal.writerow(['ImageId', 'BuildingId', 'PolygonWKT_Pix', 'PolygonWKT_Geo'])

        for geoVectorName, chipName in zip(geoJsonList, chipnameList):
            try:
                buildingList = gT.convert_wgs84geojson_to_pixgeojson(geoVectorName, '',
                                                                     image_id=chipName)
                if len(buildingList) > 0:
                    for building in buildingList:
                        imageId = os.path.basename(building['ImageId']).replace(replaceImageID,"")
                        writerTotal.writerow([imageId, building['BuildingId'],
                                              building['polyPix'], building['polyGeo']])
                else:
                    imageId = os.path.splitext(os.path.basename(chipName))[0].replace(replaceImageID,"")
                    writerTotal.writerow([imageId, -1,'"POLYGON EMPTY"', '"POLYGON EMPTY"'])
            except:
                pass

    return 1


def createCSVSummaryFromDirectory(geoJsonDirectory, rasterFileDirectoryList,
                                  aoi_num=0,
                                  aoi_name='TEST',
                                  outputDirectory='',
                                  verbose=False):
    if outputDirectory == '':
        outputDirectory == geoJsonDirectory
    outputbaseName = "AOI_{}_{}_polygons_solution".format(aoi_num, aoi_name)
    #rasterFileDirectoryList = [
    #    ['/usr/local/share/data/AOI_1_RIO/processed2/3band', '3band', '*.tif'],
    #    ['/usr/local/share/data/AOI_1_RIO/processed2/8band', '8band', '*.tif']
    #    ]

    geoJsonList = glob.glob(os.path.join(geoJsonDirectory, '*.geojson'))


    jsonList = []
    chipSummaryList8band = []
    chipSummaryList3band = []
    # AOI_2_RIO_3Band_img997.tif
    # AOI_2_RIO_img635.geojson
    chipsSummaryList = []
    for idx, rasterFile in enumerate(rasterFileDirectoryList):
        chipsSummaryList[idx] = []


    for imageId in geoJsonList:
        imageId = os.path.basename(imageId)

        for idx, rasterFile in enumerate(rasterFileDirectoryList):
            bandName = imageId.replace('.geojson', '.tif')
            bandName = bandName.replace('Geo_', rasterFile[1]+'_')
            if verbose:
                print(imageId)
                print(os.path.join(rasterFile[0], bandName))
            chipSummaryBand = {'chipName': os.path.join(rasterFile[0], bandName),
                                'geoVectorName': os.path.join(geoJsonDirectory, imageId),
                                'imageId': os.path.splitext(imageId)[0]}

            chipsSummaryList[idx].append(chipSummaryBand)


    if verbose:
        print("starting")
    for idx, rasterFile in enumerate(rasterFileDirectoryList):
        createCSVSummaryFile(chipsSummaryList[idx], os.path.join(outputDirectory,
                                                                 outputbaseName+'_'+rasterFile[1]+'.csv'),
                            replaceImageID=rasterFile[1]+'_')

    if verbose:
        print("finished")

    return 1


def createAOIName(AOI_Name, AOI_Num,
                  srcImageryListOrig,
                  srcVectorAOIFile,
                  srcVectorFileList,
                  outputDirectory,
                  clipImageryToAOI=True,
                  windowSizeMeters=200,
                  clipOverlap=0.0,
                  minpartialPerc=0.0,
                  vrtMosaic=True,
                  createPix=False,
                  createSummaryCSVChallenge=True,
                  csvLabel='All',
                  featureName='Buildings',
                  verbose=False):

    srcImageryList = []

    # Clip Imagery to the the AOI provided.  This is important for areas that are completely labeled.
    if clipImageryToAOI:


        for srcImagery in srcImageryListOrig:

            if verbose:
                print(srcImagery)

            AOI_HighResMosaicName = os.path.join(outputDirectory, 'AOI_{}_{}_{}.vrt'.format(AOI_Num, AOI_Name, srcImagery[1]))

            if vrtMosaic:
                AOI_HighResMosaicClipName = AOI_HighResMosaicName.replace('.vrt', 'clipped.vrt')
            else:
                AOI_HighResMosaicClipName = AOI_HighResMosaicName.replace('.vrt', 'clipped.TIF')

            if os.path.isfile(AOI_HighResMosaicClipName):
                os.remove(AOI_HighResMosaicClipName)

            if vrtMosaic:
                command = 'gdalwarp -of VRT ' + \
                        '-cutline ' + srcVectorAOIFile + ' ' + \
                        srcImagery[0] + ' ' + \
                        AOI_HighResMosaicClipName
            else:
                command = 'gdalwarp ' + \
                          '-cutline ' + srcVectorAOIFile + ' ' + \
                          srcImagery[0] + ' ' + \
                          AOI_HighResMosaicClipName
            print(command)
            os.system(command)
            srcImageryList.append([AOI_HighResMosaicClipName, srcImagery[1]])


    else:
        srcImageryList = srcImageryListOrig


        # rasterFileList = [['rasterLocation', 'rasterDescription']]
        # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
        #                       ['/path/to/8band_AOI_1.tif, '8band']
        #                        ]

    chipSummaryList = gT.cutChipFromMosaic(srcImageryList, srcVectorFileList, outlineSrc=srcVectorAOIFile,
                                           outputDirectory=outputDirectory, outputPrefix='',
                                           clipSizeMX=windowSizeMeters, clipSizeMY=windowSizeMeters, clipOverlap=clipOverlap,
                                           minpartialPerc=minpartialPerc, createPix=createPix,
                                           baseName='AOI_{}_{}'.format(AOI_Num, AOI_Name),
                                           imgIdStart=1,
                                           verbose=verbose)


    outputCSVSummaryName = 'AOI_{}_{}_{}_{}_solutions.csv'.format(AOI_Num, AOI_Name, csvLabel,featureName)
    createCSVSummaryFile(chipSummaryList, outputCSVSummaryName, rasterChipDirectory='', replaceImageID='',
                         createProposalsFile=False,
                         pixPrecision=2)




def pixDFToObjectLabelDict(pixGDF,
                              bboxResize=1.0,
                              objectType='building',
                              objectTypeField='',
                              objectPose='Left',
                              objectTruncatedField='',
                              objectDifficultyField=''):

    dictList = []
    # start object segment
    for row in pixGDF.iterrows():

        if objectTypeField=='':
            objectType = objectType
        else:
            objectType = row[objectTypeField]

        if objectTruncatedField=='':
            objectTruncated = 0
        else:
            objectTruncated = row[objectTruncatedField]

        if objectDifficultyField=='':
            objectDifficulty = 0
        else:
            objectDifficulty = row[objectDifficultyField]


        # .bounds returns a tuple (minX,minY, maxX maxY)

        geomBBox = box(row['geometry'].bounds)

        if bboxResize != 1.0:
            geomBBox = affinity.scale(geomBBox, xfact=bboxResize, yfact=bboxResize)

        xmin, ymin, xmax, ymax = geomBBox.bounds

        dictEntry = {'objectType': objectType,
                    'pose': objectPose,
                    'truncated': objectTruncated,
                    'difficult': objectDifficulty,
                    'bndbox': {'xmin': xmin,
                               'ymin': ymin,
                               'xmax': xmax,
                               'ymax': ymax
                               },
                    'geometry': row['geometry'].wkt,
                    }

        dictList.append(dictEntry)

    return dictList

def geoDFtoObjectDict(geoGDF,src_meta, bboxResize=1.0,
                      objectType='building',
                      objectTypeField='',
                      objectPose='Left',
                      objectTruncatedField='',
                      objectDifficultyField=''):


    pixGDF = gT.geoDFtoPixDF(geoGDF, src_meta['transform'])

    objectDictList = pixDFToObjectLabelDict(pixGDF,
                           bboxResize=bboxResize,
                           objectType=objectType,
                           objectTypeField=objectTypeField,
                           objectPose=objectPose,
                           objectTruncatedField=objectTruncatedField,
                           objectDifficultyField=objectDifficultyField)

    return objectDictList

def createRasterSummaryDict(rasterImageName, src_meta, datasetName='SpaceNet_V2',
                  annotationStyle='spaceNet'):

    dictImageDescription = {'folder': datasetName,
                                           'filename': rasterImageName,
                                           'source': {
                                               'database': datasetName,
                                               'annotation': annotationStyle
                                           },
                                           'size': {
                                               'width': src_meta['width'],
                                               'height': src_meta['height'],
                                               'depth': src_meta['count']
                                           },
                                           'segmented': '1',

                                           }


    return dictImageDescription


def geoDFtoDict(geoGDF, rasterImageName, src_meta, datasetName='SpaceNet_V2',
                  annotationStyle='spaceNet',
                  bboxResize=1.0,
                  objectType='building',
                  objectTypeField='',
                  objectPose='Left',
                  objectTruncatedField='',
                  objectDifficultyField=''
                  ):

    imageDescriptionDict = createRasterSummaryDict(rasterImageName, src_meta, datasetName=datasetName,
                  annotationStyle=annotationStyle)

    objectDictList = geoDFtoObjectDict(geoGDF,src_meta, bboxResize=bboxResize,
                           objectType=objectType,
                           objectTypeField=objectTypeField,
                           objectPose=objectPose,
                           objectTruncatedField=objectTruncatedField,
                           objectDifficultyField=objectDifficultyField)


    return imageDescriptionDict, objectDictList



def geoJsontoDict(geoJson, rasterImageName, datasetName='SpaceNet_V2',
                  annotationStyle='spaceNet',
                  bboxResize=1.0,
                  objectType='building',
                  objectTypeField='',
                  objectPose='Left',
                  objectTruncatedField='',
                  objectDifficultyField=''
                  ):

    geoGDF = gpd.read_file(geoJson)
    with rasterio.open(rasterImageName) as src:
        src_meta = src.meta.copy()

    return geoDFtoDict(geoGDF, rasterImageName, src_meta, datasetName=datasetName,
                    annotationStyle=annotationStyle,
                    bboxResize=bboxResize,
                    objectType=objectType,
                    objectTypeField=objectTypeField,
                    objectPose=objectPose,
                    objectTruncatedField=objectTruncatedField,
                    objectDifficultyField=objectDifficultyField
                    )



def convertPixDimensionToPercent(size, box):
    '''Input = image size: (w,h), box: [x0, x1, y0, y1]'''
    #TODO change box to use shapely bounding box format
    dw = 1./size[0]
    dh = 1./size[1]
    xmid = (box[0] + box[1])/2.0
    ymid = (box[2] + box[3])/2.0
    w0 = box[1] - box[0]
    h0 = box[3] - box[2]
    x = xmid*dw
    y = ymid*dh
    w = w0*dw
    h = h0*dh

    return (x,y,w,h)

def geoJsonToDARKNET(xmlFileName, geoJson, rasterImageName, im_id='',
                     dataset ='SpaceNet',
                     folder_name='spacenet',
                     annotationStyle = 'DARKNET',
                     segment=True,
                     bufferSizePix=2.5,
                     convertTo8Bit=True,
                     outputPixType='Byte',
                     outputFormat='GTiff',
                     bboxResize=1.0):
    xmlFileName = xmlFileName.replace(".xml", ".txt")
    print("creating {}".format(xmlFileName))

    buildingList = gT.convert_wgs84geojson_to_pixgeojson(geoJson, rasterImageName, image_id=[], pixelgeojson=[], only_polygons=True,
                                       breakMultiPolygonGeo=True, pixPrecision=2)
    #                        buildinglist.append({'ImageId': image_id,
                                             #'BuildingId': building_id,
                                             #'polyGeo': ogr.CreateGeometryFromWkt(geom.ExportToWkt()),
                                             #'polyPix': ogr.CreateGeometryFromWkt('POLYGON EMPTY')
                                             #})


    srcRaster = gdal.Open(rasterImageName)
    outputRaster = rasterImageName
    if convertTo8Bit:
        cmd = ['gdal_translate', '-ot', outputPixType, '-of', outputFormat, '-co', 'PHOTOMETRIC=rgb']
        scaleList = []
        for bandId in range(srcRaster.RasterCount):
            bandId = bandId+1
            band=srcRaster.GetRasterBand(bandId)
            min = band.GetMinimum()
            max = band.GetMaximum()

            # if not exist minimum and maximum values
            if min is None or max is None:
                (min, max) = band.ComputeRasterMinMax(1)
            cmd.append('-scale_{}'.format(bandId))
            cmd.append('{}'.format(0))
            cmd.append('{}'.format(max))
            cmd.append('{}'.format(0))
            cmd.append('{}'.format(255))

        cmd.append(rasterImageName)
        if outputFormat == 'JPEG':
            outputRaster = xmlFileName.replace('.txt', '.jpg')
        else:
            outputRaster = xmlFileName.replace('.txt', '.tif')

        outputRaster = outputRaster.replace('_img', '_8bit_img')
        cmd.append(outputRaster)
        print(cmd)
        subprocess.call(cmd)

    with open(xmlFileName, 'w') as f:

        for building in buildingList:


            # Get Envelope returns a tuple (minX, maxX, minY, maxY)

            boxDim = building['polyPix'].GetEnvelope()

            if bboxResize != 1.0:
                xmin = boxDim[0]
                ymin = boxDim[2]
                xmax = boxDim[1]
                ymax = boxDim[3]
                xCenter = (xmin + xmax) / 2
                yCenter = (ymin + ymax) / 2
                bboxNewHalfHeight = ((ymax - ymin) / 2) * bboxResize
                bboxNewHalfWidth = ((ymax - ymin) / 2) * bboxResize
                xmin = xCenter - bboxNewHalfWidth
                xmax = xCenter + bboxNewHalfWidth
                ymin = yCenter - bboxNewHalfHeight
                ymax = yCenter + bboxNewHalfHeight

                boxDim = [xmin, xmax, ymin, ymax]

            rasterSize = (srcRaster.RasterXSize, srcRaster.RasterYSize)

            lineOutput = convertPixDimensionToPercent(rasterSize, boxDim)
            classNum=0
            f.write('{} {} {} {} {}\n'.format(classNum,
                                             lineOutput[0],
                                             lineOutput[1],
                                             lineOutput[2],
                                             lineOutput[3])
                    )

    entry = {'rasterFileName': outputRaster,
             'geoJsonFileName': geoJson,
             'annotationName': xmlFileName,
             'width': srcRaster.RasterXSize,
             'height': srcRaster.RasterYSize,
             'depth' : srcRaster.RasterCount,
             'basename': os.path.splitext(os.path.basename(rasterImageName))[0]
             }

    return entry

def createDistanceTransform(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):

    #TODO evaluate scipy.ndimage.morphology.distance_transform_edt vs cv2.distanceTransform for speed

    # image = features.rasterize(
    #         ((g, 255) for g, v in shapes),
    #         out_shape=src.shape,
    #         transform=src.transform)
    ## open source vector file that truth data
    ## calculate pixelSize in meters i.e. GSD

    with rasterio.open(rasterSrc) as srcRas_ds:
        src_transform = srcRas_ds.transform
        src_shape = srcRas_ds.shape
        ## calculate pixel size in meters i.e. GSD
        if units == 'meters':
            geoTrans, poly, ulX, ulY, lrX, lrY = gT.getRasterExtent(srcRas_ds)
            transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(poly)
            line = LineString([(geoTrans.c, geoTrans.f),
                               (geoTrans.c + geoTrans.a, geoTrans.f + geoTrans.e)
                               ]
                              )

            line = shapely.ops.tranform(transform_WGS84_To_UTM, line)
            metersIndex = line.lenth
        else:
            metersIndex = 1

    ## Burn source layer into image
    source_layer = gpd.read_file(vectorSrc)
    shapes = ((geom,value) for geom, value in zip(source_layer.geometry, 255))
    baseImage = features.rasterize(shapes,
                               out_shape=src_shape,
                               transform=src_transform)

    ## calculate Distance between Feature point and closest background pixel
    # distance_transform_edt takes the place of ComputeProximity

    #calculate distance from Any point inside a feature to the closest background pixel
    interiorDist = morphology.distance_transform_edt(baseImage)

    #inverse Image so that background is a feature and interior feature values = 0
    inverseImage = baseImage
    inverseImage[inverseImage==0]=200
    inverseImage[inverseImage==255]=0

    # calculate distance from any point exterior a feature to the closet feature pixel
    exteriorDist = morphology.distance_transform_edt(inverseImage)


    proxTotal = interiorDist - exteriorDist
    proxTotal = proxTotal*metersIndex

    if npDistFileName != '':
        np.save(npDistFileName, proxTotal)

    return proxTotal





def convertGTiffTo8Bit(rasterImageName, outputImageName, outputFormat='GTiff',
                       minPercent=2,
                       maxPercent=98,
                       maxValue=255):

    # Other Format would be JPG


    with rasterio.open(rasterImageName) as src:
        data = src.read()
        profile = src.profile

    # resize band
    for bandId in range(profile['count']):
         # data[bandId,:,:]
        minBand = np.percentile(data[bandId,:,:], minPercent)
        maxBand = np.percentile(data[bandId,:,:], maxPercent)

        # Clip lower percent of values to Zero
        data[bandId, :, :] = data[bandId,:,:] - minBand
        data[bandId, :, :][data[bandId,:,:]<=0] = 0

        #Normalize band to between 0 and 1 saturating values above the maxPercent
        data[bandId, :, :] = data[bandId,:,:]/(maxBand-minBand)
        data[bandId, :, :][data[bandId,:,:]>1]=1

        # Rescale to max Value (255)
        data[bandId, :, :] = data[bandId,:,:]*maxValue
        data[bandId, :, :][data[bandId,:,:]>maxValue]=maxValue



    profile.update(dtype=rasterio.uint8,
                   driver=outputFormat)


    with rasterio.open(outputImageName, 'w', **profile) as dst:
        dst.write(data.astype(rasterio.uint8))


    return outputImageName

from osgeo import gdal, osr, ogr, gdalnumeric
import numpy as np
import os
import geoTools as gT
import math
import cPickle as pickle
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


def evaluateLineStringPlane(geom, label='Airplane'):
    ring = ogr.Geometry(ogr.wkbLinearRing)

    for i in range(0, geom.GetPointCount()):
        # GetPoint returns a tuple not a Geometry
        pt = geom.GetPoint(i)
        ring.AddPoint(pt[0], pt[1])
    pt = geom.GetPoint(0)
    ring.AddPoint(pt[0], pt[1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(geom)
    geom.Transform(transform_WGS84_To_UTM)
    pt0 = geom.GetPoint(0) # Tail
    pt1 = geom.GetPoint(1) # Wing
    pt2 = geom.GetPoint(2) # Nose
    pt3 = geom.GetPoint(3) # Wing
    Length = math.sqrt((pt2[0]-pt0[0])**2 + (pt2[1]-pt0[1])**2)
    Width = math.sqrt((pt3[0] - pt1[0])**2 + (pt3[1] - pt1[1])**2)
    Aspect = Length/Width
    Direction = (math.atan2(pt2[0]-pt0[0], pt2[1]-pt0[1])*180/math.pi) % 360


    geom.Transform(transform_UTM_To_WGS84)

    return [poly, Length, Width, Aspect, Direction]

def evaluateLineStringBoat(geom, label='Boat', aspectRatio=3):


    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(geom)

    geom.Transform(transform_WGS84_To_UTM)
    pt0 = geom.GetPoint(0) # Stern
    pt1 = geom.GetPoint(1) # Bow
    Length = math.sqrt((pt1[0]-pt0[0])**2 + (pt1[1]-pt0[1])**2)
    Direction = (math.atan2(pt1[0]-pt0[0], pt1[1]-pt0[1])*180/math.pi) % 360
    geom.Transform(transform_UTM_To_WGS84)

    poly, areaM, angRad, lengthM = gT.createBoxFromLine(geom, aspectRatio,
                                                              transformRequired=True,
                                                              transform_WGS84_To_UTM=transform_WGS84_To_UTM,
                                                              transform_UTM_To_WGS84=transform_UTM_To_WGS84)

    Width = Length/aspectRatio
    Aspect = aspectRatio

    return [poly, Length, Width, Aspect, Direction]


def convertLabelStringToPoly(shapeFileSrc, outGeoJSon, labelType='Airplane'):

        shapeSrc = ogr.Open(shapeFileSrc)
        source_layer = shapeSrc.GetLayer()
        source_srs = source_layer.GetSpatialRef()
        # Create the output Layer
        outDriver = ogr.GetDriverByName("geojson")
        if os.path.exists(outGeoJSon):
            outDriver.DeleteDataSource(outGeoJSon)


        outDataSource = outDriver.CreateDataSource(outGeoJSon)
        outLayer = outDataSource.CreateLayer("groundTruth", source_srs, geom_type=ogr.wkbPolygon)
        # Add input Layer Fields to the output Layer
        inLayerDefn = source_layer.GetLayerDefn()
        for i in range(0, inLayerDefn.GetFieldCount()):
            fieldDefn = inLayerDefn.GetFieldDefn(i)
            outLayer.CreateField(fieldDefn)
        outLayer.CreateField(ogr.FieldDefn("Length_m", ogr.OFTReal))
        outLayer.CreateField(ogr.FieldDefn("Width_m", ogr.OFTReal))
        outLayer.CreateField(ogr.FieldDefn("Aspect(L/W)", ogr.OFTReal))
        outLayer.CreateField(ogr.FieldDefn("compassDeg", ogr.OFTReal))

        outLayerDefn = outLayer.GetLayerDefn()
        for inFeature in source_layer:

            outFeature = ogr.Feature(outLayerDefn)

            for i in range(0, inLayerDefn.GetFieldCount()):
                outFeature.SetField(inLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))

            geom = inFeature.GetGeometryRef()
            if labelType == 'Airplane':
                poly, Length, Width, Aspect, Direction = evaluateLineStringPlane(geom, label='Airplane')
            elif labelType == 'Boat':
                poly, Length, Width, Aspect, Direction = evaluateLineStringBoat(geom, label='Boat')

            outFeature.SetGeometry(poly)
            outFeature.SetField("Length_m", Length)
            outFeature.SetField("Width_m", Width)
            outFeature.SetField("Aspect(L/W)", Aspect)
            outFeature.SetField("compassDeg", Direction)

            outLayer.CreateFeature(outFeature)


def createTruthPixelLinePickle(truthLineFile, pickleLocation=''):
    if pickleLocation=='':
        extension = os.path.splitext(truthLineFile)[1]
        pickleLocation = truthLineFile.replace(extension, 'Pixline.p')
    if truthLineFile != '':
        # get Source Line File Information
        shapef = ogr.Open(truthLineFile, 0)
        truthLayer = shapef.GetLayer()
        pt1X = []
        pt1Y = []
        pt2X = []
        pt2Y = []
        for tmpFeature in truthLayer:
            tmpGeom = tmpFeature.GetGeometryRef()
            for i in range(0, tmpGeom.GetPointCount()):
                pt = tmpGeom.GetPoint(i)

                if i == 0:
                    pt1X.append(pt[0])
                    pt1Y.append(pt[1])
                elif i == 1:
                    pt2X.append(pt[0])
                    pt2Y.append(pt[1])

        lineData = {'pt1X': np.asarray(pt1X),
                    'pt1Y': np.asarray(pt1Y),
                    'pt2X': np.asarray(pt2X),
                    'pt2Y': np.asarray(pt2Y)
                    }

        with open(pickleLocation, 'wb') as f:
            pickle.dump(lineData, f)
            # get Source Line File Information


def createTruthPixelPolyPickle(truthPoly, pickleLocation=''):
    # returns dictionary with list of minX, maxX, minY, maxY

    if pickleLocation=='':
        extension = os.path.splitext(truthPoly)[1]
        pickleLocation = truthPoly.replace(extension, 'PixPoly.p')
    if truthPoly != '':
        # get Source Line File Information
        shapef = ogr.Open(truthPoly, 0)
        truthLayer = shapef.GetLayer()
        envList = []

        for tmpFeature in truthLayer:
            tmpGeom = tmpFeature.GetGeometryRef()
            env = tmpGeom.GetEvnelope()
            envList.append(env)

        envArray = np.asarray(envList)
        envelopeData = {'minX': envArray[:,0],
                        'maxX': envArray[:,1],
                        'minY': envArray[:,2],
                        'maxY': envArray[:,3]
                        }


        with open(pickleLocation, 'wb') as f:
            pickle.dump(envelopeData, f)
            # get Source Line File Information


def createNPPixArrayDist(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):

    ## open source vector file that truth data
    source_ds = ogr.Open(vectorSrc)
    source_layer = source_ds.GetLayer()

    ## extract data from src Raster File to be emulated
    ## open raster file that is to be emulated
    srcRas_ds = gdal.Open(rasterSrc)
    cols = srcRas_ds.RasterXSize
    rows = srcRas_ds.RasterYSize
    noDataValue = 0

    if units=='meters':
        geoTrans, poly, ulX, ulY, lrX, lrY = gT.getRasterExtent(srcRas_ds)
        transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(poly)
        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(geoTrans[0], geoTrans[3])
        line.AddPoint(geoTrans[0]+geoTrans[1], geoTrans[3])

        line.Transform(transform_WGS84_To_UTM)
        metersIndex = line.Length()
    else:
        metersIndex = 1

    ## create First raster memory layer
    memdrv = gdal.GetDriverByName('MEM')
    dst_ds = memdrv.Create('', cols, rows, 1, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    dst_ds.SetProjection(srcRas_ds.GetProjection())
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(noDataValue)

    gdal.RasterizeLayer(dst_ds, [1], source_layer, burn_values=[255])
    srcBand = dst_ds.GetRasterBand(1)

    memdrv2 = gdal.GetDriverByName('MEM')
    prox_ds = memdrv2.Create('', cols, rows, 1, gdal.GDT_Int16)
    prox_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    prox_ds.SetProjection(srcRas_ds.GetProjection())
    proxBand = prox_ds.GetRasterBand(1)
    proxBand.SetNoDataValue(noDataValue)

    options = ['NODATA=0']

    gdal.ComputeProximity(srcBand, proxBand, options)

    memdrv3 = gdal.GetDriverByName('MEM')
    proxIn_ds = memdrv3.Create('', cols, rows, 1, gdal.GDT_Int16)
    proxIn_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    proxIn_ds.SetProjection(srcRas_ds.GetProjection())
    proxInBand = proxIn_ds.GetRasterBand(1)
    proxInBand.SetNoDataValue(noDataValue)
    options = ['NODATA=0', 'VALUES=0']
    gdal.ComputeProximity(srcBand, proxInBand, options)

    proxIn = gdalnumeric.BandReadAsArray(proxInBand)
    proxOut = gdalnumeric.BandReadAsArray(proxBand)

    proxTotal = proxIn.astype(float) - proxOut.astype(float)
    proxTotal = proxTotal*metersIndex

    if npDistFileName != '':
        np.save(npDistFileName, proxTotal)

    return proxTotal


def createGeoJSONFromRaster(geoJsonFileName, array2d, geom, proj,
                            layerName="BuildingID",
                            fieldName="BuildingID"):

    memdrv = gdal.GetDriverByName('MEM')
    src_ds = memdrv.Create('', array2d.shape[1], array2d.shape[0], 1)
    src_ds.SetGeoTransform(geom)
    src_ds.SetProjection(proj)
    band = src_ds.GetRasterBand(1)
    band.WriteArray(array2d)

    dst_layername = "BuildingID"
    drv = ogr.GetDriverByName("geojson")
    dst_ds = drv.CreateDataSource(geoJsonFileName)
    dst_layer = dst_ds.CreateLayer(layerName, srs=None)

    fd = ogr.FieldDefn(fieldName, ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = 1

    gdal.Polygonize(band, None, dst_layer, dst_field, [], callback=None)

    return


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



def createCSVSummaryFromDirectory(geoJsonDirectory, rasterFileDirectoryList,
                                  aoi_num=0,
                                  aoi_name='TEST',
                                  outputDirectory=''):
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
            print imageId
            print os.path.join(rasterFile[0], bandName)
            chipSummaryBand = {'chipName': os.path.join(rasterFile[0], bandName),
                                'geoVectorName': os.path.join(geoJsonDirectory, imageId),
                                'imageId': os.path.splitext(imageId)[0]}

            chipsSummaryList[idx].append(chipSummaryBand)


    print "starting"
    for idx, rasterFile in enumerate(rasterFileDirectoryList):
        createCSVSummaryFile(chipsSummaryList[idx], os.path.join(outputDirectory,
                                                                 outputbaseName+'_'+rasterFile[1]+'.csv'),
                            replaceImageID=rasterFile[1]+'_')


    print "finished"



def createRasterFromGeoJson(srcGeoJson, srcRasterFileName, outRasterFileName):
    NoData_value = 0
    source_ds = ogr.Open(srcGeoJson)
    source_layer = source_ds.GetLayer()

    srcRaster = gdal.Open(srcRasterFileName)


    # Create the destination data source
    target_ds = gdal.GetDriverByName('GTiff').Create(outRasterFileName, srcRaster.RasterXSize, srcRaster.RasterYSize, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform(srcRaster.GetGeoTransform())
    target_ds.SetProjection(srcRaster.GetProjection())
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(NoData_value)

    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])
    band.FlushCache()




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
                  featureName='Buildings'):

    srcImageryList = []
    if clipImageryToAOI:


        for srcImagery in srcImageryListOrig:

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
                                           imgIdStart=1)


    outputCSVSummaryName = 'AOI_{}_{}_{}_{}_solutions.csv'.format(AOI_Num, AOI_Name, csvLabel,featureName)
    createCSVSummaryFile(chipSummaryList, outputCSVSummaryName, rasterChipDirectory='', replaceImageID='',
                         createProposalsFile=False,
                         pixPrecision=2)




def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")



def geoJsonToPASCALVOC2012(xmlFileName, geoJson, rasterImageName, im_id='',
                           dataset ='SpaceNet',
                           folder_name='spacenet',
                           annotationStyle = 'PASCAL VOC2012',
                           segment=True,
                           bufferSizePix=2.5,
                           convertTo8Bit=True,
                           outputPixType='Byte',
                           outputFormat='GTiff',
                           bboxResize=1.0):

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
        cmd = ['gdal_translate', '-ot', outputPixType, '-of', outputFormat, '-co', '"PHOTOMETRIC=rgb"']
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
            outputRaster = xmlFileName.replace('.xml', '.jpg')
        else:
            outputRaster = xmlFileName.replace('.xml', '.tif')

        outputRaster = outputRaster.replace('_img', '_8bit_img')
        cmd.append(outputRaster)
        print(cmd)
        subprocess.call(cmd)


    if segment:
        segmented = 1  # 1=True, 0 = False
    else:
        segmented = 0

    top = Element('annotation')
    ## write header
    childFolder = SubElement(top, 'folder')
    childFolder.text = dataset
    childFilename = SubElement(top, 'filename')
    childFilename.text = rasterImageName

    # write source block
    childSource = SubElement(top, 'source')
    SubElement(childSource, 'database').text = dataset
    SubElement(childSource, 'annotation').text = annotationStyle

    # write size block
    childSize = SubElement(top, 'size')
    SubElement(childSize, 'width').text = str(srcRaster.RasterXSize)
    SubElement(childSize, 'height').text = str(srcRaster.RasterYSize)
    SubElement(childSize, 'depth').text = str(srcRaster.RasterCount)

    SubElement(top, 'segmented').text = str(segmented)

    # start object segment
    for building in buildingList:
        objectType = 'building'
        objectPose = 'Left'
        objectTruncated = 0  # 1=True, 0 = False
        objectDifficulty = 0  # 0 Easy - 3 Hard

        env = building['polyPix'].GetEnvelope()
        xmin=env[0]
        ymin=env[2]
        xmax=env[1]
        ymax=env[3]

        if bboxResize != 1.0:
            xCenter = (xmin+xmax)/2
            yCenter = (ymin+ymax)/2
            bboxNewHalfHeight = ((ymax-ymin)/2)*bboxResize
            bboxNewHalfWidth  = ((ymax - ymin) / 2)*bboxResize
            xmin = xCenter - bboxNewHalfWidth
            xmax = xCenter + bboxNewHalfWidth
            ymin = yCenter - bboxNewHalfHeight
            ymax = yCenter + bboxNewHalfHeight

        # Get Envelope returns a tuple (minX, maxX, minY, maxY)


        childObject = SubElement(top, 'object')
        SubElement(childObject, 'name').text = objectType
        SubElement(childObject, 'pose').text = objectPose
        SubElement(childObject, 'truncated').text = str(objectTruncated)
        SubElement(childObject, 'difficult').text = str(objectDifficulty)
        # write bounding box
        childBoundBox = SubElement(childObject, 'bndbox')
        SubElement(childBoundBox, 'xmin').text = str(int(round(xmin)))
        SubElement(childBoundBox, 'ymin').text = str(int(round(ymin)))
        SubElement(childBoundBox, 'xmax').text = str(int(round(xmax)))
        SubElement(childBoundBox, 'ymax').text = str(int(round(ymax)))

    with open(xmlFileName, 'w') as f:
        f.write(prettify(top))


    print('creating segmentation')
    if segment:
        NoData_value = -9999

        source_ds = ogr.Open(geoJson)
        source_layer = source_ds.GetLayer()
        srs = source_layer.GetSpatialRef()
        memDriver = ogr.GetDriverByName('MEMORY')
        outerBuffer=memDriver.CreateDataSource('outer')
        outerBufferLayer = outerBuffer.CreateLayer("test", srs, geom_type=ogr.wkbPolygon)
        innerBuffer = memDriver.CreateDataSource('inner')
        innerBufferLayer = innerBuffer.CreateLayer("test2", srs, geom_type=ogr.wkbPolygon)

        idField = ogr.FieldDefn("objid", ogr.OFTInteger)
        innerBufferLayer.CreateField(idField)

        featureDefn = innerBufferLayer.GetLayerDefn()
        bufferDist = srcRaster.GetGeoTransform()[1]*bufferSizePix
        for idx, feature in enumerate(source_layer):
            ingeom = feature.GetGeometryRef()
            geomBufferOut = ingeom.Buffer(bufferDist)
            geomBufferIn  = ingeom.Buffer(-bufferDist)
            print(geomBufferIn.ExportToWkt())
            print(geomBufferIn.IsEmpty())
            print(geomBufferIn.IsSimple())

            if geomBufferIn.GetArea()>0.0:
                outBufFeature = ogr.Feature(featureDefn)
                outBufFeature.SetGeometry(geomBufferOut)

                outerBufferLayer.CreateFeature(outBufFeature)

                inBufFeature = ogr.Feature(featureDefn)
                inBufFeature.SetGeometry(geomBufferIn)
                inBufFeature.SetField('objid', idx)
                innerBufferLayer.CreateFeature(inBufFeature)

                outBufFeature = None
                inBufFeature = None



        print('writing GTIFF sgcls')
        print('rasterToWrite = {}'.format(xmlFileName.replace('.xml', 'segcls.tif')))
        target_ds = gdal.GetDriverByName('GTiff').Create(xmlFileName.replace('.xml', 'segcls.tif'), srcRaster.RasterXSize, srcRaster.RasterYSize, 1, gdal.GDT_Byte)
        print('setTransform')
        target_ds.SetGeoTransform(srcRaster.GetGeoTransform())
        print('setProjection')
        target_ds.SetProjection(srcRaster.GetProjection())
        print('getBand')
        band = target_ds.GetRasterBand(1)
        print('setnodata')
        band.SetNoDataValue(NoData_value)

        # Rasterize
        print('rasterize outer buffer')
        gdal.RasterizeLayer(target_ds, [1], outerBufferLayer, burn_values=[255])
        print('rasterize inner buffer')
        gdal.RasterizeLayer(target_ds, [1], innerBufferLayer, burn_values=[100])
        print('writing png sgcls')
        # write to .png
        imageArray = np.array(target_ds.GetRasterBand(1).ReadAsArray())
        im = Image.fromarray(imageArray)
        im.save(xmlFileName.replace('.xml', 'segcls.png'))

        print('writing GTIFF sgobj')
        ## create objectSegment
        target_ds = gdal.GetDriverByName('GTiff').Create(xmlFileName.replace('.xml', 'segobj.tif'),
                                                         srcRaster.RasterXSize, srcRaster.RasterYSize, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform(srcRaster.GetGeoTransform())
        target_ds.SetProjection(srcRaster.GetProjection())
        band = target_ds.GetRasterBand(1)

        band.SetNoDataValue(NoData_value)

        # Rasterize
        gdal.RasterizeLayer(target_ds, [1], outerBufferLayer, burn_values=[255])
        gdal.RasterizeLayer(target_ds, [1], innerBufferLayer, burn_values=[100], options=['ATTRIBUTE=objid'])
        print('writing png sgobj')
        # write to .png
        imageArray = np.array(target_ds.GetRasterBand(1).ReadAsArray())
        im = Image.fromarray(imageArray)
        im.save(xmlFileName.replace('.xml', 'segobj.png'))

    entry = {'rasterFileName': outputRaster,
                 'geoJsonFileName': geoJson,
                 'annotationName': xmlFileName,
                 'width': srcRaster.RasterXSize,
                 'height': srcRaster.RasterYSize,
                 'depth' : srcRaster.RasterCount,
                 'basename': os.path.splitext(os.path.basename(rasterImageName))[0]
                 }

    return entry

def convertPixDimensionToPercent(size, box):
    '''Input = image size: (w,h), box: [x0, x1, y0, y1]'''
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

    ## open source vector file that truth data
    source_ds = ogr.Open(vectorSrc)
    source_layer = source_ds.GetLayer()

    ## extract data from src Raster File to be emulated
    ## open raster file that is to be emulated
    srcRas_ds = gdal.Open(rasterSrc)
    cols = srcRas_ds.RasterXSize
    rows = srcRas_ds.RasterYSize
    noDataValue = 0

    if units=='meters':
        geoTrans, poly, ulX, ulY, lrX, lrY = gT.getRasterExtent(srcRas_ds)
        transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = gT.createUTMTransform(poly)
        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(geoTrans[0], geoTrans[3])
        line.AddPoint(geoTrans[0]+geoTrans[1], geoTrans[3])

        line.Transform(transform_WGS84_To_UTM)
        metersIndex = line.Length()
    else:
        metersIndex = 1

    ## create First raster memory layer
    memdrv = gdal.GetDriverByName('MEM')
    dst_ds = memdrv.Create('', cols, rows, 1, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    dst_ds.SetProjection(srcRas_ds.GetProjection())
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(noDataValue)

    gdal.RasterizeLayer(dst_ds, [1], source_layer, burn_values=[255])
    srcBand = dst_ds.GetRasterBand(1)

    memdrv2 = gdal.GetDriverByName('MEM')
    prox_ds = memdrv2.Create('', cols, rows, 1, gdal.GDT_Int16)
    prox_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    prox_ds.SetProjection(srcRas_ds.GetProjection())
    proxBand = prox_ds.GetRasterBand(1)
    proxBand.SetNoDataValue(noDataValue)
    options = ['NODATA=0']
    
    ##compute distance to non-zero pixel values and scrBand and store in proxBand
    gdal.ComputeProximity(srcBand, proxBand, options)

    memdrv3 = gdal.GetDriverByName('MEM')
    proxIn_ds = memdrv3.Create('', cols, rows, 1, gdal.GDT_Int16)
    proxIn_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    proxIn_ds.SetProjection(srcRas_ds.GetProjection())
    proxInBand = proxIn_ds.GetRasterBand(1)
    proxInBand.SetNoDataValue(noDataValue)
    options = ['NODATA=0', 'VALUES=0']
    
    ##compute distance to zero pixel values and scrBand and store in proxInBand
    gdal.ComputeProximity(srcBand, proxInBand, options)

    proxIn = gdalnumeric.BandReadAsArray(proxInBand)
    proxOut = gdalnumeric.BandReadAsArray(proxBand)
 
    ##distance tranform is the distance to zero pixel values minus distance to non-zero pixel values
    proxTotal = proxIn.astype(float) - proxOut.astype(float)
    proxTotal = proxTotal*metersIndex

    if npDistFileName != '':
        np.save(npDistFileName, proxTotal)

    return proxTotal


def createClassSegmentation(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
    dist_trans = createDistanceTransform(rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    dist_trans[dist_trans > 0] = 1
    dist_trans[dist_trans < 0] = 0
    return dist_trans


def createClassBoundaries(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
    dist_trans = createDistanceTransform(rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    #From distance transform to boundary
    dist_trans[dist_trans > 1.0] = 255
    dist_trans[dist_trans < -1.0] = 255
    dist_trans[dist_trans != 255] = 1
    dist_trans[dist_trans == 255] = 0
    sparse_total = csr_matrix(dist_trans)
    return sparse_total.astype(np.uint8)


def createClassCategoriesPresent(vectorSrc):
    with open(vectorSrc) as my_file:
        data = json.load(my_file)
    if(len(data['features']) == 0):
       return np.array([],dtype=np.uint8)
    else:
       return np.array([1],dtype=np.uint8)

def createDistanceTransformByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
    ## open source vector file that truth data
    source_ds = ogr.Open(vectorSrc)
    source_layer = source_ds.GetLayer()
    
    #Define feature
    my_feature = source_layer[feature_index]
    
    #Spatial Reference
    srs = source_layer.GetSpatialRef()
    
    #Create feature Layer
    outDriver = ogr.GetDriverByName('MEMORY')
    outDataSource = outDriver.CreateDataSource('memData')
    Feature_Layer = outDataSource.CreateLayer("this_feature", srs, geom_type=ogr.wkbPolygon)
    
    #Add feature to layer
    Feature_Layer.CreateFeature(my_feature)
    

    ## extract data from src Raster File to be emulated
    ## open raster file that is to be emulated
    srcRas_ds = gdal.Open(rasterSrc)
    cols = srcRas_ds.RasterXSize
    rows = srcRas_ds.RasterYSize
    noDataValue = 0
    metersIndex = 1

    ## create First raster memory layer
    memdrv = gdal.GetDriverByName('MEM')
    dst_ds = memdrv.Create('', cols, rows, 1, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    dst_ds.SetProjection(srcRas_ds.GetProjection())
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(noDataValue)

    gdal.RasterizeLayer(dst_ds, [1], Feature_Layer, burn_values=[255])
    srcBand = dst_ds.GetRasterBand(1)

    memdrv2 = gdal.GetDriverByName('MEM')
    prox_ds = memdrv2.Create('', cols, rows, 1, gdal.GDT_Int16)
    prox_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    prox_ds.SetProjection(srcRas_ds.GetProjection())
    proxBand = prox_ds.GetRasterBand(1)
    proxBand.SetNoDataValue(noDataValue)

    options = ['NODATA=0']

    gdal.ComputeProximity(srcBand, proxBand, options)

    memdrv3 = gdal.GetDriverByName('MEM')
    proxIn_ds = memdrv3.Create('', cols, rows, 1, gdal.GDT_Int16)
    proxIn_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
    proxIn_ds.SetProjection(srcRas_ds.GetProjection())
    proxInBand = proxIn_ds.GetRasterBand(1)
    proxInBand.SetNoDataValue(noDataValue)
    options = ['NODATA=0', 'VALUES=0']
    gdal.ComputeProximity(srcBand, proxInBand, options)

    proxIn = gdalnumeric.BandReadAsArray(proxInBand)
    proxOut = gdalnumeric.BandReadAsArray(proxBand)

    proxTotal = proxIn.astype(float) - proxOut.astype(float)
    proxTotal = proxTotal*metersIndex

    if npDistFileName != '':
        np.save(npDistFileName, proxTotal)
        
    return proxTotal

def createSegmentationByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
    dist_trans_by_feature = createDistanceTransformByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    dist_trans_by_feature[dist_trans_by_feature > 0] = feature_index + 1
    dist_trans_by_feature[dist_trans_by_feature < 0] = 0  
    return dist_trans_by_feature.astype(np.uint8)

def createInstanceSegmentation(rasterSrc, vectorSrc):
    json_data = open(vectorSrc)
    data = json.load(json_data)
    num_features = len(data['features'])
    
    cell_array = np.zeros((num_features,), dtype=np.object)
    for i in range(num_features):
        cell_array[i] = createSegmentationByFeatureIndex(i, rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    return cell_array

def createBoundariesByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
    dist_trans_by_feature = createDistanceTransformByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    dist_trans_by_feature[dist_trans_by_feature > 1.0] = 255
    dist_trans_by_feature[dist_trans_by_feature < -1.0] = 255
    dist_trans_by_feature[dist_trans_by_feature != 255] = 1
    dist_trans_by_feature[dist_trans_by_feature == 255] = 0
    return dist_trans_by_feature.astype(np.uint8)

def createInstanceBoundaries(rasterSrc, vectorSrc):
    json_data = open(vectorSrc)
    data = json.load(json_data)
    num_features = len(data['features'])
    
    cell_array = np.zeros((num_features,), dtype=np.object)
    for i in range(num_features):
        full_boundary_matrix = createBoundariesByFeatureIndex(i, rasterSrc, vectorSrc, npDistFileName='', units='pixels')
        cell_array[i] = csr_matrix(full_boundary_matrix)
    return cell_array

def createInstanceCategories(vectorSrc):
    with open(vectorSrc) as my_file:
        data = json.load(my_file)
    if(len(data['features']) == 0):
       return np.array([],dtype=np.uint8)
    else:
       return np.ones(len(data['features']),dtype=np.uint8).reshape((len(data['features']), 1))



def geoJsonToSBD(annotationName_cls, annotationName_inst, geoJson, rasterSource):

    #Print raster file name
    my_raster_source = rasterSource
    print("Raster directory : ",my_raster_source)
    srcRaster = gdal.Open(my_raster_source)

    
    my_vector_source = geoJson
    print("Vector directory : ",my_vector_source)

    
    #Call main functions to create label datafor cls
    my_cls_segmentation = createClassSegmentation(my_raster_source, my_vector_source, npDistFileName='', units='pixels')
    my_cls_boundaries =  createClassBoundaries(my_raster_source, my_vector_source, npDistFileName='', units='pixels')
    my_cls_categories = createClassCategoriesPresent(my_vector_source)

    #Call main functions to create label datafor inst
    my_inst_segmentation = createInstanceSegmentation(my_raster_source, my_vector_source)
    my_inst_boundaries = createInstanceBoundaries(my_raster_source, my_vector_source)
    my_inst_categories = createInstanceCategories(my_vector_source)

    #Wraps for cls struct
    cls_boundaries_wrap = np.array([my_cls_boundaries])
    cls_categories_wrap = my_cls_categories

    #Wraps for inst struct
    inst_boundaries_wrap = np.array([my_inst_boundaries])
    inst_categories_wrap = my_inst_categories

    #Create a class struct
    GTcls = {'Segmentation': my_cls_segmentation , 'Boundaries': cls_boundaries_wrap, 'CategoriesPresent': cls_categories_wrap}


    #Create the instance struct
    GTinst = {'Segmentation': my_inst_segmentation , 'Boundaries': inst_boundaries_wrap, 'Categories': inst_categories_wrap}

    #Save the files
    scipy.io.savemat(annotationName_cls,{'GTcls': GTcls})
    scipy.io.savemat(annotationName_inst,{'GTinst': GTinst})

    print("Done with "+str())

    entry = {'rasterFileName': my_raster_source,
             'geoJsonFileName': geoJson,
             'annotationName' : annotationName_cls,
             'annotationName_cls': annotationName_cls,
             'annotationName_inst':annotationName_inst,
             'width': srcRaster.RasterXSize,
             'height': srcRaster.RasterYSize,
             'depth' : srcRaster.RasterCount,
             'basename': os.path.splitext(os.path.basename(my_raster_source))[0]
             }

    return entry


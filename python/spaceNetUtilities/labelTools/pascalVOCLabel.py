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

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def writePacalVocObject(objectDict, top):

    childObject = SubElement(top, 'object')
    SubElement(childObject, 'name').text = objectDict['objectType']
    SubElement(childObject, 'pose').text = objectDict['pose']
    SubElement(childObject, 'truncated').text = str(objectDict['truncated'])
    SubElement(childObject, 'difficult').text = str(objectDict['difficult'])
    # write bounding box
    childBoundBox = SubElement(childObject, 'bndbox')
    SubElement(childBoundBox, 'xmin').text = str(objectDict['truncated']['bndbox']['xmin'])
    SubElement(childBoundBox, 'ymin').text = str(objectDict['truncated']['bndbox']['ymin'])
    SubElement(childBoundBox, 'xmax').text = str(objectDict['truncated']['bndbox']['xmax'])
    SubElement(childBoundBox, 'ymax').text = str(objectDict['truncated']['bndbox']['ymax'])

    return top

def writePascalVocHeader(imageDescriptionDict, top):
    ## write header


    childFolder = SubElement(top, 'folder')
    childFolder.text = imageDescriptionDict['folder']
    childFilename = SubElement(top, 'filename')
    childFilename.text = imageDescriptionDict['filename']

    # write source block
    childSource = SubElement(top, 'source')
    SubElement(childSource, 'database').text = imageDescriptionDict['source']['database']
    SubElement(childSource, 'annotation').text = imageDescriptionDict['source']['annotation']

    # write size block
    childSize = SubElement(top, 'size')
    SubElement(childSize, 'width').text = str(imageDescriptionDict['size']['width'])
    SubElement(childSize, 'height').text = str(imageDescriptionDict['size']['height'])
    SubElement(childSize, 'depth').text = str(imageDescriptionDict['size']['depth'])

    SubElement(top, 'segmented').text = str(imageDescriptionDict['segmented'])

    return top

def writeToPascalVOCLabel(xmlFilename, imageDescriptionDict, objectDictList):



    top = Element('annotation')

    top = writePascalVocHeader(imageDescriptionDict, top)

    for objectDict in objectDictList:
        top = writePacalVocObject(objectDict, top)


    with open(xmlFilename, 'w') as f:
        f.write(prettify(top))


    return xmlFilename

def geoJsonToPASCALVOC2012Label(xmlFileName, geoJson, rasterImageName, im_id='',
                           dataset ='SpaceNet',
                           folder_name='spacenet',
                           annotationStyle = 'PASCAL VOC2012',
                           segment=True,
                           bufferSizePix=2.5,
                           convertTo8Bit=True,
                           outputPixType='Byte',
                           outputFormat='GTiff',
                           bboxResize=1.0,
                           objectType='building',
                           objectTypeField=''):

    imageDescriptionDict, objectDictList = geoJsontoDict(geoJson, rasterImageName, datasetName='SpaceNet_V2',
                  annotationStyle=annotationStyle,
                  bboxResize=bboxResize,
                  objectType=objectType,
                  objectTypeField=objectTypeField,
                  objectPose='Left',
                  objectTruncatedField='',
                  objectDifficultyField=''
                  )

    xmlFileName = writeToPascalVOCLabel(xmlFileName, imageDescriptionDict, objectDictList)



def geoJsonToPASCALVOC2012_Deprecated(xmlFileName, geoJson, rasterImageName, im_id='',
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
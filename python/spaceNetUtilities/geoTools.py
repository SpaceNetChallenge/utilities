from osgeo import gdal, osr, ogr
import numpy as np
import os
import csv
import subprocess
import math
import geopandas as gpd
import shapely
from shapely.geometry import Point
from pyproj import Proj, transform
from fiona.crs import from_epsg
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.linestring import LineString
from shapely.geometry.multilinestring import MultiLineString
try:
    import rtree
    import centerline
    import osmnx
except:
    print("rtree not installed, Will break evaluation code")


def import_summary_geojson(geojsonfilename, removeNoBuildings=True):
    # driver = ogr.GetDriverByName('geojson')
    datasource = ogr.Open(geojsonfilename, 0)

    layer = datasource.GetLayer()
    print(layer.GetFeatureCount())

    buildingList = []
    for idx, feature in enumerate(layer):

        poly = feature.GetGeometryRef()

        if poly:
            if removeNoBuildings:
                if feature.GetField('BuildingId') != -1:
                    buildingList.append({'ImageId': feature.GetField('ImageId'), 'BuildingId': feature.GetField('BuildingId'),
                                  'poly': feature.GetGeometryRef().Clone()})
            else:

                buildingList.append({'ImageId': feature.GetField('ImageId'), 'BuildingId': feature.GetField('BuildingId'),
                              'poly': feature.GetGeometryRef().Clone()})

    return buildingList


def import_chip_geojson(geojsonfilename, ImageId=''):
    # driver = ogr.GetDriverByName('geojson')
    datasource = ogr.Open(geojsonfilename, 0)
    if ImageId=='':
        ImageId = geojsonfilename
    layer = datasource.GetLayer()
    print(layer.GetFeatureCount())

    polys = []
    BuildingId = 0
    for idx, feature in enumerate(layer):

        poly = feature.GetGeometryRef()

        if poly:
            BuildingId = BuildingId + 1
            polys.append({'ImageId': ImageId, 'BuildingId': BuildingId,
                          'poly': feature.GetGeometryRef().Clone()})

    return polys


def mergePolyList(geojsonfilename):

    multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
    datasource = ogr.Open(geojsonfilename, 0)

    layer = datasource.GetLayer()
    print(layer.GetFeatureCount())


    for idx, feature in enumerate(layer):

        poly = feature.GetGeometryRef()

        if poly:

            multipolygon.AddGeometry(feature.GetGeometryRef().Clone())

    return multipolygon

def readwktcsv(csv_path,removeNoBuildings=True, groundTruthFile=True):
    #
    # csv Format Expected = ['ImageId', 'BuildingId', 'PolygonWKT_Pix', 'PolygonWKT_Geo']
    # returns list of Dictionaries {'ImageId': image_id, 'BuildingId': building_id, 'poly': poly}
    # image_id is a string,
    # BuildingId is an integer,
    # poly is a ogr.Geometry Polygon

    buildinglist = []
    with open(csv_path, 'rb') as csvfile:
        building_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        next(building_reader, None)  # skip the headers
        for row in building_reader:

            if removeNoBuildings:
                if int(row[1]) != -1:
                    polyPix = ogr.CreateGeometryFromWkt(row[2])
                    if groundTruthFile:
                        polyGeo = ogr.CreateGeometryFromWkt(row[3])
                    else:
                        polyGeo = []
                    buildinglist.append({'ImageId': row[0], 'BuildingId': int(row[1]), 'polyPix': polyPix,
                                         'polyGeo': polyGeo,
                                         })

            else:

                polyPix = ogr.CreateGeometryFromWkt(row[2])
                if groundTruthFile:
                    polyGeo = ogr.CreateGeometryFromWkt(row[3])
                else:
                    polyGeo = []
                buildinglist.append({'ImageId': row[0], 'BuildingId': int(row[1]), 'polyPix': polyPix,
                                     'polyGeo': polyGeo,
                                     })

    return buildinglist


def exporttogeojson(geojsonfilename, buildinglist):
    #
    # geojsonname should end with .geojson
    # building list should be list of dictionaries
    # list of Dictionaries {'ImageId': image_id, 'BuildingId': building_id, 'polyPix': poly,
    #                       'polyGeo': poly}
    # image_id is a string,
    # BuildingId is an integer,
    # poly is a ogr.Geometry Polygon
    #
    # returns geojsonfilename

    driver = ogr.GetDriverByName('geojson')
    if os.path.exists(geojsonfilename):
        driver.DeleteDataSource(geojsonfilename)
    datasource = driver.CreateDataSource(geojsonfilename)
    layer = datasource.CreateLayer('buildings', geom_type=ogr.wkbPolygon)
    field_name = ogr.FieldDefn("ImageId", ogr.OFTString)
    field_name.SetWidth(75)
    layer.CreateField(field_name)
    layer.CreateField(ogr.FieldDefn("BuildingId", ogr.OFTInteger))

    # loop through buildings
    for building in buildinglist:
        # create feature
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField("ImageId", building['ImageId'])
        feature.SetField("BuildingId", building['BuildingId'])
        feature.SetGeometry(building['polyPix'])

        # Create the feature in the layer (geojson)
        layer.CreateFeature(feature)
        # Destroy the feature to free resources
        feature.Destroy()

    datasource.Destroy()

    return geojsonfilename


def createmaskfrompolygons(polygons):
    pass
    ## see labelTools createRasterFromGeoJson


def latlon2pixel(lat, lon, input_raster='', targetsr='', geom_transform=''):
    # type: (object, object, object, object, object) -> object

    sourcesr = osr.SpatialReference()
    sourcesr.ImportFromEPSG(4326)

    geom = ogr.Geometry(ogr.wkbPoint)
    geom.AddPoint(lon, lat)

    if targetsr == '':
        src_raster = gdal.Open(input_raster)
        targetsr = osr.SpatialReference()
        targetsr.ImportFromWkt(src_raster.GetProjectionRef())
    coord_trans = osr.CoordinateTransformation(sourcesr, targetsr)
    if geom_transform == '':
        src_raster = gdal.Open(input_raster)
        transform = src_raster.GetGeoTransform()
    else:
        transform = geom_transform

    x_origin = transform[0]
    # print(x_origin)
    y_origin = transform[3]
    # print(y_origin)
    pixel_width = transform[1]
    # print(pixel_width)
    pixel_height = transform[5]
    # print(pixel_height)
    geom.Transform(coord_trans)
    # print(geom.GetPoint())
    x_pix = (geom.GetPoint()[0] - x_origin) / pixel_width
    y_pix = (geom.GetPoint()[1] - y_origin) / pixel_height

    return (x_pix, y_pix)


def returnBoundBox(xOff, yOff, pixDim, inputRaster, targetSR='', pixelSpace=False):
    # Returns Polygon for a specific square defined by a center Pixel and
    # number of pixels in each dimension.
    # Leave targetSR as empty string '' or specify it as a osr.SpatialReference()
    # targetSR = osr.SpatialReference()
    # targetSR.ImportFromEPSG(4326)
    if targetSR == '':
        targetSR = osr.SpatialReference()
        targetSR.ImportFromEPSG(4326)
    xCord = [xOff - pixDim / 2, xOff - pixDim / 2, xOff + pixDim / 2,
             xOff + pixDim / 2, xOff - pixDim / 2]
    yCord = [yOff - pixDim / 2, yOff + pixDim / 2, yOff + pixDim / 2,
             yOff - pixDim / 2, yOff - pixDim / 2]

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for idx in xrange(len(xCord)):
        if pixelSpace == False:
            geom = pixelToGeoCoord(xCord[idx], yCord[idx], inputRaster)
            ring.AddPoint(geom[0], geom[1], 0)
        else:
            ring.AddPoint(xCord[idx], yCord[idx], 0)

    poly = ogr.Geometry(ogr.wkbPolygon)
    if pixelSpace == False:
        poly.AssignSpatialReference(targetSR)

    poly.AddGeometry(ring)

    return poly

def createBoxFromLine(tmpGeom, ratio=1, halfWidth=-999, transformRequired=True, transform_WGS84_To_UTM='', transform_UTM_To_WGS84=''):
    # create Polygon Box Oriented with the line

    if transformRequired:
        if transform_WGS84_To_UTM == '':
            transform_WGS84_To_UTM, transform_UTM_To_WGS84 = createUTMTransform(tmpGeom)

        tmpGeom.Transform(transform_WGS84_To_UTM)


    # calculatuate Centroid
    centroidX, centroidY, centroidZ = tmpGeom.Centroid().GetPoint()
    lengthM = tmpGeom.Length()
    if halfWidth ==-999:
        halfWidth = lengthM/(2*ratio)

    envelope=tmpGeom.GetPoints()
    cX1 = envelope[0][0]
    cY1 = envelope[0][1]
    cX2 = envelope[1][0]
    cY2 = envelope[1][1]
    angRad = math.atan2(cY2-cY1,cX2-cX1)

    d_X = math.cos(angRad-math.pi/2)*halfWidth
    d_Y = math.sin(angRad-math.pi/2)*halfWidth

    e_X = math.cos(angRad+math.pi/2)*halfWidth
    e_Y = math.sin(angRad+math.pi/2)*halfWidth

    ring = ogr.Geometry(ogr.wkbLinearRing)

    ring.AddPoint(cX1+d_X, cY1+d_Y)
    ring.AddPoint(cX1+e_X, cY1+e_Y)
    ring.AddPoint(cX2+e_X, cY2+e_Y)
    ring.AddPoint(cX2+d_X, cY2+d_Y)
    ring.AddPoint(cX1+d_X, cY1+d_Y)
    polyGeom = ogr.Geometry(ogr.wkbPolygon)
    polyGeom.AddGeometry(ring)
    areaM = polyGeom.GetArea()

    if transformRequired:
        tmpGeom.Transform(transform_UTM_To_WGS84)
        polyGeom.Transform(transform_UTM_To_WGS84)


    return (polyGeom, areaM, angRad, lengthM)


def pixelToGeoCoord(xPix, yPix, inputRaster, sourceSR='', geomTransform='', targetSR=''):
    # If you want to garuntee lon lat output, specify TargetSR  otherwise, geocoords will be in image geo reference
    # targetSR = osr.SpatialReference()
    # targetSR.ImportFromEPSG(4326)
    # Transform can be performed at the polygon level instead of pixel level

    if targetSR =='':
        performReprojection=False
        targetSR = osr.SpatialReference()
        targetSR.ImportFromEPSG(4326)
    else:
        performReprojection=True

    if geomTransform=='':
        srcRaster = gdal.Open(inputRaster)
        geomTransform = srcRaster.GetGeoTransform()

        source_sr = osr.SpatialReference()
        source_sr.ImportFromWkt(srcRaster.GetProjectionRef())




    geom = ogr.Geometry(ogr.wkbPoint)
    xOrigin = geomTransform[0]
    yOrigin = geomTransform[3]
    pixelWidth = geomTransform[1]
    pixelHeight = geomTransform[5]

    xCoord = (xPix * pixelWidth) + xOrigin
    yCoord = (yPix * pixelHeight) + yOrigin
    geom.AddPoint(xCoord, yCoord)


    if performReprojection:
        if sourceSR=='':
            srcRaster = gdal.Open(inputRaster)
            sourceSR = osr.SpatialReference()
            sourceSR.ImportFromWkt(srcRaster.GetProjectionRef())
        coord_trans = osr.CoordinateTransformation(sourceSR, targetSR)
        geom.Transform(coord_trans)

    return (geom.GetX(), geom.GetY())


def geoPolygonToPixelPolygonWKT(geom, inputRaster, targetSR, geomTransform, breakMultiPolygonGeo=True,
                                pixPrecision=2):
    # Returns Pixel Coordinate List and GeoCoordinateList

    polygonPixBufferList = []
    polygonPixBufferWKTList = []
    polygonGeoWKTList = []
    if geom.GetGeometryName() == 'POLYGON':
        polygonPix = ogr.Geometry(ogr.wkbPolygon)
        for ring in geom:
            # GetPoint returns a tuple not a Geometry
            ringPix = ogr.Geometry(ogr.wkbLinearRing)

            for pIdx in xrange(ring.GetPointCount()):
                lon, lat, z = ring.GetPoint(pIdx)
                xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)

                xPix = round(xPix, pixPrecision)
                yPix = round(yPix, pixPrecision)
                ringPix.AddPoint(xPix, yPix)

            polygonPix.AddGeometry(ringPix)
        polygonPixBuffer = polygonPix.Buffer(0.0)
        polygonPixBufferList.append([polygonPixBuffer, geom])

    elif geom.GetGeometryName() == 'MULTIPOLYGON':

        for poly in geom:
            polygonPix = ogr.Geometry(ogr.wkbPolygon)
            for ring in poly:
                # GetPoint returns a tuple not a Geometry
                ringPix = ogr.Geometry(ogr.wkbLinearRing)

                for pIdx in xrange(ring.GetPointCount()):
                    lon, lat, z = ring.GetPoint(pIdx)
                    xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)

                    xPix = round(xPix, pixPrecision)
                    yPix = round(yPix, pixPrecision)
                    ringPix.AddPoint(xPix, yPix)

                polygonPix.AddGeometry(ringPix)
            polygonPixBuffer = polygonPix.Buffer(0.0)
            if breakMultiPolygonGeo:
                polygonPixBufferList.append([polygonPixBuffer, poly])
            else:
                polygonPixBufferList.append([polygonPixBuffer, geom])

    for polygonTest in polygonPixBufferList:
        if polygonTest[0].GetGeometryName() == 'POLYGON':
            polygonPixBufferWKTList.append([polygonTest[0].ExportToWkt(), polygonTest[1].ExportToWkt()])
        elif polygonTest[0].GetGeometryName() == 'MULTIPOLYGON':
            for polygonTest2 in polygonTest[0]:
                polygonPixBufferWKTList.append([polygonTest2.ExportToWkt(), polygonTest[1].ExportToWkt()])

    return polygonPixBufferWKTList

def pixelWKTToGeoWKT(geomWKT, inputRaster, targetSR='', geomTransform='', breakMultiPolygonPix=False):
    # Returns  GeoCoordinateList from PixelCoordinateList



    geomPix = ogr.CreateGeometryFromWkt(geomWKT)
    geomGeoList = pixelGeomToGeoGeom(geomPix, inputRaster, targetSR=targetSR,
                                 geomTransform=geomTransform, breakMultiPolygonPix=breakMultiPolygonPix)

    return geomGeoList

def pixelGeomToGeoGeom(geom, inputRaster, targetSR='', geomTransform='', breakMultiPolygonPix=False):


    if geomTransform=='':
        targetRaster = gdal.Open(inputRaster)
        geomTransform = targetRaster.GetGeoTransform()

    polygonGeoBufferWKTList = []
    polygonGeoBufferList = []
    if geom:
        if geom.GetGeometryName() == 'POLYGON':
            polygonGeo = ogr.Geometry(ogr.wkbPolygon)
            for ring in geom:
                # GetPoint returns a tuple not a Geometry
                ringGeo = ogr.Geometry(ogr.wkbLinearRing)

                for pIdx in xrange(ring.GetPointCount()):
                    xPix, yPix, zPix = ring.GetPoint(pIdx)
                    #xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)
                    lon, lat = pixelToGeoCoord(xPix, yPix, inputRaster=inputRaster, targetSR=targetSR, geomTransform=geomTransform)

                    ringGeo.AddPoint(lon, lat)


                polygonGeo.AddGeometry(ringGeo)
            polygonGeoBuffer = polygonGeo.Buffer(0.0)
            polygonGeoBufferList.append([polygonGeoBuffer, geom])

        elif geom.GetGeometryName() == 'MULTIPOLYGON':

            for poly in geom:
                polygonGeo = ogr.Geometry(ogr.wkbPolygon)
                for ring in poly:
                    # GetPoint returns a tuple not a Geometry
                    ringGeo = ogr.Geometry(ogr.wkbLinearRing)

                    for pIdx in xrange(ring.GetPointCount()):
                        xPix, yPix, zPix = ring.GetPoint(pIdx)
                        # xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)
                        lon, lat = pixelToGeoCoord(xPix, yPix, inputRaster=inputRaster, targetSR=targetSR,
                                                   geomTransform=geomTransform)
                        ringGeo.AddPoint(lon, lat)

                    polygonGeo.AddGeometry(ringGeo)
                polygonGeoBuffer = polygonGeo.Buffer(0.0)
                if breakMultiPolygonPix:
                    polygonGeoBufferList.append([polygonGeoBuffer, poly])
                else:
                    polygonGeoBufferList.append([polygonGeoBuffer, geom])


        elif geom.GetGeometryName() == 'LINESTRING':
            lineGeo = ogr.Geometry(ogr.wkbLineString)
            for pIdx in xrange(geom.GetPointCount()):
                xPix, yPix, zPix = geom.GetPoint(pIdx)
                lon, lat = pixelToGeoCoord(xPix, yPix, inputRaster=inputRaster, targetSR=targetSR,
                                           geomTransform=geomTransform)
                lineGeo.AddPoint(lon, lat)

            polygonGeoBufferList.append([lineGeo, geom])

        elif geom.GetGeometryName() == 'POINT':
            pointGeo = ogr.Geometry(ogr.wkbPoint)

            for pIdx in xrange(geom.GetPointCount()):
                xPix, yPix, zPix = geom.GetPoint(pIdx)
                lon, lat = pixelToGeoCoord(xPix, yPix, inputRaster=inputRaster, targetSR=targetSR,
                                           geomTransform=geomTransform)
                pointGeo.AddPoint(lon, lat)

            polygonGeoBufferList.append([pointGeo, geom])






    #for polygonTest in polygonGeoBufferList:
    #    if polygonTest[0].GetGeometryName() == 'POLYGON':
    #        polygonGeoBufferWKTList.append([polygonTest[0].ExportToWkt(), polygonTest[1].ExportToWkt()])
    #    elif polygonTest[0].GetGeometryName() == 'MULTIPOLYGON':
    #        for polygonTest2 in polygonTest[0]:
    #            polygonGeoBufferWKTList.append([polygonTest2.ExportToWkt(), polygonTest[1].ExportToWkt()])



    # [GeoWKT, PixWKT])
    return polygonGeoBufferList

def geoWKTToPixelWKT(geom, inputRaster, targetSR, geomTransform, pixPrecision=2):
    # Returns Pixel Coordinate List and GeoCoordinateList

    geom_list = []
    geom_pix_wkt_list = []
    if geom.GetGeometryName() == 'POLYGON':
        polygonPix = ogr.Geometry(ogr.wkbPolygon)
        for ring in geom:
            # GetPoint returns a tuple not a Geometry
            ringPix = ogr.Geometry(ogr.wkbLinearRing)

            for pIdx in xrange(ring.GetPointCount()):
                lon, lat, z = ring.GetPoint(pIdx)
                xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)

                xPix = round(xPix, pixPrecision)
                yPix = round(yPix, pixPrecision)
                ringPix.AddPoint(xPix, yPix)

            polygonPix.AddGeometry(ringPix)
            polygonPixBuffer = polygonPix.Buffer(0.0)
            geom_list.append([polygonPixBuffer, geom])

    elif geom.GetGeometryName() == 'MULTIPOLYGON':

        for poly in geom:
            polygonPix = ogr.Geometry(ogr.wkbPolygon)
            for ring in poly:
                # GetPoint returns a tuple not a Geometry
                ringPix = ogr.Geometry(ogr.wkbLinearRing)

                for pIdx in xrange(ring.GetPointCount()):
                    lon, lat, z = ring.GetPoint(pIdx)
                    xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)

                    xPix = round(xPix, pixPrecision)
                    yPix = round(yPix, pixPrecision)
                    ringPix.AddPoint(xPix, yPix)

                polygonPix.AddGeometry(ringPix)
                polygonPixBuffer = polygonPix.Buffer(0.0)
                geom_list.append([polygonPixBuffer, geom])
    elif geom.GetGeometryName() == 'LINESTRING':
        line = ogr.Geometry(ogr.wkbLineString)
        for pIdx in xrange(geom.GetPointCount()):
            lon, lat, z = geom.GetPoint(pIdx)
            xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)

            xPix = round(xPix, pixPrecision)
            yPix = round(yPix, pixPrecision)
            line.AddPoint(xPix, yPix)
        geom_list.append([line, geom])

    elif geom.GetGeometryName() == 'POINT':
        point = ogr.Geometry(ogr.wkbPoint)
        for pIdx in xrange(geom.GetPointCount()):
            lon, lat, z = geom.GetPoint(pIdx)
            xPix, yPix = latlon2pixel(lat, lon, inputRaster, targetSR, geomTransform)

            xPix = round(xPix, pixPrecision)
            yPix = round(yPix, pixPrecision)
            point.AddPoint(xPix, yPix)
        geom_list.append([point, geom])


    for polygonTest in geom_list:
        if polygonTest[0].GetGeometryName() == 'POLYGON' or \
                        polygonTest[0].GetGeometryName() == 'LINESTRING' or \
                        polygonTest[0].GetGeometryName() == 'POINT':
            geom_pix_wkt_list.append([polygonTest[0].ExportToWkt(), polygonTest[1].ExportToWkt()])
        elif polygonTest[0].GetGeometryName() == 'MULTIPOLYGON':
            for polygonTest2 in polygonTest[0]:
                geom_pix_wkt_list.append([polygonTest2.ExportToWkt(), polygonTest[1].ExportToWkt()])

    return geom_pix_wkt_list


def convert_wgs84geojson_to_pixgeojson(wgs84geojson, inputraster, image_id=[], pixelgeojson=[], only_polygons=True,
                                       breakMultiPolygonGeo=True, pixPrecision=2,
                                       attributeName='',
                                       objectClassDict=''):
    dataSource = ogr.Open(wgs84geojson, 0)
    layer = dataSource.GetLayer()
    #print(layer.GetFeatureCount())
    building_id = 0
    # check if geoJsonisEmpty
    feautureList = []
    if not image_id:
        image_id = inputraster.replace(".tif", "")







    if layer.GetFeatureCount() > 0:

        if len(inputraster)>0:
            if os.path.isfile(inputraster):
                srcRaster = gdal.Open(inputraster)
                targetSR = osr.SpatialReference()
                targetSR.ImportFromWkt(srcRaster.GetProjectionRef())
                geomTransform = srcRaster.GetGeoTransform()





                for feature in layer:
                    if attributeName != '':
                        featureName = feature.GetField(attributeName)
                    else:
                        featureName = 'building'

                    if len(objectClassDict) > 0:
                        try:
                            featureId = objectClassDict[featureName]['featureIdNum']
                        except:
                            print('featureName {} not recognized'.format(featureName))
                    else:
                        featureId = 1


                    geom = feature.GetGeometryRef()
                    if len(inputraster)>0:
                        ## Calculate 3 band
                        if only_polygons:
                            geom_wkt_list = geoPolygonToPixelPolygonWKT(geom, inputraster, targetSR, geomTransform,
                                                                        breakMultiPolygonGeo=breakMultiPolygonGeo,
                                                                        pixPrecision=pixPrecision)
                        else:
                            geom_wkt_list = geoWKTToPixelWKT(geom, inputraster, targetSR, geomTransform,
                                                             pixPrecision=pixPrecision)

                        for geom_wkt in geom_wkt_list:
                            building_id += 1
                            feautureList.append({'ImageId': image_id,
                                                 'BuildingId': building_id,
                                                 'polyGeo': ogr.CreateGeometryFromWkt(geom_wkt[1]),
                                                 'polyPix': ogr.CreateGeometryFromWkt(geom_wkt[0]),
                                                 'featureName': featureName,
                                                 'featureIdNum': featureId
                                                 })
                    else:
                        building_id += 1
                        feautureList.append({'ImageId': image_id,
                                             'BuildingId': building_id,
                                             'polyGeo': ogr.CreateGeometryFromWkt(geom.ExportToWkt()),
                                             'polyPix': ogr.CreateGeometryFromWkt('POLYGON EMPTY'),
                                             'featureName' : featureName,
                                             'featureIdNum': featureId
                                             })
            else:
                #print("no File exists")
                pass
        if pixelgeojson:
            exporttogeojson(pixelgeojson, buildinglist=feautureList)



    return feautureList





def convert_pixgwktList_to_wgs84wktList(inputRaster, wktPolygonPixList):
    ## returns # [[GeoWKT, PixWKT], ...]
    wgs84WKTList=[]
    if os.path.isfile(inputRaster):
        srcRaster = gdal.Open(inputRaster)
        targetSR = osr.SpatialReference()
        targetSR.ImportFromWkt(srcRaster.GetProjectionRef())
        geomTransform = srcRaster.GetGeoTransform()

    for wktPolygonPix in wktPolygonPixList:
        geom_wkt_list = pixelWKTToGeoWKT(wktPolygonPix, inputRaster, targetSR='',
                                         geomTransform=geomTransform,
                                         breakMultiPolygonPix=False)

        wgs84WKTList.extend(geom_wkt_list)

    # [[GeoWKT, PixWKT], ...]
    return wgs84WKTList

def create_rtreefromdict(buildinglist):
    # create index
    index = rtree.index.Index(interleaved=False)
    for idx, building in enumerate(buildinglist):
        index.insert(idx, building['poly'].GetEnvelope())

    return index


def create_rtree_from_poly(poly_list):
    # create index
    index = rtree.index.Index(interleaved=False)
    for idx, building in enumerate(poly_list):
        index.insert(idx, building.GetEnvelope())

    return index


def search_rtree(test_building, index):
    # input test poly ogr.Geometry  and rtree index
    if test_building.GetGeometryName() == 'POLYGON' or \
                    test_building.GetGeometryName() == 'MULTIPOLYGON':
        fidlist = index.intersection(test_building.GetEnvelope())
    else:
        fidlist = []

    return fidlist


def get_envelope(poly):
    env = poly.GetEnvelope()

    # Get Envelope returns a tuple (minX, maxX, minY, maxY)
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(env[0], env[2])
    ring.AddPoint(env[0], env[3])
    ring.AddPoint(env[1], env[3])
    ring.AddPoint(env[1], env[2])
    ring.AddPoint(env[0], env[2])

    # Create polygon
    poly1 = ogr.Geometry(ogr.wkbPolygon)
    poly1.AddGeometry(ring)

    return poly1

def utm_getZone(longitude):
    return (int(1+(longitude+180.0)/6.0))


def utm_isNorthern(latitude):
    if (latitude < 0.0):
        return 0
    else:
        return 1


def createUTMTransform(polyGeom):
    # pt = polyGeom.Boundary().GetPoint()
    utm_zone = utm_getZone(polyGeom.GetEnvelope()[0])
    is_northern = utm_isNorthern(polyGeom.GetEnvelope()[2])
    utm_cs = osr.SpatialReference()
    utm_cs.SetWellKnownGeogCS('WGS84')
    utm_cs.SetUTM(utm_zone, is_northern);
    wgs84_cs = osr.SpatialReference()
    wgs84_cs.ImportFromEPSG(4326)

    transform_WGS84_To_UTM = osr.CoordinateTransformation(wgs84_cs, utm_cs)
    transform_UTM_To_WGS84 = osr.CoordinateTransformation(utm_cs, wgs84_cs)

    return transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs


def getRasterExtent(srcImage):
    geoTrans = srcImage.GetGeoTransform()
    ulX = geoTrans[0]
    ulY = geoTrans[3]
    xDist = geoTrans[1]
    yDist = geoTrans[5]
    rtnX = geoTrans[2]
    rtnY = geoTrans[4]

    cols = srcImage.RasterXSize
    rows = srcImage.RasterYSize

    lrX = ulX + xDist * cols
    lrY = ulY + yDist * rows

    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(lrX, lrY)
    ring.AddPoint(lrX, ulY)
    ring.AddPoint(ulX, ulY)
    ring.AddPoint(ulX, lrY)
    ring.AddPoint(lrX, lrY)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    return geoTrans, poly, ulX, ulY, lrX, lrY

def createPolygonFromCenterPoint(cX,cY, radiusMeters, transform_WGS_To_UTM_Flag=True):

    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(cX, cY)

    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = createUTMTransform(point)
    if transform_WGS_To_UTM_Flag:
        point.Transform(transform_WGS84_To_UTM)

    poly = point.Buffer(radiusMeters)

    if transform_WGS_To_UTM_Flag:
        poly.Transform(transform_UTM_To_WGS84)

    return poly


def createPolygonFromCorners(lrX,lrY,ulX, ulY):
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(lrX, lrY)
    ring.AddPoint(lrX, ulY)
    ring.AddPoint(ulX, ulY)
    ring.AddPoint(ulX, lrY)
    ring.AddPoint(lrX, lrY)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    return poly


def clipShapeFile(shapeSrc, outputFileName, polyToCut, minpartialPerc=0.0, shapeLabel='Geo', debug=False):

    source_layer = shapeSrc.GetLayer()
    source_srs = source_layer.GetSpatialRef()
    # Create the output Layer

    outGeoJSon = os.path.splitext(outputFileName)[0] + '.geojson'
    if not os.path.exists(os.path.dirname(outGeoJSon)):
        os.makedirs(os.path.dirname(outGeoJSon))
    print(outGeoJSon)
    outDriver = ogr.GetDriverByName("geojson")
    if os.path.exists(outGeoJSon):
        outDriver.DeleteDataSource(outGeoJSon)

    if debug:
        outGeoJSonDebug = outputFileName.replace('.tif', 'outline.geojson')
        outDriverDebug = ogr.GetDriverByName("geojson")
        if os.path.exists(outGeoJSonDebug):
            outDriverDebug.DeleteDataSource(outGeoJSonDebug)
        outDataSourceDebug = outDriver.CreateDataSource(outGeoJSonDebug)
        outLayerDebug = outDataSourceDebug.CreateLayer("groundTruth", source_srs, geom_type=ogr.wkbPolygon)

        outFeatureDebug = ogr.Feature(source_layer.GetLayerDefn())
        outFeatureDebug.SetGeometry(polyToCut)
        outLayerDebug.CreateFeature(outFeatureDebug)


    outDataSource = outDriver.CreateDataSource(outGeoJSon)
    outLayer = outDataSource.CreateLayer("groundTruth", source_srs, geom_type=ogr.wkbPolygon)
    # Add input Layer Fields to the output Layer
    inLayerDefn = source_layer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    outLayer.CreateField(ogr.FieldDefn("partialBuilding", ogr.OFTReal))
    outLayer.CreateField(ogr.FieldDefn("partialDec", ogr.OFTReal))
    outLayerDefn = outLayer.GetLayerDefn()
    source_layer.SetSpatialFilter(polyToCut)
    for inFeature in source_layer:

        outFeature = ogr.Feature(outLayerDefn)

        for i in range (0, inLayerDefn.GetFieldCount()):
            outFeature.SetField(inLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))

        geom = inFeature.GetGeometryRef()
        geomNew = geom.Intersection(polyToCut)
        partialDec = -1
        if geomNew:

            if geomNew.GetGeometryName()=='POINT':
                outFeature.SetField("partialDec", 1)
                outFeature.SetField("partialBuilding", 1)
            else:

                if geom.GetArea() > 0:
                    partialDec = geomNew.GetArea() / geom.GetArea()
                else:
                    partialDec = 0

                outFeature.SetField("partialDec", partialDec)

                if geom.GetArea() == geomNew.GetArea():
                    outFeature.SetField("partialBuilding", 0)
                else:
                    outFeature.SetField("partialBuilding", 1)


        else:
            outFeature.SetField("partialBuilding", 1)
            outFeature.SetField("partialBuilding", 1)


        outFeature.SetGeometry(geomNew)
        if partialDec >= minpartialPerc:
            outLayer.CreateFeature(outFeature)
            #print ("AddFeature")


def cutChipFromMosaic(rasterFileList, shapeFileSrcList, outlineSrc='',outputDirectory='', outputPrefix='clip_',
                      clipSizeMX=100, clipSizeMY=100, clipOverlap=0.0, minpartialPerc=0.0, createPix=False,
                      baseName='',
                      imgIdStart=-1,
                      parrallelProcess=False,
                      noBlackSpace=False,
                      randomClip=-1):
    #rasterFileList = [['rasterLocation', 'rasterDescription']]
    # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
    #                       ['/path/to/8band_AOI_1.tif, '8band']
    #                        ]
    # open Base Image
    #print(rasterFileList[0][0])
    srcImage = gdal.Open(rasterFileList[0][0])
    geoTrans, poly, ulX, ulY, lrX, lrY = getRasterExtent(srcImage)
    # geoTrans[1] w-e pixel resolution
    # geoTrans[5] n-s pixel resolution
    if outputDirectory=="":
        outputDirectory=os.path.dirname(rasterFileList[0][0])

    rasterFileBaseList = []
    for rasterFile in rasterFileList:
        rasterFileBaseList.append(os.path.basename(rasterFile[0]))

    if not createPix:
        transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = createUTMTransform(poly)
        poly.Transform(transform_WGS84_To_UTM)

    env = poly.GetEnvelope()
    minX = env[0]
    minY = env[2]
    maxX = env[1]
    maxY = env[3]

    #return poly to WGS84
    if not createPix:
        poly.Transform(transform_UTM_To_WGS84)

    shapeSrcList = []
    for shapeFileSrc in shapeFileSrcList:
        print(shapeFileSrc[1])
        shapeSrcList.append([ogr.Open(shapeFileSrc[0],0), shapeFileSrc[1]])


    if outlineSrc == '':
        geomOutline = poly
    else:
        outline = ogr.Open(outlineSrc)
        layer = outline.GetLayer()
        featureOutLine = layer.GetFeature(0)
        geomOutlineBase = featureOutLine.GetGeometryRef()
        geomOutline = geomOutlineBase.Intersection(poly)

    chipSummaryList = []

    for rasterFile in rasterFileList:
        if not os.path.exists(os.path.join(outputDirectory, rasterFile[1])):
            os.makedirs(os.path.join(outputDirectory, rasterFile[1]))
    idx = 0
    if createPix:
        print(geoTrans)
        clipSizeMX=clipSizeMX*geoTrans[1]
        clipSizeMY=abs(clipSizeMY*geoTrans[5])

    xInterval = np.arange(minX, maxX, clipSizeMX*(1.0-clipOverlap))
    print('minY = {}'.format(minY))
    print('maxY = {}'.format(maxY))
    print('clipsizeMX ={}'.format(clipSizeMX))
    print('clipsizeMY ={}'.format(clipSizeMY))

    yInterval = np.arange(minY, maxY, clipSizeMY*(1.0-clipOverlap))
    print(xInterval)
    print(yInterval)
    for llX in xInterval:
        if parrallelProcess:
            for llY in yInterval:
                pass

        else:
            for llY in yInterval:
                uRX = llX+clipSizeMX
                uRY = llY+clipSizeMY

                # check if uRX or uRY is outside image
                if noBlackSpace:
                    if uRX > maxX:
                        uRX = maxX
                        llX = maxX - clipSizeMX
                    if uRY > maxY:
                        uRY = maxY
                        llY = maxY - clipSizeMY




                polyCut = createPolygonFromCorners(llX, llY, uRX, uRY)







                if not createPix:
                    polyCut.Transform(transform_UTM_To_WGS84)
                ## add debug line do cuts
                if (polyCut).Intersects(geomOutline):
                    print("Do it.")
                    #envCut = polyCut.GetEnvelope()
                    #minXCut = envCut[0]
                    #minYCut = envCut[2]
                    #maxXCut = envCut[1]
                    #maxYCut = envCut[3]

                    #debug for real
                    minXCut = llX
                    minYCut = llY
                    maxXCut = uRX
                    maxYCut = uRY

                    #print('minYCut = {}'.format(minYCut))
                    #print('maxYCut = {}'.format(maxYCut))
                    #print('minXCut = {}'.format(minXCut))
                    #print('maxXCut = {}'.format(maxXCut))

                    #print('clipsizeMX ={}'.format(clipSizeMX))
                    #print('clipsizeMY ={}'.format(clipSizeMY))



                    idx = idx+1
                    if imgIdStart == -1:
                        imgId = -1
                    else:
                        imgId = idx

                    chipSummary = createclip(outputDirectory, rasterFileList, shapeSrcList,
                                             maxXCut, maxYCut, minYCut, minXCut,
                                             rasterFileBaseList=rasterFileBaseList,
                                             minpartialPerc=minpartialPerc,
                                             outputPrefix=outputPrefix,
                                             createPix=createPix,
                                             rasterPolyEnvelope=poly,
                                             baseName=baseName,
                                             imgId=imgId)
                    chipSummaryList.append(chipSummary)

    return chipSummaryList

def createclip(outputDirectory, rasterFileList, shapeSrcList,
               maxXCut, maxYCut, minYCut, minXCut,
               rasterFileBaseList=[],
               minpartialPerc=0,
               outputPrefix='',
               createPix=False,
               rasterPolyEnvelope=ogr.CreateGeometryFromWkt("POLYGON EMPTY"),
               className='',
               baseName='',
               imgId=-1):

    #rasterFileList = [['rasterLocation', 'rasterDescription']]
    # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
    #                       ['/path/to/8band_AOI_1.tif, '8band']
    #                        ]

    polyCutWGS = createPolygonFromCorners(minXCut, minYCut, maxXCut, maxYCut)


    if not rasterFileBaseList:
        rasterFileBaseList = []
        for rasterFile in rasterFileList:
            rasterFileBaseList.append(os.path.basename(rasterFile[0]))

    if rasterPolyEnvelope == '':
        pass

    chipNameList = []
    for rasterFile in rasterFileList:
        if className == '':
            if imgId==-1:
                chipNameList.append(outputPrefix + rasterFile[1] +
                                    "_" + baseName + "_{}_{}.tif".format(minXCut, minYCut))
            else:
                chipNameList.append(outputPrefix + rasterFile[1] +
                                    "_" + baseName + "_img{}.tif".format(imgId))
        else:
            if imgId==-1:
                chipNameList.append(outputPrefix + className + "_" +
                                    rasterFile[1] + "_" + baseName + "_{}_{}.tif".format(minXCut, minYCut))
            else:
                chipNameList.append(outputPrefix + className + '_' +
                                rasterFile[1] + "_" + baseName + "_img{}.tif".format(imgId))

    # clip raster

    for chipName, rasterFile in zip(chipNameList, rasterFileList):
        outputFileName = os.path.join(outputDirectory, rasterFile[1], className, chipName)
        ## Clip Image
        print(rasterFile)
        print(outputFileName)
        subprocess.call(["gdalwarp", "-te", "{}".format(minXCut), "{}".format(minYCut),  "{}".format(maxXCut),
                         "{}".format(maxYCut),
                         '-co', 'PHOTOMETRIC=rgb',
                         rasterFile[0], outputFileName])

    baseLayerRasterName = os.path.join(outputDirectory, rasterFileList[0][1], className, chipNameList[0])
    outputFileName = os.path.join(outputDirectory, rasterFileList[0][1], chipNameList[0])


    ### Clip poly to cust to Raster Extent
    if rasterPolyEnvelope.GetArea() == 0:
        srcImage = gdal.Open(rasterFileList[0][0])
        geoTrans, rasterPolyEnvelope, ulX, ulY, lrX, lrY = getRasterExtent(srcImage)
        polyVectorCut = polyCutWGS.Intersection(rasterPolyEnvelope)
    else:
        polyVectorCut = polyCutWGS.Intersection(rasterPolyEnvelope)

    # Interate thorough Vector Src List
    for shapeSrc in shapeSrcList:
        if imgId == -1:
            outGeoJson = outputPrefix + shapeSrc[1] + \
                         "_" + baseName + "_{}_{}.geojson".format(minXCut, minYCut)
        else:
            outGeoJson = outputPrefix + shapeSrc[1] + \
                         "_" + baseName + "_img{}.geojson".format(imgId)

        outGeoJson = os.path.join(outputDirectory, 'geojson', shapeSrc[1], outGeoJson)

        clipShapeFile(shapeSrc[0], outGeoJson, polyVectorCut, minpartialPerc=minpartialPerc)


    chipSummary = {'rasterSource': baseLayerRasterName,
                   'chipName': chipNameList[0],
                   'geoVectorName': outGeoJson,
                   'pixVectorName': ''
                   }

    return chipSummary

def cutChipFromRasterCenter(rasterFileList, shapeFileSrc, outlineSrc='',
                            outputDirectory='', outputPrefix='clip_',
                            clipSizeMeters=50, createPix=False,
                            classFieldName = 'TYPE',
                            minpartialPerc=0.1,
                            ):
    #rasterFileList = [['rasterLocation', 'rasterDescription']]
    # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
    #                       ['/path/to/8band_AOI_1.tif, '8band']
    #                        ]
    srcImage = gdal.Open(rasterFileList[0][0])
    geoTrans, poly, ulX, ulY, lrX, lrY = getRasterExtent(srcImage)

    if outputDirectory == "":
        outputDirectory = os.path.dirname(rasterFileList[0])

    rasterFileBaseList = []
    for rasterFile in rasterFileList:
        rasterFileBaseList.append(os.path.basename(rasterFile[0]))

    transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = createUTMTransform(poly)
    poly.Transform(transform_WGS84_To_UTM)
    env = poly.GetEnvelope()

    # return poly to WGS84
    poly.Transform(transform_UTM_To_WGS84)

    shapeSrc = ogr.Open(shapeFileSrc, 0)
    if outlineSrc == '':
        geomOutline = poly
    else:
        outline = ogr.Open(outlineSrc)
        layer = outline.GetLayer()
        featureOutLine = layer.GetFeature(0)
        geomOutlineBase = featureOutLine.GetGeometryRef()
        geomOutline = geomOutlineBase.Intersection(poly)

    shapeSrcBase = ogr.Open(shapeFileSrc, 0)
    layerBase = shapeSrcBase.GetLayer()
    layerBase.SetSpatialFilter(geomOutline)
    for rasterFile in rasterFileList:
        if not os.path.exists(os.path.join(outputDirectory, rasterFile[1])):
            os.makedirs(os.path.join(outputDirectory, rasterFile[1]))
    for feature in layerBase:
        featureGeom = feature.GetGeometryRef()
        cx, cy, cz = featureGeom.Centroid().GetPoint()
        polyCut = createPolygonFromCenterPoint(cx, cy, radiusMeters=clipSizeMeters)
        print(classFieldName)
        classDescription = feature.GetField(classFieldName)
        classDescription = classDescription.replace(" ","")
        envCut = polyCut.GetEnvelope()
        minXCut = envCut[0]
        minYCut = envCut[2]
        maxXCut = envCut[1]
        maxYCut = envCut[3]
        createclip(outputDirectory, rasterFileList, shapeSrc,
                       maxXCut, maxYCut, minYCut, minXCut,
                       rasterFileBaseList=rasterFileBaseList,
                       minpartialPerc=minpartialPerc,
                       outputPrefix=outputPrefix,
                       createPix=createPix,
                       rasterPolyEnvelope=poly,
                       className=classDescription)



def rotateClip(clipFileName, sourceGeoJson, rotaionList=[0,90,180,275]):
    # will add "_{}.ext".formate(rotationList[i]
    pass



def createMaskedMosaic(input_raster, output_raster, outline_file):

    subprocess.call(["gdalwarp", "-q", "-cutline", outline_file, "-of", "GTiff", input_raster, output_raster,
                     '-wo', 'OPTIMIZE_SIZE=YES',
                     '-co', 'COMPRESS=JPEG',
                     '-co', 'PHOTOMETRIC=YCBCR',
                     '-co', 'TILED=YES'])



def explodeGeoPandasFrame(inGDF):

    #This function splits entries with MultiPolygon geometries into Polygon Geometries

    outdf = gpd.GeoDataFrame(columns=inGDF.columns)
    for idx, row in inGDF.iterrows():
        if type(row.geometry) == Polygon:
            outdf = outdf.append(row,ignore_index=True)
        if type(row.geometry) == MultiPolygon:
            multdf = gpd.GeoDataFrame(columns=inGDF.columns)
            recs = len(row.geometry)
            multdf = multdf.append([row]*recs,ignore_index=True)
            for geom in range(recs):
                multdf.loc[geom,'geometry'] = row.geometry[geom]
            multdf.head()
            outdf = outdf.append(multdf,ignore_index=True)

        if type(row.geometry) == LineString:
            outdf = outdf.append(row, ignore_index=True)

        if type(row.geometry) == MultiLineString:
            multdf = gpd.GeoDataFrame(columns=inGDF.columns)
            recs = len(row.geometry)
            multdf = multdf.append([row]*recs,ignore_index=True)
            for geom in range(recs):
                multdf.loc[geom,'geometry'] = row.geometry[geom]
            multdf.head()
            outdf = outdf.append(multdf,ignore_index=True)


    outdf.crs = inGDF.crs


    return outdf

def calculateCenterLineFromGeopandasPolygon(inGDF,
                                            centerLineDistanceInput_Meters=5,
                                            simplifyDistanceMeters=5,
                                            projectToUTM=True):

    # project To UTM for GeoSpatial Measurements
    if projectToUTM:
        tmpGDF = osmnx.project_gdf(inGDF)
    else:
        tmpGDF = inGDF

    # Explode GeoPandas
    tmpGDF1 = explodeGeoPandasFrame(tmpGDF)
    tmpGDF1.crs = tmpGDF.crs
    gdf_centerline_utm = tmpGDF1


    # Loop through Geomertries to calculate Centerline for Each Polygon
    listOfGeoms = tmpGDF1['geometry'].values
    lineStringList = []

    for geom in listOfGeoms:
        tmpGeom = centerline.Centerline(geom, centerLineDistanceInput_Meters)
        lineStringList.append(tmpGeom.createCenterline())

    gdf_centerline_utm['geometry'] = lineStringList

    lineList = gdf_centerline_utm['geometry'].values
    lineSimplifiedList = []

    for geo in lineList:


        if geo.type == 'MultiLineString':

            geoNew = shapely.ops.linemerge(geo).simplify(simplifyDistanceMeters, preserve_topology=False)

        else:

            geoNew = geo.simplify(simplifyDistanceMeters, preserve_topology=False)

        lineSimplifiedList.append(geoNew)

    simplifiedGdf_utm = gpd.GeoDataFrame({'geometry': lineSimplifiedList})
    simplifiedGdf_utm.crs = tmpGDF.crs
    print (tmpGDF.crs)

    if projectToUTM:
        gdf_simple_centerline = simplifiedGdf_utm.to_crs(inGDF.crs)
    else:
        gdf_simple_centerline = simplifiedGdf_utm


    return gdf_simple_centerline


def calculateCenterLineFromOGR(inputSrcFile, centerLineDistanceInput_Meters=5, outputShpFile=''):

    inGDF = gpd.read_file(inputSrcFile)
    outGDF = calculateCenterLineFromGeopandasPolygon(inGDF, centerLineDistanceInput_Meters=centerLineDistanceInput_Meters)

    if outputShpFile != '':
        outGDF.to_file(outputShpFile)


    return outGDF


def createBufferGeoPandas(inGDF, bufferDistanceMeters=5, bufferRoundness=1, projectToUTM=True):
    # Calculate CenterLine
    ## Define Buffer Constraints


    # Transform gdf Roadlines into UTM so that Buffer makes sense
    if projectToUTM:
        tmpGDF = osmnx.project_gdf(inGDF)
    else:
        tmpGDF = inGDF

    gdf_utm_buffer = tmpGDF

    # perform Buffer to produce polygons from Line Segments
    gdf_utm_buffer['geometry'] = tmpGDF.buffer(bufferDistanceMeters,
                                                bufferRoundness)

    gdf_utm_dissolve = gdf_utm_buffer.dissolve(by='class')
    gdf_utm_dissolve.crs = gdf_utm_buffer.crs

    if projectToUTM:
        gdf_buffer = gdf_utm_dissolve.to_crs(inGDF.crs)
    else:
        gdf_buffer = gdf_utm_dissolve


    return gdf_buffer



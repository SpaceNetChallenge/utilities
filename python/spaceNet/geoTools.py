from osgeo import gdal, osr, ogr
from pandas import pandas as pd
import os
import csv
import rtree


def importgeojson(geojsonfilename, removeNoBuildings=False):
    # driver = ogr.GetDriverByName('geojson')
    datasource = ogr.Open(geojsonfilename, 0)

    layer = datasource.GetLayer()
    print(layer.GetFeatureCount())

    polys = []
    for idx, feature in enumerate(layer):

        poly = feature.GetGeometryRef()

        if poly:
            polys.append({'ImageId': feature.GetField('ImageId'), 'BuildingId': feature.GetField('BuildingId'),
                          'poly': feature.GetGeometryRef().Clone()})

    return polys


def readwktcsv(csv_path):
    #
    # csv Format Expected = ['ImageId', 'BuildingId', 'PolygonWKT_Pix', 'PolygonWKT_Geo']
    # returns list of Dictionaries {'ImageId': image_id, 'BuildingId': building_id, 'poly': poly}
    # image_id is a string,
    # BuildingId is an integer,
    # poly is a ogr.Geometry Polygon

    # buildinglist = []
    # polys_df = pd.read_csv(csv_path)
    # image_ids = set(polys_df['ImageId'].tolist())
    # for image_id in image_ids:
    #    img_df = polys_df.loc[polys_df['ImageId'] == image_id]
    #    building_ids = set(img_df['BuildingId'].tolist())
    #    for building_id in building_ids:
    #
    #            building_df = img_df.loc[img_df['BuildingId'] == building_id]
    #            poly = ogr.CreateGeometryFromWkt(building_df.iloc[0, 2])
    #            buildinglist.append({'ImageId': image_id, 'BuildingId': building_id, 'poly': poly})
    buildinglist = []
    with open(csv_path, 'rb') as csvfile:
        building_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        next(building_reader, None)  # skip the headers
        for row in building_reader:
            poly = ogr.CreateGeometryFromWkt(row[2])
            buildinglist.append({'ImageId': row[0], 'BuildingId': int(row[1]), 'poly': poly})

    return buildinglist


def exporttogeojson(geojsonfilename, buildinglist):
    #
    # geojsonname should end with .geojson
    # building list should be list of dictionaries
    # list of Dictionaries {'ImageId': image_id, 'BuildingId': building_id, 'poly': poly}
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
        feature.SetGeometry(building['poly'])

        # Create the feature in the layer (geojson)
        layer.CreateFeature(feature)
        # Destroy the feature to free resources
        feature.Destroy()

    datasource.Destroy()

    return geojsonfilename


def createmaskfrompolygons(polygons):
    pass


def latLonToPixel(lat, lon, input_raster='', targetsr='', geomTransform=''):
    sourcesr = osr.SpatialReference()
    sourcesr.ImportFromEPSG(4326)

    geom = ogr.Geometry(ogr.wkbPoint)
    geom.AddPoint(lon, lat)

    if targetsr == '':
        src_raster = gdal.Open(input_raster)
        targetsr = osr.SpatialReference()
        targetsr.ImportFromWkt(src_raster.GetProjectionRef())
    coordTrans = osr.CoordinateTransformation(sourcesr, targetsr)
    if geomTransform == '':
        src_raster = gdal.Open(input_raster)
        transform = src_raster.GetGeoTransform()
    else:
        transform = geomTransform

    xOrigin = transform[0]
    # print xOrigin
    yOrigin = transform[3]
    # print yOrigin
    pixelWidth = transform[1]
    # print pixelWidth
    pixelHeight = transform[5]
    # print pixelHeight
    geom.Transform(coordTrans)
    # print geom.GetPoint()
    xPix = (geom.GetPoint()[0] - xOrigin) / pixelWidth
    yPix = (geom.GetPoint()[1] - yOrigin) / pixelHeight

    return (xPix, yPix)


def pixelToLatLon(xPix, yPix, inputRaster, targetSR=''):
    if targetSR == '':
        targetSR = osr.SpatialReference()
        targetSR.ImportFromEPSG(4326)

    geom = ogr.Geometry(ogr.wkbPoint)
    srcRaster = gdal.Open(inputRaster)
    sourceSR = osr.SpatialReference()
    sourceSR.ImportFromWkt(srcRaster.GetProjectionRef())
    coordTrans = osr.CoordinateTransformation(sourceSR, targetSR)

    transform = srcRaster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]

    xCoord = (xPix * pixelWidth) + xOrigin
    yCoord = (yPix * pixelHeight) + yOrigin
    geom.AddPoint(xCoord, yCoord)
    geom.Transform(coordTrans)
    return (geom.GetX(), geom.GetY())


def geoPolygonToPixelPolygonWKT(geom, inputRaster, targetSR, geomTransform):
    # Returns Pixel Coordinate List and GeoCoordinateList

    polygonPixBufferList = []
    polygonPixBufferWKTList = []
    if geom.GetGeometryName() == 'POLYGON':
        polygonPix = ogr.Geometry(ogr.wkbPolygon)
        for ring in geom:
            # GetPoint returns a tuple not a Geometry
            ringPix = ogr.Geometry(ogr.wkbLinearRing)

            for pIdx in xrange(ring.GetPointCount()):
                lon, lat, z = ring.GetPoint(pIdx)
                xPix, yPix = latLonToPixel(lat, lon, inputRaster, targetSR, geomTransform)
                ringPix.AddPoint(xPix, yPix)

            polygonPix.AddGeometry(ringPix)
            polygonPixBuffer = polygonPix.Buffer(0.0)
            polygonPixBufferList.append([polygonPixBuffer, geom])

    elif geom.GetGeometryName() == 'MULTIPOLYGON':

        for poly in geom:
            polygonPix = ogr.Geometry(ogr.wkbPolygon)
            for ring in geom:
                # GetPoint returns a tuple not a Geometry
                ringPix = ogr.Geometry(ogr.wkbLinearRing)

                for pIdx in xrange(ring.GetPointCount()):
                    lon, lat, z = ring.GetPoint(pIdx)
                    xPix, yPix = latLonToPixel(lat, lon, inputRaster, targetSR, geomTransform)
                    ringPix.AddPoint(xPix, yPix)

                polygonPix.AddGeometry(ringPix)
                polygonPixBuffer = polygonPix.Buffer(0.0)
                polygonPixBufferList.append([polygonPixBuffer, geom])

    for polygonTest in polygonPixBufferList:
        if polygonTest[0].GetGeometryName() == 'POLYGON':
            polygonPixBufferWKTList.append([polygonTest[0].ExportToWkt(), polygonTest[1].ExportToWkt()])
        elif polygonTest[0].GetGeometryName() == 'MULTIPOLYGON':
            for polygonTest2 in polygonTest[0]:
                polygonPixBufferWKTList.append([polygonTest2.ExportToWkt(), polygonTest[1].ExportToWkt()])

    return polygonPixBufferWKTList


def convert_wgs84geojson_to_pixgeojson(wgs84geojson, inputraster, image_id=[], pixelgeojson=[]):
    dataSource = ogr.Open(wgs84geojson, 0)
    layer = dataSource.GetLayer()
    print(layer.GetFeatureCount())
    building_id = 0
    # check if geoJsonisEmpty
    buildinglist = []
    if not image_id:
        image_id = inputraster.replace(".tif", "")

    if layer.GetFeatureCount() > 0:
        srcRaster = gdal.Open(inputraster)
        targetSR = osr.SpatialReference()
        targetSR.ImportFromWkt(srcRaster.GetProjectionRef())
        geomTransform = srcRaster.GetGeoTransform()

        for feature in layer:

            geom = feature.GetGeometryRef()

            ## Calculate 3 band
            polygonWKTList = geoPolygonToPixelPolygonWKT(geom, inputraster, targetSR, geomTransform)

            for polygonWKT in polygonWKTList:
                building_id += 1
                buildinglist.append({'ImageId': image_id,
                                     'BuildingId': building_id,
                                     'poly': ogr.CreateGeometryFromWkt(polygonWKT)})

    if pixelgeojson:
        exporttogeojson(pixelToLatLon, buildinglist=buildinglist)

    return buildinglist


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


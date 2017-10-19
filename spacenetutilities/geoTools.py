import numpy as np
import os
import subprocess
import math
import geopandas as gpd
import shapely
import rasterio as rio
import affine as af
import pandas as pd
from shapely.geometry import Point
import pyproj
import fiona
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.linestring import LineString
from shapely.geometry.multilinestring import MultiLineString
from shapely.geometry import shape
from tqdm import tqdm
import rtree
from functools import partial

#try:
#    import centerline
#
#except:
#    print("rtree not installed, Will break evaluation code")


def import_summary_geojson(geojsonfilename, removeNoBuildings=True):
    """read summary spacenetV2 geojson into geopandas dataFrame.

       Keyword arguments:
       geojsonfilename -- geojson to read
       removeNoBuildings -- remove all samples with BuildingId == -1 (default =True)
    """


    buildingList_df = gpd.read_file(geojsonfilename)


    if removeNoBuildings:
        buildingList_df = buildingList_df[buildingList_df['BuildingId']!=-1]

    buildingList_df['poly'] = buildingList_df.geometry

    return buildingList_df


def import_chip_geojson(geojsonfilename, ImageId=''):
    """read spacenetV2 chip geojson into geopandas dataFrame.

           Keyword arguments:
           geojsonfilename -- geojson to read
           ImageId -- Specify ImageId.  If not specified. ImageId is defined by
            os.path.splitext(os.path.basename(geojsonfilename))[0]
    """

    buildingList_df = gpd.read_file(geojsonfilename)


    if ImageId=='':
        ImageId = os.path.splitext(os.path.basename(geojsonfilename))[0]

    buildingList_df['ImageId']=ImageId
    buildingList_df['BuildingId'] = range(1, len(buildingList_df) + 1)
    buildingList_df['poly']       = buildingList_df.geometry #[shapely.wkt.loads(x) for x in buildingList_df.geometry.values]

    return buildingList_df


def mergePolyList(geojsonfilename):
    """read geoJson and return dataframe of unary_union

           Keyword arguments:
           geojsonfilename -- geojson to read

    """

    buildingList_df = gpd.read_file(geojsonfilename)

    return buildingList_df.unary_union


def readwktcsv(csv_path):
    """read spacenetV2 csv and return geopandas dataframe

               Keyword arguments:

               csv_path -- path to csv of spacenetV2 ground truth or solution submission format 
                    csv Format Expected = ['ImageId', 'BuildingId', 'PolygonWKT_Pix', 'PolygonWKT_Geo'] or
                    csv Format Expected = ['ImageId', 'BuildingId', 'PolygonWKT', 'Confidence']

            see https://community.topcoder.com/longcontest/?module=ViewProblemStatement&rd=16892&pm=14551 to 
            learn more about the spacenetV2 csv formats   
    """
    #


    df = pd.read_csv(csv_path)
    crs = {}
    if 'PolygonWKT_Geo' in df.columns:
        geometry = [shapely.wkt.loads(x) for x in df['PolygonWKT_Geo'].values]
        crs = {'init': 'epsg:4326'}
    elif 'PolygonWKT_Pix' in df.columns:
        geometry = [shapely.wkt.loads(x) for x in df['PolygonWKT_Pix'].values]
    elif 'PolygonWKT' in df.columns:
        geometry = [shapely.wkt.loads(x) for x in df['PolygonWKT'].values]

    else:
        print(
            'Eror No Geometry Column detected, column must be called "PolygonWKT_Geo", "PolygonWKT_Pix", or "PolygonWKT"')
        return -1

    geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)

    return geo_df


def exporttogeojson(geojsonfilename, geo_df):
    """Write geopandas dataframe to geo_df 

           Keyword arguments:
           geojsonfilename -- geojson to create
           geo_df          -- geopandas dataframe

    """

    #geo_df.to_file(geojsonfilename, driver='GeoJSON', crs=from_epsg(4326))
    geo_df.to_file(geojsonfilename, driver='GeoJSON')

    return geojsonfilename


def geomGeo2geomPixel(geom, affineObject=[], input_raster='', gdal_geomTransform=[]):
    # This function transforms a shapely geometry in geospatial coordinates into pixel coordinates
    # geom must be shapely geometry
    # affineObject = rasterio.open(input_raster).affine
    # gdal_geomTransform = gdal.Open(input_raster).GetGeoTransform()
    # input_raster is path to raster to gather georectifcation information
    if not affineObject:
        if input_raster != '':
            affineObject = rio.open(input_raster).affine
        else:
            affineObject = af.Affine.from_gdal(gdal_geomTransform)

    affineObjectInv = ~affineObject

    geomTransform = shapely.affinity.affine_transform(geom,
                                      [affineObjectInv.a,
                                       affineObjectInv.b,
                                       affineObjectInv.d,
                                       affineObjectInv.e,
                                       affineObjectInv.xoff,
                                       affineObjectInv.yoff]
                                      )

    return geomTransform

def geomPixel2geomGeo(geom, affineObject=[], input_raster='', gdal_geomTransform=[]):
    # This function transforms a shapely geometry in pixel coordinates into geospatial coordinates
    # geom must be shapely geometry
    # affineObject = rasterio.open(input_raster).affine
    # gdal_geomTransform = gdal.Open(input_raster).GetGeoTransform()
    # input_raster is path to raster to gather georectifcation information
    if not affineObject:
        if input_raster != '':
            affineObject = rio.open(input_raster).affine
        else:
            affineObject = af.Affine.from_gdal(gdal_geomTransform)


    geomTransform = shapely.affinity.affine_transform(geom,
                                                      [affineObject.a,
                                                       affineObject.b,
                                                       affineObject.d,
                                                       affineObject.e,
                                                       affineObject.xoff,
                                                       affineObject.yoff]
                                                      )

    return geomTransform

def geoDFtoPixDF(geoDF, affineObject=[], input_raster='', gdal_geomTransform=[]):

    geomList =[]
    for geom in geoDF.geometry.values:
        geomList.append(geomGeo2geomPixel(geom,
                                          affineObject=affineObject,
                                          input_raster=input_raster,
                                          gdal_geomTransform=gdal_geomTransform)
                        )

    pixDF = geoDF.copy()
    pixDF['geometry']=geomList


    return pixDF


def pixDFtoGeoDF(pixDF, affineObject=[], input_raster='', gdal_geomTransform=[]):
    geomList = []
    for geom in pixDF.geometry.values:
        geomList.append(geomPixel2geomGeo(geom,
                                          affineObject=affineObject,
                                          input_raster=input_raster,
                                          gdal_geomTransform=gdal_geomTransform)
                        )

    geoDF = pixDF.copy()
    geoDF['geometry'] = geomList

    return geoDF


def returnBoundBox(xCenter, yCenter, pixDim, affineObject=[], input_raster='', gdal_geomTransform=[], pixelSpace=False):

    geom = Point(xCenter, yCenter)
    geom = geom.buffer(pixDim)
    geom = geom.envelope

    if not pixelSpace:
        geom = geomPixel2geomGeo(geom, affineObject=affineObject, input_raster=input_raster, gdal_geomTransform=gdal_geomTransform)
    else:
        geom

    return geom

def createBoxFromLine(tmpGeom, ratio=1, halfWidth=-999, transformRequired=True, transform_WGS84_To_UTM='', transform_UTM_To_WGS84=''):
    # create Polygon Box Oriented with the line

    if transformRequired:
        if transform_WGS84_To_UTM == '':
            transform_WGS84_To_UTM, transform_UTM_To_WGS84 = createUTMTransform(tmpGeom)

        tmpGeom = shapely.ops.transform(transform_WGS84_To_UTM, tmpGeom)




    # calculatuate Centroid

    lengthM = tmpGeom.length
    if halfWidth ==-999:
        halfWidth = lengthM/(2*ratio)

    polyGeom = tmpGeom.buffer(halfWidth, cap_style=shapely.geometry.CAP_STYLE.flat)

    angRad = math.atan2((tmpGeom.coords[1][1]-tmpGeom.coords[0][1],
                        tmpGeom.coords[1][0] - tmpGeom.coords[0][0]))

    areaM = polyGeom.area

    if transformRequired:
        polyGeom = shapely.ops.transform(transform_UTM_To_WGS84, polyGeom)



    return (polyGeom, areaM, angRad, lengthM)


def create_rtreefromdict(buildinglist):
    # create index
    index = rtree.index.Index(interleaved=True)
    for idx, building in enumerate(buildinglist):
        index.insert(idx, building['poly'].bounds)

    return index


def create_rtree_from_poly(poly_list):
    # create index
    index = rtree.index.Index(interleaved=True)
    for idx, building in enumerate(poly_list):
        index.insert(idx, building.bounds)

    return index


def search_rtree(test_building, index):
    # input test poly ogr.Geometry  and rtree index
    if test_building.GetGeometryName() == 'POLYGON' or \
                    test_building.GetGeometryName() == 'MULTIPOLYGON':
        fidlist = index.intersection(test_building.bounds)
    else:
        fidlist = []

    return fidlist


def get_envelope(poly):

    return poly.envelope

def utm_getZone(longitude):
    return (int(1+(longitude+180.0)/6.0))


def utm_isNorthern(latitude):
    if (latitude < 0.0):
        return 0
    else:
        return 1


def createUTMandLatLonCrs(polyGeom):
    print("does not work below equator yet")
    polyCentroid = polyGeom.centroid
    utm_zone = utm_getZone(polyCentroid.x)
    is_northern = utm_isNorthern(polyCentroid.y)
    if is_northern:
        directionIndicator = '+north'
    else:
        directionIndicator = '+south'

    utm_crs = {'datum': 'NAD83',
               'ellps': 'GRS80',
               'proj': 'utm',
               'zone': utm_zone,
               'units': 'm'}

    latlong_crs = {'init': 'epsg:4326'}

    return utm_crs, latlong_crs

def createUTMTransform(polyGeom):

    polyCentroid = polyGeom.centroid
    utm_zone = utm_getZone(polyCentroid.x)
    is_northern = utm_isNorthern(polyCentroid.y)
    if is_northern:
        directionIndicator = '+north'
    else:
        directionIndicator = '+south'

    print('utm zone = {}'.format(utm_zone))
    projectTO_UTM = partial(
        pyproj.transform,
        pyproj.Proj("+proj=longlat +datum=WGS84 +no_defs"),  #Proj(proj='latlong',datum='WGS84')
        pyproj.Proj("+proj=utm +zone={} {} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(utm_zone,
                                                                                               directionIndicator))
    )


    projectTO_WGS = partial(
        pyproj.transform,
        pyproj.Proj("+proj=utm +zone={} {} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(utm_zone,
                                                                                               directionIndicator)
                    ),
        pyproj.Proj("+proj=longlat +datum=WGS84 +no_defs"),  # Proj(proj='latlong',datum='WGS84')

    )
    utm_cs = "+proj=utm +zone={} {} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(utm_zone,
                                                                                               directionIndicator)

    return projectTO_UTM,  projectTO_WGS, utm_cs

def transformGeomToUTM(geom):
    transform_WGS84_To_UTM, transform_UTM_To_WGS84 = createUTMTransform(geom)

    return shapely.ops.transform(transform_WGS84_To_UTM, geom)


def getRasterExtent(srcImage):
    'returns srcImage.transform which is an Affine Object'


    poly = Polygon(((srcImage.bounds.left, srcImage.bounds.top),
                   (srcImage.bounds.right, srcImage.bounds.top),
                   (srcImage.bounds.right, srcImage.bounds.bottom),
                   (srcImage.bounds.left, srcImage.bounds.bottom),
                    (srcImage.bounds.left, srcImage.bounds.top))
                    )



    return srcImage.transform, \
           poly, \
           srcImage.bounds.left, \
           srcImage.bounds.top, \
           srcImage.bounds.right, \
           srcImage.bounds.bottom

def createPolygonFromCenterPointXY(cX,cY, radiusMeters, transform_WGS_To_UTM_Flag=True):


    point = Point(cX, cY)

    return createPolygonFromCenterPoint(point, radiusMeters, transform_WGS_To_UTM_Flag=True)

def createPolygonFromCenterPoint(point, radiusMeters, transform_WGS_To_UTM_Flag=True):



    if transform_WGS_To_UTM_Flag:
        transform_WGS84_To_UTM, transform_UTM_To_WGS84 = createUTMTransform(point)
        point = shapely.ops.transform(transform_WGS84_To_UTM, point)

    poly = point.buffer(radiusMeters)

    if transform_WGS_To_UTM_Flag:
        poly = shapely.ops.transform(transform_UTM_To_WGS84, poly)

    return poly.envelope

def createPolygonFromCentroidGDF(gdf, radiusMeters, transform_WGS_To_UTM_Flag=True):

    #TODO needs fixing
    if transform_WGS_To_UTM_Flag:
        transform_WGS84_To_UTM, transform_UTM_To_WGS84 = createUTMTransform(gdf.centroid.values[0])
        gdf.to_crs(transform_WGS84_To_UTM)

    poly = gdf.centroids.buffer(radiusMeters)

    if transform_WGS_To_UTM_Flag:
        poly = shapely.ops.transform(transform_UTM_To_WGS84, poly)

    return poly

def createPolygonFromCorners(left,bottom,right, top):
    # Create ring
    poly = Polygon(
        (
            (left, top),
            (left, bottom),
            (right, bottom),
            (right, top),
            (left, top)
        )
    )

    return poly


def clipShapeFile(geoDF, outputFileName, polyToCut, minpartialPerc=0.0, shapeLabel='Geo', debug=False):
    # check if geoDF has origAreaField
    outGeoJSon = os.path.splitext(outputFileName)[0] + '.geojson'
    if not os.path.exists(os.path.dirname(outGeoJSon)):
        os.makedirs(os.path.dirname(outGeoJSon))
    if debug:
        print(outGeoJSon)

    if 'origarea' in geoDF.columns:
        pass
    else:
        geoDF['origarea'] = geoDF.area

    #TODO must implement different case for lines and for spatialIndex

    cutGeoDF = geoDF.copy()
    cutGeoDF.geometry=geoDF.intersection(polyToCut)
    cutGeoDF['partialDec'] = cutGeoDF.area / cutGeoDF['origarea']
    cutGeoDF = cutGeoDF.loc[cutGeoDF['partialDec'] > minpartialPerc].copy()
    #cutGeoDF = geoDF.loc[geoDF.intersection(polyToCut).area/geoDF['origarea'] > minpartialPerc].copy()
    cutGeoDF['truncated'] = (cutGeoDF['partialDec']!=1.0).astype(int)

    if cutGeoDF.empty:
        with open(outGeoJSon, 'a'):
            os.utime(outGeoJSon, None)
    else:
    ##TODO Verify this works in DockerBuild
        cutGeoDF.to_file(outGeoJSon, driver='GeoJSON')


def cutChipFromMosaic(rasterFileList, shapeFileSrcList, outlineSrc='',outputDirectory='', outputPrefix='clip_',
                      clipSizeMX=100, clipSizeMY=100, clipOverlap=0.0, minpartialPerc=0.0, createPix=False,
                      baseName='',
                      imgIdStart=-1,
                      parrallelProcess=False,
                      noBlackSpace=False,
                      randomClip=-1,
                      verbose=False):

    #rasterFileList = [['rasterLocation', 'rasterDescription']]
    # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
    #                       ['/path/to/8band_AOI_1.tif, '8band']
    #                        ]
    # open Base Image
    #print(rasterFileList[0][0])
    srcImage = rio.open(rasterFileList[0][0])
    geoTrans, poly, ulX, ulY, lrX, lrY = getRasterExtent(srcImage)
    # geoTrans.a w-e pixel resolution
    # geoTrans.e n-s pixel resolution
    if outputDirectory=="":
        outputDirectory=os.path.dirname(rasterFileList[0][0])

    rasterFileBaseList = []
    for rasterFile in rasterFileList:
        rasterFileBaseList.append(os.path.basename(rasterFile[0]))

    if not createPix:
        transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = createUTMTransform(poly)
        poly = shapely.ops.transform(transform_WGS84_To_UTM, poly)


    env = poly.bounds
    minX = env[0]
    minY = env[1]
    maxX = env[2]
    maxY = env[3]

    #return poly to WGS84
    if not createPix:
        poly = shapely.ops.transform(transform_UTM_To_WGS84, poly)

    shapeSrcList = []
    for shapeFileSrc in shapeFileSrcList:
        print(shapeFileSrc[1])
        shapeSrcList.append([gpd.read_file(shapeFileSrc[0]), shapeFileSrc[1]])


    if outlineSrc == '':
        geomOutline = poly
    else:
        with fiona.open(outlineSrc) as src:
            outline = src.next()
            geomOutlineBase = shape(outline['geometry'])
            geomOutline = geomOutlineBase.intersection(poly)

    chipSummaryList = []

    for rasterFile in rasterFileList:
        if not os.path.exists(os.path.join(outputDirectory, rasterFile[1])):
            os.makedirs(os.path.join(outputDirectory, rasterFile[1]))
    idx = 0
    if createPix:
        print(geoTrans)
        clipSizeMX=clipSizeMX*geoTrans.a
        clipSizeMY=abs(clipSizeMY*geoTrans.e)

    xInterval = np.arange(minX, maxX, clipSizeMX*(1.0-clipOverlap))
    yInterval = np.arange(minY, maxY, clipSizeMY*(1.0-clipOverlap))

    if verbose:
        print('minY = {}'.format(minY))
        print('maxY = {}'.format(maxY))
        print('clipsizeMX ={}'.format(clipSizeMX))
        print('clipsizeMY ={}'.format(clipSizeMY))
        print(xInterval)
        print(yInterval)

    pbar = tqdm(total=len(xInterval)* len(yInterval), desc='Creating Chips')




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
                    polyCut = shapely.ops.transform(transform_UTM_To_WGS84, polyCut)

                ## add debug line do cuts
                if polyCut.intersects(geomOutline):
                    if verbose:
                        print("Do it.")
                        print('minYCut = {}'.format(minYCut))
                        print('maxYCut = {}'.format(maxYCut))
                        print('minXCut = {}'.format(minXCut))
                        print('maxXCut = {}'.format(maxXCut))

                        print('clipsizeMX ={}'.format(clipSizeMX))
                        print('clipsizeMY ={}'.format(clipSizeMY))


                    minXCut = llX
                    minYCut = llY
                    maxXCut = uRX
                    maxYCut = uRY





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

                    pbar.update(1)

    pbar.close()

    return chipSummaryList

def createclip(outputDirectory, rasterFileList, shapeSrcList,
               maxXCut, maxYCut, minYCut, minXCut,
               rasterFileBaseList=[],
               minpartialPerc=0,
               outputPrefix='',
               createPix=False,
               rasterPolyEnvelope=Point(0,0),
               className='',
               baseName='',
               imgId=-1,
               s3Options=[],
               verbose=False):

    #rasterFileList = [['rasterLocation', 'rasterDescription']]
    # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
    #                       ['/path/to/8band_AOI_1.tif, '8band']
    #                        ]

    ## Create Polygon of area to Cut
    polyCutWGS = createPolygonFromCorners(minXCut, minYCut, maxXCut, maxYCut)

    # create rasterFile BaseName for later Copying
    if not rasterFileBaseList:
        rasterFileBaseList = []
        for rasterFile in rasterFileList:
            rasterFileBaseList.append(os.path.basename(rasterFile[0]))

    if rasterPolyEnvelope == '':
        #set rasterPolyEnvelope to point to indicate it has no area and rasterExtent should be used
        rasterPolyEnvelope=Point(0,0)

    # Generate chipName List depending on type of image
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
        if verbose:
            print(rasterFile)
            print(outputFileName)
        #TODO replace gdalwarp with rasterio and windowed reads

        cmd = ["gdalwarp", "-te", "{}".format(minXCut), "{}".format(minYCut),  "{}".format(maxXCut),
                         "{}".format(maxYCut),
                         '-co', 'PHOTOMETRIC=rgb',
                         rasterFile[0], outputFileName]
        cmd.extend(s3Options)
        subprocess.call(cmd)

    baseLayerRasterName = os.path.join(outputDirectory, rasterFileList[0][1], className, chipNameList[0])
    outputFileName = os.path.join(outputDirectory, rasterFileList[0][1], chipNameList[0])


    ### Clip poly to Raster Extent
    if rasterPolyEnvelope.area == 0:
        srcImage = rio.open(rasterFileList[0][0])
        geoTrans, rasterPolyEnvelope, ulX, ulY, lrX, lrY = getRasterExtent(srcImage)
        if verbose:
            print("rasterPolyEnvelope={}".format(rasterPolyEnvelope.wkt))
            print("polyCutWGS={}".format(polyCutWGS.wkt))

    else:

        polyVectorCut = polyCutWGS.intersection(rasterPolyEnvelope)

        if verbose:
            print("rasterPolyEnvelope={}".format(rasterPolyEnvelope))
            print("polyCutWGS={}".format(polyCutWGS))
            print('polyVectorCut area={}'.format(polyVectorCut.area))
            print("polyVectorCut={}".format(polyVectorCut))

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



def cutChipFromRasterCenter(rasterFileList, shapeFileSrc, shapeFileSrcList, outlineSrc='',
                            outputDirectory='', outputPrefix='clip_',
                            clipSizeMeters=50, createPix=False,
                            classFieldName='TYPE',
                            minpartialPerc=0.1,
                            preciseMatch=False,
                            verbose=False,
                            baseName=''
                            ):
    #rasterFileList = [['rasterLocation', 'rasterDescription']]
    # i.e rasterFileList = [['/path/to/3band_AOI_1.tif, '3band'],
    #                       ['/path/to/8band_AOI_1.tif, '8band']
    #                        ]
    srcImage = rio.open(rasterFileList[0][0])
    geoTrans, poly, ulX, ulY, lrX, lrY = getRasterExtent(srcImage)

    if outputDirectory == "":
        outputDirectory = os.path.dirname(rasterFileList[0])

    rasterFileBaseList = []
    for rasterFile in rasterFileList:
        rasterFileBaseList.append(os.path.basename(rasterFile[0]))
    # Check if Transform Is Neccessary
    #transform_WGS84_To_UTM, transform_UTM_To_WGS84, utm_cs = createUTMTransform(poly)
    #poly = shapely.ops.transform(transform_WGS84_To_UTM, poly)





    if outlineSrc == '':
        geomOutline = poly
        print('geomOutline poly={}'.format(geomOutline.wkt))

    else:
        with fiona.open(outlineSrc) as src:
            outline = src.next()
            geomOutlineBase = shape(outline['geometry'])
            geomOutline = geomOutlineBase.intersection(poly)

    # use sindex function of geoPandas for filter
    shapeSrcBaseGDF = shapeFileSrc
    if shapeFileSrc.crs['init'] == 'epsg:4326':
        utmConversionNeeded = True
    else:
        utmConversionNeeded = False

    print("srcFileShape = {}".format(shapeSrcBaseGDF.shape))
    shapeSrcBaseIndex = shapeSrcBaseGDF.sindex
    possible_matches_index = list(shapeSrcBaseIndex.intersection(geomOutline.bounds))
    possible_matches = shapeSrcBaseGDF.iloc[possible_matches_index]
    print("PossibleMatches = {}".format(possible_matches.shape))
    #if geomOutline is very irregular, preciseMatch will be neccessary
    if preciseMatch:
        layerBase = possible_matches[possible_matches.intersects(geomOutline)]
    else:
        layerBase = possible_matches



    for rasterFile in rasterFileList:
        if not os.path.exists(os.path.join(outputDirectory, rasterFile[1])):
            os.makedirs(os.path.join(outputDirectory, rasterFile[1]))

    shapeSrcList = []
    for shapeFileSrc in shapeFileSrcList:
        print(shapeFileSrc[1])
        shapeSrcList.append([gpd.read_file(shapeFileSrc[0]), shapeFileSrc[1]])

    chipSummaryList = []
    for idx, feature in tqdm(layerBase.iterrows(), desc="Processing Features"):
        featureGeom = feature['geometry']
        #cx = featureGeom.centroid.x
        #cy = featureGeom.centroid.y

        polyCut = createPolygonFromCenterPoint(featureGeom.centroid, radiusMeters=clipSizeMeters, transform_WGS_To_UTM_Flag=utmConversionNeeded)

        if verbose:
            print(classFieldName)
        if classFieldName!='':
            if classFieldName in feature:
                classDescription = feature[classFieldName]
            else:
                classDescription = 'unknown'

            for rasterFile in rasterFileList:
                if not os.path.exists(os.path.join(outputDirectory, rasterFile[1], classDescription)):
                    os.makedirs(os.path.join(outputDirectory, rasterFile[1], classDescription))
        else:
            classDescription=''
        #Eliminate WhiteSpace
        classDescription = classDescription.replace(" ","")

        minXCut, minYCut, maxXCut, maxYCut = polyCut.bounds

        chipSummary = createclip(outputDirectory, rasterFileList, shapeSrcList,
                       maxXCut, maxYCut, minYCut, minXCut,
                       rasterFileBaseList=rasterFileBaseList,
                       minpartialPerc=minpartialPerc,
                       outputPrefix=outputPrefix,
                       createPix=createPix,
                       rasterPolyEnvelope=poly,
                       className=classDescription,
                       baseName=baseName
                                 )

        chipSummaryList.append(chipSummary)


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

# def calculateCenterLineFromGeopandasPolygon(inGDF,
#                                             centerLineDistanceInput_Meters=5,
#                                             simplifyDistanceMeters=5,
#                                             projectToUTM=True):
#
#     # project To UTM for GeoSpatial Measurements
#     if projectToUTM:
#         tmpGDF = osmnx.project_gdf(inGDF)
#     else:
#         tmpGDF = inGDF
#
#     # Explode GeoPandas
#     tmpGDF1 = explodeGeoPandasFrame(tmpGDF)
#     tmpGDF1.crs = tmpGDF.crs
#     gdf_centerline_utm = tmpGDF1
#
#
#     # Loop through Geomertries to calculate Centerline for Each Polygon
#     listOfGeoms = tmpGDF1['geometry'].values
#     lineStringList = []
#
#     for geom in listOfGeoms:
#         tmpGeom = centerline.Centerline(geom, centerLineDistanceInput_Meters)
#         lineStringList.append(tmpGeom.createCenterline())
#
#     gdf_centerline_utm['geometry'] = lineStringList
#
#     lineList = gdf_centerline_utm['geometry'].values
#     lineSimplifiedList = []
#
#     for geo in lineList:
#
#
#         if geo.type == 'MultiLineString':
#
#             geoNew = shapely.ops.linemerge(geo).simplify(simplifyDistanceMeters, preserve_topology=False)
#
#         else:
#
#             geoNew = geo.simplify(simplifyDistanceMeters, preserve_topology=False)
#
#         lineSimplifiedList.append(geoNew)
#
#     simplifiedGdf_utm = gpd.GeoDataFrame({'geometry': lineSimplifiedList})
#     simplifiedGdf_utm.crs = tmpGDF.crs
#     print (tmpGDF.crs)
#
#     if projectToUTM:
#         gdf_simple_centerline = simplifiedGdf_utm.to_crs(inGDF.crs)
#     else:
#         gdf_simple_centerline = simplifiedGdf_utm
#
#
#     return gdf_simple_centerline
#
#
# def calculateCenterLineFromOGR(inputSrcFile, centerLineDistanceInput_Meters=5, outputShpFile=''):
#
#     inGDF = gpd.read_file(inputSrcFile)
#     outGDF = calculateCenterLineFromGeopandasPolygon(inGDF, centerLineDistanceInput_Meters=centerLineDistanceInput_Meters)
#
#     if outputShpFile != '':
#         outGDF.to_file(outputShpFile)
#
#
#     return outGDF


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





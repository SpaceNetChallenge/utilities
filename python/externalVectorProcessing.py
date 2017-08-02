from spaceNetUtilities import geoTools as gT
from spaceNetUtilities import labelTools as lT
from osgeo import gdal, ogr, osr
import argparse
import os
import subprocess
import glob



def buildTindex(rasterFolder, rasterExtention='.tif'):
    rasterList = glob.glob(os.path.join(rasterFolder, '*{}'.format(rasterExtention)))
    print(rasterList)

    print(os.path.join(rasterFolder, '*{}'.format(rasterExtention)))

    memDriver = ogr.GetDriverByName('MEMORY')
    gTindex = memDriver.CreateDataSource('gTindex')
    srcImage = gdal.Open(rasterList[0])
    spat_ref = osr.SpatialReference()
    spat_ref.SetProjection(srcImage.GetProjection())
    gTindexLayer = gTindex.CreateLayer("gtindexlayer", spat_ref, geom_type=ogr.wkbPolygon)

    # Add an ID field
    idField = ogr.FieldDefn("location", ogr.OFTString)
    gTindexLayer.CreateField(idField)

    # Create the feature and set values
    featureDefn = gTindexLayer.GetLayerDefn()



    for rasterFile in rasterList:
        srcImage = gdal.Open(rasterFile)

        geoTrans, polyToCut, ulX, ulY, lrX, lrY = gT.getRasterExtent(srcImage)

        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(polyToCut)
        feature.SetField("location", rasterFile)
        gTindexLayer.CreateFeature(feature)
        feature = None


    return gTindex, gTindexLayer


def createTiledGeoJsonFromSrc(rasterFolderLocation, vectorSrcFile, geoJsonOutputDirectory, rasterTileIndex='',
                              geoJsonPrefix='GEO', rasterFileExtenstion='.tif',
                              rasterPrefixToReplace='PAN'
                              ):
    if rasterTileIndex == '':
        gTindex, gTindexLayer = buildTindex(rasterFolderLocation, rasterExtention=rasterFileExtenstion)
    else:
        gTindex = ogr.Open(rasterTileIndex,0)
        gTindexLayer = gTindex.GetLayer()

    shapeSrc = ogr.Open(vectorSrcFile,0)
    chipSummaryList = []
    for feature in gTindexLayer:
        featureGeom = feature.GetGeometryRef()
        rasterFileName = feature.GetField('location')
        rasterFileBaseName = os.path.basename(rasterFileName)
        outGeoJson = rasterFileBaseName.replace(rasterPrefixToReplace, geoJsonPrefix)
        outGeoJson = outGeoJson.replace(rasterFileExtenstion, '.geojson')
        outGeoJson = os.path.join(geoJsonOutputDirectory, outGeoJson)

        gT.clipShapeFile(shapeSrc, outGeoJson, featureGeom, minpartialPerc=0.0, debug=False)
        imageId = rasterFileBaseName.replace(rasterPrefixToReplace+"_", "")
        chipSummary = {'chipName': rasterFileName,
                           'geoVectorName': outGeoJson,
                           'imageId': os.path.splitext(imageId)[0]}

        chipSummaryList.append(chipSummary)

    return chipSummaryList

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-imgDir", "--imgDir", type=str,
                        help="Directory of Raster Images")
    parser.add_argument("-vecSrc", "--vectorSrcFile", type=str,
                        help="Geo spatial Vector src file supported by GDAL and OGR")
    parser.add_argument("-vecPrFx", "--vectorPrefix", type=str,
                        help="Prefix to attach to image id to indicate type of geojson created",
                        default='OSM')
    parser.add_argument("-rastPrFx", "--rasterPrefix", type=str,
                        help="Prefix of raster images to replace when creating geojson of geojson created",
                        default='PAN')
    parser.add_argument("-rastExt", "--rasterExtension", type=str,
                        help="Extension of raster images to i.e. .tif, .png, .jpeg",
                        default='.tif')

    parser.add_argument("-o", "--outputCSV", type=str,
                        help="Output file name and location for truth summary CSV equivalent to SpacenetV2 competition")
    parser.add_argument("-pixPrecision", "--pixelPrecision", type=int,
                        help="Number of decimal places to include for pixel, uses round(xPix, pixPrecision)"
                             "Default = 2",
                        default=2)
    parser.add_argument("--CreateProposalFile",
                        help="Create proposals file in format approriate for SpacenetV2 competition",
                        action="store_true")




    args = parser.parse_args()



    rasterFolderLocation = args.imgDir
    vectorSrcFile = args.vectorSrcFile
    vectorPrefix = args.vectorPrefix
    rasterPrefix = args.rasterPrefix
    pixPrecision = args.pixelPrecision
    createProposalFile = args.CreateProposalFile
    rasterFileExtension = args.rasterExtension

    rasterFolderBaseName = os.path.basename(rasterFolderLocation)
    if rasterFolderBaseName == "":
        rasterFolderBaseName = os.path.basename(os.path.dirname(rasterFolderLocation))

    geoJsonOutputDirectory = os.path.join(os.path.dirname(vectorSrcFile), rasterFolderBaseName)
    chipSummaryList = createTiledGeoJsonFromSrc(rasterFolderLocation, vectorSrcFile, geoJsonOutputDirectory, rasterTileIndex='',
                              geoJsonPrefix=vectorPrefix, rasterFileExtenstion=rasterFileExtension,
                              rasterPrefixToReplace=rasterPrefix
                              )


    outputCSVFileName = geoJsonOutputDirectory+"OSM_Proposal.csv"
    lT.createCSVSummaryFile(chipSummaryList, outputCSVFileName,
                                replaceImageID=rasterPrefix+"_",
                                pixPrecision=pixPrecision,
                                createProposalsFile=createProposalFile
                                )






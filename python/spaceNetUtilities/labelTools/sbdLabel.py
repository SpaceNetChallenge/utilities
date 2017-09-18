import numpy as np
import os
from spaceNetUtilities.labelTools import coreLabelTools as clT
import scipy.io
from scipy.sparse import csr_matrix
import json
#from osgeo import gdal, osr, ogr, gdalnumeric


#TODO reimplemtation using rasterio.  Simple changed needed to BurnLayer Function see coreLabelImplementation
def createClassSegmentation(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):


    dist_trans = clT.createDistanceTransform(rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    dist_trans[dist_trans > 0] = 1
    dist_trans[dist_trans < 0] = 0
    return dist_trans


def createClassBoundaries(rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
    dist_trans = clT.createDistanceTransform(rasterSrc, vectorSrc, npDistFileName='', units='pixels')
    # From distance transform to boundary
    dist_trans[dist_trans > 1.0] = 255
    dist_trans[dist_trans < -1.0] = 255
    dist_trans[dist_trans != 255] = 1
    dist_trans[dist_trans == 255] = 0
    sparse_total = csr_matrix(dist_trans)
    return sparse_total.astype(np.uint8)


def createClassCategoriesPresent(vectorSrc):
    with open(vectorSrc) as my_file:
        data = json.load(my_file)
    if (len(data['features']) == 0):
        return np.array([], dtype=np.uint8)
    else:
        return np.array([1], dtype=np.uint8)
#
#
# def createDistanceTransformByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
#     ## open source vector file that truth data
#     source_ds = ogr.Open(vectorSrc)
#     source_layer = source_ds.GetLayer()
#
#     # Define feature
#     my_feature = source_layer[feature_index]
#
#     # Spatial Reference
#     srs = source_layer.GetSpatialRef()
#
#     # Create feature Layer
#     outDriver = ogr.GetDriverByName('MEMORY')
#     outDataSource = outDriver.CreateDataSource('memData')
#     Feature_Layer = outDataSource.CreateLayer("this_feature", srs, geom_type=ogr.wkbPolygon)
#
#     # Add feature to layer
#     Feature_Layer.CreateFeature(my_feature)
#
#     ## extract data from src Raster File to be emulated
#     ## open raster file that is to be emulated
#     srcRas_ds = gdal.Open(rasterSrc)
#     cols = srcRas_ds.RasterXSize
#     rows = srcRas_ds.RasterYSize
#     noDataValue = 0
#     metersIndex = 1
#
#     ## create First raster memory layer
#     memdrv = gdal.GetDriverByName('MEM')
#     dst_ds = memdrv.Create('', cols, rows, 1, gdal.GDT_Byte)
#     dst_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
#     dst_ds.SetProjection(srcRas_ds.GetProjection())
#     band = dst_ds.GetRasterBand(1)
#     band.SetNoDataValue(noDataValue)
#
#     gdal.RasterizeLayer(dst_ds, [1], Feature_Layer, burn_values=[255])
#     srcBand = dst_ds.GetRasterBand(1)
#
#     memdrv2 = gdal.GetDriverByName('MEM')
#     prox_ds = memdrv2.Create('', cols, rows, 1, gdal.GDT_Int16)
#     prox_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
#     prox_ds.SetProjection(srcRas_ds.GetProjection())
#     proxBand = prox_ds.GetRasterBand(1)
#     proxBand.SetNoDataValue(noDataValue)
#
#     options = ['NODATA=0']
#
#     gdal.ComputeProximity(srcBand, proxBand, options)
#
#     memdrv3 = gdal.GetDriverByName('MEM')
#     proxIn_ds = memdrv3.Create('', cols, rows, 1, gdal.GDT_Int16)
#     proxIn_ds.SetGeoTransform(srcRas_ds.GetGeoTransform())
#     proxIn_ds.SetProjection(srcRas_ds.GetProjection())
#     proxInBand = proxIn_ds.GetRasterBand(1)
#     proxInBand.SetNoDataValue(noDataValue)
#     options = ['NODATA=0', 'VALUES=0']
#     gdal.ComputeProximity(srcBand, proxInBand, options)
#
#     proxIn = gdalnumeric.BandReadAsArray(proxInBand)
#     proxOut = gdalnumeric.BandReadAsArray(proxBand)
#
#     proxTotal = proxIn.astype(float) - proxOut.astype(float)
#     proxTotal = proxTotal * metersIndex
#
#     if npDistFileName != '':
#         np.save(npDistFileName, proxTotal)
#
#     return proxTotal

#
# def createSegmentationByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
#     dist_trans_by_feature = createDistanceTransformByFeatureIndex(feature_index, rasterSrc, vectorSrc,
#                                                                   npDistFileName='', units='pixels')
#     dist_trans_by_feature[dist_trans_by_feature > 0] = feature_index + 1
#     dist_trans_by_feature[dist_trans_by_feature < 0] = 0
#     return dist_trans_by_feature.astype(np.uint8)
#
#
# def createInstanceSegmentation(rasterSrc, vectorSrc):
#     json_data = open(vectorSrc)
#     data = json.load(json_data)
#     num_features = len(data['features'])
#
#     cell_array = np.zeros((num_features,), dtype=np.object)
#     for i in range(num_features):
#         cell_array[i] = createSegmentationByFeatureIndex(i, rasterSrc, vectorSrc, npDistFileName='', units='pixels')
#     return cell_array
#
#
# def createBoundariesByFeatureIndex(feature_index, rasterSrc, vectorSrc, npDistFileName='', units='pixels'):
#     dist_trans_by_feature = createDistanceTransformByFeatureIndex(feature_index, rasterSrc, vectorSrc,
#                                                                   npDistFileName='', units='pixels')
#     dist_trans_by_feature[dist_trans_by_feature > 1.0] = 255
#     dist_trans_by_feature[dist_trans_by_feature < -1.0] = 255
#     dist_trans_by_feature[dist_trans_by_feature != 255] = 1
#     dist_trans_by_feature[dist_trans_by_feature == 255] = 0
#     return dist_trans_by_feature.astype(np.uint8)
#
#
# def createInstanceBoundaries(rasterSrc, vectorSrc):
#     json_data = open(vectorSrc)
#     data = json.load(json_data)
#     num_features = len(data['features'])
#
#     cell_array = np.zeros((num_features,), dtype=np.object)
#     for i in range(num_features):
#         full_boundary_matrix = createBoundariesByFeatureIndex(i, rasterSrc, vectorSrc, npDistFileName='',
#                                                               units='pixels')
#         cell_array[i] = csr_matrix(full_boundary_matrix)
#     return cell_array
#
#
# def createInstanceCategories(vectorSrc):
#     with open(vectorSrc) as my_file:
#         data = json.load(my_file)
#     if (len(data['features']) == 0):
#         return np.array([], dtype=np.uint8)
#     else:
#         return np.ones(len(data['features']), dtype=np.uint8).reshape((len(data['features']), 1))
#
#
# def geoJsonToSBD(annotationName_cls, annotationName_inst, geoJson, rasterSource):
#     # Print raster file name
#     my_raster_source = rasterSource
#     print("Raster directory : ", my_raster_source)
#     srcRaster = gdal.Open(my_raster_source)
#
#     my_vector_source = geoJson
#     print("Vector directory : ", my_vector_source)
#
#     # Call main functions to create label datafor cls
#     my_cls_segmentation = createClassSegmentation(my_raster_source, my_vector_source, npDistFileName='', units='pixels')
#     my_cls_boundaries = createClassBoundaries(my_raster_source, my_vector_source, npDistFileName='', units='pixels')
#     my_cls_categories = createClassCategoriesPresent(my_vector_source)
#
#     # Call main functions to create label datafor inst
#     my_inst_segmentation = createInstanceSegmentation(my_raster_source, my_vector_source)
#     my_inst_boundaries = createInstanceBoundaries(my_raster_source, my_vector_source)
#     my_inst_categories = createInstanceCategories(my_vector_source)
#
#     # Wraps for cls struct
#     cls_boundaries_wrap = np.array([my_cls_boundaries])
#     cls_categories_wrap = my_cls_categories
#
#     # Wraps for inst struct
#     inst_boundaries_wrap = np.array([my_inst_boundaries])
#     inst_categories_wrap = my_inst_categories
#
#     # Create a class struct
#     GTcls = {'Segmentation': my_cls_segmentation, 'Boundaries': cls_boundaries_wrap,
#              'CategoriesPresent': cls_categories_wrap}
#
#     # Create the instance struct
#     GTinst = {'Segmentation': my_inst_segmentation, 'Boundaries': inst_boundaries_wrap,
#               'Categories': inst_categories_wrap}
#
#     # Save the files
#     scipy.io.savemat(annotationName_cls, {'GTcls': GTcls})
#     scipy.io.savemat(annotationName_inst, {'GTinst': GTinst})
#
#     print("Done with " + str())
#
#     entry = {'rasterFileName': my_raster_source,
#              'geoJsonFileName': geoJson,
#              'annotationName': annotationName_cls,
#              'annotationName_cls': annotationName_cls,
#              'annotationName_inst': annotationName_inst,
#              'width': srcRaster.RasterXSize,
#              'height': srcRaster.RasterYSize,
#              'depth': srcRaster.RasterCount,
#              'basename': os.path.splitext(os.path.basename(my_raster_source))[0]
#              }
#
#     return entry
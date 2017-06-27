import numpy as np
import geoTools as gT
from osgeo import ogr
import os

def iou(test_poly, truth_polys, truth_index=[]):
    fidlistArray = []
    iou_list = []
    if truth_index:
        fidlist = gT.search_rtree(test_poly, truth_index)

        for fid in fidlist:
            if not test_poly.IsValid():
                test_poly = test_poly.Buffer(0.0)

            intersection_result = test_poly.Intersection(truth_polys[fid])
            fidlistArray.append(fid)

            if intersection_result.GetGeometryName() == 'POLYGON' or \
                            intersection_result.GetGeometryName() == 'MULTIPOLYGON':
                intersection_area = intersection_result.GetArea()
                union_area = test_poly.Union(truth_polys[fid]).GetArea()
                iou_list.append(intersection_area / union_area)

            else:
                iou_list.append(0)

    else:

        for idx, truth_poly in enumerate(truth_polys):
            intersection_result = test_poly.Intersection(truth_poly)
            intersection_result.GetGeometryName()

            if intersection_result.GetGeometryName() == 'POLYGON' or \
                            intersection_result.GetGeometryName() == 'MULTIPOLYGON':
                intersection_area = intersection_result.GetArea()
                union_area = test_poly.Union(truth_poly).GetArea()
                iou_list.append(intersection_area / union_area)

            else:

                iou_list.append(0)

    return iou_list, fidlistArray


def score(test_polys, truth_polys, threshold=0.5, truth_index=[],
          resultGeoJsonName = [],
          imageId = []):

    # Define internal functions

    # Find detections using threshold/argmax/IoU for test polygons
    true_pos_count = 0
    false_pos_count = 0
    truth_poly_count = len(truth_polys)

    if resultGeoJsonName:
        if not imageId:
            imageId = os.path.basename(os.path.splitext(resultGeoJsonName)[0])

        driver = ogr.GetDriverByName('geojson')
        if os.path.exists(resultGeoJsonName):
            driver.DeleteDataSource(resultGeoJsonName)
        datasource = driver.CreateDataSource(resultGeoJsonName)
        layer = datasource.CreateLayer('buildings', geom_type=ogr.wkbPolygon)
        field_name = ogr.FieldDefn("ImageId", ogr.OFTString)
        field_name.SetWidth(75)
        layer.CreateField(field_name)
        layer.CreateField(ogr.FieldDefn("BuildingId", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("IOUScore", ogr.OFTReal))


    for test_poly in test_polys:
        if truth_polys:
            iou_list, fidlist = iou(test_poly, truth_polys, truth_index)
            if not iou_list:
                maxiou = 0
            else:
                maxiou = np.max(iou_list)

            if maxiou >= threshold:
                true_pos_count += 1
                truth_index.delete(fidlist[np.argmax(iou_list)], truth_polys[fidlist[np.argmax(iou_list)]].GetEnvelope())
                #del truth_polys[fidlist[np.argmax(iou_list)]]
                if resultGeoJsonName:
                    feature = ogr.Feature(layer.GetLayerDefn())
                    feature.SetField('ImageId', imageId)
                    feature.SetField('BuildingId', fidlist[np.argmax(iou_list)])
                    feature.SetField('IOUScore', maxiou)
                    feature.SetGeometry(test_poly)


            else:
                false_pos_count += 1

                if resultGeoJsonName:
                    feature = ogr.Feature(layer.GetLayerDefn())
                    feature.SetField('ImageId', imageId)
                    feature.SetField('BuildingId', -1)
                    feature.SetField('IOUScore', maxiou)
                    feature.SetGeometry(test_poly)






        else:
            false_pos_count += 1

            if resultGeoJsonName:
                feature = ogr.Feature(layer.GetLayerDefn())
                feature.SetField('ImageId', imageId)
                feature.SetField('BuildingId', 0)
                feature.SetField('IOUScore', 0)
                feature.SetGeometry(test_poly)

        if resultGeoJsonName:
            layer.CreateFeature(feature)
            feature.Destroy()

    if resultGeoJsonName:
        datasource.Destroy()


    false_neg_count = truth_poly_count - true_pos_count

    return true_pos_count, false_pos_count, false_neg_count


def evalfunction((image_id, test_polys, truth_polys, truth_index),
                 resultGeoJsonName=[],
                 threshold = 0.5):


    if len(truth_polys)==0:
        true_pos_count = 0
        false_pos_count = len(test_polys)
        false_neg_count = 0
    else:
        true_pos_count, false_pos_count, false_neg_count = score(test_polys, truth_polys.tolist(),
                                                                 truth_index=truth_index,
                                                                 resultGeoJsonName=resultGeoJsonName,
                                                                 imageId=image_id,
                                                                 threshold=threshold
                                                                 )


    if (true_pos_count > 0):

        precision = float(true_pos_count) / (float(true_pos_count) + float(false_pos_count))
        recall = float(true_pos_count) / (float(true_pos_count) + float(false_neg_count))
        F1score = 2.0 * precision * recall / (precision + recall)
    else:
        F1score = 0
    return ((F1score, true_pos_count, false_pos_count, false_neg_count), image_id)


def  create_eval_function_input((image_ids, (prop_polysIdList, prop_polysPoly), (sol_polysIdsList, sol_polysPoly))):

    evalFunctionInput = []


    for image_id in image_ids:
        test_polys = prop_polysPoly[np.argwhere(prop_polysIdList == image_id).flatten()]
        truth_polys = sol_polysPoly[np.argwhere(sol_polysIdsList == image_id).flatten()]
        truth_index = gT.create_rtree_from_poly(truth_polys)
        evalFunctionInput.append([image_id, test_polys, truth_polys, truth_index])

    return evalFunctionInput



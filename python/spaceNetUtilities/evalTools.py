import numpy as np
import geoTools as gT
from shapely.geometry import mapping, Polygon
import fiona
from tqdm import tqdm
import os

def iou(test_poly, truth_polys, truth_index=[]):
    fidlistArray = []
    iou_list = []
    if truth_index:
        fidlist = gT.search_rtree(test_poly, truth_index)
        #print(test_poly)
        for fid in fidlist:
            if not test_poly.is_valid:
                test_poly = test_poly.buffer(0.0)

            intersection_result = test_poly.intersection(truth_polys[fid])
            fidlistArray.append(fid)

            if intersection_result.geom_type == 'Polygon' or \
                            intersection_result.geom_type == 'MultiPolygon':
                intersection_area = intersection_result.area
                union_area = test_poly.union(truth_polys[fid]).area
                iou_list.append(intersection_area / union_area)

            else:
                iou_list.append(0)


    return iou_list, fidlistArray

def write_geojson(geojson_name,
                  feature_list,
                  output_crs={'init': 'epsg:4326'},
                  output_schema={'geometry': 'Polygon',
                                 'properties': {'ImageId': 'str',
                                                'IOUScore': 'float:15.5',
                                                'BuildingId': 'int'}
                                 },
                  output_driver='GeoJSON'
                  ):
    with fiona.open(geojson_name,'w',
                    driver=output_driver,
                    crs=output_crs,
                    schema=output_schema) as sink:

        for feature in feature_list:
            sink.write(feature)


def score(test_polys, truth_polys, threshold=0.5, truth_index=[],
          resultGeoJsonName = [],
          imageId = []):

    # Define internal functions
    #
    output_schema = {'geometry': 'Polygon',
                    'properties': {'ImageId': 'str',
                                   'IOUScore': 'float:15.5',
                                   'BuildingId': 'int'}
                    }
    output_crs = {'init': 'epsg:4326'}

    # Find detections using threshold/argmax/IoU for test polygons
    true_pos_count = 0
    false_pos_count = 0
    truth_poly_count = len(truth_polys)

    if resultGeoJsonName:
        if not imageId:
            imageId = os.path.basename(os.path.splitext(resultGeoJsonName)[0])

    feature_list = []

    for test_poly in tqdm(test_polys):
        if truth_polys:
            iou_list, fidlist = iou(test_poly, truth_polys, truth_index)
            if not iou_list:
                maxiou = 0
            else:
                maxiou = np.max(iou_list)

            if maxiou >= threshold:
                true_pos_count += 1
                truth_index.delete(fidlist[np.argmax(iou_list)], truth_polys[fidlist[np.argmax(iou_list)]].bounds)
                #del truth_polys[fidlist[np.argmax(iou_list)]]

                feature = {'geometry': mapping(test_poly),
                           'properties': {'ImageId': imageId,
                                          'IOUScore': maxiou,
                                          'BuildingId': fidlist[np.argmax(iou_list)]
                                          }
                           }

            else:
                false_pos_count += 1

                feature = {'geometry': mapping(test_poly),
                           'properties': {'ImageId': imageId,
                                          'IOUScore': maxiou,
                                          'BuildingId': -1
                                          }
                           }

        else:
            false_pos_count += 1
            feature = {'geometry': mapping(test_poly),
                       'properties': {'ImageId': imageId,
                                      'IOUScore': 0,
                                      'BuildingId': 0
                                      }
                       }

        feature_list.append(feature)

    if resultGeoJsonName:
        write_geojson(resultGeoJsonName, feature_list)


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



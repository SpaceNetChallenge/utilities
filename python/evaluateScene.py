from spaceNet import evalTools as eT
from spaceNet import geoTools as gT
import numpy as np
import sys
import multiprocessing
import time

if __name__ == "__main__":

    # load Truth and Test File Locations
    if len(sys.argv) > 1:
        truth_fp = sys.argv[1]
        test_fp = sys.argv[2]
    else:
        test_fp = '../testData/public_polygons_solution_3Band_envelope.geojson'
        truth_fp = '../testData/public_polygons_solution_3Band.geojson'
    # check for cores available
    if len(sys.argv) > 3:
        max_cpu = int(sys.argv[3])
    else:
        max_cpu = multiprocessing.cpu_count()
    parallel=False

    # initialize scene counts
    true_pos_counts = []
    false_pos_counts = []
    false_neg_counts = []

    t0 = time.time()
    # Start Ingest Of Truth and Test Case
    sol_polys = gT.importgeojson(truth_fp, removeNoBuildings=True)
    prop_polys = gT.importgeojson(test_fp)

    t1 = time.time()
    total = t1 - t0
    print('time of ingest: ', total)

    # Speed up search by preprocessing ImageId and polygonIds

    test_image_ids = set([item['ImageId'] for item in prop_polys if item['ImageId'] > 0])
    prop_polysIdList = np.asarray([item['ImageId'] for item in prop_polys if item["ImageId"] > 0 and \
                                   item['BuildingId']!=-1])
    prop_polysPoly = np.asarray([item['poly'] for item in prop_polys if item["ImageId"] > 0 and \
                                   item['BuildingId']!=-1])

    sol_polysIdsList = np.asarray([item['ImageId'] for item in sol_polys if item["ImageId"] > 0 and \
                                   item['BuildingId']!=-1])
    sol_polysPoly = np.asarray([item['poly'] for item in sol_polys if item["ImageId"] > 0 and \
                                   item['BuildingId']!=-1])
    bad_count = 0
    F1ScoreList = []
    cpu_count = min(multiprocessing.cpu_count(), max_cpu)
    print('{}'.format(max_cpu))
    p = multiprocessing.Pool(processes=cpu_count)
    ResultList = []

    eval_function_input_list = eT.create_eval_function_input((test_image_ids,
                                                         (prop_polysIdList, prop_polysPoly),
                                                         (sol_polysIdsList, sol_polysPoly)))
    # evalFunctionInput =  creatEevalFunctionInput((test_image_ids,
    #                                               (prop_polysIdList, prop_polysPoly),
    #                                               (sol_polysIdsList, sol_polysPoly)))
    # Calculate Values
    t3 = time.time()
    print('time For DataCreation {}s'.format(t3-t1))
    #result_list = p.map(eT.evalfunction, eval_function_input_list)
    if parallel==False:
        result_list = []
        for eval_input in eval_function_input_list:

            result_list.append(eT.evalfunction(eval_input))
    else:
        result_list = p.map(eT.evalfunction, eval_function_input_list)

    result_sum = np.sum(result_list, axis=0)
    true_pos_total = result_sum[1]
    false_pos_total = result_sum[2]
    false_neg_total = result_sum[3]
    print('True_Pos_Total', true_pos_total)
    print('False_Pos_Total', false_pos_total)
    print('False_Neg_Total', false_neg_total)
    precision = float(true_pos_total) / (float(true_pos_total) + float(false_pos_total))
    recall = float(true_pos_total) / (float(true_pos_total) + float(false_neg_total))
    F1ScoreTotal = 2.0 * precision*recall / (precision + recall)
    print('F1Total', F1ScoreTotal)

    t2 = time.time()
    total = t2-t0
    print('time of evaluation: {}'.format(t2-t1))
    print('time of evaluation {}s/imageId'.format((t2-t1)/len(result_list)))
    print('Total Time {}s'.format(total))
    print(result_list)
    print(np.mean(result_list))

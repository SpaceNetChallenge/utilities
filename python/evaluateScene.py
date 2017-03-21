from spaceNetUtilities import evalTools as eT
from spaceNetUtilities import geoTools as gT
import numpy as np
import csv
import multiprocessing
import time
import argparse



def evaluateSpaceNetSolution(summaryTruthFile, summaryProposalFile, resultsOutputFile='', processgeoJson=False,
                             useParallelProcessing=False):

    truth_fp = summaryTruthFile
    test_fp = summaryProposalFile
    # check for cores available
    if useParallelProcessing:

        max_cpu = multiprocessing.cpu_count()
        parallel = True
    else:
        max_cpu = 1
        parallel = False

    # initialize scene counts
    true_pos_counts = []
    false_pos_counts = []
    false_neg_counts = []

    t0 = time.time()
    # Start Ingest Of Truth and Test Case
    if args.geoJson:
        sol_polys = gT.import_summary_geojson(truth_fp, removeNoBuildings=True)
        prop_polys = gT.import_summary_geojson(test_fp)
        polyFlag = 'poly'
    else:
        sol_polys = gT.readwktcsv(truth_fp, removeNoBuildings=True)
        prop_polys = gT.readwktcsv(test_fp)
        polyFlag = 'polyPix'

    t1 = time.time()
    total = t1 - t0
    print('time of ingest: ', total)

    # Speed up search by preprocessing ImageId and polygonIds
    polyFlag = ''
    test_image_ids = set([item['ImageId'] for item in prop_polys if item['ImageId'] > 0])
    prop_polysIdList = np.asarray([item['ImageId'] for item in prop_polys if item["ImageId"] > 0 and \
                                   item['BuildingId'] != -1])
    prop_polysPoly = np.asarray([item[polyFlag] for item in prop_polys if item["ImageId"] > 0 and \
                                 item['BuildingId'] != -1])

    sol_polysIdsList = np.asarray([item['ImageId'] for item in sol_polys if item["ImageId"] > 0 and \
                                   item['BuildingId'] != -1])
    sol_polysPoly = np.asarray([item[polyFlag] for item in sol_polys if item["ImageId"] > 0 and \
                                item['BuildingId'] != -1])
    bad_count = 0
    F1ScoreList = []
    cpu_count = min(multiprocessing.cpu_count(), max_cpu)
    print('{}'.format(max_cpu))
    p = multiprocessing.Pool(processes=cpu_count)
    ResultList = []

    eval_function_input_list = eT.create_eval_function_input((test_image_ids,
                                                              (prop_polysIdList, prop_polysPoly),
                                                              (sol_polysIdsList, sol_polysPoly)))

    # Calculate Values
    t3 = time.time()
    print('time For DataCreation {}s'.format(t3 - t1))

    # result_list = p.map(eT.evalfunction, eval_function_input_list)
    if parallel == False:
        result_list = []
        for eval_input in eval_function_input_list:
            result_list.append(eT.evalfunction(eval_input))
    else:
        result_list = p.map(eT.evalfunction, eval_function_input_list)

    result_listNP = np.asarray([item[0] for item in result_list])

    result_sum = np.sum(result_listNP, axis=0)
    true_pos_total = result_sum[1]
    false_pos_total = result_sum[2]
    false_neg_total = result_sum[3]
    print('True_Pos_Total', true_pos_total)
    print('False_Pos_Total', false_pos_total)
    print('False_Neg_Total', false_neg_total)
    precision = float(true_pos_total) / (float(true_pos_total) + float(false_pos_total))
    recall = float(true_pos_total) / (float(true_pos_total) + float(false_neg_total))
    F1ScoreTotal = 2.0 * precision * recall / (precision + recall)
    print('F1Total', F1ScoreTotal)

    t2 = time.time()
    total = t2 - t0
    print('time of evaluation: {}'.format(t2 - t1))
    print('time of evaluation {}s/imageId'.format((t2 - t1) / len(result_list)))
    print('Total Time {}s'.format(total))
    print(result_list)
    print(np.mean(result_listNP))

    if resultsOutputFile != '':
        with open(resultsOutputFile, 'w') as csvFile:
            csvwriter = csv.writer(csvFile, delimiter=',')
            csvwriter.writerow['TruthFile', truth_fp]
            csvwriter.writerow['ProposalFile', test_fp]
            csvwriter.writerow['Summary Results']
            csvwriter.writerow['F1Score Total', F1ScoreTotal]
            csvwriter.writerow['Precision', precision]
            csvwriter.writerow['Recall', recall]
            csvwriter.writerow['True Positive Total', true_pos_total]
            csvwriter.writerow['False Positive Total', false_pos_total]
            csvwriter.writerow['False Negative Total', false_neg_total]
            csvwriter.writerow['']
            csvwriter.writerow['Per Image Stats']
            csvwriter.writerow['ImageId', 'F1Score', 'True Positive Count', 'False Positive Count']
            for result in result_list:
                tmpList = [result[1]]
                tmpList.extend(result[0])
                csvwriter.writerow(tmpList)


    resultsDict = {'TruthFile': truth_fp,
                   'ProposalFile': test_fp,
                   'F1ScoreTotal': F1ScoreTotal,
                   'PrecisionTotal': precision,
                   'RecalTotal': recall,
                   'TruePositiveTotal': true_pos_total,
                   'FalsePositiveTotal': false_pos_total,
                   'FalseNegativeTotal': false_neg_total,
                   'PerImageStatsResultList': result_list,
                   'OutputSummaryFile': resultsOutputFile}

    return resultsDict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Evaluate Score for SpaceNet')
    parser.add_argument("summaryTruthFile",
                        help="The Location of Summary Ground Truth csv File"
                             "Format is {},{},{},{}.format(ImageId, BuildingId, polygonPixWKT, polygonGeoPix),"
                             "unless --geoJson flag is set"
                             "See spaceNet competition details for more information about file format"
                        )
    parser.add_argument("summaryProposalFile",
                        help="The Location of summary Propsal csv File"
                             "Format is {},{},{},{}.format(ImageId, BuildingId, polygonPixWKT, Confidence)"
                             "unless --geoJson flag is set"
                        )
    parser.add_argument("--resultsOutputFile",
                        help="If you would like summary data outwritten to a file, specify the file",
                        default='')
    parser.add_argument("--geoJson",
                        help='Convert Image from Native format to 8bit',
                        action='store_true')

    parser.add_argument("--useParallelProcessing",
                        help='Convert Image from Native format to 8bit',
                        action='store_true')

    args = parser.parse_args()
    # load Truth and Test File Locations

    summaryDict = evaluateSpaceNetSolution(args.summaryTruthFile,
                                           args.summaryProposalFile,
                                           resultsOutputFile=args.resultsOutputFile,
                                           processgeoJson=args.geoJson,
                                           useParallelProcessing=args.useParallelProcessing)




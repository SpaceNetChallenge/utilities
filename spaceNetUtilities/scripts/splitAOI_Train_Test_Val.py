from spaceNetUtilities import labelTools as lT
import glob
import argparse
import random
import os
import errno
import shutil


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-imgDir", "--imgDir", action='append', type=str,
                        help="Directory of Raster Images BaseLayer")
    parser.add_argument("-geoDir", "--geojsonDir", type=str,
                        help="Directory of geojson files")

    parser.add_argument("-o", "--outputCSV", type=str,
                        help="Output File Name and Location for CSV")
    parser.add_argument("-pixPrecision", "--pixelPrecision", type=int,
                        help="Number of decimal places to include for pixel, uses round(xPix, pixPrecision)"
                             "Default = 2",
                        default=2)
    parser.add_argument("-outputDir", "--outputDir", type=str,
                        help="baseLocation PlaceOutproduct")

    parser.add_argument("--CreateProposalFile", help="Create ProposalsFile",
                        action="store_true")
    parser.add_argument("-strip", "--stripOutFromGeoJson", type=str,
                        help="string delimited")
    parser.add_argument("--DontstripFirstUnderScore", action="store_false")
    args = parser.parse_args()

    rasterDirectory = args.imgDir
    geoJsonDirectory = args.geojsonDir
    outputCSVFileName = args.outputCSV
    createProposalFile = args.CreateProposalFile
    if args.stripOutFromGeoJson:
        stripList = args.stripOutFromGeoJson.split(' ')
    else:
        stripList = []

    imageDirList = args.imgDir

    jsonList = []
    chipSummaryList = []

    rasterDirList = []
    secondaryImagePrefix = []
    for imageDir in imageDirList:
        rasListTemp = glob.glob(os.path.join(imageDir, '*.tif'))
        rasterDirList.append([imageDir, os.path.basename(rasListTemp[0]).split("_")[0], rasListTemp
                              ])

    geoJsonList = glob.glob(os.path.join(geoJsonDirectory, '*.geojson'))
    for imageId in geoJsonList:
        imageId = os.path.basename(imageId)
        rasterName = imageId.replace('.geojson', '.tif')

        for stripItem in stripList:
            rasterName = rasterName.replace(stripItem, '')

        chipNameList = []
        for rastDir in rasterDirList:
            if args.DontstripFirstUnderScore:

                rasterName = rastDir[1] + "_" + rasterName.split('_', 1)[1]

            else:
                rasterName = rastDir[1] + "_" + rasterName
            chipName = [rastDir[1], os.path.join(rastDir[0], rasterName)]
            chipNameList.append(chipName)
        print(imageId)

        chipSummary = {'chipName': chipNameList,
                       'geoVectorName': os.path.join(geoJsonDirectory, imageId),
                       'imageId': os.path.splitext(imageId)[0]}

        chipSummaryList.append(chipSummary)

    random.shuffle(chipSummaryList)
    trainSplitPoint = int(round(len(chipSummaryList)*0.6))
    valSplitPoint = int(round(len(chipSummaryList)*0.8))

    splitInformationList = [['train', chipSummaryList[0:trainSplitPoint]],
                            ['test',  chipSummaryList[trainSplitPoint+1:valSplitPoint]],
                            ['validate', chipSummaryList[valSplitPoint+1:]]
                            ]

    outputDir = args.outputDir

    # Perform Split
    for splitInformation in splitInformationList:

        for rastDir in rasterDirList:
            print(outputDir)
            print(splitInformationList)

            mkdir_p(os.path.join(outputDir, splitInformation[0], rastDir[1]))

        mkdir_p(os.path.join(outputDir, splitInformation[0], 'geojson', 'buildings'))

        for chipSummary in splitInformation[1]:
            for chip in chipSummary['chipName']:
                shutil.copy(chip[1], os.path.join(outputDir, splitInformation[0], chip[0]))

            shutil.copy(chipSummary['geoVectorName'], os.path.join(outputDir,
                                                                   splitInformation[0],
                                                                   'geojson',
                                                                   'buildings'))

        rasterPrefix = 'PAN'
        chipSummaryListTmp = []

        for chipSummary in splitInformation[1]:

            chipSummaryTmp = chipSummary
            chipSummaryTmp['chipName'] = chipSummary['chipName'][0][1]
            chipSummaryListTmp.append(chipSummaryTmp)

        print("starting")
        outputCSVFileName = os.path.join(outputDir, args.outputCSV + splitInformation[0] + ".csv")

        lT.createCSVSummaryFile(chipSummaryListTmp, outputCSVFileName,
                                replaceImageID=rasterPrefix+"_",
                                pixPrecision=args.pixelPrecision)

        print("finished")

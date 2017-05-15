import os
import glob
import random
from spaceNetUtilities import labelTools as lT
from spaceNetUtilities import geoTools as gT
import argparse


def processRasterChip(rasterImage, rasterDescription, geojson, geojsonDescription, outputDirectory='',
                      imagePixSize=-1, clipOverlap=0.0, randomClip=False,
                      minpartialPerc=0.0,
                      outputPrefix=''):

    # cut Image to Size
    chipSummaryList=[]
    if imagePixSize>0:

        rasterFileList = [[rasterImage, rasterDescription]]
        shapeFileSrcList = [[geojson, geojsonDescription]]
        # cut image to size
        print(rasterFileList)
        chipSummaryList = gT.cutChipFromMosaic(rasterFileList,
                                           shapeFileSrcList,
                                           outputDirectory=outputDirectory,
                                           outputPrefix=outputPrefix,
                                           clipSizeMX=imagePixSize,
                                           clipSizeMY=imagePixSize,
                                           minpartialPerc=minpartialPerc,
                                           createPix=True,
                                           clipOverlap=clipOverlap,
                                           noBlackSpace=True,
                                           randomClip=-1
                                           )





    else:
        chipSummary =  {'rasterSource': rasterImage,
                   'chipName': rasterImage,
                   'geoVectorName': geojson,
                   'pixVectorName': ''
                   }

        chipSummaryList.append(chipSummary)


    return chipSummaryList


def processChipSummaryList(chipSummaryList, outputDirectory='', annotationType='PASCALVOC2012', outputFormat='GTiff',
                           outputPixType='',
                           datasetName='spacenetV2',
                           folder_name='folder_name',
                           bboxResize=1.0
                           ):

    if outputPixType == '':
        convertTo8Bit = False
    else:
        convertTo8Bit = True

    entryList = []
    for chipSummary in chipSummaryList:



        annotationName = os.path.basename(chipSummary['rasterSource'])
        annotationName = annotationName.replace('.tif', '.xml')
        annotationName = os.path.join(outputDirectory, annotationName)



        if annotationType=='PASCALVOC2012':
            entry = lT.geoJsonToPASCALVOC2012(annotationName, chipSummary['geoVectorName'], chipSummary['rasterSource'],
                                              dataset='spacenetV2',
                                              folder_name='spacenetV2',
                                              annotationStyle=annotationType,
                                              segment=True,
                                              bufferSizePix=2.5,
                                              convertTo8Bit=convertTo8Bit,
                                              outputPixType=outputPixType,
                                              outputFormat=outputFormat,
                                              bboxResize=bboxResize
                                              )
        elif annotationType=='DARKNET':
            entry = lT.geoJsonToDARKNET(annotationName, chipSummary['geoVectorName'], chipSummary['rasterSource'],
                                        dataset='spacenetV2',
                                        folder_name='spacenetV2',
                                        annotationStyle=annotationType,
                                        convertTo8Bit=convertTo8Bit,
                                        outputPixType=outputPixType,
                                        outputFormat=outputFormat,
                                        bboxResize=bboxResize
                                        )

        elif annotationType=='SBD':
            basename = os.path.basename(chipSummary['rasterSource'])
            annotationName = basename.replace('.tif', '.mat')
            annotationName_cls = os.path.join(outputDirectory,'cls', annotationName)
            annotationName_inst = os.path.join(outputDirectory,'inst', annotationName)

            #Check to make sure output directories exist, if not make it.
            if not os.path.exists(os.path.join(outputDirectory,'cls')):
              os.makedirs(os.path.join(outputDirectory,'cls'))
            if not os.path.exists(os.path.join(outputDirectory,'inst')):
              os.makedirs(os.path.join(outputDirectory,'inst'))
            
            entry = lT.geoJsonToSBD(annotationName_cls, annotationName_inst, chipSummary['geoVectorName'], chipSummary['rasterSource'])

        else:

            print("Annotation Type = {} is not supported yet".format(annotationType))
            break



        entryList.append(entry)

    return entryList

def createTrainTestSplitSummary(entryList, trainTestSplit=0.8,
                                outputDirectory='',
                                annotationSummaryPrefix='',
                                annotationType='PASCALVOC2012',
                                shuffleList=True,
                           ):

    if shuffleList:
        random.shuffle(entryList)


    splitPoint=int(round(len(entryList)*trainTestSplit))
    trainvalList = entryList[0:splitPoint]
    testList     = entryList[splitPoint+1:]


    trainValFileName = os.path.join(outputDirectory, annotationSummaryPrefix+'trainval.txt')
    print('creating trainval.txt {} entries'.format(len(trainvalList)))
    print('Writing to TrainVal List to file: {}'.format(trainValFileName))
    with open(trainValFileName, 'w') as f:
        for entry in trainvalList:
            if annotationType=='SBD':
                f.write('{} {} {}\n'.format(entry['rasterFileName'], entry['annotationName_cls'],
                                            entry['annotationName_inst']))
            else:
                f.write('{} {}\n'.format(entry['rasterFileName'], entry['annotationName']))


    testFileName = os.path.join(outputDirectory, annotationSummaryPrefix+'test.txt')
    testNameSizeFileName = os.path.join(outputDirectory, annotationSummaryPrefix + 'test_name_size.txt')
    print('creating test.txt {} entries'.format(len(testList)))
    print('Writing to Test List to file: {}'.format(testFileName))
    with open(testFileName, 'w') as f, \
            open(testNameSizeFileName, 'w') as fname:
        for entry in testList:
            if annotationType == 'SBD':
                f.write('{} {} {}\n'.format(entry['rasterFileName'], entry['annotationName_cls'],
                                            entry['annotationName_inst']))
                fname.write('{} {} {}\n'.format(entry['basename'], entry['width'], entry['height']))
            else:
                f.write('{} {}\n'.format(entry['rasterFileName'], entry['annotationName']))
                fname.write('{} {} {}\n'.format(entry['basename'], entry['width'], entry['height']))


    return (trainValFileName, testFileName, testNameSizeFileName)

if __name__ == '__main__':

    #python createDataSpaceNet.py /data/spacenet_sample/AOI_2_Vegas_Train/ RGB-PanSharpen \
    #                             --outputDirectory /data/spacenet_sample/annotations/ \
    #                             --imgSizePix 416

    # python createDataSpaceNet.py /data/spacenet_sample/AOI_2_Vegas_Train/
    # RGB-PanSharpen --outputDirectory /data/spacenet_sample/annotations/ --imgSizePix 416
    # --annotationType PASCALVOC2012 --convertTo8Bit


    # python createDataSpaceNet.py /data/spacenet_sample/AOI_2_Vegas_Train/
    # RGB-PanSharpen --outputDirectory /data/spacenet_sample/annotations/ --imgSizePix 416 --annotationType DARKNET
    #  --convertTo8Bit
    parser = argparse.ArgumentParser(description='Process SrcData for Region ComputerVision Dataset')
    parser.add_argument("srcSpaceNetFolder",
                        help="location of Spacenet AOI Data i.e. '/path/to/AOI_2_Vegas")
    parser.add_argument("--srcImageryDirectory",
                        help="The Folder to look for imagery in"
                             "i.e. RGB-PanSharpen",
                        default='RGB-PanSharpen')
    parser.add_argument("--geoJsonDirectory",
                        help="name of geojson folder typedirectory to use"
                             "i.e. 'buildings'",
                        default='buildings')
    parser.add_argument("--outputDirectory",
                        help="Location To place processed Files"
                             "If not used output directory will be"
                             "os.path.join(srcSpacenetFolder, 'annotations')",
                        default='annotations')
    parser.add_argument("--annotationSummaryPrefix",
                        help="Prefix to attach to trainval.txt and test.txt",
                        default='')
    parser.add_argument("--imgSizePix",
                        help="set the dimensions of the square image in pixels"
                             "Default is -1, images are not modified",
                        type=int,
                        default=-1)
    parser.add_argument("--annotationType",
                        help="Set the annotationType.  Currently Supported are:"
                             "1.  PASCALVOC2012"
                             "2.  DARKNET"
                             "3.  SBD"
                             "default is PASCALVOC2012",
                        default='PASCALVOC2012')
    parser.add_argument("--outputFileType",
                        help="What type of image type would you like to output to currently supported are:"
                             "1. GTiff"
                             "2. JPEG",
                        default='GTiff')
    parser.add_argument("--convertTo8Bit",
                        help='Convert Image from Native format to 8bit',
                        action='store_true')
    parser.add_argument("--featureName",
                        help='Type of feature to be summarized by csv (i.e. Building)'
                             'Currently in SpaceNet V2 Building is only label',
                        type=str,
                        default='Buildings')
    parser.add_argument("--spacenetVersion",
                        help='Spacenet Version to process,  '
                             'Version 1 supports AOI_1_RIO, '
                             'Version 2 is AOI_2_Vegas-AOI_5_Khartoum',
                        type=int,
                        default=2)
    parser.add_argument("--trainTestSplit",
                        help='Decimal of data to use for training i.e. 0.8 = 80% of data for Training',
                        type=float,
                        default=0.8)
    parser.add_argument("--boundingBoxResize",
                        help='Decimal Resize Annotation, i.e 0.8 = shrink annotation by 20 percent'
                             'default is 1.0',
                        type=float,
                        default=1.0)

    args = parser.parse_args()

    entryList = []
    srcSpaceNetDirectory = args.srcSpaceNetFolder

    #listOfAOIs = [subdir for subdir in os.listdir(spaceNetDirectory) if
    #              os.path.isdir(os.path.join(spaceNetDirectory, subdir))]

    listOfAOIs = [srcSpaceNetDirectory]
    srcImageryDirectory = args.srcImageryDirectory  # 'PAN', 'MUL, 'MUL-PanSharpen', 'RGB-PanSharpen'
    if args.spacenetVersion == 2:
        geojsonDirectory = os.path.join('geojson', args.geoJsonDirectory) # 'geojson/buildings/'
    elif args.spacenetVersion == 1:
        geojsonDirectory = os.path.join('vectordata','geojson')
        # 'vectordata/geojson'
    else:
        print('Bad Spacenet Version Submitted,  Version {} is not supported'.foramt(args.spacenetVersion))

    if args.convertTo8Bit:

        outputDataType = 'Byte'
        outputFileType = args.outputFileType

    else:
        outputDataType = ''
        outputFileType = ''

    if args.outputDirectory == 'annotations':
        fullPathAnnotationsDirectory = os.path.join(srcSpaceNetDirectory, 'annotations')
    else:
        fullPathAnnotationsDirectory = args.outputDirectory

    for aoiSubDir in listOfAOIs:
        fullPathSubDir = os.path.join(srcSpaceNetDirectory, aoiSubDir)

        ## Create Annotations directory
        #fullPathAnnotationsDirectory = os.path.join(fullPathSubDir, annotationsDirectory)
        if not os.path.exists(fullPathAnnotationsDirectory):
            os.makedirs(fullPathAnnotationsDirectory)
        if not os.path.exists(os.path.join(fullPathAnnotationsDirectory, 'annotations')):
            os.makedirs(os.path.join(fullPathAnnotationsDirectory, 'annotations'))

        fullPathImageDirectory = os.path.join(fullPathSubDir, srcImageryDirectory)
        fullPathGeoJsonDirectory = os.path.join(fullPathSubDir, geojsonDirectory)

        listofRaster = sorted(glob.glob(os.path.join(fullPathImageDirectory, '*.tif')))
        listofgeojson = sorted(glob.glob(os.path.join(fullPathGeoJsonDirectory, '*.geojson')))

        print('fullpathImageDirectory = {}'.format(fullPathImageDirectory))
        print('fullpathGeoJsonDirectory = {}'.format(fullPathGeoJsonDirectory))
        if len(listofRaster) != len(listofgeojson):
            print('Error lists do not match fix source errors')

            break

        else:

            for rasterImage, geoJson in zip(listofRaster, listofgeojson):

                chipSummaryList = processRasterChip(rasterImage, srcImageryDirectory,
                                                    geoJson, args.geoJsonDirectory,
                                                    outputDirectory=fullPathAnnotationsDirectory,
                                                    imagePixSize=args.imgSizePix, clipOverlap=0.0, randomClip=False,
                                                    minpartialPerc=0.0,
                                                    outputPrefix=''
                                                    )

                entryListTmp = processChipSummaryList(chipSummaryList,
                                                      outputDirectory=os.path.join(fullPathAnnotationsDirectory, 'annotations'),
                                                      annotationType=args.annotationType,
                                                      outputFormat=outputFileType,
                                                      outputPixType=outputDataType,
                                                      datasetName='spacenetV2',
                                                      folder_name='folder_name',
                                                      bboxResize= args.boundingBoxResize
                                       )
                print(entryListTmp)
                entryList.extend(entryListTmp)

    createTrainTestSplitSummary(entryList,
                                trainTestSplit=args.trainTestSplit,
                                outputDirectory=fullPathAnnotationsDirectory,
                                annotationSummaryPrefix=args.annotationSummaryPrefix,
                                annotationType=args.annotationType,
                                shuffleList=True
                                )


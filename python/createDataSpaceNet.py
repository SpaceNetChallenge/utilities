import os
import sys
import glob
from osgeo import gdal
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
        print rasterFileList
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


def processChipSummaryList(chipSummaryList, outputDirectory='', annotationType='PASCAL VOC2012', outputFormat='GTiff',
                           outputPixType='',
                           datasetName='spacenetV2',
                           folder_name='folder_name'
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



        if annotationType=='PASCAL VOC2012':
            entry = lT.geoJsonToPascalVOC(annotationName, chipSummary['geoVectorName'], chipSummary['rasterSource'],
                                          dataset='spacenetV2',
                                          folder_name='spacenetV2',
                                          annotationStyle=annotationType,
                                          segment=True,
                                          bufferSizePix=2.5,
                                          convertTo8Bit=convertTo8Bit,
                                          outputPixType=outputPixType
                                          )
        elif annotationType=='YOLO':
            entry = lT.geoJsonToYolo(annotationName, chipSummary['geoVectorName'], chipSummary['rasterSource'],
                                           dataset='spacenetV2',
                                           folder_name='spacenetV2',
                                           annotationStyle=annotationType,
                                           convertTo8Bit=convertTo8Bit,
                                           outputPixType=outputPixType,
                                           outputFormat=outputFormat
                                  )

        elif annotationType=='MNC':
            print('MNC is not supported yet')
            return -1
        else:
            print ("Annotation Type = {} is not supported yet".format(annotationType))
            return -1



        entryList.append(entry)

        return entryList



if __name__ == '__main__':

    #python createDataSpaceNet.py /data/spacenet_sample/AOI_2_Vegas_Train/ RGB-PanSharpen \
    #                             --outputDirectory /data/spacenet_sample/annotations/ \
    #                             --imgSizePix 416
    parser = argparse.ArgumentParser(description='Process SrcData for Region ComputerVision Dataset')
    parser.add_argument("srcSpaceNetFolder", help="location of Spacenet AOI Data i.e. '/path/to/AOI_2_Vegas")
    parser.add_argument("srcImageryDirectory", help="folder to look for imagery in i.e. 'RGB-PanSharpen'")

    parser.add_argument("--geoJsonDirectory", help="name of geojson folder typedirectory to use"
                                            "i.e. 'buildings'",
                        default='buildings')
    parser.add_argument("--outputDirectory", help="Location To place processed Files"
                                                 "If not used output directory will be"
                                                 "os.path.join(srcSpacenetFolder, 'annotations')",
                                                default='annotations')
    parser.add_argument("--imgSizePix",
                        help="set the dimensions of the square image in pixels"
                             "Default is -1, images are not modified",
                        type=int,
                        default=-1)
    parser.add_argument("--annotationType",
                        help="Set the annotationType.  Currently Supported is YOLO and 'PASCAL VOC'"
                             "default is 'PASCAL VOC'",
                        default='PASCAL VOC2012')
    parser.add_argument("--convertTo8Bit",
                        help='Convert Image from Native format to 8bit',
                        action='store_true')

    parser.add_argument("--featureName",
                        help='Type of feature to be summarized by csv (i.e. Building)',
                        type=str,
                        default='Buildings')

    args = parser.parse_args()

    entryList = []
    srcSpaceNetDirectory = args.srcSpaceNetFolder

    #listOfAOIs = [subdir for subdir in os.listdir(spaceNetDirectory) if
    #              os.path.isdir(os.path.join(spaceNetDirectory, subdir))]

    listOfAOIs = [srcSpaceNetDirectory]
    srcImageryDirectory = args.srcImageryDirectory  # 'PAN', 'MUL, 'MUL-PanSharpen', 'RGB-PanSharpen'
    geojsonDirectory = os.path.join('geojson', args.geoJsonDirectory) # 'geojson/buildings/'


    if args.convertTo8Bit:
        outputDataType = 'Byte'
        outputFileType = 'JPEG'
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
                                                    outputPrefix='')

                entryListTmp = processChipSummaryList(chipSummaryList,
                                                      outputDirectory=os.path.join(fullPathAnnotationsDirectory, 'annotations'),
                                                      annotationType=args.annotationType,
                                                      outputFormat=outputFileType,
                                                      outputPixType=outputDataType,
                                                      datasetName='spacenetV2',
                                                      folder_name='folder_name'
                                       )
                print entryListTmp
                entryList.extend(entryListTmp)


    random.shuffle(entryList)

    splitPoint=int(round(len(entryList)*0.8))
    trainvalList = entryList[0:splitPoint]
    testList     = entryList[splitPoint+1:]


    print('creating trainval.txt {} entries'.format(len(trainvalList)))
    with open(os.path.join(fullPathAnnotationsDirectory, 'trainval.txt'), 'w') as f:
        for entry in trainvalList:
            f.write('{} {}\n'.format(entry['rasterFileName'], entry['annotationName']))

    print('creating test.txt {} entries'.format(len(testList)))
    with open(os.path.join(fullPathAnnotationsDirectory, 'test.txt'), 'w') as f, open(os.path.join(fullPathAnnotationsDirectory, 'test_name_size.txt'), 'w') as fname:
        for entry in testList:
            f.write('{} {}\n'.format(entry['rasterFileName'], entry['annotationName']))
            fname.write('{} {} {}\n'.format(entry['basename'], entry['width'], entry['height']))






















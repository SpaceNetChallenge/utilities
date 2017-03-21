import geoTools as gT
import numpy as np
import os
from osgeo import ogr
import cv2
from math import ceil
import math



def chipImage(rasterFileName, shapeFile, outputDirectory='', numPixelWidth=30,
              finalImageSize = 256, rotationList=np.array([]), outputFormat='.png', rotateNorth=False,
              numPixelSetByLine='', windowSize = 'static',
              truthFile='', sortbyCount='Yes', outputFile=''
              ):
    baseName, fileExt = os.path.splitext(rasterFileName)
    baseDir = os.path.dirname(rasterFileName)
    baseName = os.path.basename(baseName)

    # % Create Directories
    if outputDirectory == '':
        outputDirectory = os.path.join(baseDir, baseName)
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)

    originalDirectory = os.path.join(outputDirectory, 'noRotation')
    if not os.path.isdir(originalDirectory):
        os.mkdir(originalDirectory)
    if rotateNorth:
        rotateNorthDirectory = os.path.join(outputDirectory, 'rotateNorth')
        if not os.path.isdir(rotateNorthDirectory):
            os.mkdir(rotateNorthDirectory)
    if rotationList.size != 0:
        rotateListDirectory = []
        for idx, rotateNum in enumerate(rotationList):
            rotateListDirectory.append(os.path.join(outputDirectory, 'rot_' + str(int(math.ceil(rotateNum)))))
            if not os.path.isdir(rotateListDirectory[idx]):
                os.mkdir(rotateListDirectory[idx])

    fullimg = cv2.imread(rasterFileName)

    # Import Source Shape File
    shapef = ogr.Open(shapeFile, 0)
    chipLayer = shapef.GetLayer()

    # create Truth File
    if truthFile != '':
        shapeT = ogr.Open(truthFile)
        truthLayer = shapeT.GetLayer()
    else:
        # create an output datasource in memory
        # create an output datasource in memory
        outdriver = ogr.GetDriverByName('MEMORY')
        source = outdriver.CreateDataSource('memData')

        # open the memory datasource with write access
        tmp = outdriver.Open('memData', 1)
        pipes_mem = source.CopyLayer(chipLayer, 'truthLayer', ['OVERWRITE=YES'])
        truthLayer = source.GetLayer('truthLayer')
        chipLayer.ResetReading()

    idx = 0
    imgShape = fullimg.shape
    for tmpFeature in chipLayer:
        tmpGeom = tmpFeature.GetGeometryRef()
        angFromNor = tmpFeature.GetField('compassDeg')
        # label   = tmpFeature.GetField("id")
        idx = idx + 1
        # if label is None:
        #    label = 0
        # Convert the layer extent to image pixel coordinates
        centroidX, centroidY, centroidZ = tmpGeom.Centroid().GetPoint()
        envelope = tmpGeom.GetEnvelope() #minX, maxX, minY, maxY

        cX1, cY1 = gT.latlon2pixel(envelope[0], envelope[2], rasterFileName)
        cX2, cY2 = gT.latlon2pixel(envelope[1], envelope[3], rasterFileName)

        if windowSize == 'static':
            numPixelWidth = numPixelWidth
        elif windowSize == 'step4':
            numPixelWidth = ceil((math.sqrt(abs(cY2 - cY1) ** 2 + abs(cX2 - cX1) ** 2) * 1.3) / 4) * 4

        elif windowSize == 'step8':
            numPixelWidth = ceil((math.sqrt(abs(cY2 - cY1) ** 2 + abs(cX2 - cX1) ** 2) * 1.3) / 8) * 8
        elif windowSize == 'adjust':
            numPixelWidth = ceil(math.sqrt(abs(cY2 - cY1) ** 2 + abs(cX2 - cX1) ** 2) * 1.3)

        if numPixelWidth == 0.0:
            continue
        cX, cY = gT.latlon2pixel(centroidY, centroidX, rasterFileName)


        if ((cX > numPixelWidth and cX < imgShape[1] - numPixelWidth) and
                (cY > numPixelWidth and cY < imgShape[0] - numPixelWidth)):

            ulX = cX + ceil(numPixelWidth)
            lrX = cX - ceil(numPixelWidth)
            ulY = cY + ceil(numPixelWidth)
            lrY = cY - ceil(numPixelWidth)

            # geoUtils.pixelToLonLat(xPi)
            polyGeom = gT.returnBoundBox(cX, cY, numPixelWidth,
                                         rasterFileName, targetSR='')
            # calculate # of features within object

            truthLayer.SetSpatialFilter(polyGeom)
            countInBox = math.ceil(truthLayer.GetFeatureCount())
            label = countInBox

            # Calculate the pixel size of the new image
            print(idx)
            res = fullimg[lrY:ulY, lrX:ulX, :]
            resImgShape = res.shape
            ulXdst = resImgShape[1] / 2 + ceil(numPixelWidth / 2)
            lrXdst = resImgShape[1] / 2 - ceil(numPixelWidth / 2)
            ulYdst = resImgShape[0] / 2 + ceil(numPixelWidth / 2)
            lrYdst = resImgShape[0] / 2 - ceil(numPixelWidth / 2)

            # % Save Original Chip
            if finalImageSize == -1:
                dstFinal = res[lrYdst:ulYdst, lrXdst:ulXdst, :]

            else:
                print(finalImageSize)
                print(numPixelWidth)

                dstFinal = cv2.resize(res[lrYdst:ulYdst, lrXdst:ulXdst, :], (finalImageSize, finalImageSize),
                                      interpolation=cv2.INTER_CUBIC)


            fileName = os.path.join(originalDirectory,
                                    baseName + str(int(ulX)) + '_' +
                                    str(int(lrY)) +
                                    '_ROT_' + str(round(angFromNor, 1)) +
                                    '_label_' + str(int(label)) +
                                    outputFormat)
            cv2.imwrite(fileName, dstFinal)
            idstring = baseName + str(int(ulX)) + '_' + str(int(lrY)) + '_ROT_' + str(
                round(angFromNor, 1)) + '_label_' + str(int(label)) + outputFormat


            if rotateNorth:

                M = cv2.getRotationMatrix2D((resImgShape[1] / 2, resImgShape[0] / 2), angFromNor, 1)
                dst = cv2.warpAffine(res, M, (resImgShape[1], resImgShape[0]))
                dstRotateNorth = dst[lrYdst:ulYdst, lrXdst:ulXdst, :]

                if finalImageSize == -1:
                    dstFinal = dstRotateNorth
                else:
                    dstFinal = cv2.resize(dstRotateNorth, (finalImageSize, finalImageSize),
                                          interpolation=cv2.INTER_CUBIC)

                fileName = os.path.join(rotateNorthDirectory,
                                        baseName + str(int(ulX)) + '_' +
                                        str(int(lrY)) +
                                        '_ROT_' + str(int(angFromNor)) +
                                        '_label_' + str(int(label)) +
                                        outputFormat)
                cv2.imwrite(fileName, dstFinal)

            if rotationList.size != 0:
                resImgShape = res.shape
                for idx, rotationNum in enumerate(rotationList):
                    if rotateNorth:
                        M = cv2.getRotationMatrix2D((resImgShape[1] / 2, resImgShape[0] / 2), -rotationNum, 1)
                        dstRotate = cv2.warpAffine(dst, M, (resImgShape[1], resImgShape[0]))
                        dstRotate = dstRotate[lrYdst:ulYdst, lrXdst:ulXdst, :]

                    else:
                        M = cv2.getRotationMatrix2D((resImgShape[1] / 2, resImgShape[0] / 2), -rotationNum, 1)
                        dstRotate = cv2.warpAffine(res, M, (resImgShape[1], resImgShape[0]))
                        dstRotate = dstRotate[lrYdst:ulYdst, lrXdst:ulXdst, :]

                    if finalImageSize == -1:
                        dstFinal = dstRotate
                    else:
                        dstFinal = cv2.resize(dstRotate, (finalImageSize, finalImageSize),
                                              interpolation=cv2.INTER_CUBIC)

                    fileName = os.path.join(rotateListDirectory[idx],
                                            baseName + str(int(ulX)) + '_' +
                                            str(int(lrY)) +
                                            '_ROT_' + str(int(rotationNum)) +
                                            '_label_' + str(int(label)) +
                                            outputFormat)
                    cv2.imwrite(fileName, dstFinal)

    return outputDirectory



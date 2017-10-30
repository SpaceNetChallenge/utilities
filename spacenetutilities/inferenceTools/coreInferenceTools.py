import numpy as np
import rasterio
from rasterio import Affine
from rasterio.warp import reproject, Resampling
import types
import tqdm
from rasterio.windows import Window
from rasterio.coords import BoundingBox
import shapely
import pyproj
from functools import partial
from spacenetutilities import geoTools as gT


def resampleImage(array, spatialScaleFactor, src_meta=[], src_transform=[], src_crs=[], dst_transform=[],
                  newarr=np.array([])):
    # with rasterio.open('path/to/file.tif') as src:
    #   src_meta= src.meta
    # or
    # src_transform = src.affine
    # src_crs = src.crs
    if newarr.size == 0:
        newarr = np.empty(shape=(array.shape[0],  # same number of bands
                                 round(array.shape[1] * spatialScaleFactor),  # 150% resolution
                                 round(array.shape[2] * spatialScaleFactor)))

    # adjust the new affine transform to the 150% smaller cell size
    if src_meta:
        src_transform = src_meta['transform']
        src_crs = src_meta['crs']
    else:
        if not src_transform or not src_crs:
            print('Error src transform and src_crs must be defined if src file is not defined')
            return -1


    if not dst_transform:
        dst_transform = Affine(src_transform.a / spatialScaleFactor, src_transform.b, src_transform.c,
                               src_transform.d, src_transform.e / spatialScaleFactor, src_transform.f)

    reproject(
        array, newarr,
        src_transform=src_transform,
        dst_transform=dst_transform,
        src_crs=src_crs,
        dst_crs=src_crs,
        resample=Resampling.bilinear)

    return newarr, dst_transform, src_crs, src_transform, src_crs


def imageCombiner(listOfImages, image, strideInPixels, windowSize=(256, 256), debug=False):
    if debug:
        print(image.shape)
    newImage = np.empty(shape=(1,
                               image.shape[1],
                               image.shape[2]))

    newImageCount = np.empty(shape=(1,
                                    image.shape[1],
                                    image.shape[2]))
    if isinstance(listOfImages, types.GeneratorType):
        listOfImages = list(listOfImages)
    idx = 0
    for y in range(0, image.shape[1], strideInPixels):
        for x in range(0, image.shape[2], strideInPixels):
            if y + windowSize[1] > image.shape[1]:
                y = image.shape[1] - windowSize[1]

            if x + windowSize[0] > image.shape[2]:
                x = image.shape[2] - windowSize[0]

            newImage[:, y:y + windowSize[1], x:x + windowSize[0]] += listOfImages[idx]

            newImageCount[:, y:y + windowSize[1], x:x + windowSize[0]] += 1

            idx += 1

    return newImage, newImageCount


def sceneTilerGeneratorCount(datapathPSMUL, sceneSize, percentOverlap=1.0, debug=False):
    readWindowList = []
    with rasterio.open(datapathPSMUL, 'r') as src:

        src_meta = src.meta
        srcShape = src.shape

        strideInPixels = round(percentOverlap * sceneSize[0])
        for y in range(0, srcShape[0], strideInPixels):
            for x in range(0, srcShape[1], strideInPixels):
                # ((row_start, row_stop), (col_start, col_stop))
                readWindow = ((y, y+sceneSize[1]), (x, x+sceneSize[0]))
                readWindowList.append(readWindow)

    return readWindowList


def sceneTilerGenerator(datapathPSMUL, sceneSize, percentOverlap=1.0, debug=False):

    with rasterio.open(datapathPSMUL, 'r') as src:

        src_meta = src.meta
        srcShape = src.shape

        strideInPixels = round(percentOverlap * sceneSize[0])
        for y in range(0, srcShape[0], strideInPixels):
            for x in range(0, srcShape[1], strideInPixels):
                # ((row_start, row_stop), (col_start, col_stop))
                readWindow = ((y, y+sceneSize[1]), (x, x+sceneSize[0]))

                # yield the current window
                array = src.read(window=readWindow)
                if debug:
                    print(array.shape)

                tmp_meta = src_meta
                tmp_meta['transform'] = src.window_transform(readWindow)
                tmp_meta['bounds'] = BoundingBox(*src.window_bounds(readWindow))

                yield (array, tmp_meta, readWindow)

def createImageBoundsList(imageShape,
                   strideInPixels,
                   windowSize=(256, 256),
                   debug=False
                   ):
    # imageShape = rastrio.open().src.shape
    # slide a window across the image
    # returns list of tuples of (minx, miny, maxx, maxy)
    imageBoundsList = []
    for y in range(0, imageShape[1], strideInPixels):
        for x in range(0, imageShape[2], strideInPixels):

            if y + windowSize[1] > imageShape[1]:
                y = imageShape[1] - windowSize[1]

            if x + windowSize[0] > imageShape[2]:
                x = imageShape[2] - windowSize[0]

            imageBoundsList.append((x, y, x + windowSize[0], y + windowSize[1]))

    if debug:
        print('Processing {} windows'.format(len(imageBoundsList)))


    return imageBoundsList

def returnImgArrayForTensorFlowFromRasterio(dataArray,
                                            bandsToInclude=[0,1,2],
                                            minPercentToClip=0,
                                            maxPercentToClip=100,
                                            maxBandValue=255,
                                            dataFormat=np.uint8,
                                            transposeOrder=[1,2,0]
                                            ):

    return returnImgArrayFromArray(dataArray,
                                   bandsToInclude=bandsToInclude,
                                   minPercentToClip=minPercentToClip,
                                   maxPercentToClip=maxPercentToClip,
                                   maxBandValue=maxBandValue,
                                   dataFormat=dataFormat,
                                   transposeOrder=transposeOrder)

def returnImgArrayFromArray(dataArray, bandsToInclude=[], maxBandValue=-1,
                            minPercentToClip=0,
                            maxPercentToClip=100,
                            verbose=False,
                            dataFormat=np.uint8,
                            transposeOrder=[1,2,0]):


    ## Set maxBandValue to -1 to perform zero normalization
    maxImgValue = np.percentile(dataArray, maxPercentToClip)
    if bandsToInclude:
        img = dataArray[bandsToInclude]
    else:
        img = dataArray

    scale_factor = 1
    if maxBandValue>-1:
        #maxImgValue = np.percentile(img, maxPercentToClip)
        scale_factor = maxBandValue/maxImgValue

    img = scale_factor*img
    img = np.clip(img, 0, maxBandValue)



    # transpose images (height, width, bandNum)
    imgTranspose = img.transpose(transposeOrder)
    imgTranspose = imgTranspose.astype(dataFormat)

    return imgTranspose




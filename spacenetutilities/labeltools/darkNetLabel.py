import os
from spacenetutilities.labeltools import coreLabelTools as clT
import rasterio
from shapely.geometry import box

def convertPixDimensionToPercent(size, shpBox):
    '''Input = image size: (w,h), box: [x0, y0, x1, y1]'''
    #TODO change box to use shapely bounding box format

    # calculate normalized Spaceing
    dw = 1./size[0]
    dh = 1./size[1]

    #calculate center of bounding box
    xmid = (box[0] + box[2])/2.0
    ymid = (box[1] + box[3])/2.0

    #calculate width
    w0 = box[2] - box[0]
    h0 = box[3] - box[1]

    # calculate final values
    x = xmid*dw
    y = ymid*dh
    w = w0*dw
    h = h0*dh

    return (x,y,w,h)


def writeTODARKNETLabel(labelFileName, imageDescriptionDict, objectDictList):

    #TODO implement multiclass
    imageSize = (imageDescriptionDict['size']['width'], imageDescriptionDict['size']['height'])
    with open(labelFileName, 'w') as f:
        for objectDict in objectDictList:

            lineOutput = convertPixDimensionToPercent(imageSize,
                                                      [objectDict['bndBox']['xmin'],
                                                                   objectDict['bndBox']['ymin'],
                                                                   objectDict['bndBox']['xmax'],
                                                                   objectDict['bndBox']['ymax']]
                                                      )

            classNum = 0
            f.write('{} {} {} {} {}\n'.format(classNum,
                                          lineOutput[0],
                                          lineOutput[1],
                                          lineOutput[2],
                                          lineOutput[3])
                    )

    return labelFileName

def geoJsonToDARKNETLabel(xmlFileName, geoJson, rasterImageName, im_id='',
                           dataset ='SpaceNet',
                           folder_name='spacenet',
                           annotationStyle = 'PASCAL VOC2012',
                           segment=True,
                           bufferSizePix=2.5,
                           convertTo8Bit=True,
                           outputPixType='Byte',
                           outputFormat='GTiff',
                           bboxResize=1.0,
                           objectType='building',
                           objectTypeField=''):

    imageDescriptionDict, objectDictList = clT.geoJsontoDict(geoJson, rasterImageName, datasetName=dataset,
                  annotationStyle=annotationStyle,
                  bboxResize=bboxResize,
                  objectType=objectType,
                  objectTypeField=objectTypeField,
                  objectPose='Left',
                  objectTruncatedField='',
                  objectDifficultyField=''
                  )

    xmlFileName = writeTODARKNETLabel(xmlFileName, imageDescriptionDict, objectDictList)

    return xmlFileName


def geoJsonToDARKNET(xmlFileName, geoJson, rasterImageName, im_id='',
                     dataset='SpaceNet',
                     folder_name='spacenet',
                     annotationStyle='DARKNET',
                     segment=True,
                     bufferSizePix=2.5,
                     convertTo8Bit=True,
                     outputPixType='Byte',
                     outputFormat='JPEG',
                     bboxResize=1.0,
                     objectType='building',
                     objectTypeField='',
                     outputRasterName=''):



    geoJsonToDARKNETLabel(xmlFileName, geoJson, rasterImageName, im_id='',
                                dataset=dataset,
                                folder_name=folder_name,
                                annotationStyle=annotationStyle,
                                segment=segment,
                                bufferSizePix=bufferSizePix,
                                convertTo8Bit=convertTo8Bit,
                                outputPixType=outputPixType,
                                outputFormat=outputFormat,
                                bboxResize=bboxResize,
                                objectType=objectType,
                                objectTypeField=objectTypeField)


    if convertTo8Bit:
        if outputRasterName=='':
            if outputFormat == 'JPEG':
                outputImageName = xmlFileName.replace('.xml', '_8bit.jpg')
            else:
                outputImageName = xmlFileName.replace('.xml', '_8bit.tif')

            clT.convertGTiffTo8Bit(rasterImageName, outputImageName, outputFormat=outputFormat)
    else:
        outputImageName = rasterImageName


    with rasterio.open(rasterImageName) as src:
        src_meta = src.meta.copy()
        src_profile = src.profile.copy()


    entry = {'rasterFileName': outputImageName,
             'geoJsonFileName': geoJson,
             'annotationName': xmlFileName,
             'width': src_meta['width'],
             'height': src_meta['height'],
             'depth': src_meta['count'],
             'src_meta': src_meta,
             'basename': os.path.splitext(os.path.basename(rasterImageName))[0]
             }

    return entry



# SpaceNet Utilities

This repository has two python packages, geoTools and evalTools. The geoTools packages is intended to assist in the preprocessing of [SpaceNet](https://spacenetchallenge.github.io/) satellite imagery data corpus hosted on [SpaceNet on AWS](https://aws.amazon.com/public-datasets/spacenet/) to a format that is consumable by machine learning algorithms. The evalTools package is used to evaluate the effectiveness of object detection algorithms using ground truth.
This is version 2.0 and has been updated with more capabilities to 
## Download Instructions

Further download instructions for the [SpaceNet Dataset](https://github.com/SpaceNetChallenge/utilities/tree/master/content/download_instructions) can be found [here](https://github.com/SpaceNetChallenge/utilities/tree/master/content/download_instructions)


## Evaluation Metric
The evaluation metric for this competition is an F1 score with the matching algorithm inspired by Algorithm 2 in the [ILSVRC paper applied to the detection of building footprints](https://arxiv.org/pdf/1409.0575v3.pdf). For each building there is a geospatially defined polygon label to represent the footprint of the building. A SpaceNet entry will generate polygons to represent proposed building footprints.  Each proposed building footprint is either a “true positive” or a “false positive”.

* The proposed footprint is a “true positive” if the proposal is the closest (measured by the IoU) proposal to a labeled polygon AND the IoU between the proposal and the label is about the prescribed threshold of 0.5.
* Otherwise, the proposed footprint is a “false positive”.

There is at most one “true positive” per labeled polygon.
The measure of proximity between labeled polygons and proposed polygons is the Jaccard similarity or the “Intersection over Union (IoU)”, defined as:

![alt text](https://github.com/SpaceNetChallenge/utilities/blob/master/content/IoU.jpg "IoU")

The value of IoU is between 0 and 1, where closer polygons have higher IoU values.

The F1 score is the harmonic mean of precision and recall, combining the accuracy in the precision measure and the completeness in the recall measure. For this competition, the number of true positives and false positives are aggregated over all of the test imagery and the F1 score is computed from the aggregated counts.

For example, suppose there are N polygon labels for building footprints that are considered ground truth and suppose there are M proposed polygons by an entry in the SpaceNet competition.  Let tp denote the number of true positives of the M proposed polygons.  The F1 score is calculated as follows:

![alt text](https://github.com/SpaceNetChallenge/utilities/blob/master/content/F1.jpg "IoU")

The F1 score is between 0 and 1, where larger numbers are better scores.

Hints:
* The images provided could contain anywhere from zero to multiple buildings.
* All proposed polygons should be legitimate (they should have an area, they should have points that at least make a triangle instead of a point or a line, etc).
* Use the [metric implementation code](https://github.com/SpaceNetChallenge/utilities/blob/master/python/evaluateScene.py) to self evaluate.
To run the metric you can use the following command.
```
python python/evaluateScene.py /path/to/SpaceNetTruthFile.csv \
                               /path/to/SpaceNetProposalFile.csv \
                               --resultsOutputFile /path/to/SpaceNetResults.csv
```

## Data Transformation Code

To make the Spacenet dataset easier to use we have created a tool createDataSpaceNet.py
This tool currently supports the creation of datasets with annotation to support 3 Formats
1. [PASCAL VOC2012](http://host.robots.ox.ac.uk/pascal/VOC/)
2. [Darknet](https://pjreddie.com/darknet/yolo/)
3. [Segmenation Boundaries Dataset (SBD)](http://home.bharathh.info/pubs/codes/SBD/download.html)

It will create the appropriate annotation files and a summary trainval.txt and test.txt in the outputDirectory

### Create an PASCAL VOC2012 Compatiable Dataset
The final product will have image dimensions of 420 pixels
```
python python/createDataSpaceNet.py /path/to/spacenet_sample/AOI_2_Vegas_Train/ \
           --srcImageryDirectory RGB-PanSharpen
           --outputDirectory /path/to/spacenet_sample/annotations/ \
           --annotationType PASCALVOC2012 \
           --imgSizePix 400

```
### Changing the raster format
Some GIS Images have 16-bit pixel values which openCV has trouble with.  createDataSpaceNet.py can convert the 16bit GeoTiff to an 8bit GeoTiff or 8bit JPEG 

To create the 8bit GeoTiff
```
python python/createDataSpaceNet.py /path/to/spacenet_sample/AOI_2_Vegas_Train/ \
           --srcImageryDirectory RGB-PanSharpen
           --outputDirectory /path/to/spacenet_sample/annotations/ \
           --annotationType PASCALVOC2012 \
           --convertTo8Bit \
           --outputFileType GTiff \
           --imgSizePix 400
    
```

To create the 8bit JPEG
```
python python/createDataSpaceNet.py /path/to/spacenet_sample/AOI_2_Vegas_Train/ \
           --srcImageryDirectory RGB-PanSharpen
           --outputDirectory /path/to/spacenet_sample/annotations/ \
           --annotationType PASCALVOC2012 \
           --convertTo8Bit \
           --outputFileType JPEG \
           --imgSizePix 400

```

For more Features
```
python python/createDataSpaceNet.py -h

```



## Use our Docker Container
We have created two Docker files at /docker/standalone/cpu and /docker/standalone/gpu
These Dockerfiles will build a docker container with all packages neccessary to run the package

More documenation to follow


## Dependencies
All dependencies can be found in [requirements.txt](./python/requirements.txt)

## License
See [LICENSE](./LICENSE).

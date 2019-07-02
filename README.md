# This repository is no longer being updated. Future development of code tools for geospatial machine learning analysis will be done at https://github.com/cosmiq/solaris.

# SpaceNet Utilities

This repository has three python packages, geoTools and evalTools and labelTools. The geoTools packages is intended to assist in the preprocessing of [SpaceNet](https://spacenetchallenge.github.io/) satellite imagery data corpus hosted on [SpaceNet on AWS](https://aws.amazon.com/public-datasets/spacenet/) to a format that is consumable by machine learning algorithms. 
The evalTools package is used to evaluate the effectiveness of object detection algorithms using ground truth.
The labelTools package assists in transfering geoJson labels into common label schemes for machine learning frameworks
This is version 3.0 and has been updated with more capabilities to allow for computer vision applications using remote sensing data

## Download Instructions
Further download instructions for the [SpaceNet Dataset](https://github.com/SpaceNetChallenge/utilities/tree/master/content/download_instructions) can be found [here](https://github.com/SpaceNetChallenge/utilities/tree/master/content/download_instructions)



## Installation Instructions
Several packages require binaries to be installed before pip installing the other packages.  Conda is a simple way to install everything and their dependencies

* Install GDAL binaries and scripts 
```commandline
conda install -c conda-forge gdal
```
* Install [Rtree](http://toblerity.org/rtree/install.html) 
```commandline 
conda install -c conda-forge rtree
```

* Install [pyproj](https://pypi.python.org/pypi/pyproj)
```commandline 
conda install -c conda-forge pyproj
```
* Install [geopandas](https://pypi.python.org/pypi/geopandas)
```commandline 
conda install -c conda-forge geopandas
```

* Install [shapely](https://pypi.python.org/pypi/shapely)
```commandline 
conda install -c conda-forge shapely
```

* Install [rasterio](https://pypi.python.org/pypi/rasterio)
```commandline 
conda install -c conda-forge rasterio
```


* Pip Install from github 
```commandline
    git clone -b spacenetV3 https://github.com/SpaceNetChallenge/utilities.git
    cd utilities
    pip install -e .
    
```

or 
```commandline
    pip install --upgrade git+https://github.com/SpaceNetChallenge/utilities.git
```




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

## Using SpaceNet Utilities to Process Imagery and Vector Data
The SpaceNet imagery provided for this challenge must be processed and transformed into a deep-learning compatible format.  SpaceNet utilites helps to achieve this transformation.  A traditional implementation strategy may look similar to this:

1. Chipping and clipping SpaceNet imagery into smaller areas to allow for deep learning consumption (create_spacenet_AOI.py)

2. Split imagery and vector datasets (such as building or road labels) into training, testing, and validation datasets randomly (splitAOI_Train_Test_Val.py)

3. Easily add or update your vector datasets to seamlessly match existing SpaceNet imagery chips. (externalVectorProcessing.py)

4. Translate SpaceNet image chips and vector data into various machine learning and deep learning consumable formats such as PASCAL VOC2012, DarkNet, or Semantic Boundaries Dataset (SBD). (createDataSpaceNet.py)

5. Evaluate your deep learning outputs against validation datasets to determine your results' accuracy and estimate the amount of comission and omission errors ocurring. (evaluateScene.py)

6.  Various other maintenance utility scripts to enhance ease of use.


## Chipping Imagery Code
The script create_spacenet_AOI.py is used to create a SpaceNet competition dataset (chips) from a larger imagery dataset.  Its base function is to create an N x N meters (or N x N pxiels) chip with associated object labels (such as buildings or roads).  The script will only create chips in the area where labeled items exist, thus saving space and reducing computational intensity.

The script requires a few pre-processing steps, a recommended process for this would be

1. Build a VRT file to point to the source imagery.  A VRT is essentailly a virtual mosaic that links all the files, but does not build an entirely new (and monsterous) mosaic.  http://www.gdal.org/gdalbuildvrt.html is one of the best ways to do this easily.  A seperate VRT should be built for each type of imagery data you plan to use (Ex: Pan, Multi-spectral, etc..)
    
2.  Build two CSV pointer files that point to the specific location of both your imagery VRT's and the labeled vector data (buildings, roads, etc..).  This file will have two columns with NO headers. 

```     
    Example raster CSV:
   
    Column A:   Column B:
    PAN         C:/SpaceNet/Imagery/Vegas_PAN.vrt
    MUL         C:/SpaceNet/Imagery/Vegas_MUL.vrt
    MUL-PS      C:/SpaceNet/Imagery/Vegas_MUL-PS.vrt
 
    Example vector CSV:
    
    Column A:   Column B:
    Buildings   C:/SpaceNet/Vector/Vegas_BuildingLabels.geojson
```

    
Script Inputs:
1. CSV of raster imagery VRT locations
2. CSV of vector labels (geojson)
3. SRC_Outline- The outline of where labelling is ocurring in a geojson format
4.  Other optional inputs are also availble, more information on these can be gleaned by looking into the raw code itself or using      the -h help feature.
    
   
    
The script will then chip and clip the source SpaceNet imagery and vector labels.  An example prompt of running the script is as follows:
```commandline 
python create_spacenet_AOI.py --srcOutline /data/vectorData/Shanghai_AOI_Fixed.geojson --outputDirectory /data/output --AOI_Name Shanghai --AOI_Num 4 --createSummaryCSV --featureName Building /data/AOI_4_Shanghai_srcRasterList.csv /data/AOI_4_Shanghai_srcVectorList.csv 
```

This will output chipped imagery into your outputDirectory folder for further usage.



## Data Transformation Code

To make the Spacenet dataset easier to use we have created a tool createDataSpaceNet.py
This tool currently supports the creation of datasets with annotation to support 3 Formats
1. [PASCAL VOC2012](http://host.robots.ox.ac.uk/pascal/VOC/)
2. [Darknet](https://pjreddie.com/darknet/yolo/)
3. [Segmenation Boundaries Dataset (SBD)](http://home.bharathh.info/pubs/codes/SBD/download.html)

It will create the appropriate annotation files and a summary trainval.txt and test.txt in the outputDirectory

### Create an PASCAL VOC2012 Compatiable Dataset
The final product will have image dimensions of 400 pixels
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
All dependencies can be found in the docker file [Dockerfile](./docker/standalone/gpu/Dockerfile)

## License
See [LICENSE](./LICENSE).

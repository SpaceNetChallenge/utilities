# SpaceNet Utilities

This repository has two python packages, geoTools and evalTools. The geoTools packages is intended to assist in the preprocessing of [SpaceNet](https://aws.amazon.com/public-data-sets/spacenet/) satellite imagery data corpus to a format that is consumable by machine learning algorithms. The evalTools package is used to evaluate the effectiveness of object detection algorithms using ground truth.

## Dependencies
All dependencies can be found in [requirements.txt](./python/requirements.txt)

## License
See [LICENSE](./LICENSE).

## Metric
The evaluation metric for this competition is an F1 score with the matching algorithm inspired by Algorithm 2 in the [ILSVRC paper applied to the detection of building footprints](https://arxiv.org/pdf/1409.0575v3.pdf). For each building there is a geospatially defined polygon label to represent the footprint of the building. A SpaceNet entry will generate polygons to represent proposed building footprints.  Each proposed building footprint is either a “true positive” or a “false positive”.

* The proposed footprint is a “true positive” if the proposal is the closest (measured by the IoU) proposal to a labeled polygon AND the IoU between the proposal and the label is about the prescribed threshold of 0.5.
* Otherwise, the proposed footprint is a “false positive”.

There is at most one “true positive” per labeled polygon.
The measure of proximity between labeled polygons and proposed polygons is the Jaccard similarity or the “Intersection over Union (IoU)”, defined as:

![alt text](https://github.com/SpaceNetChallenge/utilities/blob/master/content/IoU.jpg "IoU")

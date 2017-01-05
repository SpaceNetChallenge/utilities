# SpaceNet Utilities

Instruction for the download of [SpaceNet](https://aws.amazon.com/public-data-sets/spacenet/) satellite imagery data corpus to a format that is consumable by machine learning algorithms.
 
 

## Dependencies
Required that AWS CLI is installed and that an active AWS account with credit card is associated with the aws cli

Configure the AWS CLI using aws configure

SpaceNet AWS Structure
```
s3://spacenet-dataset/
-- AOI_1_Rio
    |-- processedData
    |   -- processedBuildingLabels.tar.gz  # Compressed 3band and 8band 200m x 200m tiles with associated building foot print labels
                                           # This dataset is the Training Dataset for the first [Top Coder Competition](https://community.topcoder.com/longcontest/?module=ViewProblemStatement&rd=16835&pm=14439)
    `-- srcData
        |-- rasterData
        |   |-- 3-Band.tar.gz # 3band (RGB) Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by [WorldView-2](http://satimagingcorp.s3.amazonaws.com/site/pdf/WorldView-2_datasheet.pdf)
        |    -- 8-Band.tar.gz # 8band Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by [WorldView-2](http://satimagingcorp.s3.amazonaws.com/site/pdf/WorldView-2_datasheet.pdf)
         -- vectorData
            |-- Rio_BuildingLabels.tar.gz # Source Dataset that contains Building the building foot prints traced from the Mosaic
            |-- Rio_HGIS_Metro.gdb.tar.gz # Source Point of Interest Dataset in GeoDatabase Format.  Best if Used with ESRI
             -- Rio_HGIS_Metro_extract.tar # Source Point of Interest Dataset in GeoJSON with associated .jpg.  Easy to Use without ESRI toolset
```

#To download the Imagery

##To download processed 200mx200m Tiles with associated building foot prints for building foot print extraction tests do the following
```
## Warning this file is 3.4 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/processedData/processedBuildingLabels.tar.gz --request-payer requester processedBuildingLabels.tar.gz
```


##To download the Source Imagery Mosaic
```
## Warning this file is 2.3 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/rasterData/3-Band.tar.gz --request-payer requester 3-Band.tar.gz
## Warning this file is 6.5 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/rasterData/8-Band.tar.gz --request-payer requester 8-Band.tar.gz
```

##To download the Source Vector Data for the Building Extraction Challenge
```
## Warning this file is 0.18 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/vectorData/Rio_BuildingLabels.tar.gz --request-payer requester Rio_BuildingLabels.tar.gz

```

##To download the Rio Point of Interest Dataset in ESRI GeoDatabase Form
```
## Warning this file is 31 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/vectorData/Rio_HGIS_Metro.gdb.tar.gz --request-payer requester Rio_HGIS_Metro.gdb.tar.gz

```

##To download the Rio Point of Interest Dataset Extracted into GeoJSONs with associated .jpg
```
## Warning this file is 29 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/vectorData/Rio_HGIS_Metro_extract.tar --request-payer requester Rio_HGIS_Metro_extract.tar

```




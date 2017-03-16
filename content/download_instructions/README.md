## Hosting
[SpaceNet](https://aws.amazon.com/public-datasets/spacenet/) is a corpus of commercial satellite imagery and labeled training data to use for machine learning research. The dataset is currently hosted as an [Amazon Web Services (AWS) Public Dataset](https://aws.amazon.com/public-datasets/).

## Catalog
1. Area of Interest 1 (AOI 1) - Location: Rio de Janeiro. 50cm imagery collected from DigitalGlobe’s [WorldView-2 satellite](http://satimagingcorp.s3.amazonaws.com/site/pdf/WorldView-2_datasheet.pdf). The dataset includes building footprints and 8-band multispectral data.
2. Area of Interest 2 (AOI 2) - Location: Vegas. 30cm imagery collected from DigitalGlobe’s [WorldView-3 satellite](https://www.spaceimagingme.com/downloads/sensors/datasheets/DG_WorldView3_DS_2014.pdf). The dataset includes building footprints and 8-band multispectral data.
3. Area of Interest 3 (AOI 3) - Location: Paris. 30cm imagery collected from DigitalGlobe’s [WorldView-3 satellite](https://www.spaceimagingme.com/downloads/sensors/datasheets/DG_WorldView3_DS_2014.pdf). The dataset includes building footprints and 8-band multispectral data.
4. Area of Interest 4 (AOI 4) - Location: Shanghai. 30cm imagery collected from DigitalGlobe’s [WorldView-3 satellite](https://www.spaceimagingme.com/downloads/sensors/datasheets/DG_WorldView3_DS_2014.pdf). The dataset includes building footprints and 8-band multispectral data.
5. Area of Interest 5 (AOI 5) - Location: Khartoum. 30cm imagery collected from DigitalGlobe’s [WorldView-3 satellite](https://www.spaceimagingme.com/downloads/sensors/datasheets/DG_WorldView3_DS_2014.pdf). The dataset includes building footprints and 8-band multispectral data.
6. Point of Interest (POI) Dataset- Location: Rio de Janeiro. The dataset includes POIs.

## Dependencies
The [AWS Command Line Interface (CLI)](https://aws.amazon.com/cli/) must be installed with an active AWS account. Configure the AWS CLI using 'aws configure'


## SpaceNet Simple Storage Service (S3) Directory Structure (AOI 1)
```
s3://spacenet-dataset/
-- AOI_1_Rio
    |-- processedData
    |   -- processedBuildingLabels.tar.gz  # Compressed 3band and 8band 200m x 200m tiles with associated building foot print labels                                 # This dataset is the Training Dataset for the first Top Coder Competition
    `-- srcData
        |-- rasterData
        |   |-- 3-Band.tar.gz # 3band (RGB) Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by WorldView-2
        |    -- 8-Band.tar.gz # 8band Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by WorldView-2
         -- vectorData
            |-- Rio_BuildingLabels.tar.gz # Source Dataset that contains Building the building foot prints traced from the Mosaic
            |-- Rio_HGIS_Metro.gdb.tar.gz  # Source Point of Interest Dataset in GeoDatabase Format.  Best if Used with ESRI
             -- Rio_HGIS_Metro_extract.tar # Source Point of Interest Dataset in GeoJSON with associated .jpg.  Easy to Use without ESRI toolset
-- AOI_1_Rio
    |-- processedData
    |   -- processedBuildingLabels.tar.gz  # Compressed 3band and 8band 200m x 200m tiles with associated building foot print labels                                 # This dataset is the Training Dataset for the first Top Coder Competition
    `-- srcData
        |-- rasterData
        |   |-- 3-Band.tar.gz # 3band (RGB) Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by WorldView-2
        |    -- 8-Band.tar.gz # 8band Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by WorldView-2
         -- vectorData
            |-- Rio_BuildingLabels.tar.gz # Source Dataset that contains Building the building foot prints traced from the Mosaic
            |-- Rio_HGIS_Metro.gdb.tar.gz # Source Point of Interest Dataset in GeoDatabase Format.  Best if Used with ESRI
             -- Rio_HGIS_Metro_extract.tar # Source Point of Interest Dataset in GeoJSON with associated .jpg.  Easy to Use without ESRI toolset
-- AOI_1_Rio
    |-- processedData
    |   -- processedBuildingLabels.tar.gz  # Compressed 3band and 8band 200m x 200m tiles with associated building foot print labels                                 # This dataset is the Training Dataset for the first Top Coder Competition
    `-- srcData
        |-- rasterData
        |   |-- 3-Band.tar.gz # 3band (RGB) Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by WorldView-2
        |    -- 8-Band.tar.gz # 8band Raster Mosaic for Rio De Jenairo area (2784 sq KM) collected by WorldView-2
         -- vectorData
            |-- Rio_BuildingLabels.tar.gz # Source Dataset that contains Building the building foot prints traced from the Mosaic
            |-- Rio_HGIS_Metro.gdb.tar.gz # Source Point of Interest Dataset in GeoDatabase Format.  Best if Used with ESRI
             -- Rio_HGIS_Metro_extract.tar # Source Point of Interest Dataset in GeoJSON with associated .jpg.  Easy to Use without ESRI toolset
```

## SpaceNet Simple Storage Service (S3) Directory Structure (AOI 2-5)
```
├── AOI_[Num]_[City]_Train
│   ├── geojson
│   │   └── buildings  # Contains GeoJson labels of buildings for each tile
│   ├── MUL            # Contains Tiles of 8-Band Multi-Spectral raster data from WorldView-3
│   ├── MUL-PanSharpen # Contains Tiles of 8-Band Multi-Spectral raster data pansharpened to 0.3m
│   ├── PAN            # Contains Tiles of Panchromatic raster data from Worldview-3
│   ├── RGB-PanSharpen # Contains Tiles of RGB raster data from Worldview-3
│   └── summaryData    # Contains CSV with pixel based labels for each building in the Tile Set.
```

## Download instructions

### AOI 1 - Rio de Janeiro
To download processed 200mx200m tiles of AOI 1 (3.4 GB) with associated building footprints do the following:
```
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/processedData/processedBuildingLabels.tar.gz --request-payer requester processedBuildingLabels.tar.gz
```
To download the Source Imagery Mosaic (3-band = 2.3 GB and 8-band = 6.5 GB):
```
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/rasterData/3-Band.tar.gz --request-payer requester 3-Band.tar.gz
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/rasterData/8-Band.tar.gz --request-payer requester 8-Band.tar.gz
```
To download the Source Vector Data (0.18 GB):
```
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/vectorData/Rio_BuildingLabels.tar.gz --request-payer requester Rio_BuildingLabels.tar.gz
```

### AOI 2 - Vegas
To download processed 200mx200m tiles of AOI 2 (23 GB) with associated building footprints do the following:
```
aws s3api get-object --bucket spacenet-dataset --key AOI_2_Vegas/AOI_2_Vegas_Train.tar.gz --request-payer requester AOI_2_Vegas_Train.tar.gz
```

### AOI 3 - Paris
To download processed 200mx200m tiles of AOI 3 (5 GB) with associated building footprints do the following:
```
## Warning this file is 5 GB
aws s3api get-object --bucket spacenet-dataset --key AOI_3_Paris/AOI_3_Paris_Train.tar.gz --request-payer requester AOI_3_Paris_Train.tar.gz
```

### AOI 4 - Shanghai
To download processed 200mx200m tiles of AOI 4 (23 GB) with associated building footprints do the following:
```
aws s3api get-object --bucket spacenet-dataset --key AOI_4_Shanghai/AOI_4_Shanghai_Train.tar.gz --request-payer requester AOI_4_Shanghai_Train.tar.gz
```

### AOI 5 - Khartoum
To download processed 200mx200m tiles of AOI 5 (4 GB) with associated building footprints do the following:
```
aws s3api get-object --bucket spacenet-dataset --key AOI_5_Khartoum/AOI_5_Khartoum_Train.tar.gz --request-payer requester AOI_5_Khartoum_Train.tar.gz
```

### Point of Interest Dataset in ESRI GeoDatabase Form (31 GB)
```
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/vectorData/Rio_HGIS_Metro.gdb.tar.gz --request-payer requester Rio_HGIS_Metro.gdb.tar.gz
```

### Point of Interest Dataset Extracted into GeoJSONs with associated .jpg (29 GB)
```
aws s3api get-object --bucket spacenet-dataset --key AOI_1_Rio/srcData/vectorData/Rio_HGIS_Metro_extract.tar --request-payer requester Rio_HGIS_Metro_extract.tar
```




from spaceNetUtilities import labelTools as lT
import argparse
import os
import csv

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process SrcData for Region into small chips')
    parser.add_argument("srcRasterList", help="csv file with path/to/raster.vrt, rasterDescripiton "
                                              "i.e, path/to/AOI_#_Num_3band.vrt, 3band ")
    parser.add_argument("geoJsonList", help="csv file with path/to/vector_buildings.geojson, vectorDescription"
                                              "i.e, path/to/AOI_#_buildings.geojson, buildings")

    parser.add_argument("--srcOutline", help="Vector Area describing extent of labeled area"
                                            "If not specified, All of raster will be assumed to be labeled",
                        default='')

    parser.add_argument("--outputDirectory", help="Location To place processed Files"
                                                 "If not used output directory will be"
                                                 "os.path.join(os.getcwd(), processed)",
                                                default='')
    parser.add_argument("--imgSizeM",
                        help="set the dimensions of the square image in meters  "
                             "Default is 200m",
                        type=int,
                        default=200)
    parser.add_argument("--AOI_Name",
                        help="AOI City Name i.e. RIO",
                        default='TEST')
    parser.add_argument("--AOI_Num",
                        help='AOI Number i.e 3',
                        type=int,
                        default='0')
    parser.add_argument("--mosaicGTIFF",
                        help='By default, all mosaic actions are done using virtual raster tiles.  This is a low memory'
                             'and is generally faster.  However no physical geotiff mosaic is made.   '
                             'Enable this flag to create a GTIFF of the mosaiced area before tiling',
                        action="store_false")
    parser.add_argument("--createPix",
                        help='Use imageSize in Pixels as opposed to Meters',
                        action="store_true")
    parser.add_argument("--createSummaryCSV",
                        help='Create Summary CSV used for grading of SpaceNet Challenge V2',
                        action="store_true")
    parser.add_argument("--csvLabel",
                        help='Type of csv to be created (i.e. All, Train, Test, Validate',
                        type=str,
                        default='All')
    parser.add_argument("--featureName",
                        help='Type of feature to be summarized by csv (i.e. Building)',
                        type=str,
                        default='Buildings')

    # geoJSON AOI boundary
    args = parser.parse_args()



    AOI_Name = 'RIO'
    AOI_Num = 3
    # outputDirectory Base Location


    AOI_Name = args.AOI_Name
    AOI_Num = args.AOI_Num
    srcOutline = args.srcOutline
    clipImageryToAOI=True
    if srcOutline=='':
        clipImageryToAOI=False

    outputDirectory = args.outputDirectory

    if outputDirectory=='':
    # outputDirectory Base Location
        outputDirectory = os.path.join(os.getcwd(), 'processed')

    if not os.path.isdir(outputDirectory):
        os.makedirs(outputDirectory)
    # location of imagery

    srcImageryList = []
    with open(args.srcRasterList, 'rb') as csvfile:
        print(args.srcRasterList)
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            print(row)
            if row:
                if not row[0].startswith("#"):
                    srcImageryList.append([x.strip() for x in row])

    srcVectorFileList = []
    with open(args.geoJsonList, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            if row:
                if not row[0].startswith("#"):
                    srcVectorFileList.append([x.strip() for x in row])

    lT.createAOIName(AOI_Name, AOI_Num,
                     srcImageryList,
                     srcOutline,
                     srcVectorFileList,
                     outputDirectory=outputDirectory,
                     clipImageryToAOI=clipImageryToAOI,
                     windowSizeMeters=args.imgSizeM,
                     vrtMosaic=args.mosaicGTIFF,
                     createPix=args.createPix,
                     createSummaryCSVChallenge=args.createSummaryCSV,
                     csvLabel=args.csvLabel,
                     featurName=args.featureName
                     )
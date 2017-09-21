from spacenetutilities.labeltools import coreLabelTools
import json
import glob
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Removes Id field from geojson which causes errors with certain tools')
    parser.add_argument("geojsonFolder", help="Folder to process geojsons")
    parser.add_argument("--geojsonIDModifierer", help="modifier to extension of geojson",
                        default='noid')

    args = parser.parse_args()
    coreLabelTools.removeIdinGeoJSONFolder(args.geojsonFolder,
                                           args.geojsonIDModifierer)




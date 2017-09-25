import osmnx as ox
from collections import OrderedDict
from itertools import groupby
from dateutil import parser as date_parser
from shapely.geometry import Point
from shapely.geometry import LineString, LinearRing
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.ops import unary_union
import geopandas as gpd
import pandas as pd

def createGDF_from_polygon(polygon, infrastructureList=['way["highway"]']):


    for infrastructure in infrastructureList:
        yield create_gdf_from_responseJSON(
            ox.core.osm_net_download(polygon=polygon, network_type='all_private',
                    timeout=180, memory=None,
                    max_query_area_size=50*1000*50*1000,
                    infrastructure=infrastructure)
        )


def createGDF_from_polygon_power(polygon, infrastructureList=['node["power"]',
                                                          'way["power"]']):

    edgeList = []
    nodeList = []
    for infrastructure in infrastructureList:

        tmpEdge, tmpNode = create_gdf_from_responseJSON(
            ox.core.osm_net_download(polygon=polygon, network_type='all_private',
                    timeout=180, memory=None,
                    max_query_area_size=50*1000*50*1000,
                    infrastructure=infrastructure)
        )

        edgeList.append(tmpEdge)
        nodeList.append(tmpNode)


    gdf_edges = pd.concat(edgeList)
    gdf_nodes = pd.concat(nodeList)

    return gdf_edges, gdf_nodes





def create_gdf_from_responseJSON(responseJsons):

    nodesList = []
    edgesList = []

    for element in responseJsons[0]['elements']:
        tmpElement = element
        if 'tags' in tmpElement:
            tmpElement.update(element['tags'])
        if tmpElement['type'] == 'node':
            tmpElement.update({'geometry': Point(tmpElement['lon'], tmpElement['lat'])})
            nodesList.append(tmpElement)
        elif tmpElement['type'] == 'way':
            edgesList.append(tmpElement)

    gdf_nodes = gpd.GeoDataFrame(nodesList)
    edgesList2 = []
    for edge in edgesList:
        gdflist = gdf_nodes['id'].map(lambda x: x in edge['nodes'])
        print(len(edge))
        print(edge)
        if edge['nodes'][0] == edge['nodes'][-1]:
            edge.update({'geometry': Polygon(LineString(list(gdf_nodes[gdflist].geometry.values)))
                         })
        else:
            edge.update({'geometry': LineString(list(gdf_nodes[gdflist].geometry.values))
                     })
        edgesList2.append(edge)
    gdf_edges = gpd.GeoDataFrame(edgesList2)

    return gdf_edges, gdf_nodes


if __name__ == '__main__':
    from spacenetutilities import datasets
    from spacenetutilities.osmtools import coreosmtools
    from shapely.geometry import box

    gdfOutline = datasets.load_tindex_dataset_to_gdf('AOI_2_Vegas')
    gdf_edges, gdf_nodes = coreosmtools.createGDF_from_polygon_power(box(*gdfOutline.unary_union.bounds))


import osmnx as ox


def createGDF_from_polygon(polygon, infrastructureList=['way["highway"]']):


    for infrastructure in infrastructureList:
        yield ox.core.graph_from_polygon(polygon, network_type='all_private', simplify=True,
                       retain_all=False, truncate_by_edge=False, name='unnamed',
                       timeout=180, memory=None,
                       max_query_area_size=50*1000*50*1000,
                       clean_periphery=True, infrastructure=infrastructure)


def createGDF_from_polygon_power(polygon, infrastructureList=['node["power"]',
                                                          'way["power"]']):


    for infrastructure in infrastructureList:

        yield ox.core.graph_from_polygon(polygon, network_type='all_private', simplify=True,
                       retain_all=False, truncate_by_edge=False, name='unnamed',
                       timeout=180, memory=None,
                       max_query_area_size=50*1000*50*1000,
                       clean_periphery=True, infrastructure=infrastructure)



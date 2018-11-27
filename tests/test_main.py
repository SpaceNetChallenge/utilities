"""tests ml_export.tile_generator.base"""
import os


def test_import_geoTools():
    """ Tests import of geoTools and labeltools"""

    from spacenetutilities import geoTools

    return 0

def test_import_evalTools():

    from spacenetutilities import evalTools

    return 0

def test_import_dataTools():
    from spacenetutilities import dataTools

    return 0

def test_import_osmttools():

    from spacenetutilities.osmtools import coreosmtools

    return 0

def test_import_labeltools():

    from spacenetutilities.labeltools import coreLabelTools
    from spacenetutilities.labeltools import darkNetLabel
    from spacenetutilities.labeltools import geojsonPrepTools
    from spacenetutilities.labeltools import pascalVOCLabel
    from spacenetutilities.labeltools import sbdLabel
    from spacenetutilities.labeltools import tfRecordLabel

    return 0

def test_import_inferencetools():

    from spacenetutilities.inferenceTools import coreInferenceTools


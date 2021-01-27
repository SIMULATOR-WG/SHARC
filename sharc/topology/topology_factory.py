# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:41:25 2017

@author: edgar
@modified: Luciano Camilo Tue Nov 17 09:26:25 2020
"""
import sys

from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.topology.topology_hotspot import TopologyHotspot
from sharc.topology.topology_indoor import TopologyIndoor
from sharc.topology.topology_single_base_station import TopologySingleBaseStation
from sharc.topology.topology_hibs import TopologyHIBS
from sharc.parameters.parameters import Parameters


class TopologyFactory(object):

    @staticmethod
    def createtopology(parameters: Parameters) -> Topology:
        if parameters.imt.topology == "SINGLE_BS":
            return TopologySingleBaseStation(parameters.imt.intersite_distance*2/3, parameters.imt.num_clusters)
        elif parameters.imt.topology == "MACROCELL":
            return TopologyMacrocell(parameters.imt.intersite_distance, parameters.imt.num_clusters)
        elif parameters.imt.topology == "HOTSPOT":
            return TopologyHotspot(parameters.hotspot, parameters.imt.intersite_distance, parameters.imt.num_clusters)
        elif parameters.imt.topology == "INDOOR":
            return TopologyIndoor(parameters.indoor)
        elif parameters.imt.topology == "HIBS":
            return TopologyHIBS(parameters.hibs.intersite_distance, parameters.hibs.cell_radius,
                                parameters.hibs.num_clusters, parameters.hibs.num_sectors, parameters.hibs.bs_height,
                                parameters.hibs.azimuth3, parameters.hibs.azimuth7, parameters.hibs.azimuth19,
                                parameters.hibs.elevation3, parameters.hibs.elevation7, parameters.hibs.elevation19)
        else:
            sys.stderr.write("ERROR\nInvalid topology: " + parameters.imt.topology)
            sys.exit(1)

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:41:25 2017

@author: edgar
"""
import sys

from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.topology.topology_hotspot import TopologyHotspot
from sharc.topology.topology_indoor import TopologyIndoor
from sharc.topology.topology_ntn import TopologyNTN
from sharc.topology.topology_single_base_station import TopologySingleBaseStation
from sharc.parameters.parameters import Parameters


class TopologyFactory(object):

    @staticmethod
    def createTopology(parameters: Parameters) -> Topology:
        if parameters.imt.topology.type == "SINGLE_BS":
            return TopologySingleBaseStation(
                parameters.imt.topology.single_bs.cell_radius,
                parameters.imt.topology.single_bs.num_clusters
            )
        elif parameters.imt.topology.type == "MACROCELL":
            return TopologyMacrocell(
                parameters.imt.topology.macrocell.intersite_distance,
                parameters.imt.topology.macrocell.num_clusters
            )
        elif parameters.imt.topology.type == "HOTSPOT":
            return TopologyHotspot(
                parameters.imt.topology.hotspot,
                parameters.imt.topology.hotspot.intersite_distance,
                parameters.imt.topology.hotspot.num_clusters
            )
        elif parameters.imt.topology.type == "INDOOR":
            return TopologyIndoor(parameters.imt.topology.indoor)
        elif parameters.imt.topology.type == "NTN":
            return TopologyNTN(
                parameters.imt.topology.ntn.intersite_distance,
                parameters.imt.topology.ntn.cell_radius,
                parameters.imt.topology.ntn.bs_height,
                parameters.imt.topology.ntn.bs_azimuth,
                parameters.imt.topology.ntn.bs_elevation,
                parameters.imt.topology.ntn.num_sectors,
            )
        else:
            sys.stderr.write(
                "ERROR\nInvalid topology: " +
                parameters.imt.topology,
            )
            sys.exit(1)

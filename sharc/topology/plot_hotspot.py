# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:22:58 2017

@author: edgar
"""
import math
import matplotlib.pyplot as plt
import numpy as np

from sharc.parameters.parameters_hotspot import ParametersHotspot
from sharc.topology.topology_hotspot import TopologyHotspot


def plot_topology(topology: TopologyHotspot):
    psi = np.radians([60, 120, 240, 300])

    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
    ax = fig.gca()
    
    r = topology.topology_macrocell.intersite_distance/3
    for x, y in zip(topology.topology_macrocell.x, topology.topology_macrocell.y):
        se = list([[x,y]])
        se.extend([[se[-1][0] + r, se[-1][1]]])
        se.extend([[se[-1][0] + r*math.cos(psi[0]), se[-1][1] + r*math.sin(psi[0])]])
        se.extend([[se[-1][0] + r*math.cos(psi[1]), se[-1][1] + r*math.sin(psi[1])]])
        se.extend([[se[-1][0] - r, se[-1][1]]])
        se.extend([[se[-1][0] + r*math.cos(psi[2]), se[-1][1] + r*math.sin(psi[2])]])
        sector = plt.Polygon(se, fill=None, edgecolor='k')
        ax.add_patch(sector)

        se = list([[x,y]])
        se.extend([[se[-1][0] + r*math.cos(psi[1]), se[-1][1] + r*math.sin(psi[1])]])
        se.extend([[se[-1][0] - r, se[-1][1]]])
        se.extend([[se[-1][0] + r*math.cos(psi[2]), se[-1][1] + r*math.sin(psi[2])]])
        se.extend([[se[-1][0] + r*math.cos(psi[3]), se[-1][1] + r*math.sin(psi[3])]])
        se.extend([[se[-1][0] + r, se[-1][1]]])
        sector = plt.Polygon(se, fill=None, edgecolor='k')
        ax.add_patch(sector)
        
        se = list([[x,y]])
        se.extend([[se[-1][0] + r, se[-1][1]]])
        se.extend([[se[-1][0] + r*math.cos(psi[3]), se[-1][1] + r*math.sin(psi[3])]])
        se.extend([[se[-1][0] + r*math.cos(psi[2]), se[-1][1] + r*math.sin(psi[2])]])
        se.extend([[se[-1][0] - r, se[-1][1]]])
        se.extend([[se[-1][0] + r*math.cos(psi[1]), se[-1][1] + r*math.sin(psi[1])]])
        sector = plt.Polygon(se, fill=None, edgecolor='k')
        ax.add_patch(sector)       
        

    # plot hotspot centers
    #plt.scatter(topology.hotspot_x, topology.hotspot_y, color='k', edgecolor="w", linewidth=0.5)
    
    # plot small cells
    plt.scatter(topology.x, topology.y, color='r', edgecolor="w", linewidth=0.5, label="Small cell")
    
    # plot hotspots coverage area
    for hx, hy in zip(topology.hotspot_x, topology.hotspot_y):
        circ = plt.Circle((hx, hy), radius=50, color='g', fill=False, linewidth=0.5)
        ax.add_patch(circ)
    
    # macro cell base stations
    plt.scatter(topology.topology_macrocell.x, topology.topology_macrocell.y, color='k', edgecolor="k", linewidth=4, label="Macro cell")

    # sector centers
    #plt.scatter(-sector_y, sector_x, color='g', edgecolor="g")
    
    # plot macro cell coverage area
    #ax = fig.gca()
#    for mx, my in zip(topology.topology_macrocell.x, topology.topology_macrocell.y):
#        circ = plt.Circle((mx, my), radius=666.667*math.sqrt(3)/2-70, color='b', fill=False, linewidth=0.5)
#        ax.add_patch(circ)    
    
    plt.axis('image') 
    plt.title("Macro cell topology with hotspots")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    #plt.xlim((-3000, 3000))
    #plt.ylim((-3000, 3000))                
    plt.legend(loc="upper left", scatterpoints=1)
    plt.tight_layout()    
    plt.show()


if __name__ == '__main__':
    param = ParametersHotspot()
    intersite_distance = 600
    num_clusters = 1
    topology = TopologyHotspot(param, intersite_distance, num_clusters)
    topology.calculate_coordinates()
    
    plot_topology(topology)
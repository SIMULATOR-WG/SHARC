# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:52:28 2017

@author: Calil
"""

from area import area
from numpy import cos, sin, tan, arctan, deg2rad, rad2deg, arccos, pi, linspace, arcsin, isnan, vstack
import matplotlib.pyplot as plt

class Footprint(object):
    """
    Defines a satellite footprint region and calculates its area.
    
    Construction:
        FootprintArea(bore_lat_deg, bore_subsat_long_deg, beam)
            bore_lat_deg (float): latitude of boresight point
            bore_subsat_long_deg (float): longitude of boresight with respect
                to sub-satellite point, taken positive when to the west of the
                sub-satellite point
            beam_deg (float): half of beam width in degrees
    """
    def __init__(self,bore_lat_deg: float, bore_subsat_long_deg: float, beam_deg:float):
        # Initialize attributes
        self.bore_lat_deg = bore_lat_deg
        self.bore_subsat_long_deg = bore_subsat_long_deg
        self.beam_width_deg = beam_deg
        
        # Convert to radians
        self.bore_lat_rad = deg2rad(bore_lat_deg)
        self.bore_subsat_long_rad = deg2rad(bore_subsat_long_deg)
        self.beam_width_rad = deg2rad(beam_deg)
        
        # Calculate tilt
        self.beta = arccos(cos(self.bore_lat_rad)*\
                           cos(self.bore_subsat_long_rad))
        self.bore_tilt = arctan(sin(self.beta)/(6.6235 - cos(self.beta)))
        
    def calc_footprint(self, n: int):
        """
        Defines footprint polygonal approximation
        
        Input:
            n (int): number of vertices on polygonal
            
        Outputs:
            pt_long (np.array): longitude of verices in deg
            pt_lat (np.array): latiture of verices in deg
        """
        # Projection angles
        phi = linspace(0,2*pi,num = n)
        
        cos_gamma_n = cos(self.bore_tilt)*cos(self.beam_width_rad) + \
                      sin(self.bore_tilt)*sin(self.beam_width_rad)*\
                      cos(phi)
        cot_phi_n = (sin(self.bore_tilt)*self.cot(self.beam_width_rad) - \
                     cos(self.bore_tilt)*cos(phi))/sin(phi)
        
        gamma_n = arccos(cos_gamma_n)
        phi_n = self.arccot(cot_phi_n)
        
        eps_n = arctan(sin(self.bore_subsat_long_rad)/tan(self.bore_lat_rad)) + \
                phi_n
        eps_n[isnan(eps_n)] = arctan(1) + phi_n[isnan(eps_n)]
                
        beta_n = arcsin(6.6235*sin(gamma_n)) - gamma_n
        
        pt_lat  = arcsin(sin(beta_n)*cos(eps_n))
        pt_long = arctan(tan(beta_n)*sin(eps_n))
        
        return rad2deg(pt_long), rad2deg(pt_lat)
    
    def calc_area(self, n:int):
        """
        Returns footprint area in km^2
        
        Input:
            n (int): number of vertices on polygonal approximation
        Output:
            a (float): footprint area
        """
        long, lat = self.calc_footprint(n)
        
        long_lat = vstack((long, lat)).T
        
        obj = {'type':'Polygon',
               'coordinates':[long_lat.tolist()]}
        
        return area(obj)*1e-6
        
    def cot(self,angle):
        return tan(pi/2 - angle)
    
    def arccot(self,x):
        return pi/2 - arctan(x)
        
if __name__ == '__main__':
    # Earth  [km]
    R = 6371
    
    # Create object
    fprint = Footprint(0,0,0.325)
    
    # Define coordinates
    long, lat = fprint.calc_footprint(1000)
    plt.plot(long,lat)
    plt.show()
        
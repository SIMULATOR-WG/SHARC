# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:52:28 2017

@author: Calil
"""

from area import area
from numpy import cos, sin, tan, arctan, deg2rad, rad2deg, arccos, pi, linspace, arcsin, vstack, arctan2, where
import matplotlib.pyplot as plt

class Footprint(object):
    """
    Defines a satellite footprint region and calculates its area.
    Method for generating footprints (Siocos,1973) is found in the book 
    "Satellite Communication Systems" by M. Richharia ISBN 0-07-134208-7
    
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
        self.bore_tilt = arctan2(sin(self.beta),(6.6235 - cos(self.beta)))
        
        # Maximum tilt and latitute coverage
        self.max_gamma_rad = deg2rad(8.6833)
        self.max_beta_rad = deg2rad(81.3164)
        
    def calc_footprint(self, n: int):
        """
        Defines footprint polygonal approximation
        
        Input:
            n (int): number of vertices on polygonal
            
        Outputs:
            pt_long (np.array): longitude of vertices in deg
            pt_lat (np.array): latiture of vertices in deg
        """
        # Projection angles
        phi = linspace(0,2*pi,num = n)
        
        cos_gamma_n = cos(self.bore_tilt)*cos(self.beam_width_rad) + \
                      sin(self.bore_tilt)*sin(self.beam_width_rad)*\
                      cos(phi)
#        tan_phi_n = sin(phi)/(sin(self.bore_tilt)*self.cot(self.beam_width_rad) - \
#                     cos(self.bore_tilt)*cos(phi))
        
        gamma_n = arccos(cos_gamma_n)
#        phi_n = arctan(tan_phi_n) 
        phi_n = arctan2(sin(phi),(sin(self.bore_tilt)*self.cot(self.beam_width_rad) - \
                     cos(self.bore_tilt)*cos(phi))) 
        
        eps_n = arctan2(sin(self.bore_subsat_long_rad),tan(self.bore_lat_rad)) + \
                phi_n
                
        beta_n = arcsin(6.6235*sin(gamma_n)) - gamma_n
        beta_n[where(gamma_n >  self.max_gamma_rad)] = self.max_beta_rad
        
        pt_lat  = arcsin(sin(beta_n)*cos(eps_n))
        pt_long = arctan(tan(beta_n)*sin(eps_n))
        
        return rad2deg(pt_long), rad2deg(pt_lat)
    
    def calc_area(self, n:int):
        """
        Returns footprint area in km^2
        
        Input:
            n (int): number of vertices on polygonal approximation
        Output:
            a (float): footprint area in km^2
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
    fprint = Footprint(45,61,0.325)
    
    # Define coordinates
    long, lat = fprint.calc_footprint(100)
    col = linspace(-10,10,num=len(long))
#    plt.scatter(long,lat,c=col,cmap="inferno")
    plt.plot(long,lat)
#    plt.xlim([-2, +2])
#    plt.ylim([-2, +2])
    plt.show()
        
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
            beam_deg (float): half of beam width in degrees
            elevation_deg (float): optional. Satellite elevation at 
            boresight bore_lat_deg (float): optional, default = 0. 
                Latitude of boresight point. If elevation is given this 
                parameter is not used. Default = 0
            bore_subsat_long_deg (float): longitude of boresight with respect
                to sub-satellite point, taken positive when to the west of the
                sub-satellite point. If elevation is given this 
                parameter is not used. Default = 0
    """
    def __init__(self,beam_deg:float,**kwargs):
        # Initialize attributes
        if 'elevation_deg' in kwargs.keys():
            self.elevation_deg = kwargs['elevation_deg']
            self.bore_lat_deg = 0.0
            self.bore_subsat_long_deg = self.calc_beta(self.elevation_deg)
        else:
            self.bore_lat_deg = 0.0
            self.bore_subsat_long_deg = 0.0
            if 'bore_lat_deg' in kwargs.keys():
                self.bore_lat_deg = kwargs['bore_lat_deg']
            if 'bore_subsat_long_deg' in kwargs.keys():
                self.bore_subsat_long_deg = kwargs['bore_subsat_long_deg']
            self.elevation_deg = \
                self.calc_elevation(self.bore_lat_deg,self.bore_subsat_long_deg)
        
        self.beam_width_deg = beam_deg
        
        # Convert to radians
        self.elevation_rad = deg2rad(self.elevation_deg)
        self.bore_lat_rad = deg2rad(self.bore_lat_deg)
        self.bore_subsat_long_rad = deg2rad(self.bore_subsat_long_deg)
        self.beam_width_rad = deg2rad(self.beam_width_deg)
        
        # Calculate tilt
        self.beta = arccos(cos(self.bore_lat_rad)*\
                           cos(self.bore_subsat_long_rad))
        self.bore_tilt = arctan2(sin(self.beta),(6.6235 - cos(self.beta)))
        
        # Maximum tilt and latitute coverage
        self.max_gamma_rad = deg2rad(8.6833)
        self.max_beta_rad = deg2rad(81.3164)
        
    def calc_beta(self,elev_deg: float):
        """
        Calculates elevation angle based on given elevation. Beta is the 
        subpoint to earth station great-circle distance
        
        Input:
            elev_deg (float): elevation in degrees
            
        Output:
            beta (float): beta angle in degrees  
        """
        elev_rad = deg2rad(elev_deg)
        beta = 90 - elev_deg - rad2deg(arcsin(cos(elev_rad)/6.6235))
        return beta
    
    def calc_elevation(self,lat_deg: float, long_deg: float):
        """
        Calculates elevation for given latitude of boresight point and 
        longitude of boresight with respect to sub-satellite point.
        
        Inputs:
            lat_deg (float): latitude of boresight point in degrees
            long_deg (float): longitude of boresight with respect
                to sub-satellite point, taken positive when to the west of the
                sub-satellite point, in degrees
        
        Output:
            elev (float): elevation in degrees
        """
        lat_rad = deg2rad(lat_deg)
        long_rad = deg2rad(long_deg)
        beta = arccos(cos(lat_rad)*cos(long_rad))
        elev = arctan2((cos(beta) - 0.1510),sin(beta))
        
        return rad2deg(elev)
    
    def set_elevation(self,elev: float):
        """
        Resets elevation angle to given value
        """
        self.elevation_deg = elev
        self.bore_lat_deg = 0.0
        self.bore_subsat_long_deg = self.calc_beta(self.elevation_deg)
        
        # Convert to radians
        self.elevation_rad = deg2rad(self.elevation_deg)
        self.bore_lat_rad = deg2rad(self.bore_lat_deg)
        self.bore_subsat_long_rad = deg2rad(self.bore_subsat_long_deg)
        
        # Calculate tilt
        self.beta = arccos(cos(self.bore_lat_rad)*\
                           cos(self.bore_subsat_long_rad))
        self.bore_tilt = arctan2(sin(self.beta),(6.6235 - cos(self.beta)))
        
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
        
        gamma_n = arccos(cos_gamma_n) 
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
    fprint90 = Footprint(0.325,elevation_deg=90)
    fprint45 = Footprint(0.325,elevation_deg=45)
    fprint30 = Footprint(0.325,elevation_deg=30)
    fprint20 = Footprint(0.325,elevation_deg=20)
    
    # Plot coordinates
    plt.figure(figsize=(15,2))
    long, lat = fprint90.calc_footprint(100)
    plt.plot(long,lat,'k',label='$90^o$')
    long, lat = fprint45.calc_footprint(100)
    plt.plot(long,lat,'b',label='$45^o$')
    long, lat = fprint30.calc_footprint(100)
    plt.plot(long,lat,'r',label='$30^o$')
    long, lat = fprint20.calc_footprint(100)
    plt.plot(long,lat,'g',label='$20^o$')
    plt.legend(loc='upper right')
    plt.xlabel('Longitude [deg]')
    plt.ylabel('Latitude [deg]')
    plt.xlim([-5, 80])
    plt.grid()
    plt.show()
    
    # Print areas
    n = 1000
    print("Sat elevation 90 deg: area = {}".format(fprint90.calc_area(n)))
    print("Sat elevation 45 deg: area = {}".format(fprint45.calc_area(n)))
    print("Sat elevation 30 deg: area = {}".format(fprint30.calc_area(n)))
    print("Sat elevation 20 deg: area = {}".format(fprint20.calc_area(n)))
    
    
        
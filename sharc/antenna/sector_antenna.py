# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 14:50:27 2017

@author: Calil
"""

import numpy as np

from sharc.antenna.antenna import Antenna

class SectorAntenna(Antenna):
    
    def __init__(self,azimuth,position):
        """
        Angles in DEGREES.
        """
        self.azimuth = azimuth
        self.position = position
        
    def antenna_gain(self,ue_pos):
        
        hBs = self.position[2]
        point_vec = ue_pos - self.position
        
        distUe = np.linalg.norm(point_vec)
        
        # antenna model
        tilt_deg = 4 # down tilt in degrees
        gain_dBi = 18
        HPBW_deg = 65 # half-power beamwidth in degrees
        HPBWVert_deg = 6 # vertical half-power beamwidth in degrees
        FBR_dB = 30 # front-to-back ratio in dB
        SLL_dB = 18 # side-lobe level in -dB
        
        # add antennaGain, model from
        # Gunnarson, "Downtilted Base Station Antennas – A Simulation Model Proposal and
        # Impact on HSPA and LTE Performance, VTC-2008
        horAngle = self.azimuth - np.arctan2(point_vec[1],point_vec[0])
        
        horAngle_deg = horAngle * 180 / np.pi
        horAngle_deg  = ( horAngle_deg + 180 ) % 360 - 180 # angle between -180 and 180 degrees
        horGain = -np.minimum(12 * (horAngle_deg / HPBW_deg) ** 2, FBR_dB) + gain_dBi
        vertGain = np.maximum(-12 * ((np.arctan(hBs / distUe) * 180 / np.pi - tilt_deg) / HPBWVert_deg) ** 2, -SLL_dB)
        
        self.gain = horGain + vertGain
        
        return self.gain
    
    def get_antenna_obj(self,ue_pos):
        return Antenna(self.antenna_gain(ue_pos))
        
        
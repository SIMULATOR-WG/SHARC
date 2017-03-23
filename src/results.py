# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:47:46 2017

@author: edgar
"""


class Results(object):
    
    def __init__(self):
        self.__coupling_loss_dl = list()
        self.__tx_power_dl = list()
        self.__sinr_dl = list()
        self.__snr_dl = list()
        self.__throughoput_dl = list()
        
    def add_coupling_loss_dl(self, sample):
        self.coupling_loss_dl.extend(sample)
        
    def add_tx_power_dl(self, sample):
        self.tx_power_dl.extend(sample)
        
    def add_sinr_dl(self, sample):
        self.sinr_dl.extend(sample)
        
    def add_snr_dl(self, sample):
        self.snr_dl.extend(sample)        
        
    @property
    def coupling_loss_dl(self):
        return self.__coupling_loss_dl
        
    @property
    def tx_power_dl(self):
        return self.__tx_power_dl

    @property
    def sinr_dl(self):
        return self.__sinr_dl

    @property
    def snr_dl(self):
        return self.__snr_dl

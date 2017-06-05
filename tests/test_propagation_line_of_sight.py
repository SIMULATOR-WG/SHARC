# -*- coding: utf-8 -*-


import unittest
import numpy as np
import matplotlib.pyplot as plt
import numpy.testing as npt

from sharc.propagation.P452.propagation_line_of_sight import PropagationLineOfSight

class PropagationLineOfSightTest(unittest.TestCase):
    
    def setUp(self):
        self.__LineOfSight = PropagationLineOfSight()
        
    def test_loss(self):
        d = 10000
        f = 27000
        Ph = 1013 
        T = 288
        ro = 7.5
        
        npt.assert_allclose(142.226, 
                         self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)
        
        d = 20000
        f = 40000
        
        npt.assert_allclose(153.152, 
                         self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)
        
        
        #Grafico da perda de difraçao em funçao da distancia e da frequencia
        data1 = []
        data2 = []
        data3 = []
        data4 = []
        data5 = []
        eixo_x = [] 
        
        for n in range (1,6,1):
    
            if (n==1):
                f = 5000
                    
                for d in range(1, 100100,1000):
           
                    loss = self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
        
                    data1.append(loss)
                    eixo_x.append(d/1000)
                    
            if (n==2):
                f = 10000
                    
                for d in range(1, 100100,1000):
                     loss = self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
                     data2.append(loss)

            if (n==3):
                f = 20000
                    
                for d in range(1, 100100,1000):
                     loss = self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
                     data3.append(loss)

            if (n==4):
                f = 30000
                    
                for d in range(1, 100100,1000):
                     loss = self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
                     data4.append(loss)


            if (n==5):
                f = 40000
                    
                for d in range(1, 100100,1000):
                     loss = self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
                     data5.append(loss)
                 
        fig = plt.figure(5)
        f  = ['5 GHz','10 GHz','20 GHz','30 GHz','40 GHz']
        ax = fig.add_subplot(111)
        ax.plot(eixo_x, data1)
        ax.plot(eixo_x, data2)
        ax.plot(eixo_x, data3)
        ax.plot(eixo_x, data4)
        ax.plot(eixo_x, data5)
    #    ax.plot(eixo_x, data6)  
        
        # Add legend, title and axis labels
        lgd = ax.legend( [ 'f = ' + str(lag) for lag in f], loc='upper center', bbox_to_anchor=(0.8, 0.45))
        ax.set_title('Line Of Sight attenuation')
        ax.set_xlabel('Distance (Km)')
        ax.set_ylabel('Attenuation (dB)')
        ax.set_xlim([0,100])
        ax.set_ylim([100,180])
        ax.grid(True)
        fig.savefig('line_of_sight.png', dpi=350, format='png') 
                          
    
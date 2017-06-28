# -*- coding: utf-8 -*-
"""
Created on Mon Mar  13 15:14:34 2017

@author: edgar
"""

import unittest
import numpy as np
import matplotlib.pyplot as plt

from sharc.propagation.propagation_free_space import PropagationFreeSpace

class PropagationFreeSpaceTest(unittest.TestCase):
    
    def setUp(self):
        self.__freeSpace = PropagationFreeSpace()
        
    def test_loss(self):
        d = 10
        f = 10
        self.assertEqual(12.45, 
                         self.__freeSpace.get_loss(distance=d, frequency=f))

        d = [ 10, 100 ]
        f = [ 10, 100 ]
        self.assertTrue(np.all(np.equal([12.45, 52.45], 
                         self.__freeSpace.get_loss(distance=d, frequency=f))))

        d = [ 10, 100, 1000 ]
        f = [ 10, 100, 1000 ]
        self.assertTrue(np.all(np.equal([12.45, 52.45, 92.45], 
                         self.__freeSpace.get_loss(distance=d, frequency=f))))

        d = [[10, 20, 30],[40, 50, 60]]
        f = [ 100 ]
        ref_loss = [[ 32.45,  38.47,  41.99],
                    [ 44.49,  46.42,  48.01]]
        loss = self.__freeSpace.get_loss(distance=d, frequency=f)
        self.assertTrue(np.all(np.isclose(ref_loss, loss, atol=1e-2)))         
        
        
        #Grafico da perda de difraçao em funçao da distancia (km) e da frequencia (GHz)
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
           
                    loss = self.__freeSpace.get_loss(distance=d, frequency=f)
                    data1.append(loss)
                    eixo_x.append(d/1000)
                    
            if (n==2):
                f = 10000
                    
                for d in range(1, 100100,1000):
                     loss = self.__freeSpace.get_loss(distance=d, frequency=f)
                     data2.append(loss)

            if (n==3):
                f = 20000
                    
                for d in range(1, 100100,1000):
                     loss = self.__freeSpace.get_loss(distance=d, frequency=f)
                     data3.append(loss)

            if (n==4):
                f = 30000
                    
                for d in range(1, 100100,1000):
                     loss = self.__freeSpace.get_loss(distance=d, frequency=f)
                     data4.append(loss)


            if (n==5):
                f = 40000
                    
                for d in range(1, 100100,1000):
                     loss = self.__freeSpace.get_loss(distance=d, frequency=f)
                     data5.append(loss)
                 
        fig = plt.figure(6)
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
        ax.set_title('Free space attenuation')
        ax.set_xlabel('Distance (Km)')
        ax.set_ylabel('Attenuation (dB)')
        ax.set_xlim([0,100])
        ax.set_ylim([100,170])
        ax.grid(True)
        fig.savefig('free_space.png', dpi=350, format='png') 
               
  
        
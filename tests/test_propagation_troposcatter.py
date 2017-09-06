
import unittest

import matplotlib.pyplot as plt
import numpy.testing as npt

from sharc.propagation.P452.propagation_troposcatter import PropagationTropScatter

class PropagationTropScatterTest(unittest.TestCase):
    
    def setUp(self):
        self.__gasAtt = PropagationTropScatter()
        
    def test_loss(self):
       

        d = 300000       #distance in m
        f = 27000       #Frequency in GHz
        p = 40          #percentage p
        Gt = 30         #Transmition antena gain
        Gr = 30         #Reception antena gain
        thetaT = 10     #Transmit horizon elevation (mrad)
        thetaR = 10     #Receive horizon elevation (mrad)
        N0 = 355        #Sea-level surface refractivity (use the map)
        deltaN = 60     #Average radio-refractive (use the map)
        Ph = 1013 
        T = 288
        ro = 3
       
        npt.assert_allclose(262.273, 
                         self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p),atol=1e-3)

        d = 600000       #distance in m
        f = 40000       #Frequency in GHz
        npt.assert_allclose(320.272, 
                         self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p),atol=1e-3)

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
           
                    loss = self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p)
                      
                    data1.append(loss)
                    eixo_x.append(d/1000)
                    
            if (n==2):
                f = 10000
                    
                for d in range(1, 100100,1000):
                     loss = self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p)
                     data2.append(loss)

            if (n==3):
                f = 20000
                    
                for d in range(1, 100100,1000):
                     loss = self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p)
                     data3.append(loss)

            if (n==4):
                f = 30000
                    
                for d in range(1, 100100,1000):
                     loss = self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p)
                     data4.append(loss)


            if (n==5):
                f = 40000
                    
                for d in range(1, 100100,1000):
                     loss = self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p)
                     data5.append(loss)
                 
        fig = plt.figure(4)
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
        ax.set_title('Tropospheric scatter attenuation')
        ax.set_xlabel('Distance (Km)')
        ax.set_ylabel('Attenuation (dB)')
        ax.set_xlim([0,100])
        ax.set_ylim([160,240])
        ax.grid(True)
        fig.savefig('tropospheric_scatter.png', dpi=350, format='png') 
                          
    

import unittest

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

if __name__ == '__main__':
    unittest.main()
       

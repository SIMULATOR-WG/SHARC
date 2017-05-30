
import unittest

import numpy.testing as npt

from sharc.propagation.P452.propagation_duting_reflection import PropagationDutingReflection

class PropagationDutingReflectionTest(unittest.TestCase):
    
    def setUp(self):
        self.__dutingRef = PropagationDutingReflection()
        
    def test_loss(self):
       

        f = 27000          #Frequency in GHz
        p = 40          #percentage p
        d = 10000         #distance in km
        Dlt = 30        #distance from the transmit antennas to their respective horizons (km)
        Dlr = 10        #distance from the receive antennas to their respective horizons (km)
        Dct = 10        #Distance over land from the transmit and receive antennas to the coast (km)
        Dcr = 10        #Distance over land from the transmit and receive antennas to the coast (km)
        Hts = 280       #Antenna centre height above mean sea level (m)
        Hrs = 244       #Antenna centre height above mean sea level (m)
        Hte = 50        #Effective height of interfering antenna (m)
        Hre = 50        #Effective height of interfered-with antenna (m)
        thetaT = 10     #Transmit horizon elevation (mrad)
        thetaR = 10     #Receive horizon elevation (mrad)
        N0 = 355        #Sea-level surface refractivity (use the map)
        deltaN = 60     #Average radio-refractive (use the map)
        omega = 0       #Fraction of the total path over water, 0 for totally overland paths
        phi = 60        #path centre latitude (degrees).
        dtm = 0.8       #longest continuous land (inland + coastal) section of the great-circle path (km)
        dlm = 0.8       #longest continuous inland section of the great-circle path (km)                                    
        epsilon = 3.5                                     
        hm = 15             # ? nao sei o que é, pag. 18 eq (56)
        Ph = 1013 
        T = 288
        ro = 3
        
        npt.assert_allclose(307.170, 
                         self.__dutingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm),atol=1e-3)

        f = 40000
        d = 20000
        npt.assert_allclose(332.832, 
                         self.__dutingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm),atol=1e-3)

#        d = [ 10000, 2000]
#        f = [ 27000, 40000 ]
#        self.assertTrue(np.all(np.equal([307.403, 333.586], 
#                         self.__dutingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm))))

if __name__ == '__main__':
    unittest.main()
       

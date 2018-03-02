
import unittest
import matplotlib.pyplot as plt
import numpy.testing as npt
import numpy as np

from sharc.propagation.propagation_ducting_reflection import PropagationDuctingReflection
from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation

class PropagationDuctingReflectionTest(unittest.TestCase):

    def setUp(self):

        random_number_gen = np.random.RandomState()
        self.__ductingRef = PropagationDuctingReflection(random_number_gen,
                                                         PropagationGasesAttenuation(random_number_gen))

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
                         self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm),atol=1e-3)

        f = 40000
        d = 20000
        npt.assert_allclose(332.832,
                         self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm),atol=1e-3)

#        d = [ 10000, 2000]
#        f = [ 27000, 40000 ]
#        self.assertTrue(np.all(np.equal([307.403, 333.586],
#                         self.__dutingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm))))

   #Grafico da perda de difraçao em funçao da distancia e da frequencia
#        data1 = []
#        data2 = []
#        data3 = []
#        data4 = []
#        data5 = []
#        eixo_x = []
#
#
#
#        for n in range (1,6,1):
#
#            if (n==1):
#                f = 5000
#
#                for d in range(1000, 100100,1000):
#                    loss = self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)
#
#                    data1.append(loss)
#                    eixo_x.append(d/1000)
#
#            if (n==2):
#                f = 10000
#
#                for d in range(1000, 100100,1000):
#                     loss = self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)
#                     data2.append(loss)
#
#            if (n==3):
#                f = 20000
#
#                for d in range(1000, 100100,1000):
#                     loss = self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)
#                     data3.append(loss)
#
#            if (n==4):
#                f = 30000
#
#                for d in range(1000, 100100,1000):
#                     loss = self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)
#                     data4.append(loss)
#
#
#            if (n==5):
#                f = 40000
#
#                for d in range(1000, 100100,1000):
#                     loss = self.__ductingRef.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)
#                     data5.append(loss)

#        fig = plt.figure(3)
#        f  = ['5 GHz','10 GHz','20 GHz','30 GHz','40 GHz']
#        ax = fig.add_subplot(111)
#        ax.plot(eixo_x, data1)
#        ax.plot(eixo_x, data2)
#        ax.plot(eixo_x, data3)
#        ax.plot(eixo_x, data4)
#        ax.plot(eixo_x, data5)
#
#        # Add legend, title and axis labels
#        lgd = ax.legend( [ 'f = ' + str(lag) for lag in f], loc='upper center', bbox_to_anchor=(0.8, 0.45))
#        ax.set_title('Ducting/reflection attenuation')
#        ax.set_xlabel('Distance (Km)')
#        ax.set_ylabel('Attenuation (dB)')
#        ax.set_xlim([0,100])
#        ax.set_ylim([200,380])
#        ax.grid(True)
#        fig.savefig('ducting_ref_att.png', dpi=350, format='png')

if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-


import unittest
import numpy as np
import matplotlib.pyplot as plt
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy.testing as npt

from sharc.propagation.propagation_clear_air_452 import PropagationClearAir

class PropagationClearAirTest(unittest.TestCase):

    def setUp(self):
        self.__ClearAir = PropagationClearAir(np.random.RandomState())

    def test_loss(self):

        params = ParametersFssEs()

        d = 10000
        f = 27000
        params.atmospheric_pressure = 1013
        params.air_temperature = 288
        params.water_vapour = 7.5
        params.Dlt = 30
        params.Dlr = 10
        params.Dct = 10
        params.Dcr = 10

        params.Hts = 244
        params.Hrs = 280
        params.Hte = 50
        params.Hre = 50

        params.theta_tx = 20
        params.theta_rx = 20

        params.N0 = 355
        params.delta_N = 60
        params.percentage_p = 40

        params.omega = 0
        params.phi = 60
        params.dtm = .8
        params.dlm = .8

        params.epsilon = 3.5

        params.hm = 15
        params.Hsr = 45
        params.Hst = 48

        params.H0 = 15
        params.Hn = 17

        params.thetaJ = 0.3
        params.par_ep = 0.8

        params.clutter_loss = False

        Gt = 10
        Gr = 10

        di = [1,1,1]
        hi = [2,4,6]


#        Ld50, Ldbeta, Ldb = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte, Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p)
#
#        npt.assert_allclose(158.491,Ldb,atol=1e-3)

        #Grafico da perda de difraçao em funçao da distancia e da frequencia
#        data1 = []
#        data2 = []
#        data3 = []
#        data4 = []
#        data5 = []
#        eixo_x = []
#
#        for n in range (1,6,1):
#
#            params.eta = n
#
#            if params.eta==1:
#                f = 10000
#
#                for d in range(1000, 100100,1000):
#
#                    d = np.array(d, ndmin=2)
#                    Loss = self.__ClearAir.get_loss(distance_3D=d, frequency=f, es_params=params,
#                                                    tx_gain=Gt, rx_gain=Gr, di=di, hi=hi)
#
#                    data1.append(Loss)
#
#                    eixo_x.append(d/1000)
#
#            if params.eta==2:
#                f = 20000
#
#                for d in range(1000, 100100,1000):
#                    d = np.array(d, ndmin=2)
#
#                    Loss = self.__ClearAir.get_loss(distance_3D=d, frequency=f,es_params=params,
#                                                    tx_gain=Gt, rx_gain=Gr, di=di, hi=hi)
#                    data2.append(Loss)
#
#            if params.eta==3:
#                f = 30000
#
#                for d in range(1000, 100100,1000):
#
#                    d = np.array(d, ndmin=2)
#                    Loss = self.__ClearAir.get_loss(distance_3D=d, frequency=f,es_params=params,
#                                                    tx_gain=Gt, rx_gain=Gr, di=di, hi=hi)
#                    data3.append(Loss)
#
#            if params.eta==4:
#                f = 40000
#
#                for d in range(1000, 100100,1000):
#                    d = np.array(d, ndmin=2)
#                    Loss = self.__ClearAir.get_loss(distance_3D=d, frequency=f,es_params=params,
#                                                    tx_gain=Gt, rx_gain=Gr, di=di, hi=hi)
#                    data4.append(Loss)
#
#
#            if params.eta==5:
#                f = 50000
#
#                for d in range(1000, 100100,1000):
#                    d = np.array(d, ndmin=2)
#                    Loss = self.__ClearAir.get_loss(distance_3D=d, frequency=f,es_params=params,
#                                                    tx_gain=Gt, rx_gain=Gr, di=di, hi=hi)
#                    data5.append(Loss)

    #     fig = plt.figure(2)
    #     f = ['10 GHz','20 GHz','30 GHz','40 GHz','50 GHz','60 GHz']
    #     ax = fig.add_subplot(111)
    #     ax.plot(eixo_x, data1)
    #     ax.plot(eixo_x, data2)
    #     ax.plot(eixo_x, data3)
    #     ax.plot(eixo_x, data4)
    #     ax.plot(eixo_x, data5)
    # #    ax.plot(eixo_x, data6)
    #
    #     # Add legend, title and axis labels
    #     lgd = ax.legend( [ 'f = ' + str(lag) for lag in f], loc='upper center', bbox_to_anchor=(0.16, 1))
    #     ax.set_title('Overall prediction attenuation')
    #     ax.set_xlabel('Distance (Km)')
    #     ax.set_ylabel('Attenuation (dB)')
    #     ax.set_xlim([0,100])
    #     ax.set_ylim([0,60])
    #     ax.grid(True)
    #     fig.savefig('clear_air_att.png', dpi=350, format='png')


if __name__ == '__main__':
    unittest.main()

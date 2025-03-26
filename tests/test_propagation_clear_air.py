# -*- coding: utf-8 -*-


import unittest
import numpy as np
import numpy.testing as npt
from sharc.parameters.parameters_p452 import ParametersP452
from sharc.propagation.propagation_clear_air_452 import PropagationClearAir


class PropagationClearAirTest(unittest.TestCase):

    def setUp(self):
        param_p452 = ParametersP452()
        # param_p452.atmospheric_pressure = 1013
        # param_p452.air_temperature = 288
        # param_p452.Dct = 10
        # param_p452.Dcr = 10

        # param_p452.Hte = 50
        # param_p452.Hre = 50

        # param_p452.N0 = 355
        # param_p452.delta_N = 60
        # param_p452.percentage_p = 40

        # param_p452.clutter_loss = False

        self.prop_clear_air = PropagationClearAir(
            np.random.RandomState(), param_p452)

    def test_loss(self):

        # distance between stations in meters
        distances = np.ones((1, 1), dtype=float) * 1000
        frequencies = np.ones((1, 1), dtype=float) * 27000  # frequency in MHz
        indoor_stations = np.zeros((1, 1), dtype=bool)
        # elevation between stations in degrees
        elevations = np.zeros((1, 1), dtype=float)
        tx_gain = np.ones((1, 1), dtype=float) * 0
        rx_gain = np.ones((1, 1), dtype=float) * 0

        loss = self.prop_clear_air.get_loss(
            distances,
            frequencies,
            indoor_stations,
            elevations,
            tx_gain,
            rx_gain)
        # npt.assert_allclose(158.491, loss, atol=1e-3)

    #    Ld50, Ldbeta, Ldb = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph,
    # air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte, Hre=hre, Hsr=hsr, Hst=hst, H0=h0,
    # Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p)

    #    npt.assert_allclose(158.491,Ldb,atol=1e-3)

        # Grafico da perda de difraçao em funçao da distancia e da frequencia
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

# -*- coding: utf-8 -*-


import unittest
import numpy as np
import matplotlib.pyplot as plt

import numpy.testing as npt

from sharc.propagation.propagation_diffraction import PropagationDiffraction
from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation

class PropagationDiffractionTest(unittest.TestCase):

    def setUp(self):
        random_number_gen = np.random.RandomState()
        self.__Diffraction = PropagationDiffraction(random_number_gen,
                                                    PropagationGasesAttenuation(random_number_gen))

    def test_loss(self):
        d = 10000
        f = 27000
        Ph = 1013
        T = 288
        ro = 7.5
        deltaN = 60
        hrs = 280
        hts = 244
        hte = 50
        hre = 50
        hsr = 45
        hst = 48
        h0 = 15
        hn = 17
        di = [1,1,1]
        hi = [2,4,6]
        omega = 0
        dlt = 30
        dlr = 10
        p = 40
        Beta = 60

        C0 = 2.515516698
        C1 = 0.802853
        C2 = 0.010328
        D1 = 1.432788
        D2 = 0.189269
        D3 = 0.001308


        Ld50, Ldbeta,Ldp, Ldb  = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte,
                                                             Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p,C0=C0,C1=C1,C2=C2,D1=D1,D2=D2,D3=D3)

        npt.assert_allclose(158.491,Ldb,atol=1e-3)





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
#            if (n==1):
#                f = 10000
#
#                for d in range(1000, 100100,1000):
#
#                    Ld50, Ldbeta,Ldp, Ldb  = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte,
#                                                                         Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p,C0=C0,C1=C1,C2=C2,D1=D1,D2=D2,D3=D3)
#                    data1.append(Ldb)
#                    eixo_x.append(d/1000)
#
#
#            if (n==2):
#                f = 20000
#
#                for d in range(1000, 100100,1000):
#                     Ld50, Ldbeta,Ldp, Ldb  = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte,
#                                                                         Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p,C0=C0,C1=C1,C2=C2,D1=D1,D2=D2,D3=D3)
#                     data2.append(Ldb)
#
#            if (n==3):
#                f = 30000
#
#                for d in range(1000, 100100,1000):
#                     Ld50, Ldbeta,Ldp, Ldb  = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte,
#                                                                         Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p,C0=C0,C1=C1,C2=C2,D1=D1,D2=D2,D3=D3)
#                     data3.append(Ldb)
#
#            if (n==4):
#                f = 40000
#
#                for d in range(1000, 100100,1000):
#                     Ld50, Ldbeta,Ldp, Ldb  = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte,
#                                                                         Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p,C0=C0,C1=C1,C2=C2,D1=D1,D2=D2,D3=D3)
#                     data4.append(Ldb)
#
#
#            if (n==5):
#                f = 50000
#
#                for d in range(1000, 100100,1000):
#                     Ld50, Ldbeta,Ldp, Ldb  = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte,
#                                                                         Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p,C0=C0,C1=C1,C2=C2,D1=D1,D2=D2,D3=D3)
#                     data5.append(Ldb)

#        fig = plt.figure(1)
#        f  = ['10 GHz','20 GHz','30 GHz','40 GHz','50 GHz','60 GHz']
#        ax = fig.add_subplot(111)
#        ax.plot(eixo_x, data1)
#        ax.plot(eixo_x, data2)
#        ax.plot(eixo_x, data3)
#        ax.plot(eixo_x, data4)
#        ax.plot(eixo_x, data5)
#    #    ax.plot(eixo_x, data6)
#
#        # Add legend, title and axis labels
#        lgd = ax.legend( [ 'f = ' + str(lag) for lag in f], loc='upper center', bbox_to_anchor=(0.16, 1))
#        ax.set_title('Diffraction attenuation')
#        ax.set_xlabel('Distance (Km)')
#        ax.set_ylabel('Attenuation (dB)')
#        ax.set_xlim([0,100])
#        ax.set_ylim([130,220])
#        ax.grid(True)
#        fig.savefig('dif_att.png', dpi=350, format='png')


if __name__ == '__main__':
    unittest.main()
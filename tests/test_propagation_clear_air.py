# -*- coding: utf-8 -*-


import unittest
import numpy as np
import matplotlib.pyplot as plt

import numpy.testing as npt

from sharc.propagation.P452.propagation_clear_air_452 import PropagationClearAir

class PropagationClearAirTest(unittest.TestCase):
    
    def setUp(self):
        self.__ClearAir = PropagationClearAir()
        
    def test_loss(self):
        d = 10000
        f = 27000
        Ph = 1013 
        T = 288
        ro = 7.5
        Dlt = 30
        Dlr = 10
        Dct = 10
        Dcr = 10
        
        Hts = 244
        Hrs = 280
        Hte = 50
        Hre = 50
        thetaT = 20
        thetaR = 20
        N0 = 355
        deltaN = 60 
        p = 40
        
        omega = 0
        phi = 60
        dtm = 0.8
        dlm = 0.8                         
        epsilon = 3.5                            
        hm = 15
        Gt = 10
        Gr = 10

        Hsr = 45
        Hst = 48
        H0 = 15
        Hn = 17
        thetaJ = 0.3 
        ep = 0.8
        dsw = 20  
        k = 0.5
        di = [1,1,1]
        hi = [2,4,6]
        n = 2.5
        Aht = 0
        Ahr = 0
              
        C0 = 2.515516698
        C1 = 0.802853
        C2 = 0.010328
        D1 = 1.432788
        D2 = 0.189269 
        D3 = 0.001308
    
#        Ld50, Ldbeta, Ldb = self.__Diffraction.get_loss(beta = Beta, distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro, delta_N=deltaN, Hrs=hrs, Hts=hts, Hte=hte, Hre=hre, Hsr=hsr, Hst=hst, H0=h0, Hn=hn, dist_di=di, hight_hi=hi, omega=omega, Dlt=dlt ,Dlr=dlr, percentage_p=p)
#        
#        npt.assert_allclose(158.491,Ldb,atol=1e-3)
 
        #Grafico da perda de difraçao em funçao da distancia e da frequencia
        data1 = []
        data2 = []
        data3 = []
        data4 = []
        data5 = []
        eixo_x = [] 
        
        for n in range (1,6,1):
    
            if (n==1):
                f = 10000
                    
                for d in range(1000, 100100,1000):
           
                    Loss = self.__ClearAir.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro,Dlt=Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr,
                                                                    Hts=Hts,Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx=thetaT, theta_rx=thetaR, N0=N0, delta_N=deltaN, percentage_p=p,
                                                                    omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm, tx_gain=Gt, rx_gain=Gr, Hsr=Hsr, Hst=Hst, H0=H0, Hn=Hn,
                                                                    thetaJ=thetaJ, par_ep=ep, dsw=dsw, k=k, dist_di=di, hight_hi=hi, eta=n, Aht=Aht, Ahr=Ahr, C0=C0,C1=C1, C2=C2, D1=D1, D2=D2, D3=D3)
                                                                    
                    data1.append(Loss)
                    
                    eixo_x.append(d/1000)
                    
            if (n==2):
                f = 20000
                    
                for d in range(1000, 100100,1000):
                    Loss = self.__ClearAir.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro,Dlt=Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr,
                                                                    Hts=Hts,Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx=thetaT, theta_rx=thetaR, N0=N0, delta_N=deltaN, percentage_p=p,
                                                                    omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm, tx_gain=Gt, rx_gain=Gr, Hsr=Hsr, Hst=Hst, H0=H0, Hn=Hn,
                                                                    thetaJ=thetaJ, par_ep=ep, dsw=dsw, k=k, dist_di=di, hight_hi=hi, eta=n, Aht=Aht, Ahr=Ahr, C0=C0,C1=C1, C2=C2, D1=D1, D2=D2, D3=D3)
                    data2.append(Loss) 

            if (n==3):
                f = 30000
                    
                for d in range(1000, 100100,1000):
                    Loss = self.__ClearAir.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro,Dlt=Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr,
                                                                    Hts=Hts,Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx=thetaT, theta_rx=thetaR, N0=N0, delta_N=deltaN, percentage_p=p,
                                                                    omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm, tx_gain=Gt, rx_gain=Gr, Hsr=Hsr, Hst=Hst, H0=H0, Hn=Hn,
                                                                    thetaJ=thetaJ, par_ep=ep, dsw=dsw, k=k, dist_di=di, hight_hi=hi, eta=n, Aht=Aht, Ahr=Ahr, C0=C0,C1=C1, C2=C2, D1=D1, D2=D2, D3=D3)
                    data3.append(Loss)  

            if (n==4):
                f = 40000
                    
                for d in range(1000, 100100,1000):
                    Loss = self.__ClearAir.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro,Dlt=Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr,
                                                                    Hts=Hts,Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx=thetaT, theta_rx=thetaR, N0=N0, delta_N=deltaN, percentage_p=p,
                                                                    omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm, tx_gain=Gt, rx_gain=Gr, Hsr=Hsr, Hst=Hst, H0=H0, Hn=Hn,
                                                                    thetaJ=thetaJ, par_ep=ep, dsw=dsw, k=k, dist_di=di, hight_hi=hi, eta=n, Aht=Aht, Ahr=Ahr, C0=C0,C1=C1, C2=C2, D1=D1, D2=D2, D3=D3)
                    data4.append(Loss)  


            if (n==5):
                f = 50000
                    
                for d in range(1000, 100100,1000):
                    Loss = self.__ClearAir.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro,Dlt=Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr,
                                                                    Hts=Hts,Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx=thetaT, theta_rx=thetaR, N0=N0, delta_N=deltaN, percentage_p=p,
                                                                    omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm, tx_gain=Gt, rx_gain=Gr, Hsr=Hsr, Hst=Hst, H0=H0, Hn=Hn,
                                                                    thetaJ=thetaJ, par_ep=ep, dsw=dsw, k=k, dist_di=di, hight_hi=hi, eta=n, Aht=Aht, Ahr=Ahr, C0=C0,C1=C1, C2=C2, D1=D1, D2=D2, D3=D3)
                    data5.append(Loss)  

        fig = plt.figure(2)
        f  = ['10 GHz','20 GHz','30 GHz','40 GHz','50 GHz','60 GHz']
        ax = fig.add_subplot(111)
        ax.plot(eixo_x, data1)
        ax.plot(eixo_x, data2)
        ax.plot(eixo_x, data3)
        ax.plot(eixo_x, data4)
        ax.plot(eixo_x, data5)
    #    ax.plot(eixo_x, data6)  
        
        # Add legend, title and axis labels
        lgd = ax.legend( [ 'f = ' + str(lag) for lag in f], loc='upper center', bbox_to_anchor=(0.16, 1))
        ax.set_title('Overall prediction attenuation')
        ax.set_xlabel('Distance (Km)')
        ax.set_ylabel('Attenuation (dB)')
        ax.set_xlim([0,100])
        ax.set_ylim([0,60])
        ax.grid(True)
        fig.savefig('clear_air_att.png', dpi=350, format='png') 
           
          
  
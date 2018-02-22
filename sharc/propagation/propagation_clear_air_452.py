# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:10:11 2017

@author: LeticiaValle_Mac
"""
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation
from sharc.propagation.propagation_ducting_reflection import PropagationDuctingReflection
from sharc.propagation.propagation_troposcatter import PropagationTropScatter
from sharc.propagation.propagation_diffraction import PropagationDiffraction
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.support.enumerations import StationType

import numpy as np

class PropagationClearAir(Propagation):
    """
    Basic transmission loss due to free-space propagation and attenuation by atmospheric gases
    """
    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)

        self.propagationAg = PropagationGasesAttenuation(random_number_gen)
        self.propagationDucting = PropagationDuctingReflection(random_number_gen, self.propagationAg)
        self.propagationTropoScatter = PropagationTropScatter(random_number_gen, self.propagationAg)
        self.propagationDiffraction = PropagationDiffraction(random_number_gen, self.propagationAg)
        self.clutter = PropagationClutterLoss(random_number_gen)

        self.building_loss = 20


        ###########################################################################
        # Model coeficients
        self.coeff_r = 12.64
        self.coeff_s = 3.72
        self.coeff_t = 0.96
        self.coeff_u = 9.6
        self.coeff_v = 2
        self.coeff_w = 9.1
        self.coeff_x = -3
        self.coeff_y = 4.5
        self.coeff_z = -2

        ###########################################################################
        # Coeficients for the Approximation to the inverse cumulative normal distribution function
        self.C0 = 2.515516698
        self.C1 = 0.802853
        self.C2 = 0.010328
        self.D1 = 1.432788
        self.D2 = 0.189269
        self.D3 = 0.001308

        self.dsw = 20
        self.k = 0.5

        ###########################################################################
        # The distance of the i-th profile point
        self.dist_di = [0,0,0]
        #self.dist_di = [.1,.2,.5]


        ###########################################################################
        #Height of the i-th profile point
        self.height_hi = [2.2,4.3,6.4]
        #self.height_hi = [22, 100, 64]

    def get_loss(self, *args, **kwargs) -> np.array:

        d_km = np.asarray(kwargs["distance_3D"])*(1e-3)   #Km
        f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
        number_of_sectors = kwargs.pop("number_of_sectors",1)
        indoor_stations = kwargs.pop("indoor_stations",1)

        f = np.unique(f)
        if len(f) > 1:
            error_message = "different frequencies not supported in P619"
            raise ValueError(error_message)

        es_params =kwargs["es_params"]
        Ph = np.asarray(es_params.atmospheric_pressure)
        T = np.asarray(es_params.air_temperature)
        ro = np.asarray(es_params.water_vapour)

        Dlt = np.asarray(es_params.Dlt)
        Dlr = np.asarray(es_params.Dlr)
        Dct = np.asarray(es_params.Dct)
        Dcr = np.asarray(es_params.Dcr)

        Hts = np.asarray(es_params.Hts)
        Hrs = np.asarray(es_params.Hrs)
        Hte = np.asarray(es_params.Hte)
        Hre = np.asarray(es_params.Hre)
        thetaT = np.asarray(es_params.theta_tx)
        thetaR = np.asarray(es_params.theta_rx)
        N0 = np.asarray(es_params.N0)
        deltaN = np.asarray(es_params.delta_N)
        p = np.asarray(es_params.percentage_p)

        omega = np.asarray(es_params.omega)
        phi = np.asarray(es_params.phi)
        dtm = np.asarray(es_params.dtm)
        dlm = np.asarray(es_params.dlm)
        epsilon = np.asarray(es_params.epsilon)
        hm = np.asarray(es_params.hm)

        Gt = np.asarray(kwargs["tx_gain"])
        Gr = np.asarray(kwargs["rx_gain"])

        Hsr = np.asarray(es_params.Hsr)
        Hst = np.asarray(es_params.Hst)
        H0 = np.asarray(es_params.H0)
        Hn = np.asarray(es_params.Hn)
        thetaJ = np.asarray(es_params.thetaJ)
        ep = np.asarray(es_params.par_ep)
        dsw = np.asarray(self.dsw)
        k = np.asarray(self.k)
        di = kwargs.pop("di",self.dist_di)
        hi = kwargs.pop("hi",self.height_hi)
        n = np.asarray(es_params.eta)

        Stim = -np.inf * np.ones(d_km.size)

        #To find beta value
        tau = 1 - np.exp(-(4.12*(10**-4)*dlm**2.41))
        mu1 = (10**(-dtm/(16 - 6.6*tau)) + (10**-(0.496 + 0.354*tau))**5)**0.2

        if (abs(phi)<=70):
           mu4 = 10**(-0.935 + 0.0176*abs(phi)*np.log10(mu1))
        else:
           mu4 = 10**(0.3*np.log10(mu1))

        if (abs(phi)<=70):
           beta0 = 10**(-0.015*abs(phi)+ 1.67)*mu1*mu4
        else:
           beta0 = (4.17*mu1*mu4)


        k50 = 157/(157 - deltaN)
        ae = 6371*k50
        Ce = 1/ae

        val_hi =len(hi)

        # for i in range(0, val_hi,1):
        #     Stim_old = (hi[i] + 500*Ce*di[i]*(d_km - di[i]) - Hts)/di[i]
        #
        #     if Stim_old>Stim:
        #        Stim = Stim_old

        for k in range(d_km.size):
            for i in range(0, val_hi,1):
                Stim_old = (hi[i] + 500*Ce*di[i]*(d_km[0,k] - di[i]) - Hts)/di[i]

                if Stim_old>Stim[k]:
                   Stim[k] = Stim_old

        Str = (Hrs - Hts)/d_km

        #Approximation to the inverse cumulative normal distribution function
        Tx1 = np.sqrt((-2*np.log(p/100)))
        Ex1 = (((self.C2*Tx1+self.C1)*Tx1) +self.C0)/(((self.D3*Tx1 + self.D2)*Tx1 + self.D1)*Tx1 + 1)

        Tx2 = np.sqrt((-2*np.log(beta0/100)))
        Ex2 = (((self.C2*Tx2+self.C1)*Tx2) +self.C0)/(((self.D3*Tx2 + self.D2)*Tx2 + self.D1)*Tx2 + 1)

        I1 = Ex1 - Tx1
        I2 = Ex2 - Tx2
        Fi = I1/ I2


        Fj = 1.0 - 0.5*(1.0 + np.tanh(3.0*ep*(Stim - Str)/thetaJ))
        Fk = 1.0 - 0.5*(1.0 + np.tanh(3.0*k*(d_km - dsw)/dsw))

        Ag = self.propagationAg.get_loss(distance=d_km, frequency=f,
                                               atmospheric_pressure=Ph,
                                               air_temperature=T,
                                               water_vapour=ro)
        Lbfsg = 92.5 + 20 * np.log10(f) + 20*np.log10(d_km) + Ag
        # Ducting/layer reflection
        Lba =  self.propagationDucting.get_loss(distance=d_km*1000,frequency=f*1000,
                                                atmospheric_pressure=Ph, air_temperature= T,
                                                water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct,
                                                Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre,
                                                theta_tx = thetaT, theta_rx = thetaR, N0 = N0,
                                                delta_N = deltaN, percentage_p = p, omega=omega,
                                                phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)
        Lbs = self.propagationTropoScatter.get_loss(distance=d_km*1000, frequency=f*1000,
                                                    atmospheric_pressure=Ph, air_temperature= T,
                                                    water_vapour=ro, tx_gain = Gt, rx_gain = Gr,
                                                    theta_tx=thetaT, theta_rx=thetaR, N0=N0,
                                                    delta_N=deltaN, percentage_p=p,
                                                    number_of_sectors=number_of_sectors)       #Tropospheric scatter


        Esp = 2.6*(1 - np.exp(-0.1*(Dlt + Dlr)))*np.log10(p/50)
        Esbeta = 2.6*(1 - np.exp(-0.1*(Dlt + Dlr)))*np.log10(beta0/50)
        Lb0p = Lbfsg + Esp
        Lb0beta = Lbfsg + Esbeta


        Ld50 = np.empty(d_km.size)
        Ldbeta = np.empty(d_km.size)
        Ldp = np.empty(d_km.size)
        Ldb = np.empty(d_km.size)

        for index in range(d_km.size):
            Ld50[index], Ldbeta[index],Ldp[index], Ldb[index]  = \
                self.propagationDiffraction.get_loss(beta = beta0, distance=d_km[0,index]*1000,
                                                                   frequency=f*1000,
                                                                   atmospheric_pressure=Ph,
                                                                   air_temperature=T,
                                                                   water_vapour=ro,
                                                                   delta_N=deltaN,
                                                                   Hrs=Hrs, Hts=Hts,
                                                                   Hte=Hte, Hre=Hre,
                                                                   Hsr=Hsr, Hst=Hst,
                                                                   H0=H0, Hn=Hn,
                                                                   dist_di=di, hight_hi=hi,
                                                                   omega=omega, Dlt=Dlt ,Dlr=Dlr,
                                                                   percentage_p=p,
                                                                   C0=self.C0,C1=self.C1,C2=self.C2,
                                                                   D1=self.D1,D2=self.D2,D3=self.D3)

        Lbd50 = Lbfsg + Ld50
        Lbd = Lb0p + Ldp

        if (p < beta0):
            Lminb0p = Lb0p + (1 - omega)*Ldp
        if (p >= beta0):
            Lminb0p = Lbd50 + (Lb0beta + (1 - omega)*Ldp - Lbd50)*Fi


        Lminbap = n*np.log(np.exp(Lba/n) + np.exp(Lb0p/n))


        Lbda = Lbd
        index = np.where(Lminbap[0,:] <= Lbd[0,:])
        if index[0].size:
            Lbda[index] = Lminbap + (Lbd[index] - Lminbap)*Fk

        Lbam = Lbda + (Lminb0p - Lbda)*Fj

        if es_params.clutter_loss:
            clutter_loss = self.clutter.get_loss(frequency=f * 1000,
                                                 distance=d_km * 1000,
                                                 station_type=StationType.FSS_ES)
        else:
            clutter_loss = np.zeros(d_km.shape)

        building_loss = self.building_loss*indoor_stations

        if number_of_sectors > 1:
            Lbam = np.repeat(Lbam, number_of_sectors, 1)
            clutter_loss = np.repeat(clutter_loss, number_of_sectors, 1)
            building_loss = np.repeat(building_loss, number_of_sectors, 1)


        Lb = -5*np.log10(10**(-0.2*Lbs) + 10**(-0.2*Lbam))  + clutter_loss + building_loss




        #print(Lb)
        return Lb


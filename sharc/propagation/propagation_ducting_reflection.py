"""
Created on Tue Mai 08 12:05:38 2017

@author: LeticiaValle_Mac
"""
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation

import numpy as np

class PropagationDuctingReflection(Propagation):
    """
    Basic transmission loss due to free-space propagation and attenuation by atmospheric gases
    """

    def __init__(self, random_number_gen: np.random.RandomState, propagation_ag: PropagationGasesAttenuation):
        super().__init__(random_number_gen)

        self.propagation = propagation_ag

    def get_loss(self, *args, **kwargs) -> np.array:

        d = np.asarray(kwargs["distance"])*(1e-3) #Km
        f = np.asarray(kwargs["frequency"])*(1e-3) #GHz
        Ph = np.asarray(kwargs["atmospheric_pressure"])
        T = np.asarray(kwargs["air_temperature"])
        ro = np.asarray(kwargs["water_vapour"])

        Dlt = np.asarray(kwargs["Dlt"])
        Dlr = np.asarray(kwargs["Dlr"])
        Dct = np.asarray(kwargs["Dct"])
        Dcr = np.asarray(kwargs["Dcr"])

        Hts = np.asarray(kwargs["Hts"])
        Hrs = np.asarray(kwargs["Hrs"])
        Hte = np.asarray(kwargs["Hte"])
        Hre = np.asarray(kwargs["Hre"])
        thetaT = np.asarray(kwargs["theta_tx"])
        thetaR = np.asarray(kwargs["theta_rx"])
        deltaN = np.asarray(kwargs["delta_N"])
        p = np.asarray(kwargs["percentage_p"])

        omega = np.asarray(kwargs["omega"])
        phi = np.asarray(kwargs["phi"])
        dtm = np.asarray(kwargs["dtm"])
        dlm = np.asarray(kwargs["dlm"])
        epsilon = np.asarray(kwargs["epsilon"])
        hm = np.asarray(kwargs["hm"])


        #β0 (%), the time percentage for which refractive index lapse-rates exceeding
        #100 N-units/km can be expected in the first 100 m of the lower atmosphere,
        tau = 1 - np.exp(-(4.12*(10**-4)*dlm**2.41))
        alpha = -0.6 - epsilon*(10**-9)*(d**3.1)*tau
        dI = np.minimum(d - Dlt - Dlr, 40)

        k50 = 157/(157 - deltaN)
        Ae = 6371*k50               #Effective earth radius

        mu1 = (10**(-dtm/(16 - 6.6*tau)) + (10**-(0.496 + 0.354*tau))**5)**0.2

        mu2 = ((500/Ae)*((d**2)/((np.sqrt(Hte)+np.sqrt(Hre))**2)))**alpha



        if (hm <= 10):
            mu3 = 1
        else:
            mu3 = np.exp(-4.6*(10**-5)*(hm-10)*(43 + 6*dI))


        if (abs(phi)<=70):
            mu4 = 10**(-0.935 + 0.0176*abs(phi)*np.log10(mu1))
        else:
            mu4 = 10**(0.3*np.log10(mu1))


        if (abs(phi)<=70):
            B0 = (10**(-0.015*abs(phi)+ 1.67)*mu1*mu4)
        else:
            B0 = (4.17*mu1*mu4)

        if (f<0.5):
            Alf = 45.375 - 137*f + 92.5*f**2
        else:
            Alf = 0

        thetaT_line =  thetaT - 0.1*Dlt
        thetaR_line =  thetaR - 0.1*Dlr

        if (thetaT_line > 0):
            Ast = 20*np.log10(1 + 0.361*thetaT_line*(f*Dlt)**(1/2)) + 0.264*thetaT_line*f**(1/3)
        else:
            Ast = 0

        if (thetaR_line > 0):
            Asr = 20*np.log10(1 + 0.361*thetaR_line*(f*Dlr)**(1/2)) + 0.264*thetaR_line*f**(1/3)
        else:
            Asr = 0

        #over-sea surface duct coupling corrections for the interfering and interferedwith
        #stations respectively

        if ((omega >= 0.75) and (Dct >= 5) and (Dct > Dlt)):

            Act = -3*np.exp(-0.25*Dct**2)*(1 + np.tanh(0.07*(50 - Hts)))
        else:
            Act = 0

        if ((omega >= 0.75) and (Dcr >= 5) and (Dcr > Dlr)):

            Acr = -3*np.exp(-0.25*Dcr**2)*(1 + np.tanh(0.07*(50 - Hrs)))
        else:
            Acr = 0

        #Total of fixed coupling losses (dB)
        Af = 102.45 + 20*np.log10(f) + 20*np.log10(Dlt + Dlr) + Alf + Ast + Asr + Act +Acr


        Yd = 5*10**-5*Ae*f**(1/3)   #specific attenuation


        if (thetaT <= 0.1*Dlt):
            thetaT_oneline = thetaT
        else:
            thetaT_oneline = 0.1*Dlt

        if (thetaR <= 0.1*Dlr):
            thetaR_oneline = thetaR
        else:
            thetaR_oneline = 0.1*Dlr

        #Definition of the angular distance (mrad)
        teta_line = d*(10**3)/Ae +thetaT_oneline + thetaR_oneline

        beta = (B0*mu2*mu3)
        Gama = (1.076/((2.0058-np.log10(beta))**1.012))*np.exp(-((9.51-4.8*np.log10(beta))+ 0.198*(np.log10(beta))**2)*10**-6*d**1.13)

        Ap = -12 + (1.2 + 3.7*(10**-3)*d)*np.log10(p/beta) + 12*(p/beta)**Gama
        #Time percentage and time percentage and angular-distance dependent losses
        Ad = Yd*teta_line + Ap

        #Atmospheric gases attenuation
        loss_Ag = self.propagation.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)


        loss = Af + Ad + loss_Ag


        return loss

# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:10:11 2017
"""
from sharc.propagation.propagation import Propagation

from sharc.propagation.clear_air_452_aux import p676_ga
from sharc.propagation.clear_air_452_aux import inv_cum_norm
from sharc.support.enumerations import StationType
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss

import numpy as np

class PropagationClearAir(Propagation):
    """
    Basic transmission loss due to free-space propagation and attenuation by atmospheric gases
    """
    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)

        self.clutter = PropagationClutterLoss(random_number_gen)
        self.building_entry = PropagationBuildingEntryLoss(self.random_number_gen)

        self.building_loss = 20


    @staticmethod
    def closs_corr(f, d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r):
        # closs clutter loss correction according to P.452 - 16
        index1 = 1
        index2 = d.size

        htgc = htg
        hrgc = hrg

        Aht = 0
        Ahr = 0

        ha = ha_t
        dk = dk_t

        if ha > htg:

            Ffc = 0.25 + 0.375 * (1 + np.tanh(7.5 * (f - 0.5))) # (57a)
            Aht = 10.25 * Ffc * np.exp(-dk) * (1 - np.tanh(6 * (htg / ha - 0.625))) - 0.33 # (57)

            flagAht = 1

            kk = np.nonzero(d >= dk)

            if kk.size:
                index1 = kk[0]
            else:
                index1 = d.size

            htgc = ha_t

        ha = ha_r
        dk = dk_r

        if ha > hrg:
            Ffc = 0.25 + 0.375 * (1 + np.tanh(7.5 * (f - 0.5))) # (57a)
            Ahr = 10.25 * Ffc * np.exp(-dk) * (1 - np.tanh(6 * (hrg / ha - 0.625))) - 0.33 # (57)

            flagAhr = 1

            kk = np.nonzero(d <= d[-1] - dk)
            if kk.size:
                index2 = kk[-1]
            else:
                index2 = 1

            hrgc = ha_r

        # Modify the path

        if (index2 - index1 < 3): # at least two points between the clutter at Tx and Rx sides
            error_message = "tl_p452: closs_corr: the sum of clutter nominal distances is larger than the path length."
            raise ValueError(error_message)

        dc = d[index1-1:index2] - d[index1-1]
        hc = h[index1-1:index2]
        zonec = zone[index1-1:index2]

        return dc, hc, zonec,htgc, hrgc, Aht, Ahr

    @staticmethod
    def longest_cont_dist(d, zone, zone_r):
        dm = 0

        if zone_r == 12:
            aux = (zone == 1) + (zone == 2)
        else:
            aux = zone == zone_r

        aux = np.append(0,np.append(aux,0))
        aux = np.diff(aux)
        start = np.where(aux==1)[0]
        stop = np.where(aux==-1)[0] - 1

        start = np.atleast_1d(start)
        stop = np.atleast_1d(stop)
        n = start.size

        for i in range(n):
            delta = 0
            if (d[stop[i]] < d[-1]):
                delta = delta + (d[stop[i] + 1] - d[stop[i]]) / 2.0

            if (d[start[i]] > 0):
                delta = delta + (d[stop[i]] - d[stop[i] - 1]) / 2.0

            dm = max(d[stop[i]] - d[start[i]] + delta, dm)

        return dm

    @staticmethod
    def beta0(phi, dtm, dlm):

        tau = 1 - np.exp(-(4.12 * 1e-4 * dlm ** 2.41)) # (3a)

        mu1 = ( 10 ** (-dtm / (16 - 6.6 * tau)) + 10 ** (-5 * (0.496 + 0.354 * tau)))** 0.2

        indices = np.nonzero(mu1 > 1)
        mu1[indices] = 1

        if abs(phi) <= 70:
            mu4 = 10 ** ((-0.935 + 0.0176 * abs(phi)) * np.log10(mu1))
            b0 = 10 ** (-0.015 * abs(phi) + 1.67) * mu1 * mu4
        else:
            mu4 = 10 ^ (0.3 * np.log10(mu1))
            b0 = 4.17 * mu1 * mu4

        return b0

    @staticmethod
    def earth_rad_eff(DN):
        k50 = 157 / (157 - DN)
        ae = 6371 * k50

        kbeta = 3

        ab = 6371 * kbeta

        return ae, ab

    @staticmethod
    def smooth_earth_heights(d, h, htg, hrg, ae, f):
        n = d.size
        dtot = d[-1]

        # Tx and Rx antenna heights above mean sea level amsl(m)
        hts = h[0] + htg
        hrs = h[-1] + hrg

        # Section 5.1.6.2
        v1 = 0
        for ii in range (1,n):
            v1 = v1 + (d[ii] - d[ii - 1]) * (h[ii] + h[ii - 1]) # Eq(161)

        v2 = 0
        for ii in range(2, n):
            v2 = v2 + (d[ii] - d[ii - 1]) * (h[ii] * (2 * d[ii] + d[ii - 1]) + h[ii - 1] * (d[ii] + 2 * d[ii - 1])) # Eq(162)

        hst = (2 * v1 * dtot - v2) / dtot ** 2 # Eq(163)
        hsr = (v2 - v1 * dtot) / dtot ** 2 # Eq(164)

        # Section 5.1.6.3
        HH = h - (hts * (dtot - d) + hrs * d) / dtot # Eq(165d)
        hobs = max(HH[1:n - 1]) # Eq(165a)

        alpha_obt = max(HH[1:n - 1]/ d[1:n - 1]) # Eq(165b)
        alpha_obr = max(HH[1:n - 1]/ (dtot - d[1:n - 1])) # Eq(165c)

        # Calculate provisional values for the Tx and Rx smooth surface heights
        gt = alpha_obt / (alpha_obt + alpha_obr) # Eq(166e)
        gr = alpha_obr / (alpha_obt + alpha_obr) # Eq(166f)

        if hobs <= 0:
            hstp = hst
            hsrp = hsr
        else:
            hstp = hst - hobs * gt
            hsrp = hsr - hobs * gr

        # calculate the final values as required by the diffraction model
        if hstp >= h[0]:
            hstd = h[0]
        else:
            hstd = hstp

        if hsrp > h[-1]:
            hsrd = h[-1]
        else:
            hsrd = hsrp

        # Interfering antenna horizon elevation angle and distance
        ii = np.arange(1, n - 1)

        theta = 1000 * np.arctan((h[ii] - hts)/ (1000 * d[ii]) - d[ii] / (2 * ae))

        #theta(theta < 0) = 0; % condition below equation(152)

        theta_t = max(theta)

        theta_td = 1000 * np.arctan((hrs - hts)/ (1000 * dtot) - dtot / (2 * ae))
        theta_rd = 1000 * np.arctan((hts - hrs)/ (1000 * dtot) - dtot / (2 * ae))

        if theta_t > theta_td:
            pathtype = 2 # transhorizon
        else:
            pathtype = 1 # los

        kindex = np.nonzero(theta == theta_t)

        lt = kindex[0] + 1

        dlt = d[lt]

        # Interfered-with antenna horizon elevation angle and distance

        theta = 1000 * np.arctan((h[ii] - hrs)/(1000 * (dtot - d[ii])) - (dtot - d[ii])/(2 * ae))

        # theta(theta < 0) = 0;

        theta_r = max(theta)

        kindex = np.nonzero(np.ravel(theta) == theta_r)
        lr = kindex[-1] + 1

        dlr = dtot - d[lr]

        if pathtype == 1:
            theta_t = theta_td
            theta_r = theta_rd

            ii = np.arange(1,n - 1)

            lamb = 0.3 / f
            Ce = 1 / ae

            nu = (h[ii] + 500 * Ce * d[ii] * (dtot-d[ii])- (hts * (dtot- d[ii]) + hrs * d[ii]) / dtot)* \
                 np.sqrt(0.002 * dtot/( lamb * d[ii]*(dtot-d[ii])))
            numax = max(nu)

            kindex = np.nonzero(nu == numax)
            lt = kindex[-1] + 1
            dlt = d[lt]
            dlr = dtot - dlt
            kindex = np.nonzero(dlr <= dtot -d[ii])
            lr = kindex[0][-1] + 1

        # Angular distance

        theta_tot = 1e3 * dtot / ae + theta_t + theta_r

        #Section 5.1.6.4 Ducting / layer-reflection model

        # Calculate the smooth-Earth heights at transmitter and receiver as
        # required for the roughness factor

        hst = min(hst, h[0])
        hsr = min(hsr, h[-1])

        # Slope of the smooth - Earth surface
        m = (hsr - hst) / dtot

        #The terminal effective heigts for the ducting / layer - reflection model
        hte = htg + h[0] - hst
        hre = hrg + h[-1] - hsr

        ii = np.arange(lt,lr+1)
        hm = max(h[ii] - (hst + m * d[ii]))

        return hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype

    @staticmethod
    def path_fraction(d, zone, zone_r):
        dm = 0

        aux = np.nonzero(zone == zone_r)
        start = aux[0] # actually find_intervals
        stop = aux[-1]
        start = np.atleast_1d(start)
        stop = np.atleast_1d(stop)

        n = start.size

        for i in range(n):
            delta = 0
            if (d(stop[1]) < d[-1]):
                delta = delta + (d(stop[i] + 1) - d(stop[i])) / 2.0

            if (d(start[i]) > 0):
                delta = delta + (d(stop[i]) - d(stop[i] - 1)) / 2.0

            dm = dm + d(stop[i]) - d(start[i]) + delta

        omega = dm / (d[-1] - d[0])

        return omega

    @staticmethod
    def pl_los(d, f, p, b0, w, T, press, dlt, dlr):

        # water vapor density
        rho = 7.5 + 2.5 * w

        # compute specific attenuation due to dry air and water vapor:
        [g_0, g_w] = p676_ga(f, press, rho, T, True)

        Ag = (g_0 + g_w) * d

        # Basic transmission loss due to free - space propagation and attenuation
        # by atmospheric gases
        Lbfsg = 92.5 + 20.0 * np.log10(f) + 20.0 * np.log10(d) + Ag

        # Corrections for multipath and focusing effects at p and b0
        Esp = 2.6 * (1 - np.exp(-0.1 * (dlt + dlr))) * np.log10(p / 50)
        Esb = 2.6 * (1 - np.exp(-0.1 * (dlt + dlr))) * np.log10(b0 / 50)

        # Basic transmission loss not exceeded for time percentage p % due to
        # LoS propagation
        Lb0p = Lbfsg + Esp

        # Basic transmission loss not exceeded for time percentage b0 % due to
        # LoS propagation
        Lb0b = Lbfsg + Esb

        return Lbfsg, Lb0p, Lb0b

    @staticmethod
    def tl_tropo(dtot, theta, f, p, T, press, N0, Gt, Gr ):

        # Frequency dependent loss

        Lf = 25 * np.log10(f) - 2.5 * (np.log10(f / 2)) ** 2

        # aperture to  medium coupling loss(dB)
        Lc = 0.051 * np.exp(0.055 * (Gt + Gr))

        # gaseous absorbtion derived from equation (9) using rho = 3, g / m ^ 3 for the
        # whole path length

        rho = 3

        # compute specific attenuation due to dry air and water vapor:
        [g_0, g_w] = p676_ga(f, press, rho, T, True)

        Ag = (g_0 + g_w) * dtot

        # the basic transmission loss due to troposcatter not exceeded for any time
        # percentage p, below 50
        # is given

        Lbs = 190 + Lf + 20 * np.log10(dtot) + 0.573 * theta - 0.15 * N0 + Lc + Ag - 10.1 * (-np.log10(p / 50)) ** (0.7)
        return Lbs

    @staticmethod
    def tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, T, press,
                     omega, ae, b0):
        Alf = 0

        if f < 0.5:
            Alf = 45.375 - 137.0 * f + 92.5 * f * f

        # site - shielding diffraction losses for the interfering and interfered -with
        # stations(48)

        theta_t1 = theta_t - 0.1 * dlt
        theta_r1 = theta_r - 0.1 * dlr

        Ast = 0
        Asr = 0

        if theta_t1 > 0:
            Ast = 20 * np.log10(1 + 0.361 * theta_t1 * np.sqrt(f * dlt)) + 0.264 * theta_t1 * f ** (1 / 3)

        if theta_r1 > 0:
            Asr = 20 * np.log10(1 + 0.361 * theta_r1 * np.sqrt(f * dlr)) + 0.264 * theta_r1 * f ** (1 / 3)

        # over - sea surface duct coupling correction for the interfering and
        # interfered-with stations(49) and (49a)
        Act = 0
        Acr = 0

        if dct <= 5:
            if dct <= dlt:
                if omega >= 0.75:
                    Act = -3 * np.exp(-0.25 * dct * dct) * (1 + np.tanh(0.07 * (50 - hts)))

        if dcr <= 5:
            if dcr <= dlr:
                if omega >= 0.75:
                    Acr = -3 * np.exp(-0.25 * dcr * dcr) * (1 + np.tanh(0.07 * (50 - hrs)))

        # specific attenuation(51)
        gamma_d = 5e-5 * ae * f ** (1 / 3)

        # angular distance(corrected where appropriate) (52 - 52a)

        theta_t1 = theta_t
        theta_r1 = theta_r

        if theta_t > 0.1 * dlt:
            theta_t1 = 0.1 * dlt

        if theta_r > 0.1 * dlr:
            theta_r1 = 0.1 * dlr

        theta1 = 1e3 * dtot / ae + theta_t1 + theta_r1

        dI = min(dtot - dlt - dlr, 40)

        mu3 = 1

        if hm > 10:
            mu3 = np.exp(-4.6e-5 * (hm - 10) * (43 + 6 * dI))

        tau = 1 - np.exp(-(4.12e-4 * dlm ** 2.41))

        epsilon = 3.5

        alpha = -0.6 - epsilon * 1e-9 * dtot ** (3.1) * tau

        if alpha < -3.4:
            alpha = -3.4
        # correction for path geometry:
        mu2 = (500 / ae * dtot ** 2 / (np.sqrt(hte) + np.sqrt(hre)) ** 2) ** alpha

        if mu2 > 1:
            mu2 = 1

        beta = b0 * mu2 * mu3

        #beta = max(beta, eps); % to avoid division by zero

        Gamma = 1.076 / (2.0058 - np.log10(beta)) ** 1.012 * \
                np.exp(-(9.51 - 4.8 * np.log10(beta) + 0.198 * (np.log10(beta)) ** 2) * 1e-6 * dtot ** (1.13))

        # time percentage variablity(cumulative distribution):
        Ap = -12 + (1.2 + 3.7e-3 * dtot) * np.log10(p / beta) + 12 * (p / beta) ** Gamma

        # time percentage and angular - distance dependent losses within the
        # anomalous propagation mechanism
        Adp = gamma_d * theta1 + Ap

        # gaseous absorbtion derived from equation (9) using rho = 3 g / m ^ 3 for the
        # whole path length

        # water vapor density
        rho = 7.5 + 2.5 * omega

        # compute specific attenuation due to dry air and water vapor:
        [g_0, g_w] = p676_ga(f, press, rho, T, True)

        Ag = (g_0 + g_w) * dtot

        # total of fixed coupling losses(except for local clutter losses) between
        # the antennas and the anomalous propagation structure within the atmosphere (47)
        Af = 102.45 + 20 * np.log10(f) + 20 * np.log10(dlt + dlr) + Alf + Ast + Asr + Act + Acr;

        # total basic transmission loss occuring during periods of anomalaous
        # propagation

        Lba = Af + Adp + Ag;

        return Lba

    @staticmethod

    def dl_bull(d, h, hts, hrs, ap, f):

        # Effective Earth curvature Ce(km ^ -1)
        Ce = 1 / ap

        # Wavelength in meters
        lamb = 0.3 / f

        # Complete path length
        dtot = d[-1]-d[0]

        # Find the intermediate profile point with the highest slope of the line
        # from the transmitter to the point

        di = d[1: -1]
        hi = h[1:- 1]

        Stim = np.max((hi + 500 * Ce * di * (dtot - di) - hts) / di)

        # Calculate the slope of the line from transmitter to receiver assuming a
        # LoS path
        Str = (hrs - hts) / dtot

        if Stim < Str: #Case 1, Path is LoS

            # Find the intermediate profile point with the highest diffraction parameter nu:
            numax = np.max(
                        (hi + 500 * Ce * di* (dtot - di) - (hts * (dtot - di) + hrs * di) / dtot)*
                        np.sqrt(0.002 * dtot/ (lamb *di * (dtot - di))))

            Luc = 0
            if numax > -0.78:
                Luc = 6.9 + 20 * np.log10(np.sqrt((numax - 0.1) ** 2 + 1) + numax - 0.1)
        else:
            # Path is transhorizon
            # Find the intermediate profile pointwith the highest slope of the
            # line from the receiver to the point
            Srim = np.max((hi + 500 * Ce * di * (dtot - di) - hrs) / (dtot - di))

            # Calculate the distance of the Bullington point from the transmitter:
            dbp = (hrs - hts + Srim * dtot) / (Stim + Srim)

            # Calculate the diffraction parameter, nub, for the Bullington point
            nub = (hts + Stim * dbp - (hts * (dtot - dbp) + hrs * dbp) / dtot) * \
                   np.sqrt(0.002 * dtot / (lamb *dbp*(dtot - dbp)))

            # The knife - edge loss for the Bullington point is given by
            Luc = 0
            if nub > -0.78:
                Luc = 6.9 + 20 * np.log10(np.sqrt((nub - 0.1) ** 2 + 1) + nub - 0.1)

        # For Luc calculated using either(17) or (21), Bullington diffraction loss
        # for the path is given by
        Lbull = Luc + (1 - np.exp(-Luc / 6.0)) * (10 + 0.02 * dtot)
        return Lbull

    @staticmethod
    def dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f):
        # Normalized factor for surface admittance for horizontal (1) and vertical
        # (2) polarizations
        K =  np.empty(2)
        K[0] = 0.036 * (adft * f)** (-1/3) * ((epsr - 1) ** 2 + (18 * sigma / f)** 2)** (-1 / 4)
        K[1] = K[0] * (epsr** 2 + (18 * sigma / f)** 2)** (1/2)

        # Earth ground / polarization parameter
        beta_dft = (1 + 1.6 * K** 2 + 0.67 * K**4)/(1 + 4.5 * K** 2 + 1.53 * K** 4)

        # Normalized distance
        X = 21.88 * beta_dft * (f/ adft ** 2)** (1 / 3) * d

        # Normalized transmitter and receiver heights
        Yt = 0.9575 * beta_dft * (f** 2 / adft) ** (1 / 3) * hte
        Yr = 0.9575 * beta_dft * (f ** 2 / adft) ** (1 / 3) * hre

        # Calculate the distance term given by:
        Fx = np.empty(2)
        for ii in range(2):
            if X[ii] >= 1.6:
                Fx[ii] = 11 + 10 * np.log10(X[ii]) - 17.6 * X[ii]
            else:
                Fx[ii] = -20 * np.log10(X[ii]) - 5.6488 * (X[ii])** 1.425

        Bt = beta_dft * Yt
        Br = beta_dft * Yr

        GYt = np.empty(2)
        GYr = np.empty(2)

        for ii in range(2):
            if Bt[ii] > 2:
                GYt[ii] = 17.6 * (Bt[ii] - 1.1) ** 0.5 - 5 * np.log10(Bt[ii] - 1.1) - 8
            else:
                GYt[ii] = 20 * np.log10(Bt[ii] + 0.1 * Bt[ii] ** 3)

            if Br[ii] > 2:
                GYr[ii] = 17.6 * (Br[ii] - 1.1)** 0.5 - 5 * np.log10(Br[ii] - 1.1) - 8
            else:
                GYr[ii] = 20 * np.log10(Br[ii] + 0.1 * Br[ii] ** 3)

            if GYr[ii] < 2 + 20 * np.log10(K[ii]):
                GYr[ii] = 2 + 20 * np.log10(K[ii])

            if GYt[ii] < 2 + 20 * np.log10(K[ii]):
                GYt[ii] = 2 + 20 * np.log10(K[ii])

        Ldft = -Fx - GYt - GYr

        return Ldft

    @staticmethod
    def dl_se_ft(d, hte, hre, adft, f, omega):
        # First - term part of the spherical - Earth diffraction loss over land
        epsr = 22
        sigma = 0.003

        Ldft_land = PropagationClearAir.dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f)

        # First - term part of the spherical - Earth diffraction loss over sea
        epsr = 80
        sigma = 5

        Ldft_sea = PropagationClearAir.dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f)

        # First - term spherical diffraction loss
        Ldft = omega * Ldft_sea + (1 - omega) * Ldft_land

        return Ldft


    @staticmethod
    def dl_se(d, hte, hre, ap, f, omega):
        # Wavelength in meters
        lamb = 0.3 / f

        # Calculate the marginal LoS distance for a smooth path
        dlos = np.sqrt(2 * ap) * (np.sqrt(0.001 * hte) + np.sqrt(0.001 * hre))

        if d >= dlos:
            #calculate diffraction loss Ldft using the method in Sec.4.2.2.1 for
            # adft = ap and set Ldsph to Ldft
            Ldsph = PropagationClearAir.dl_se_ft(d, hte, hre, ap, f, omega)
        else:
            #calculate the smallest clearance between the curved - Earth path and
            #the ray between the antennas, hse
            c = (hte - hre) / (hte + hre)
            m = 250 * d * d / (ap * (hte + hre))

            b = 2 * np.sqrt((m + 1) / (3*m)) * np.cos(np.pi / 3 + 1 / 3 * np.arccos(3*c / 2 * np.sqrt(3*m/(m+1)** 3)))

            dse1 = d / 2 * (1 + b)
            dse2 = d - dse1

            hse = (hte - 500 * dse1 * dse1 / ap) * dse2 + (hre - 500 * dse2 * dse2 / ap) * dse1
            hse = hse / d

            # Calculate the required clearance for zero diffraction loss
            hreq = 17.456 * np.sqrt(dse1 * dse2 * lamb / d)

            if hse > hreq:
                Ldsph =np.array([0,0])
            else:
                # calculate the modified effective Earth radius aem, which gives
                # marginal LoS at distance d
                aem = 500 * (d / (np.sqrt(hte) + np.sqrt(hre)))**2

                # Use the method in Sec.4.2.2.1 for adft ) aem to obtain Ldft
                Ldft = PropagationClearAir.dl_se_ft(d, hte, hre, aem, f, omega)

                if (Ldft < 0).any():
                    Ldsph =np.array([0,0])
                else:
                    Ldsph = (1 - hse / hreq) * Ldft

        return Ldsph

    @staticmethod
    def dl_delta_bull(d, h, hts, hrs, hstd, hsrd, ap, f, omega):
        # Use the method in 4.2 .1 for the actual terrain profile and antenna
        # heights.Set the resulting Bullington diffraction loss for the actual
        # path to Lbulla
        Lbulla = PropagationClearAir.dl_bull(d, h, hts, hrs, ap, f)

        # Use the method in 4.2.1 for a second time, with all profile heights hi
        # set to zero and modified antenna heights given by
        hts1 = hts - hstd
        hrs1 = hrs - hsrd
        h1 = np.zeros(h.size)

        # where hstd and hsrd are given in 5.1 .6.3 of Attachment 2. Set the
        # resulting Bullington diffraction loss for this smooth path to Lbulls
        Lbulls = PropagationClearAir.dl_bull(d, h1, hts1, hrs1, ap, f)

        # Use the method in 4.2.2 to calculate the spherical - Earth diffraction loss
        # for the actual path length (dtot) with
        hte = hts1
        hre = hrs1
        dtot = d[-1] - d[0]

        Ldsph = PropagationClearAir.dl_se(dtot, hte, hre, ap, f, omega)

        # Diffraction loss for the general paht is now given by
        Ld = np.empty(2)
        Ld[0] = Lbulla + max(Ldsph[0] - Lbulls, 0)
        Ld[1] = Lbulla + max(Ldsph[1] - Lbulls, 0)

        return Ld

    @staticmethod
    def dl_p( d, h, hts, hrs, hstd, hsrd, f, omega, p, b0, DN ):
        # Use the method in 4.2.3 to calculate diffraction loss Ld for effective
        # Earth radius ap = ae as given by equation(6a). Set median diffractino
        # loss to Ldp50

        [ae, ab] = PropagationClearAir.earth_rad_eff(DN)

        ap = ae

        Ld50 = PropagationClearAir.dl_delta_bull(d, h, hts, hrs, hstd, hsrd, ap, f, omega)

        if p == 50:
            Ldp = Ld50
        elif p < 50:
            # Use the method in 4.2.3 to calculate diffraction loss Ld for effective
            # Earth radius ap = abeta, as given in equation(6b).Set diffraction loss
            # not exceeded for beta0 % time Ldb = Ld
            ap = ab

            Ldb = PropagationClearAir.dl_delta_bull(d, h, hts, hrs, hstd, hsrd, ap, f, omega);

            # Compute the interpolation factor Fi
            if p > b0:
                Fi = inv_cum_norm(p / 100) / inv_cum_norm(b0 / 100)
            else:
                Fi = 1

            # The diffraction loss Ldp not exceeded for p of time is now given by
            Ldp = Ld50 + Fi * (Ldb - Ld50)

        return Ldp, Ld50


    def get_loss(self, *args, **kwargs) -> np.array:

        d_km = np.asarray(kwargs["distance_3D"])*(1e-3)   #Km
        f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
        number_of_sectors = kwargs.pop("number_of_sectors",1)
        indoor_stations = kwargs.pop("indoor_stations",1)
        elevation = kwargs["elevation"]

        f = np.unique(f)
        if len(f) > 1:
            error_message = "different frequencies not supported in P619"
            raise ValueError(error_message)

        es_params =kwargs["es_params"]
        Ph = np.asarray(es_params.atmospheric_pressure)
        T = np.asarray(es_params.air_temperature)
        Dct = np.asarray(es_params.Dct)
        Dcr = np.asarray(es_params.Dcr)
        Hte = np.asarray(es_params.Hte)
        Hre = np.asarray(es_params.Hre)
        N0 = np.asarray(es_params.N0)
        deltaN = np.asarray(es_params.delta_N)
        p = np.asarray(es_params.percentage_p)

        tx_lat = es_params.tx_lat
        rx_lat = es_params.rx_lat

        Gt = np.ravel(np.asarray(kwargs["tx_gain"]))
        Gr = np.ravel(np.asarray(kwargs["rx_gain"]))

        # Modify the path according to Section 4.5.4, Step 1  and compute clutter losses
        # consider no obstacles profile
        profile_length = 100
        num_dists = d_km.size
        d = np.empty([num_dists, profile_length])
        for ii in range(num_dists):
            d[ii, :] = np.linspace(0,d_km[0][ii],profile_length)

        h = np.zeros(d.shape)

        ha_t = []
        ha_r = []
        dk_t = []
        dk_r = []

        # Compute the path profile parameters
        # Path center latitude
        phi_path = (tx_lat + rx_lat) / 2

        # Compute dtm - the longest continuous land(inland + coastal) section of the great - circle path(km)
        # Compute  dlm - the longest continuous inland section of the great-circle path (km)

        dtm = np.empty(num_dists)
        dlm = np.empty(num_dists)

        zone = np.ones(profile_length) * 2
        for index in range(num_dists):

            zone_r = 12
            dtm[index] = self.longest_cont_dist(d[index,:], zone, zone_r)

            zone_r = 2
            dlm[index] = self.longest_cont_dist(d[index,:], zone, zone_r)

        #compute beta0
        b0 = self.beta0(phi_path, dtm, dlm)
        [ae, ab] = self.earth_rad_eff(deltaN)

        # Compute the path fraction over sea
        omega = self.path_fraction(d.transpose(), zone, 3)

        # Modify the path according to Section 4.5.4, Step 1 and compute clutter losses
        # only if not isempty ha_t and ha_r
        #[dc, hc, zonec, htgc, hrgc, Aht, Ahr] = self.closs_corr(f, d, h, zone, Hte, Hre, ha_t, ha_r, dk_t, dk_r)

        Lb = np.empty([1,num_dists])

        # Effective Earth curvature Ce(km ^ -1)
        Ce = 1 / ae

        # Wavelength in meters
        lamb = 0.3 / f

        # Calculate an interpolation factor Fj to take account of the path angular
        # distance(58)
        THETA = 0.3
        KSI = 0.8

        for ii in range(num_dists):
            [dc, hc, zonec, htg, hrg, Aht, Ahr] = self.closs_corr(f, d[ii,:], h[ii,:], zone, Hte, Hre, ha_t, ha_r, dk_t, dk_r)
            d[ii,:] = dc
            h[ii,:] = hc

            [hst, hsr, hstd, hsrd, hte,hre, hm, dlt,
             dlr, theta_t, theta_r, theta, pathtype] = self.smooth_earth_heights(d[ii,:], h[ii,:], htg, hrg, ae, f)

            dtot = d[ii,-1] - d[ii,0]

            # Tx and Rx antenna heights above mean sea level amsl(m)
            hts = hc[0] + htg
            hrs = hc[-1] + hrg

            # Find the intermediate profile point with the highest slope of the line
            # from the transmitter to the point

            if len(d[ii]) < 4:
                error_message = "tl_p452: path profile requires at least 4 points."
                raise ValueError(error_message)

            di = d[ii,1: -1]
            hi = h[ii,1: -1]

            Stim = max((hi + 500 * Ce * di * (dtot - di) - hts) / di)

            # Calculate the slope of the line from transmitter to receiver assuming a
            # LoS path
            Str = (hrs - hts) / dtot

            # changed the definition for Fj on 15DEC16.
            # Fj = 1.0 - 0.5 * (1.0 + tanh(3.0 * KSI * (theta - THETA) / THETA))
            Fj = 1.0 - 0.5 * (1.0 + np.tanh(3.0 * KSI * (Stim - Str) / THETA))

            # Calculate an interpolation factor, Fk, to take account of the great
            # circle path distance:
            dsw = 20
            kappa = 0.5

            Fk = 1.0 - 0.5 * (1.0 + np.tanh(3.0 * kappa * (dtot - dsw) / dsw))

            [Lbfsg, Lb0p, Lb0b] = self.pl_los(dtot, f, p, b0[ii], omega[ii], T, Ph, dlt, dlr)

            [Ldp, Ld50] = self.dl_p(d[ii], h[ii], hts, hrs, hstd, hsrd, f, omega[ii], p, b0[ii], deltaN)

            # The median basic transmission loss associated with diffraction Eq (43)
            Lbd50 = Lbfsg + Ld50

            # The basic tranmission loss associated with diffraction not exceeded for p % time Eq(44)
            Lbd = Lb0p + Ldp

            # A notional minimum basic transmission loss associated with LoS
            # propagation and over-sea sub-path diffraction
            Lminb0p = Lb0p + (1 - omega[ii]) * Ldp

            if p >= b0[ii]:
                Fi = inv_cum_norm(p / 100) / inv_cum_norm(b0[ii] / 100)
                Lminb0p = Lbd50 + (Lb0b + (1 - omega[ii]) * Ldp - Lbd50) * Fi

            # Calculate a notional minimum basic transmission loss associated with LoS
            # and transhorizon signal enhancements
            eta = 2.5

            Lba = self.tl_anomalous(dtot, dlt, dlr, Dct, Dcr, dlm[ii], hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, T, Ph,
                                    omega[ii], ae, b0[ii])

            Lminbap = eta * np.log(np.exp(Lba / eta) + np.exp(Lb0p / eta))

            # Calculate a notional basic transmission loss associated with diffraction
            # and LoS or ducting / layer reflection enhancements
            Lbda = Lbd
            if (Lbd >= Lminbap).any():
                Lbda = Lminbap + (Lbd - Lminbap) * Fk

            # Calculate a modified basic transmission loss, which takes diffraction and
            # LoS or ducting / layer - reflection enhancements into account
            Lbam = Lbda + (Lminb0p - Lbda) * Fj

            # Calculate the basic transmission loss due to troposcatter not exceeded
            # for any time percantage p
            Lbs = self.tl_tropo(dtot, theta, f, p, T, Ph, N0, Gt[ii], Gr[ii])

            # Calculate the final transmission loss not exceeded for p % time
            Lb_pol = -5 * np.log10(10 ** (-0.2 * Lbs) + 10** (-0.2 * Lbam)) + Aht + Ahr

            if (es_params.polarization).lower() == "horizontal":
                Lb[0,ii] = Lb_pol[0]
            elif (es_params.polarization).lower() == "vertical":
                Lb[0,ii] = Lb_pol[1]
            else:
                error_message = "invalid polarization"
                raise ValueError(error_message)

        if es_params.clutter_loss:
            clutter_loss = self.clutter.get_loss(frequency=f * 1000,
                                                 distance=d_km * 1000,
                                                 station_type=StationType.FSS_ES)
        else:
            clutter_loss = np.zeros(d_km.shape)

#        building_loss = self.building_loss * indoor_stations
        b_loss = np.transpose(self.building_entry.get_loss(f, elevation))
        building_loss = b_loss * indoor_stations

        if number_of_sectors > 1:
            Lb = np.repeat(Lb, number_of_sectors, 1)
            clutter_loss = np.repeat(clutter_loss, number_of_sectors, 1)
            building_loss = np.repeat(building_loss, number_of_sectors, 1)

        Lb_new = Lb + clutter_loss + building_loss

        return Lb_new


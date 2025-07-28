# -*- coding: utf-8 -*-

from sharc.support.enumerations import StationType
from sharc.mask.spectral_mask import SpectralMask

import numpy as np
import math
import matplotlib.pyplot as plt


class SpectralMaskMSS(SpectralMask):
    """
    Implements spectral mask for all MSS Space Stations and some MSS Earth Stations
    according to REC ITU-R SM.1541. This is a generic mask and should only be used in case
    another isn't given.

    Spurious boundary has a default of 200% of the bandwidth from band edge.

    Attributes:
        spurious_emissions (float): level of power emissions at spurious
            domain [dBm/MHz].
        delta_f_lin (np.array): mask delta f breaking limits in MHz. Delta f
            values for which the spectral mask changes value. In this context, delta f is the frequency distance to
            the transmission's edge frequencies
        freq_lim (no.array): frequency values for which the spectral mask
            changes emission value
        freq_mhz (float): center frequency of station in MHz
        band_mhs (float): transmitting bandwidth of station in MHz
        p_tx (float): station's transmit power in dBm/MHz
        mask_dbm (np.array): spectral mask emission values in dBm
    """

    ALREADY_WARNED_AGAINST_LONG_CALCULATIONS = False

    def __init__(
        self,
        freq_mhz: float,
        band_mhz: float,
        spurious_emissions: float,
    ):
        """
        Class constructor.

        Parameters:
            freq_mhz (float): center frequency of station in MHz
            band_mhs (float): transmitting bandwidth of station in MHz
            spurious_emissions (float): level of spurious emissions [dBm/MHz].
        """
        # Spurious domain limits [dBm/MHz]
        self.spurious_emissions = spurious_emissions

        if freq_mhz < 15000:
            if band_mhz > 20 and not self.ALREADY_WARNED_AGAINST_LONG_CALCULATIONS:
                self.ALREADY_WARNED_AGAINST_LONG_CALCULATIONS = True
                print("WARNING: SpectralMaskMSS may take noticeably long to calculate. Consider changing its integral step.")
            self.reference_bandwidth = 0.004
        else:
            self.reference_bandwidth = 1
        self.spurious_boundary = 2 * band_mhz

        self.delta_f_lim = np.linspace(
            0., self.spurious_boundary - self.reference_bandwidth,
            math.ceil(self.spurious_boundary / self.reference_bandwidth)
        )
        self.delta_f_lim = np.concatenate((self.delta_f_lim, [self.spurious_boundary]))

        # Attributes
        self.band_mhz = band_mhz
        self.freq_mhz = freq_mhz

        self.freq_lim = np.concatenate((
            (freq_mhz - band_mhz / 2) - self.delta_f_lim[::-1],
            (freq_mhz + band_mhz / 2) + self.delta_f_lim,
        ))

    def set_mask(self, p_tx):
        """
        Set the spectral mask (mask_dbm attribute) based on station type, operating frequency and transmit power.

        Parameters:
            p_tx (float): station transmit power.
        """
        # dBm/MHz
        # this should work for the document's dBsd definition
        # when we have a uniform PSD in assigned band
        self.p_tx = p_tx - 10 * np.log10(self.band_mhz) + 30

        # attenuation mask
        mask_dbsd = 40 * np.log10(
            # both in MHz, percentage is just division
            100 * (self.delta_f_lim[:-1] / self.band_mhz) / 50 + 1
            # 100 * ((self.delta_f_lim[:-1] + self.reference_bandwidth/2) / self.band_mhz) / 50 + 1
        )

        # functionally same as 0, but won't be ignored at spectral mask calculation
        # TODO: something better than this...
        mask_dbsd[mask_dbsd == 0.] = 1e-14

        mask_dbm = self.p_tx - mask_dbsd

        mask_dbm = np.concatenate((mask_dbm, [self.spurious_emissions]))

        self.mask_dbm = np.concatenate((
            mask_dbm[::-1],
            np.array([self.p_tx]),
            mask_dbm,
        ))


if __name__ == '__main__':
    # Initialize variables
    p_tx = 34.061799739838875
    freq = 2100
    band = 5
    spurious_emissions_dbm_mhz = -30

    # Create mask
    msk = SpectralMaskMSS(freq, band, spurious_emissions_dbm_mhz)
    msk.set_mask(p_tx)

    # Frequencies
    freqs = np.linspace(-15, 15, num=1000) + freq

    # Mask values
    mask_val = np.ones_like(freqs) * msk.mask_dbm[0]
    for k in range(len(msk.freq_lim) - 1, -1, -1):
        mask_val[np.where(freqs < msk.freq_lim[k])] = msk.mask_dbm[k]
        # set as p_tx instead of 0 for plotting
        if k == len(msk.delta_f_lim):
            mask_val[np.where(freqs < msk.freq_lim[k])] = msk.p_tx

    # Plot
    # plt.plot(freqs, 10**(mask_val/10))
    plt.plot(freqs, mask_val)
    plt.xlim([freqs[0], freqs[-1]])
    plt.xlabel(r"f [MHz]")
    plt.ylabel("Spectral Mask [dBm]")
    plt.grid()
    plt.show()

    oob_idxs = np.where((freqs > freq - msk.spurious_boundary) & (freqs < freq + msk.spurious_boundary))[0]

    # Plot
    plt.plot(freqs[oob_idxs], mask_val[oob_idxs])
    # plt.plot(freqs[oob_idxs], 10 **(mask_val[oob_idxs]/10))
    plt.xlim([freqs[oob_idxs][0], freqs[oob_idxs][-1]])
    plt.xlabel(r"f [MHz]")
    plt.ylabel("Spectral Mask before spurious region [mW]")
    plt.grid()
    plt.show()

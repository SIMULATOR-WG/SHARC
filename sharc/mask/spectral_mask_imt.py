# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from sharc.support.enumerations import StationType
from sharc.mask.spectral_mask import SpectralMask

import numpy as np
import matplotlib.pyplot as plt


class SpectralMaskImt(SpectralMask):
    """
    Implements spectral masks for IMT-2020 according to Document 5-1/36-E.
        The masks are in the document's tables 1 to 8.
    Implements alternative spectral masks for IMT-2020 for cases Document 5-1/36-E does not implement
        according to Documents ITU-R SM.1541-6, ITU-R SM.1539-1 and ETSI TS 138 104 V16.6.0.
        Uses alternative when: outdoor BS's with freq < 24.25GHz

    Attributes:
        spurious_emissions (float): level of power emissions at spurious
            domain [dBm/MHz].
        delta_f_lin (np.array): mask delta f breaking limits in MHz. Delta f
            values for which the spectral mask changes value. In this context, delta f is the frequency distance to
            the transmission's edge frequencies
        freq_lim (no.array): frequency values for which the spectral mask
            changes emission value
        sta_type (StationType): type of station to which consider the spectral
            mask. Possible values are StationType.IMT_BS and StationType.IMT_UE
        freq_mhz (float): center frequency of station in MHz
        band_mhs (float): transmitting bandwidth of station in MHz
        scenario (str): INDOOR or OUTDOOR scenario
        p_tx (float): station's transmit power in dBm/MHz
        mask_dbm (np.array): spectral mask emission values in dBm
        alternative_mask_used (bool): represents whether the alternative mask should be used or not
        ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE (int): A hardcoded value that specifies how many samples of a diagonal
            line may be taken. Is needed because oob_power is calculated expecting rectangles, so we approximate
            the diagonal with 'SAMPLESIZE' rectangles
    """
    ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE: int = 40

    def __init__(
        self,
        sta_type: StationType,
        freq_mhz: float,
        band_mhz: float,
        spurious_emissions: float,
        scenario="OUTDOOR",
    ):
        """
        Class constructor.

        Parameters:
            sta_type (StationType): type of station to which consider the spectral
                mask. Possible values are StationType.IMT_BS and StationType.
                IMT_UE
            freq_mhz (float): center frequency of station in MHz
            band_mhs (float): transmitting bandwidth of station in MHz
            spurious_emissions (float): level of spurious emissions [dBm/MHz].
            scenario (str): INDOOR or OUTDOOR scenario
        """
        # Spurious domain limits [dBm/MHz]
        self.spurious_emissions = spurious_emissions

        # conditions to use alternative mask
        self.alternative_mask_used = freq_mhz < 24250                 \
            and scenario != "INDOOR"             \
            and sta_type == StationType.IMT_BS   \
            and spurious_emissions in [-13, -30]

        # Mask delta f breaking limits [MHz]
        if self.alternative_mask_used:
            self.delta_f_lim = self.get_alternative_mask_delta_f_lim(
                freq_mhz, band_mhz,
            )
        else:
            # use value from 5-1/36-E
            self.delta_f_lim = np.array([0, 20, 400])

        # Attributes
        self.sta_type = sta_type
        self.scenario = scenario
        self.band_mhz = band_mhz
        self.freq_mhz = freq_mhz

        self.freq_lim = np.concatenate((
            (freq_mhz - band_mhz / 2) - self.delta_f_lim[::-1],
            (freq_mhz + band_mhz / 2) + self.delta_f_lim,
        ))

    def set_mask(self, p_tx=0):
        """
        Sets the spectral mask (mask_dbm attribute) based on station type,
        operating frequency and transmit power.

        Parameters:
            p_tx (float): station transmit power. Default = 0
        """
        if self.alternative_mask_used:
            self.mask_dbm = self.get_alternative_mask_mask_dbm(p_tx)
            return

        self.p_tx = p_tx - 10 * np.log10(self.band_mhz)

        # Set new transmit power value
        if self.sta_type is StationType.IMT_UE:
            # Table 8
            mask_dbm = np.array([-5, -13, self.spurious_emissions])

        elif self.sta_type is StationType.IMT_BS and self.scenario == "INDOOR":
            # Table 1
            mask_dbm = np.array([-5, -13, self.spurious_emissions])

        else:

            if (self.freq_mhz > 24250 and self.freq_mhz < 33400):
                if p_tx >= 34.5:
                    # Table 2
                    mask_dbm = np.array([-5, -13, self.spurious_emissions])
                else:
                    # Table 3
                    mask_dbm = np.array([
                        -5, np.max((p_tx - 47.5, -20)),
                        self.spurious_emissions,
                    ])
            elif (self.freq_mhz > 37000 and self.freq_mhz < 52600):
                if p_tx >= 32.5:
                    # Table 4
                    mask_dbm = np.array([-5, -13, self.spurious_emissions])
                else:
                    # Table 5
                    mask_dbm = np.array([
                        -5, np.max((p_tx - 45.5, -20)),
                        self.spurious_emissions,
                    ])
            elif (self.freq_mhz > 66000 and self.freq_mhz < 86000):
                if p_tx >= 30.5:
                    # Table 6
                    mask_dbm = np.array([-5, -13, self.spurious_emissions])
                else:
                    # Table 7
                    mask_dbm = np.array([
                        -5, np.max((p_tx - 43.5, -20)),
                        self.spurious_emissions,
                    ])
            else:
                # this will only be reached when spurious emission has been manually set to something invalid and
                # alternative mask should be used
                raise ValueError(
                    "SpectralMaskIMT cannot be used with current parameters. \
                    You may have set spurious emission to a value not in [-13,-30]",
                )

        self.mask_dbm = np.concatenate((
            mask_dbm[::-1], np.array([self.p_tx]),
            mask_dbm,
        ))

    def get_alternative_mask_delta_f_lim(self, freq_mhz: float, band_mhz: float) -> np.array:
        """
            Implements spectral masks for IMT-2020 outdoor BS's when freq < 26GHz,
                according to Documents ITU-R SM.1541-6, ITU-R SM.1539-1 and ETSI TS 138 104 V16.6.0.
            Reference tables are:
                - Table 1 in ITU-R SM.1541-6
                - Table 2 in ITU-R SM.1539-1
                - Table 6.6.4.2.1-2   in ETSI TS 138 104 V16.6.0
                - Table 6.6.4.2.2.1-2 in ETSI TS 138 104 V16.6.0
                - Table 6.6.5.2.1-1 in ETSI TS 138 104 V16.6.0 (to choose spurious emission, not an impementation table)
                - Table 6.6.5.2.1-2 in ETSI TS 138 104 V16.6.0 (to choose spurious emission, not an impementation table)
        """
        # ITU-R SM.1539-1 Table 2
        if (freq_mhz > 0.009 and freq_mhz <= 0.15):
            B_L = 0.00025
            B_U = 0.01
        elif (freq_mhz > 0.15 and freq_mhz <= 30):
            B_L = 0.004
            B_U = 0.1
        elif (freq_mhz > 30 and freq_mhz <= 1000):
            B_L = 0.025
            B_U = 10
        elif (freq_mhz > 1000 and freq_mhz <= 3000):
            B_L = 0.1
            B_U = 50
        elif (freq_mhz > 3000 and freq_mhz <= 10000):
            B_L = 0.1
            B_U = 100
        elif (freq_mhz > 10000 and freq_mhz <= 15000):
            B_L = 0.3
            B_U = 250
        elif (freq_mhz > 15000 and freq_mhz <= 26000):
            B_L = 0.5
            B_U = 500
        else:
            raise ValueError(f"Invalid frequency value {freq_mhz} for mask ITU-R SM.1539-1")

        # ITU-R SM.1541-6 Table 1 (same as using only ITU-R SM.1539-1 Table 2, but with less hardcoded values)
        B_L_separation = 2.5 * B_L
        B_U_separation = 1.5 * band_mhz + B_U
        B_N_separation = 2.5 * band_mhz

        if band_mhz < B_L:
            delta_f_spurious = B_L_separation
        elif band_mhz > B_U:
            delta_f_spurious = B_U_separation
        else:
            delta_f_spurious = B_N_separation

        diagonal = np.array([
            i * 5 / self.ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE for i in range(
                self.ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE,
            )
        ])

        # band/2 is subtracted from delta_f_spurious beacuse that specific interval is from frequency center
        rest_of_oob_and_spurious = np.array(
            [5, 10.0, delta_f_spurious - band_mhz / 2],
        )

        return np.concatenate((diagonal, rest_of_oob_and_spurious))

    def get_alternative_mask_mask_dbm(self, power: float = 0) -> np.array:
        """
            Implements spectral masks for IMT-2020 outdoor BS's when freq < 26GHz,
                according to Documents ITU-R SM.1541-6, ITU-R SM.1539-1 and ETSI TS 138 104 V16.6.0.
            Reference tables are:
                - Table 1 in ITU-R SM.1541-6
                - Table 2 in ITU-R SM.1539-1
                - Table 6.6.4.2.1-2   in ETSI TS 138 104 V16.6.0
                - Table 6.6.4.2.2.1-2 in ETSI TS 138 104 V16.6.0
                - Table 6.6.5.2.1-1   in ETSI TS 138 104 V16.6.0 (to choose spurious emission, not an impementation table)
                - Table 6.6.5.2.1-2   in ETSI TS 138 104 V16.6.0 (to choose spurious emission, not an impementation table)
        """
        self.p_tx = power - 10 * np.log10(self.band_mhz)

        if self.spurious_emissions == -13:
            # use document ETSI TS 138 104 V16.6.0 Table 6.6.4.2.1-2
            diagonal_samples = np.array([
                -7 - 7 / 5 *
                (i * 5 / self.ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE)
                for i in range(self.ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE)
            ])
            rest_of_oob_and_spurious = np.array([
                -14, -13, self.spurious_emissions,
            ])
            mask_dbm = np.concatenate(
                (diagonal_samples, rest_of_oob_and_spurious),
            )
        elif self.spurious_emissions == -30:
            # use document ETSI TS 138 104 V16.6.0 Table 6.6.4.2.2.1-2
            diagonal_samples = np.array([
                -7 - 7 / 5 *
                (i * 5 / self.ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE)
                for i in range(self.ALTERNATIVE_MASK_DIAGONAL_SAMPLESIZE)
            ])
            rest_of_oob_and_spurious = np.array(
                [-14, -15, self.spurious_emissions],
            )
            mask_dbm = np.concatenate(
                (diagonal_samples, rest_of_oob_and_spurious),
            )
        else:
            raise ValueError(
                "Alternative mask should only be used for spurious emissions -13 and -30",
            )

        return np.concatenate((
            mask_dbm[::-1], np.array([self.p_tx]),
            mask_dbm,
        ))


if __name__ == '__main__':
    # Initialize variables
    sta_type = StationType.IMT_BS
    p_tx = 34.061799739838875
    freq = 9000
    band = 100
    spurious_emissions_dbm_mhz = -30

    # Create mask
    msk = SpectralMaskImt(sta_type, freq, band, spurious_emissions_dbm_mhz)
    msk.set_mask(p_tx)

    # Frequencies
    freqs = np.linspace(-600, 600, num=1000) + freq

    # Mask values
    mask_val = np.ones_like(freqs) * msk.mask_dbm[0]
    for k in range(len(msk.freq_lim) - 1, -1, -1):
        mask_val[np.where(freqs < msk.freq_lim[k])] = msk.mask_dbm[k]

    # Plot
    plt.plot(freqs, mask_val)
    plt.xlim([freqs[0], freqs[-1]])
    plt.xlabel(r"$\Delta$f [MHz]")
    plt.ylabel("Spectral Mask [dBc]")
    plt.grid()
    plt.show()

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:05:58 2017

@author: edgar
"""

class ParametersImt(object):

    __instance = None

    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersImt.__instance is None:
            ParametersImt.__instance = object.__new__(cls)
        return ParametersImt.__instance

    ###########################################################################
    # Network topology. Possible values are "MACROCELL", "HOTSPOT", "SINGLE_BS"
    topology = "MACROCELL"

    ###########################################################################
    # Number of macrocell sites per cluster (must set to 19 in macrocell network)
    num_macrocell_sites = 19

    ###########################################################################
    # Number of clusters in macro cell topology
    num_clusters = 1

    ###########################################################################
    # Inter-site distance in macrocell network topology
    intersite_distance = 150

    ###########################################################################
    # Minimum 2D separation distance from BS to UE
    minimum_separation_distance_bs_ue = 10
    
    ###########################################################################
    # Defines if IMT service is the interferer or interfered-with service
    interfered_with = False

    ###########################################################################
    # IMT center frequency [MHz]
    frequency = 26000

    ###########################################################################
    # IMT bandwidth [MHz]
    bandwidth = 200

    ###########################################################################
    # IMT resource block bandwidth [MHz]
    rb_bandwidth = 0.200

    ###########################################################################
    # Amount of guard band wrt total bandwidth. Setting this parameter to 0.1
    # means that 10% of the total bandwidth will be used as guard band: 5% in
    # the lower
    guard_band_ratio = 0.1

    ###########################################################################
    # The load probability (or activity factor) models the statistical
    # variation of the network load by defining the number of fully loaded
    # base stations that are simultaneously transmitting
    bs_load_probability = .5

    ###########################################################################
    # Maximum base station transmit power [dBm]
    bs_tx_power = 46

    ###########################################################################
    # Base station height [m]
    bs_height = 6

    ###########################################################################
    # Adjacent channel leakage power Ratio of the base station [dB]
    bs_aclr = 40

    ###########################################################################
    # Adjacent channel selectivity of the base station [dB]
    bs_acs = 30

    ###########################################################################
    # Base station noise figure [dB]
    bs_noise_figure = 10

    ###########################################################################
    # User equipment noise temperature [K]
    bs_noise_temperature = 290

    ###########################################################################
    # Base station feed loss [dB]
    bs_feed_loss = 3

    ###########################################################################
    # Uplink attenuation factor used in link-to-system mapping
    ul_attenuation_factor = 0.4
    
    ###########################################################################
    # Uplink minimum SINR of the code set [dB]
    ul_sinr_min = -10

    ###########################################################################
    # Uplink maximum SINR of the code set [dB]
    ul_sinr_max = 22
    
    ###########################################################################
    # Number of UE that is allocated to each cell within to handover margin.
    # Remenber that in macrocell network each base station has 3 cells (sectors)
    ue_k = 3

    ###########################################################################
    # Multiplication factor that is used to ensure that the sufficient number
    # of UE's will distributed throughout ths system area such that the number
    # of K users is allocated to each cell. Normally, this values varies
    # between 2 and 10 according to the user drop method
    ue_k_m = 2

    ###########################################################################
    # Percentage of indoor UE's
    ue_indoor_percent = 0.05

    ###########################################################################
    # Power control algorithm
    # ue_tx_power_control = "ON",power control On
    # ue_tx_power_control = "OFF",power control Off
    ue_tx_power_control = "ON"

    ###########################################################################
    # Power per RB used as target value in dBm
    ue_tx_power_target = -95

    ###########################################################################
    # Alfa is the balacing factor for UEs with bad channel
    # and UEs with good channel
    ue_tx_power_alfa = 1

    ###########################################################################
    # UE maximum transmit power [dBm]
    ue_tx_power = 22
    
    ###########################################################################
    # UE height [m]
    ue_height = 1.5

    ###########################################################################
    # Adjacent channel leakage power Ratio of the user equipment [dB]
    ue_aclr = 35

    ###########################################################################
    # Adjacent channel selectivity of the user equipment [dB]
    ue_acs = 25

    ###########################################################################
    # User equipment noise figure [dB]
    ue_noise_figure = 10

    ###########################################################################
    # User equipment feed loss [dB]
    ue_feed_loss = 3

    ###########################################################################
    # User equipment body loss [dB]
    ue_body_loss = 4

    ###########################################################################
    # Downlink attenuation factor used in link-to-system mapping
    dl_attenuation_factor = 0.6
    
    ###########################################################################
    # Downlink minimum SINR of the code set [dB]
    dl_sinr_min = -10

    ###########################################################################
    # Downlink maximum SINR of the code set [dB]
    dl_sinr_max = 30

    ###########################################################################
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "CI" (close-in FS reference distance)
    #                                    "UMa" (Urban Macro - 3GPP)
    #                                    "UMi" (Urban Micro - 3GPP)
    #                                    "ABG" (Alpha-Beta-Gamma)
    channel_model = "ABG"
    
    ###########################################################################
    # Probability of line-of-sight (CI)
    line_of_sight_prob = 0.95

    ###########################################################################
    # If shadowing should be applied or not
    shadowing = False

    ###########################################################################
    # System receive noise temperature [K]
    noise_temperature = 290

    BOLTZMANN_CONSTANT = 1.38064852e-23


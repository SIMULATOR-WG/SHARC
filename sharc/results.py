# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:47:46 2017

@author: edgar
"""

from sharc.plot import Plot

import numpy as np
import os

class Results(object):
    
    def __init__(self):
        self.imt_ul_tx_power_density = list()
        self.imt_ul_tx_power = list()
        self.imt_ul_sinr_ext = list()
        self.imt_ul_sinr = list()
        self.imt_ul_snr = list()
        self.imt_ul_tput_ext = list()
        self.imt_ul_tput = list()

        self.imt_path_loss = list()
        self.imt_coupling_loss = list()
        self.imt_bs_antenna_gain = list()
        self.imt_ue_antenna_gain = list()
        
        self.system_imt_antenna_gain = list()
        self.imt_system_antenna_gain = list()

        self.imt_dl_tx_power_density = list()
        self.imt_dl_tx_power = list()
        self.imt_dl_sinr_ext = list()
        self.imt_dl_sinr = list()
        self.imt_dl_snr = list()
        self.imt_dl_tput_ext = list()
        self.imt_dl_tput = list()
        
        self.system_ul_coupling_loss = list()
        self.system_ul_interf_power = list()

        self.system_dl_coupling_loss = list()
        self.system_dl_interf_power = list()

        self.system_inr = list()
        self.system_inr_scaled = list()
        self.output_directory = "output"

        
    def generate_plot_list(self, n_bins):
        self.plot_list = list()
        if len(self.system_imt_antenna_gain) > 0:
            values, base = np.histogram(self.system_imt_antenna_gain, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Antenna gain [dBi]"
            y_label = "Probability of antenna gain < $X$"
            title = "[SYS] CDF of system antenna gain towards IMT stations"
            file_name = title
            #x_limits = (0, 25)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))        
        if len(self.imt_system_antenna_gain) > 0:
            values, base = np.histogram(self.imt_system_antenna_gain, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Antenna gain [dBi]"
            y_label = "Probability of antenna gain < $X$"
            title = "[IMT] CDF of IMT station antenna gain towards system"
            file_name = title
            #x_limits = (0, 25)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))        
        if len(self.imt_bs_antenna_gain) > 0:
            values, base = np.histogram(self.imt_bs_antenna_gain, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Antenna gain [dBi]"
            y_label = "Probability of antenna gain < $X$"
            title = "[IMT] CDF of BS antenna gain towards the UE"
            file_name = title
            x_limits = (0, 25)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))        
        if len(self.imt_ue_antenna_gain) > 0:
            values, base = np.histogram(self.imt_ue_antenna_gain, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Antenna gain [dBi]"
            y_label = "Probability of antenna gain < $X$"
            title = "[IMT] CDF of UE antenna gain towards the BS"
            file_name = title
            x_limits = (0, 25)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))        
        if len(self.imt_ul_tx_power_density) > 0:
            values, base = np.histogram(self.imt_ul_tx_power_density, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Transmit power density [dBm/Hz]"
            y_label = "Probability of transmit power density < $X$"
            title = "[IMT] CDF of UE transmit power density"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_ul_tx_power) > 0:
            values, base = np.histogram(self.imt_ul_tx_power, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Transmit power [dBm]"
            y_label = "Probability of transmit power < $X$"
            title = "[IMT] CDF of UE transmit power"
            file_name = title
            x_limits = (-40, 30)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_ul_sinr_ext) > 0:
            values, base = np.histogram(self.imt_ul_sinr_ext, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL SINR [dB]"
            y_label = "Probability of SINR < $X$"
            title = "[IMT] CDF of UL SINR with external interference"
            file_name = title
            x_limits = (-15, 20)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_ul_sinr) > 0:
            values, base = np.histogram(self.imt_ul_sinr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL SINR [dB]"
            y_label = "Probability of SINR < $X$"
            title = "[IMT] CDF of UL SINR"
            file_name = title
            x_limits = (-15, 20)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_ul_snr) > 0:
            values, base = np.histogram(self.imt_ul_snr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL SNR"
            x_label = "UL SNR [dB]"
            y_label = "Probability of SNR < $X$"
            file_name = title
            x_limits = (-15, 20)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_ul_tput_ext) > 0:
            values, base = np.histogram(self.imt_ul_tput_ext, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL throughput with external interference"
            x_label = "UL throughput [bits/s/Hz]"
            y_label = "Probability of UL throughput < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))            
        if len(self.imt_ul_tput) > 0:
            values, base = np.histogram(self.imt_ul_tput, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL throughput"
            x_label = "UL throughput [bits/s/Hz]"
            y_label = "Probability of UL throughput < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_path_loss) > 0:
            values, base = np.histogram(self.imt_path_loss, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of path loss"
            x_label = "Path loss [dB]"
            y_label = "Probability of path loss < $X$"
            file_name = title
            x_limits = (60, 160)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_coupling_loss) > 0:
            values, base = np.histogram(self.imt_coupling_loss, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of coupling loss"
            x_label = "Coupling loss [dB]"
            y_label = "Probability of coupling loss < $X$"
            file_name = title
            x_limits = (60, 160)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_dl_tx_power) > 0:
            values, base = np.histogram(self.imt_dl_tx_power, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "Transmit power [dBm]"
            y_label = "Probability of transmit power < $X$"
            title = "[IMT] CDF of DL transmit power"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_sinr_ext) > 0:
            values, base = np.histogram(self.imt_dl_sinr_ext, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "SINR [dB]"
            y_label = "Probability of SINR < $X$"
            title = "[IMT] CDF of DL SINR with external interference"
            file_name = title
            x_limits = (-20, 80)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_dl_sinr) > 0:
            values, base = np.histogram(self.imt_dl_sinr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "SINR [dB]"
            y_label = "Probability of SINR < $X$"
            title = "[IMT] CDF of DL SINR"
            file_name = title
            x_limits = (-20, 80)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_dl_snr) > 0:
            values, base = np.histogram(self.imt_dl_snr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of DL SNR"
            x_label = "SNR [dB]"
            y_label = "Probability of SNR < $X$"
            file_name = title
            x_limits = (-20, 80)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))
        if len(self.imt_dl_tput_ext) > 0:
            values, base = np.histogram(self.imt_dl_tput_ext, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of DL throughput with external interference"
            x_label = "Throughput [bits/s/Hz]"
            y_label = "Probability of throughput < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_tput) > 0:
            values, base = np.histogram(self.imt_dl_tput, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of DL throughput"
            x_label = "Throughput [bits/s/Hz]"
            y_label = "Probability of throughput < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.system_inr_scaled) > 0:
            values, base = np.histogram(self.system_inr_scaled, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[SYS] CDF of scaled system INR"
            x_label = "INR [dB]"
            y_label = "Probability of INR < $X$"
            file_name = title
            x_limits = (-80, -20)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))            
        if len(self.system_inr) > 0:
            values, base = np.histogram(self.system_inr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[SYS] CDF of system INR"
            x_label = "INR [dB]"
            y_label = "Probability of INR < $X$"
            file_name = title
            x_limits = (-80, -20)
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, x_lim=x_limits, y_lim=y_limits))            
            ###################################################################
            # now we plot INR samples
            x = np.arange(len(self.system_inr))
            y = np.array(self.system_inr)
            title = "[SYS] INR samples"
            x_label = "Number of samples"
            y_label = "INR [dB]"
            file_name = title
            #x_limits = (0, 800)
            #y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))            
            
    def write_files(self, snapshot_number: int):
        n_bins = 200
        file_extension = ".txt"
        header_text = "Results collected after " + str(snapshot_number) + " snapshots."
        self.generate_plot_list(n_bins)
            
        for plot in self.plot_list:
            np.savetxt(os.path.join(self.output_directory, plot.file_name + file_extension), 
                       np.transpose([plot.x, plot.y]), 
                       fmt="%.5f", delimiter="\t", header=header_text, 
                       newline=os.linesep)

      
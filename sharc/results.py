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
        self.tx_power_ul = list()
        self.interf_power_ul = list()
        self.sinr_ul = list()
        self.snr_ul = list()
        self.throughoput_ul = list()
        self.inr = list()
        self.coupling_loss_ul = list()
        self.coupling_loss_ue_sat = list()
        self.coupling_loss_bs_sat = list()
        
        self.imt_ul_tx_power_density = list()
        self.imt_ul_tx_power = list()
        self.imt_ul_sinr = list()
        self.imt_ul_snr = list()
        self.imt_ul_tput = list()
        self.imt_ul_coupling_loss = list()

        self.imt_dl_tx_power = list()
        self.imt_dl_sinr = list()
        self.imt_dl_snr = list()
        self.imt_dl_tput = list()
        self.imt_dl_coupling_loss = list()
        
        self.system_ul_coupling_loss = list()
        self.system_ul_interf_power = list()

        self.system_dl_coupling_loss = list()
        self.system_dl_interf_power = list()

        self.system_inr = list()
        self.output_directory = "output"

    def add_interf_power_ul(self, sample):
        self.interf_power_ul.extend(sample)
        
        self.plot_list = list()
        
    def generate_plot_list(self, n_bins):
        self.plot_list = list()
        if len(self.imt_ul_tx_power_density) > 0:
            values, base = np.histogram(self.imt_ul_tx_power_density, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL transmit power density [dBm/Hz]"
            y_label = "Probability of UL transmit power density < $X$"
            title = "[IMT] CDF of UL transmit power density"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_ul_tx_power) > 0:
            values, base = np.histogram(self.imt_ul_tx_power, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL transmit power [dBm]"
            y_label = "Probability of UL transmit power < $X$"
            title = "[IMT] CDF of UL transmit power"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_ul_sinr) > 0:
            values, base = np.histogram(self.imt_ul_sinr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL SINR [dB]"
            y_label = "Probability of SINR < $X$"
            title = "[IMT] CDF of UL SINR"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_ul_snr) > 0:
            values, base = np.histogram(self.imt_ul_snr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL SNR"
            x_label = "UL SNR [dB]"
            y_label = "Probability of SNR < $X$"
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
        if len(self.imt_ul_coupling_loss) > 0:
            values, base = np.histogram(self.imt_ul_coupling_loss, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL coupling loss"
            x_label = "UL coupling loss [dB]"
            y_label = "Probability of UL coupling loss < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_tx_power) > 0:
            values, base = np.histogram(self.imt_dl_tx_power, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL transmit power [dBm]"
            y_label = "Probability of UL transmit power < $X$"
            title = "[IMT] CDF of UL transmit power"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_sinr) > 0:
            values, base = np.histogram(self.imt_dl_sinr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            x_label = "UL SINR [dB]"
            y_label = "Probability of SINR < $X$"
            title = "[IMT] CDF of UL SINR"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_snr) > 0:
            values, base = np.histogram(self.imt_dl_snr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL SNR"
            x_label = "UL SNR [dB]"
            y_label = "Probability of SNR < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_tput) > 0:
            values, base = np.histogram(self.imt_dl_tput, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL throughput"
            x_label = "UL throughput [bits/s/Hz]"
            y_label = "Probability of UL throughput < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.imt_dl_coupling_loss) > 0:
            values, base = np.histogram(self.imt_dl_coupling_loss, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[IMT] CDF of UL coupling loss"
            x_label = "UL coupling loss [dB]"
            y_label = "Probability of UL coupling loss < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))
        if len(self.system_inr) > 0:
            values, base = np.histogram(self.system_inr, bins=n_bins)
            cumulative = np.cumsum(values)
            x = base[:-1]
            y = cumulative / cumulative[-1]
            title = "[FSS] CDF of FSS INR"
            x_label = "INR [dB]"
            y_label = "Probability of INR < $X$"
            file_name = title
            y_limits = (0, 1)
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name, y_lim=y_limits))            
            
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

      
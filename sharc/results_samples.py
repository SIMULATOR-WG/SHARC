# -*- coding: utf-8 -*-
"""
Created on Fri May 11 08:39:08 2018

@author: Calil
"""

from sharc.plot import Plot
from sharc.results import Results

import numpy as np
from os.path import join
from os import fdopen, remove, linesep
from tempfile import mkstemp
from shutil import move

class ResultsSamples(Results):
    
    def __init__(self, out_dir):
        super().__init__(out_dir)
        
    def generate_plot_list(self, *args):
        self.plot_list = list()
        if len(self.system_imt_antenna_gain) > 0:
            x = np.arange(len(self.system_imt_antenna_gain))
            y = np.array(self.system_imt_antenna_gain)
            title = "[SYS] Samples of system antenna gain towards IMT stations"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.system_imt_bs_antenna_gain) > 0:
            x = np.arange(len(self.system_imt_bs_antenna_gain))
            y = np.array(self.system_imt_bs_antenna_gain)
            title = "[SYS] Samples of system antenna gain towards IMT BS stations"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.system_imt_ue_antenna_gain) > 0:
            x = np.arange(len(self.system_imt_ue_antenna_gain))
            y = np.array(self.system_imt_ue_antenna_gain)
            title = "[SYS] Samples of system antenna gain towards IMT UE stations"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_system_antenna_gain) > 0:
            x = np.arange(len(self.imt_system_antenna_gain))
            y = np.array(self.imt_system_antenna_gain)
            title = "[IMT] Samples of IMT station antenna gain towards system"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_bs_system_antenna_gain) > 0:
            x = np.arange(len(self.imt_bs_system_antenna_gain))
            y = np.array(self.imt_bs_system_antenna_gain)
            title = "[IMT] Samples of IMT BS antenna gain towards system"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ue_system_antenna_gain) > 0:
            x = np.arange(len(self.imt_ue_system_antenna_gain))
            y = np.array(self.imt_ue_system_antenna_gain)
            title = "[IMT] Samples of IMT UE antenna gain towards system"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_bs_antenna_gain) > 0:
            x = np.arange(len(self.imt_bs_antenna_gain))
            y = np.array(self.imt_bs_antenna_gain)
            title = "[IMT] Samples of BS antenna gain towards the UE"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_bs_bs_antenna_gain) > 0:
            x = np.arange(len(self.imt_bs_bs_antenna_gain))
            y = np.array(self.imt_bs_bs_antenna_gain)
            title = "[IMT] Samples of BS antenna gain towards BS"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ue_antenna_gain) > 0:
            x = np.arange(len(self.imt_ue_antenna_gain))
            y = np.array(self.imt_ue_antenna_gain)
            title = "[IMT] Samples of UE antenna gain towards the BS"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ue_ue_antenna_gain) > 0:
            x = np.arange(len(self.imt_ue_ue_antenna_gain))
            y = np.array(self.imt_ue_ue_antenna_gain)
            title = "[IMT] Samples of UE antenna gain towards the UE"
            x_label = "Number of samples"
            y_label = "Antenna gain [dBi]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_tx_power_density) > 0:
            x = np.arange(len(self.imt_ul_tx_power_density))
            y = np.array(self.imt_ul_tx_power_density)
            title = "[IMT] Samples of UE transmit power density"
            x_label = "Number of samples"
            y_label = "Transmit power density [dBm/Hz]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_tx_power) > 0:
            x = np.arange(len(self.imt_ul_tx_power))
            y = np.array(self.imt_ul_tx_power)
            title = "[IMT] Samples of UE transmit power"
            x_label = "Number of samples"
            y_label = "Transmit power [dBm]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_rx_power) > 0:
            x = np.arange(len(self.imt_ul_rx_power))
            y = np.array(self.imt_ul_rx_power)
            title = "[IMT] Samples of UE receive power"
            x_label = "Number of samples"
            y_label = "Receive power [dBm]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_sinr_ext) > 0:
            x = np.arange(len(self.imt_ul_sinr_ext))
            y = np.array(self.imt_ul_sinr_ext)
            title = "[IMT] Samples of UL SINR with external interference"
            x_label = "Number of samples"
            y_label = "UL SINR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_sinr) > 0:
            x = np.arange(len(self.imt_ul_sinr))
            y = np.array(self.imt_ul_sinr)
            title = "[IMT] Samples of UL SINR"
            x_label = "Number of samples"
            y_label = "UL SINR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_bs_interf) > 0:
            x = np.arange(len(self.imt_ul_bs_interf))
            y = np.array(self.imt_ul_bs_interf)
            title = "[IMT] Samples of UL Interference"
            x_label = "Number of samples"
            y_label = "UL Interference [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_snr) > 0:
            x = np.arange(len(self.imt_ul_snr))
            y = np.array(self.imt_ul_snr)
            title = "[IMT] Samples of UL SNR"
            x_label = "Number of samples"
            y_label = "UL SNR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_inr) > 0:
            x = np.arange(len(self.imt_ul_inr))
            y = np.array(self.imt_ul_inr)
            title = "[IMT] Samples of UL interference-to-noise ratio"
            x_label = "Number of samples"
            y_label = "$I/N$ [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ul_tput_ext) > 0:
            x = np.arange(len(self.imt_ul_tput_ext))
            y = np.array(self.imt_ul_tput_ext)
            title = "[IMT] Samples of UL throughput with external interference"
            x_label = "Number of samples"
            y_label = "UL throughput [Mbits/s]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))        
        if len(self.imt_ul_tput) > 0:
            x = np.arange(len(self.imt_ul_tput))
            y = np.array(self.imt_ul_tput)
            title = "[IMT] Samples of UL throughput"
            x_label = "Number of samples"
            y_label = "UL throughput [Mbits/s]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_total_tput) > 0:
            x = np.arange(len(self.imt_total_tput))
            y = np.array(self.imt_total_tput)
            title = "[IMT] Samples of IMT total throughput"
            x_label = "Number of samples"
            y_label = "Throughput [Mbits/s]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_path_loss) > 0:
            x = np.arange(len(self.imt_path_loss))
            y = np.array(self.imt_path_loss)
            title = "[IMT] Samples of path loss"
            x_label = "Number of samples"
            y_label = "Path loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_bs_bs_path_loss) > 0:
            x = np.arange(len(self.imt_bs_bs_path_loss))
            y = np.array(self.imt_bs_bs_path_loss)
            title = "[IMT] Samples of BS-BS path loss"
            x_label = "Number of samples"
            y_label = "Path loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_ue_ue_path_loss) > 0:
            x = np.arange(len(self.imt_ue_ue_path_loss))
            y = np.array(self.imt_ue_ue_path_loss)
            title = "[IMT] Samples of UE-UE path loss"
            x_label = "Number of samples"
            y_label = "Path loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_coupling_loss_bs_bs) > 0:
            x = np.arange(len(self.imt_coupling_loss_bs_bs))
            y = np.array(self.imt_coupling_loss_bs_bs)
            title = "[IMT] Samples of BS to BS coupling loss"
            x_label = "Number of samples"
            y_label = "Coupling loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.system_dl_coupling_loss) > 0:
            x = np.arange(len(self.system_dl_coupling_loss))
            y = np.array(self.system_dl_coupling_loss)
            title = "[SYS] Samples of System to BS coupling loss"
            x_label = "Number of samples"
            y_label = "Coupling loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.system_ul_coupling_loss) > 0:
            x = np.arange(len(self.system_ul_coupling_loss))
            y = np.array(self.system_ul_coupling_loss)
            title = "[SYS] Samples of System to UE coupling loss"
            x_label = "Number of samples"
            y_label = "Coupling loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_coupling_loss_ue_ue) > 0:
            x = np.arange(len(self.imt_coupling_loss_ue_ue))
            y = np.array(self.imt_coupling_loss_ue_ue)
            title = "[IMT] Samples of UE to UE coupling loss"
            x_label = "Number of samples"
            y_label = "Coupling loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_coupling_loss) > 0:
            x = np.arange(len(self.imt_coupling_loss))
            y = np.array(self.imt_coupling_loss)
            title = "[IMT] Samples of coupling loss"
            x_label = "Number of samples"
            y_label = "Coupling loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_coupling_loss_all) > 0:
            x = np.arange(len(self.imt_coupling_loss_all))
            y = np.array(self.imt_coupling_loss_all)
            title = "[IMT] Samples of coupling loss between all BSs and UEs"
            x_label = "Number of samples"
            y_label = "Coupling loss [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_tx_power) > 0:
            x = np.arange(len(self.imt_dl_tx_power))
            y = np.array(self.imt_dl_tx_power)
            title = "[IMT] Samples of BS transmit power"
            x_label = "Number of samples"
            y_label = "Transmit power [dBm]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_rx_power) > 0:
            x = np.arange(len(self.imt_dl_rx_power))
            y = np.array(self.imt_dl_rx_power)
            title = "[IMT] Samples of BS received power"
            x_label = "Number of samples"
            y_label = "Receive power [dBm]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_sinr_ext) > 0:
            x = np.arange(len(self.imt_dl_sinr_ext))
            y = np.array(self.imt_dl_sinr_ext)
            title = "[IMT] Samples of DL SINR with external interference"
            x_label = "Number of samples"
            y_label = "SINR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_sinr) > 0:
            x = np.arange(len(self.imt_dl_sinr))
            y = np.array(self.imt_dl_sinr)
            title = "[IMT] Samples of DL SINR"
            x_label = "Number of samples"
            y_label = "SINR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_ue_interf) > 0:
            x = np.arange(len(self.imt_dl_ue_interf))
            y = np.array(self.imt_dl_ue_interf)
            title = "[IMT] Samples of DL Interference"
            x_label = "Number of samples"
            y_label = "Interference [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_snr) > 0:
            x = np.arange(len(self.imt_dl_snr))
            y = np.array(self.imt_dl_snr)
            title = "[IMT] Samples of DL SNR"
            x_label = "Number of samples"
            y_label = "SNR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_inr) > 0:
            x = np.arange(len(self.imt_dl_inr))
            y = np.array(self.imt_dl_inr)
            title = "[IMT] Samples of DL interference-to-noise ratio"
            x_label = "Number of samples"
            y_label = "$I/N$ [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_tput_ext) > 0:
            x = np.arange(len(self.imt_dl_tput_ext))
            y = np.array(self.imt_dl_tput_ext)
            title = "[IMT] Samples of DL throughput with external interference"
            x_label = "Number of samples"
            y_label = "Throughput [Mbits/s]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))
        if len(self.imt_dl_tput) > 0:
            x = np.arange(len(self.imt_dl_tput))
            y = np.array(self.imt_dl_tput)
            title = "[IMT] Samples of DL throughput"
            x_label = "Number of samples"
            y_label = "Throughput [Mbits/s]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))  
        if len(self.system_inr_scaled) > 0:
            x = np.arange(len(self.system_inr_scaled))
            y = np.array(self.system_inr_scaled)
            title = "[SYS] Samples of scaled system INR"
            x_label = "Number of samples"
            y_label = "INR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))  
        if len(self.system_ul_inr_scaled) > 0:
            x = np.arange(len(self.system_ul_inr_scaled))
            y = np.array(self.system_ul_inr_scaled)
            title = "[SYS] Samples of scaled system UL INR"
            x_label = "Number of samples"
            y_label = "INR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))  
        if len(self.system_dl_inr_scaled) > 0:
            x = np.arange(len(self.system_dl_inr_scaled))
            y = np.array(self.system_dl_inr_scaled)
            title = "[SYS]Samples of scaled system DL INR"
            x_label = "Number of samples"
            y_label = "INR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))               
        if len(self.system_inr) > 0:
            x = np.arange(len(self.system_inr))
            y = np.array(self.system_inr)
            title = "[SYS] INR samples"
            x_label = "Number of samples"
            y_label = "INR [dB]"
            file_name = title
            self.plot_list.append(Plot(x, y, x_label, y_label, title, file_name))  
                      
            
    def write_files(self, snapshot_number: int):
        file_extension = ".txt"
        new_header_text = "Results collected after " + str(snapshot_number) + " snapshots."
        old_header_text = "Results collected after " + str(snapshot_number) + " snapshots."
        self.generate_plot_list()
            
        for plot in self.plot_list:
            
            file_name = join(self.output_directory, plot.file_name + file_extension)
            self.replace(file_name,old_header_text,new_header_text)
            
            with open(file_name,'a') as f:
                np.savetxt(f,np.transpose([plot.x, plot.y]),
                           fmt="%.5f", delimiter="\t", 
                           newline=linesep)
            
    def replace(self,file_path, pattern, subst):
        #Create temp file
        fh, abs_path = mkstemp()
        with fdopen(fh,'w') as new_file:
            with open(file_path) as old_file:                
                for line in old_file:
                    new_file.write(line.replace(pattern, subst))
        #Remove original file
        remove(file_path)
        #Move new file
        move(abs_path, file_path)

      
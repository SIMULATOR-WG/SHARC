# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:47:46 2017

@author: edgar
"""

import glob
import os
import datetime
import re
import pathlib
import pandas as pd
from shutil import copy


class SampleList(list):
    """
    This class only exists so that no list property can be confused with a SampleList
    """

    pass


class Results(object):
    """Handle the output of the simulator"""

    # This should always be true for 1st samples flush
    overwrite_sample_files = True

    def __init__(self):
        # Transmit power density [dBm/Hz]
        self.imt_ul_tx_power_density = SampleList()
        self.imt_ul_tx_power = SampleList()
        # SINR [dB]
        self.imt_ul_sinr_ext = SampleList()
        # SINR [dB]
        self.imt_ul_sinr = SampleList()
        # SNR [dB]
        self.imt_ul_snr = SampleList()
        self.imt_ul_inr = SampleList()
        # Throughput [bits/s/Hz]
        self.imt_ul_tput_ext = SampleList()
        # Throughput [bits/s/Hz]
        self.imt_ul_tput = SampleList()

        self.imt_path_loss = SampleList()
        self.imt_coupling_loss = SampleList()
        # Antenna gain [dBi]
        self.imt_bs_antenna_gain = SampleList()
        # Antenna gain [dBi]
        self.imt_ue_antenna_gain = SampleList()

        # Antenna gain [dBi]
        self.system_imt_antenna_gain = SampleList()
        # Antenna gain [dBi]
        self.imt_system_antenna_gain = SampleList()
        # Path Loss [dB]
        self.imt_system_path_loss = SampleList()
        # Building entry loss [dB]
        self.imt_system_build_entry_loss = SampleList()
        # System diffraction loss [dB]
        self.imt_system_diffraction_loss = SampleList()

        self.imt_dl_tx_power_density = SampleList()
        # Transmit power [dBm]
        self.imt_dl_tx_power = SampleList()
        # SINR [dB]
        self.imt_dl_sinr_ext = SampleList()
        # SINR [dB]
        self.imt_dl_sinr = SampleList()
        # SNR [dB]
        self.imt_dl_snr = SampleList()
        # I/N [dB]
        self.imt_dl_inr = SampleList()
        # Throughput [bits/s/Hz]
        self.imt_dl_tput_ext = SampleList()
        # Throughput [bits/s/Hz]
        self.imt_dl_tput = SampleList()

        self.system_ul_coupling_loss = SampleList()
        self.system_ul_interf_power = SampleList()
        # Interference Power [dBm]

        self.system_dl_coupling_loss = SampleList()
        self.system_dl_interf_power = SampleList()
        # Interference Power [dBm]

        self.system_inr = SampleList()
        self.system_pfd = SampleList()
        self.system_rx_interf = SampleList()

        self.__sharc_dir = pathlib.Path(__file__).parent.resolve()

    def prepare_to_write(
        self,
        parameters_filename: str,
        overwrite_output: bool,
        output_dir="output",
        output_dir_prefix="output",
    ):
        self.output_dir_parent = output_dir

        if not overwrite_output:
            today = datetime.date.today()

            results_number = 1
            results_dir_head = (
                output_dir_prefix + "_" + today.isoformat() + "_" + "{:02n}"
            )
            self.create_dir(results_number, results_dir_head)
            copy(parameters_filename, self.output_directory)
        else:
            self.output_directory = self.__sharc_dir / self.output_dir_parent
            try:
                os.makedirs(self.output_directory)
            except FileExistsError:
                pass

        return self

    def create_dir(self, results_number: int, dir_head: str):
        """Creates the output directory if it doesn't exist.

        Parameters
        ----------
        results_number : int
            Increment used in directory name
        dir_head : str
            Directory name prefix

        Returns
        -------
        str
            output directory name
        """

        dir_head_complete = (
            self.__sharc_dir / self.output_dir_parent / dir_head.format(results_number)
        )

        try:
            os.makedirs(dir_head_complete)
            self.output_directory = dir_head_complete
        except FileExistsError:
            self.create_dir(results_number + 1, dir_head)

    def get_relevant_attributes(self):
        """
        Returns the attributes that are used for storing samples
        """
        self_dict = self.__dict__

        results_relevant_attr_names = list(
            filter(lambda x: isinstance(getattr(self, x), SampleList), self_dict)
        )

        return results_relevant_attr_names

    def write_files(self, snapshot_number: int):
        """Writes the sample data to the output file

        Parameters
        ----------
        snapshot_number : int
            Current snapshot number
        """
        results_relevant_attr_names = self.get_relevant_attributes()
        for attr_name in results_relevant_attr_names:
            file_path = os.path.join(
                self.output_directory,
                attr_name + ".csv",
            )
            samples = getattr(self, attr_name)
            if len(samples) == 0:
                continue
            df = pd.DataFrame({"samples": samples})
            if self.overwrite_sample_files:
                df.to_csv(file_path, mode="w", index=False)
            else:
                df.to_csv(file_path, mode="a", index=False, header=False)
            setattr(self, attr_name, SampleList())

        if self.overwrite_sample_files:
            self.overwrite_sample_files = False

    @staticmethod
    def load_many_from_dir(root_dir: str, *, only_latest=True) -> list["Results"]:
        output_dirs = list(glob.glob(f"{root_dir}/output_*"))

        if len(output_dirs) == 0:
            print("[WARNING]: Results.load_many_from_dir did not find any results")

        if only_latest:
            output_dirs = Results.get_most_recent_outputs_for_each_prefix(output_dirs)

        all_res = []
        for output_dir in output_dirs:
            res = Results()
            res.load_from_dir(output_dir)
            all_res.append(res)

        return all_res

    def load_from_dir(self, abs_path: str) -> "Results":
        self.output_directory = abs_path

        self_dict = self.__dict__
        results_relevant_attr_names = filter(
            lambda x: isinstance(getattr(self, x), SampleList), self_dict
        )

        for attr_name in results_relevant_attr_names:
            file_path = os.path.join(abs_path, f"{attr_name}.csv")
            if os.path.exists(file_path):
                try:
                    # Try reading the .csv file using pandas with different delimiters
                    try:
                        data = pd.read_csv(file_path, delimiter=",")
                    except pd.errors.ParserError:
                        data = pd.read_csv(file_path, delimiter=";")

                    # Ensure the data has exactly one column
                    if data.shape[1] != 1:
                        raise Exception(
                            f"The file with samples of {attr_name} should have a single column.",
                        )

                    # Remove rows that do not contain valid numeric values
                    data = data.apply(pd.to_numeric, errors="coerce").dropna()

                    # Ignore if there is no data
                    if data.empty:
                        continue
                    # Check if there is enough data to load results from.

                    setattr(self, attr_name, SampleList(data.to_numpy()[:, 0]))

                except Exception as e:
                    print(e)
                    raise Exception(
                        f"Error processing the sample file ({attr_name}.csv) for {attr_name}: {e}"
                    )

        return self

    @staticmethod
    def get_most_recent_outputs_for_each_prefix(dirnames: list[str]) -> list[str]:
        """
        Input:
            A list of output directories.
        Returns:
            A list containing the most recent output dirname for each output_prefix.
            Note that if full paths are provided, full paths are returned
        """
        res = {}

        for dirname in dirnames:
            prefix, date, id = Results.get_prefix_date_and_id(dirname)
            res.setdefault(prefix, {"date": date, "id": id, "dirname": dirname})
            if date > res[prefix]["date"]:
                res[prefix]["date"] = date
                res[prefix]["id"] = id
                res[prefix]["dirname"] = dirname
            if date == res[prefix]["date"] and id > res[prefix]["id"]:
                res[prefix]["id"] = id
                res[prefix]["dirname"] = dirname

        return list(map(lambda x: x["dirname"], res.values()))

    @staticmethod
    def get_prefix_date_and_id(dirname: str) -> (str, str, str):
        mtch = re.search("(.*)(20[2-9][0-9]-[0-1][0-9]-[0-3][0-9])_([0-9]{2})", dirname)
        prefix, date, id = mtch.group(1), mtch.group(2), mtch.group(3)
        return prefix, date, id

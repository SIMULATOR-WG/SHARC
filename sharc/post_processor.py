from sharc.results import Results

from dataclasses import dataclass, field
import plotly.graph_objects as go
import os
import numpy as np
import scipy
import typing


class FieldStatistics:
    field_name: str
    median: float
    mean: float
    variance: float
    confidence_interval: (float, float)
    standard_deviation: float

    def load_from_sample(
        self, field_name: str, sample: list[float], *, confidence=0.95
    ) -> "FieldStatistics":
        self.field_name = field_name
        self.median = np.median(sample)
        self.mean = np.mean(sample)
        self.variance = np.var(sample)
        self.standard_deviation = np.std(sample)
        # @important TODO: check if using t distribution here is correct
        self.confidence_interval = scipy.stats.norm.interval(
            confidence, loc=self.mean, scale=self.standard_deviation
        )
        return self

    def __str__(self):
        attr_names = filter(
            lambda x: x != "field_name"
            and not x.startswith("__")
            and not callable(getattr(self, x)),
            dir(self),
        )
        readable_attrs = "\n".join(
            list(
                map(
                    lambda attr_name: f"\t{attr_name}: {getattr(self, attr_name)}",
                    attr_names,
                )
            )
        )
        return f"""{self.field_name}:
{readable_attrs}"""


class ResultsStatistics:
    fields_statistics: list[FieldStatistics]
    results_output_dir: str = "default_output"

    def load_from_results(self, result: Results) -> "ResultsStatistics":
        """
        Loads all relevant attributes from result and generates their statistics
        """
        self.results_output_dir = result.output_directory
        self.fields_statistics = []
        attr_names = result.get_relevant_attributes()
        for attr_name in attr_names:
            samples = getattr(result, attr_name)
            if len(samples) == 0:
                continue
            self.fields_statistics.append(
                FieldStatistics().load_from_sample(attr_name, samples)
            )

        return self

    def write_to_results_dir(self, filename="stats.txt") -> "ResultsStatistics":
        """
        Writes statistics file to the same directory of the results loaded into this class
        """
        with open(os.path.join(self.results_output_dir, filename), "w") as f:
            f.write(str(self))

        return self

    def get_stat_by_name(self, field_name: str) -> typing.Union[None, FieldStatistics]:
        """
        Gets a single field's statistics by its name.
        E.g.: get_stat_by_name("system_dl_interf_power")
        Returns
            None if not found
            FieldStatistics if found only one match
        """
        stats_found = filter(lambda field_stat: field_stat.field_name == field_name, self.fields_statistics)

        if len(stats_found) > 1:
            raise Exception(
                f"ResultsStatistics.get_stat_by_name found more than one statistic by the field name '{field_name}'\n"
                + "You probably loaded more than one result to the same ResultsStatistics object"
            )

        if len(stats_found) == 0:
            return None

        return stats_found[0]

    def __str__(self):
        return f"[{self.results_output_dir}]\n{'\n'.join(list(map(str, self.fields_statistics)))}"


@dataclass
class PostProcessor:
    IGNORE_FIELD = {
        "title": "ANTES NAO PLOTAVAMOS ISSO, ENTÃO CONTINUA SEM PLOTAR",
        "x_label": "",
    }
    # TODO: move units in the result to the Results class instead of Plot Info?
    # TODO: rename x_label to something else when making plots other than cdf
    RESULT_FIELDNAME_TO_PLOT_INFO = {
        "imt_ul_tx_power_density": {
            "x_label": "Transmit power density [dBm/Hz]",
            "title": "[IMT] UE transmit power density",
        },
        "imt_ul_tx_power": {
            "x_label": "Transmit power [dBm]",
            "title": "[IMT] UE transmit power",
        },
        "imt_ul_sinr_ext": {
            "x_label": "SINR [dB]",
            "title": "[IMT] UL SINR with external interference",
        },
        "imt_ul_snr": {
            "title": "[IMT] UL SNR",
            "x_label": "SNR [dB]",
        },
        "imt_ul_inr": {
            "title": "[IMT] UL interference-to-noise ratio",
            "x_label": "$I/N$ [dB]",
        },
        "imt_ul_sinr": {
            "x_label": "SINR [dB]",
            "title": "[IMT] UL SINR",
        },
        "imt_system_build_entry_loss": {
            "x_label": "Building entry loss [dB]",
            "title": "[SYS] IMT to system building entry loss",
        },
        "imt_ul_tput_ext": {
            "title": "[IMT] UL throughput with external interference",
            "x_label": "Throughput [bits/s/Hz]",
        },
        "imt_ul_tput": {
            "title": "[IMT] UL throughput",
            "x_label": "Throughput [bits/s/Hz]",
        },
        "imt_path_loss": {
            "title": "[IMT] path loss",
            "x_label": "Path loss [dB]",
        },
        "imt_coupling_loss": {
            "title": "[IMT] coupling loss",
            "x_label": "Coupling loss [dB]",
        },
        "imt_bs_antenna_gain": {
            "x_label": "Antenna gain [dBi]",
            "title": "[IMT] BS antenna gain towards the UE",
        },
        "imt_ue_antenna_gain": {
            "x_label": "Antenna gain [dBi]",
            "title": "[IMT] UE antenna gain towards the BS",
        },
        "system_imt_antenna_gain": {
            "x_label": "Antenna gain [dBi]",
            "title": "[SYS] system antenna gain towards IMT stations",
        },
        "imt_system_antenna_gain": {
            "x_label": "Antenna gain [dBi]",
            "title": "[IMT] IMT station antenna gain towards system",
        },
        "imt_system_path_loss": {
            "x_label": "Path Loss [dB]",
            "title": "[SYS] IMT to system path loss",
        },
        "system_dl_interf_power": {
            "x_label": "Interference Power [dBm/MHz]",
            "title": "[SYS] system interference power from IMT DL",
        },
        "imt_system_diffraction_loss": {
            "x_label": "Building entry loss [dB]",
            "title": "[SYS] IMT to system diffraction loss",
        },
        "imt_dl_sinr_ext": {
            "x_label": "SINR [dB]",
            "title": "[IMT] DL SINR with external interference",
        },
        "imt_dl_sinr": {
            "x_label": "SINR [dB]",
            "title": "[IMT] DL SINR",
        },
        "imt_dl_snr": {
            "title": "[IMT] DL SNR",
            "x_label": "SNR [dB]",
        },
        "imt_dl_inr": {
            "title": "[IMT] DL interference-to-noise ratio",
            "x_label": "$I/N$ [dB]",
        },
        "imt_dl_tput_ext": {
            "title": "[IMT] DL throughput with external interference",
            "x_label": "Throughput [bits/s/Hz]",
        },
        "imt_dl_tput": {
            "title": "[IMT] DL throughput",
            "x_label": "Throughput [bits/s/Hz]",
        },
        "system_ul_interf_power": {
            "title": "[SYS] system interference power from IMT UL",
            "x_label": "Interference Power [dBm]",
        },
        "system_inr": {
            "title": "[SYS] system INR",
            "x_label": "INR [dB]",
        },
        "system_pfd": {
            "title": "[SYS] system PFD",
            "x_label": "PFD [dBm/m^2]",
        },
        "imt_dl_tx_power": {
            "x_label": "Transmit power [dBm]",
            "title": "[IMT] DL transmit power",
        },
        # these ones were not plotted already, so will continue to not be plotted:
        "imt_dl_tx_power_density": IGNORE_FIELD,
        "system_ul_coupling_loss": IGNORE_FIELD,
        "system_dl_coupling_loss": IGNORE_FIELD,
        "system_rx_interf": IGNORE_FIELD,
    }

    plot_legend_patterns: list = field(default_factory=list)

    plots: list[go.Figure] = field(default_factory=list)
    results: list[Results] = field(default_factory=list)

    def add_plot_legend_pattern(
        self, *, dir_name_contains: str, legend: str
    ) -> "PostProcessor":
        self.plot_legend_patterns.append(
            {"dir_name_contains": dir_name_contains, "legend": legend}
        )
        self.plot_legend_patterns.sort(key=lambda p: -len(p["dir_name_contains"]))

        return self

    def generate_cdf_plots_from_results(
        self, results: list[Results], *, n_bins=200
    ) -> list[go.Figure]:
        figs: dict[str, list[go.Figure]] = {}

        for res in results:
            possible_legends_mapping = list(
                filter(
                    lambda pl: pl["dir_name_contains"]
                    in os.path.basename(res.output_directory),
                    self.plot_legend_patterns,
                )
            )

            if len(possible_legends_mapping):
                legend = possible_legends_mapping[0]["legend"]
            else:
                legend = res.output_directory

            attr_names = res.get_relevant_attributes()

            for attr_name in attr_names:
                attr_val = getattr(res, attr_name)
                if not len(attr_val):
                    continue
                if attr_name not in PostProcessor.RESULT_FIELDNAME_TO_PLOT_INFO:
                    print(
                        f"[WARNING]: {attr_name} is not a plottable field, because it does not have a configuration set on PostProcessor."
                    )
                    continue
                attr_plot_info = PostProcessor.RESULT_FIELDNAME_TO_PLOT_INFO[attr_name]
                if attr_plot_info == PostProcessor.IGNORE_FIELD:
                    print(
                        f"[WARNING]: {attr_name} is currently being ignored on plots."
                    )
                    continue
                if attr_name not in figs:
                    figs[attr_name] = go.Figure()
                    figs[attr_name].update_layout(
                        title=f'CDF Plot for {attr_plot_info["title"]}',
                        xaxis_title=attr_plot_info["x_label"],
                        yaxis_title="CDF",
                        yaxis=dict(tickmode="array", tickvals=[0, 0.25, 0.5, 0.75, 1]),
                        xaxis=dict(tickmode="linear", dtick=5),
                        legend_title="Labels",
                        meta={"related_results_attribute": attr_name},
                    )

                # TODO: take this fn as argument, to plot more than only cdf's
                x, y = PostProcessor.cdf_from(attr_val, n_bins=n_bins)

                fig = figs[attr_name]

                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=f"{legend}",
                    ),
                )

        return figs.values()

    def add_plots(self, plots: list[go.Figure]) -> None:
        self.plots.extend(plots)

    def add_results(self, results: list[Results]) -> None:
        self.results.extend(results)

    def get_results_by_output_dir(self, dir_name_contains: str, *, single_result=True):
        filtered_results = list(
            filter(
                lambda res: dir_name_contains in os.path.basename(res.output_directory),
                self.results,
            )
        )
        if len(filtered_results) == 0:
            raise ValueError(
                f"Could not find result that contains '{dir_name_contains}'"
            )

        if len(filtered_results) > 1:
            raise ValueError(
                f"There is more than one possible result with pattern '{dir_name_contains}'"
            )

        return filtered_results[0]

    def get_plot_by_results_attribute_name(self, attr_name: str) -> go.Figure:
        """
        You can get a plot using an attribute name from Results.
        See Results class to check what attributes exist.
        """
        filtered = list(
            filter(
                lambda x: x.layout.meta["related_results_attribute"] == attr_name,
                self.plots,
            )
        )

        if 0 == len(filtered):
            return None

        return filtered[0]

    @staticmethod
    def aggregate_results(
        *,
        dl_samples: list[float],
        ul_samples: list[float],
        ul_tdd_factor: float,
        n_bs_sim: int,
        n_bs_actual: int,
        random_number_gen=np.random.RandomState(31),
    ):
        """
        The method was adapted from document 'TG51_201805_E07_FSS_Uplink_ study 48GHz_GSMA_v1.5.pdf',
        a document created for Task Group 5/1.
        This is used to aggregate both uplink and downlink interference towards another system
            into a result that makes more sense to the case study.
        Inputs:
            downlink/uplink_result: list[float]
                Samples that should be aggregated.
            ul_tdd_factor: float
                The tdd ratio that uplink is activated for.
            n_bs_sim: int
                Number of simulated base stations.
                Should probably be 7 * 19 * 3 * 3 or 1 * 19 * 3 * 3
            n_bs_actual: int
                The number of base stations the study wants to have conclusions for.
            random_number_gen: np.random.RandomState
                Since this methods uses another montecarlo to aggregate results,
                it needs a random number generator
        """
        if ul_tdd_factor > 1 or ul_tdd_factor < 0:
            raise ValueError(
                "PostProcessor.aggregate_results() was called with invalid ul_tdd_factor parameter."
                + f"ul_tdd_factor must be in interval [0, 1], but is {ul_tdd_factor}"
            )

        segment_factor = round(n_bs_actual / n_bs_sim)

        dl_tdd_factor = 1 - ul_tdd_factor

        if ul_tdd_factor == 0:
            n_aggregate = len(dl_samples)
        elif dl_tdd_factor == 0:
            n_aggregate = len(ul_samples)
        else:
            n_aggregate = min(len(ul_samples), len(dl_samples))

        aggregate_samples = np.empty(n_aggregate)

        for i in range(n_aggregate):
            # choose S random samples
            ul_random_indexes = np.floor(
                random_number_gen.random(size=segment_factor)
                * len(ul_samples)
            )
            dl_random_indexes = np.floor(
                random_number_gen.random(size=segment_factor)
                * len(dl_samples)
            )
            aggregate_samples[i] = 0

            if ul_tdd_factor:
                for j in ul_random_indexes:  # random samples
                    aggregate_samples[i] += (
                        np.power(10, ul_samples[int(j)] / 10) * ul_tdd_factor
                    )

            if dl_tdd_factor:
                for j in dl_random_indexes:  # random samples
                    aggregate_samples[i] += (
                        np.power(10, dl_samples[int(j)] / 10) * dl_tdd_factor
                    )

            # convert back to dB or dBm (as was previously)
            aggregate_samples[i] = 10 * np.log10(
                aggregate_samples[i]
            )

        return aggregate_samples

    @staticmethod
    def cdf_from(data: list[float], *, n_bins=200) -> (list[float], list[float]):
        """
        Takes a dataset and returns both axis of a cdf (x, y)
        """
        values, base = np.histogram(
            data,
            bins=n_bins,
        )
        cumulative = np.cumsum(values)
        x = base[:-1]
        y = cumulative / cumulative[-1]

        return (x, y)

    @staticmethod
    def generate_statistics(result: Results) -> ResultsStatistics:
        return ResultsStatistics().load_from_results(result)

    @staticmethod
    def generate_sample_statistics(fieldname: str, sample: list[float]) -> ResultsStatistics:
        return FieldStatistics().load_from_sample(fieldname, sample)

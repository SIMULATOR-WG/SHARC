import os
from pathlib import Path
from sharc.results import Results
from sharc.post_processor import PostProcessor
import argparse

# Command line argument parser
parser = argparse.ArgumentParser(description="Generate and plot results.")
parser.add_argument("--auto_open", action="store_true", default=False, help="Set this flag to open plots in a browser.")
parser.add_argument("--scenario", type=int, choices=[0, 1], required=True,
                    help="Scenario parameter: 0 or 1. 0 for MSS-D2D to IMT-UL/DL,"
                    "1 for MSS-D2D to IMT-UL/DL with varying latitude.")
args = parser.parse_args()
scenario = args.scenario
auto_open = args.auto_open

local_dir = os.path.dirname(os.path.abspath(__file__))
campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
if scenario == 0:
    many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"),
                                              filter_fn=lambda x: "_co_channel_system_A" in x,
                                              only_latest=True)
    post_processor\
        .add_plot_legend_pattern(
            dir_name_contains="_mss_d2d_to_imt_ul_co_channel_system_A",
            legend="MSS-D2D to IMT-UL"
        )

    post_processor\
        .add_plot_legend_pattern(
            dir_name_contains="_mss_d2d_to_imt_dl_co_channel_system_A",
            legend="MSS-D2D to IMT-DL"
        )
elif scenario == 1:
    many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"),
                                              filter_fn=lambda x: "_lat_" in x,
                                              only_latest=True)
    for link in ["ul", "dl"]:
        for i in range(0, 70, 10):
            post_processor.add_plot_legend_pattern(
                dir_name_contains="_lat_" + link + "_" + str(i) + "_deg",
                legend="IMT-Link=" + link.upper() + " latitude=" + str(i) + "deg"
            )
else:
    raise ValueError("Invalid scenario. Choose 0 or 1.")

# ^: typing.List[Results]

post_processor.add_results(many_results)

plots = post_processor.generate_cdf_plots_from_results(
    many_results
)

post_processor.add_plots(plots)

plots = post_processor.generate_ccdf_plots_from_results(
    many_results
)

post_processor.add_plots(plots)

# Add a protection criteria line:
protection_criteria = -6
imt_dl_inr = post_processor.get_plot_by_results_attribute_name("imt_dl_inr", plot_type="ccdf")
imt_dl_inr.add_vline(protection_criteria, line_dash="dash", annotation=dict(
    text="Protection Criteria: " + str(protection_criteria) + " dB",
    xref="x", yref="y",
    x=protection_criteria + 0.5, y=0.8,
    font=dict(size=12, color="red")
))
imt_dl_inr.update_layout(template="plotly_white")
imt_ul_inr = post_processor.get_plot_by_results_attribute_name("imt_ul_inr", plot_type="ccdf")
imt_ul_inr.add_vline(protection_criteria, line_dash="dash")

# Combine INR plots into one:

for trace in imt_ul_inr.data:
    imt_dl_inr.add_trace(trace)

# Update layout if needed
imt_dl_inr.update_layout(title_text="CCDF Plot for IMT Downlink and Uplink INR",
                         xaxis_title="IMT INR [dB]",
                         yaxis_title="Cumulative Probability",
                         legend_title="Legend")

file = os.path.join(campaign_base_dir, "output", "imt_dl_ul_inr.html")
imt_dl_inr.write_html(file=file, include_plotlyjs="cdn", auto_open=auto_open)

file = os.path.join(campaign_base_dir, "output", "imt_system_antenna_gain.html")
imt_system_antenna_gain = post_processor.get_plot_by_results_attribute_name("imt_system_antenna_gain")
imt_system_antenna_gain.write_html(file=file, include_plotlyjs="cdn", auto_open=auto_open)

file = os.path.join(campaign_base_dir, "output", "system_imt_antenna_gain.html")
system_imt_antenna_gain = post_processor.get_plot_by_results_attribute_name("system_imt_antenna_gain")
system_imt_antenna_gain.write_html(file=file, include_plotlyjs="cdn", auto_open=auto_open)

file = os.path.join(campaign_base_dir, "output", "sys_to_imt_coupling_loss.html")
imt_system_path_loss = post_processor.get_plot_by_results_attribute_name("imt_system_path_loss")
imt_system_path_loss.write_html(file=file, include_plotlyjs="cdn", auto_open=auto_open)

for result in many_results:
    # This generates the mean, median, variance, etc
    stats = PostProcessor.generate_statistics(
        result=result
    ).write_to_results_dir()

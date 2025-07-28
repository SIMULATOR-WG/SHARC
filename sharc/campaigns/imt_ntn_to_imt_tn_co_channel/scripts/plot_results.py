import os
import numpy as np
from pathlib import Path
from sharc.results import Results
# import plotly.graph_objects as go
from sharc.post_processor import PostProcessor

post_processor = PostProcessor()

# Distance from topology boarders in meters
border_distances_array = np.array(
    [0, 10e3, 20e3, 30e3, 40e3, 50e3, 60e3, 70e3, 80e3, 90e3, 100e3])

# Add a legend to results in folder that match the pattern
for dist in border_distances_array:
    post_processor.add_plot_legend_pattern(
        dir_name_contains=f"_imt_ntn_to_imt_tn_co_channel_sep_{dist / 1e3}_km",
        legend=f"separation {dist / 1e3} Km"
    )

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"), only_latest=True)
# ^: typing.List[Results]

post_processor.add_results(many_results)

plots = post_processor.generate_cdf_plots_from_results(
    many_results
)

post_processor.add_plots(plots)

# Add a protection criteria line:
protection_criteria = -6
post_processor\
    .get_plot_by_results_attribute_name("imt_dl_inr")\
    .add_vline(protection_criteria, line_dash="dash")

# Plot every plot
plots_dir = Path(campaign_base_dir) / "output" / "plots"
plots_dir.mkdir(parents=True, exist_ok=True)
for plot in plots:
    plot.update_layout(legend_traceorder="normal")
    plot.write_html(
        plots_dir / f"{plot.layout.meta['related_results_attribute']}.html",
        include_plotlyjs="cdn",
        auto_open=False,
        config={"displayModeBar": False}
    )

for result in many_results:
    # This generates the mean, median, variance, etc
    stats = PostProcessor.generate_statistics(
        result=result
    ).write_to_results_dir()

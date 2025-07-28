import os
import pandas as pd
import plotly.graph_objects as go


def plot_cdf(
    base_dir, file_prefix, passo_xticks=5, xaxis_title='Value', legends=None, subfolders=None, save_file=True,
    show_plot=False,
):

    label = 'CDF'
    # Ensure at least one of save_file or show_plot is true
    if not save_file and not show_plot:
        raise ValueError("Either save_file or show_plot must be True.")

    # Define the base directory dynamically based on user input
    workfolder = os.path.dirname(os.path.abspath(__file__))
    csv_folder = os.path.abspath(
        os.path.join(
            workfolder, '..', "campaigns", base_dir, "output",
        ),
    )
    figs_folder = os.path.abspath(
        os.path.join(
            workfolder, '..', "campaigns", base_dir, "output", "figs",
        ),
    )

    # Check if the figs output folder exists, if not, create it
    if not os.path.exists(figs_folder):
        os.makedirs(figs_folder)

    # List all subfolders in the base directory or only those specified by the user
    if subfolders:
        subdirs = [
            os.path.join(csv_folder, d) for d in subfolders if os.path.isdir(
                os.path.join(csv_folder, d),
            )
        ]
    else:
        subdirs = [
            os.path.join(csv_folder, d) for d in os.listdir(csv_folder)
            if os.path.isdir(os.path.join(csv_folder, d)) and d.startswith(f"output_{base_dir}_")
        ]

    # Validate the number of legends
    if legends and len(legends) != len(subdirs):
        raise ValueError(
            "The number of provided legends does not match the number of found subfolders.",
        )

    # Initialize global min and max values
    global_min = float('inf')
    global_max = float('-inf')

    # First, calculate the global min and max values
    for subdir in subdirs:
        all_files = [
            f for f in os.listdir(subdir) if f.endswith('.csv',) and label in f and file_prefix in f
        ]

        for file_name in all_files:
            file_path = os.path.join(subdir, file_name)
            if os.path.exists(file_path):
                try:
                    # Try reading the .csv file using pandas with different delimiters
                    try:
                        data = pd.read_csv(
                            file_path, delimiter=',', skiprows=1,
                        )
                    except pd.errors.ParserError:
                        data = pd.read_csv(
                            file_path, delimiter=';', skiprows=1,
                        )

                    # Ensure the data has at least two columns
                    if data.shape[1] < 2:
                        print(
                            f"The file {file_name} does not have enough columns to plot.",
                        )
                        continue

                    # Remove rows that do not contain valid numeric values
                    data = data.apply(pd.to_numeric, errors='coerce').dropna()

                    if not data.empty:
                        global_min = min(global_min, data.iloc[:, 0].min())
                        global_max = max(global_max, data.iloc[:, 0].max())
                except Exception as e:
                    print(f"Error processing the file {file_name}: {e}")

    # If no valid data was found, set reasonable defaults for global_min and global_max
    if global_min == float('inf') or global_max == float('-inf'):
        global_min, global_max = 0, 1

    # Plot the graphs adjusting the axes
    fig = go.Figure()
    for idx, subdir in enumerate(subdirs):
        all_files = [
            f for f in os.listdir(subdir) if f.endswith('.csv',) and label in f and file_prefix in f
        ]
        legenda = legends[idx] if legends else os.path.basename(
            subdir,
        ).split(f"output_{base_dir}_")[1]

        for file_name in all_files:
            file_path = os.path.join(subdir, file_name)
            if os.path.exists(file_path):
                try:
                    # Try reading the .csv file using pandas with different delimiters
                    try:
                        data = pd.read_csv(
                            file_path, delimiter=',', skiprows=1,
                        )
                    except pd.errors.ParserError:
                        data = pd.read_csv(
                            file_path, delimiter=';', skiprows=1,
                        )

                    # Ensure the data has at least two columns
                    if data.shape[1] < 2:
                        print(
                            f"The file {file_name} does not have enough columns to plot.",
                        )
                        continue

                    # Remove rows that do not contain valid numeric values
                    data = data.apply(pd.to_numeric, errors='coerce').dropna()

                    # Check if there are enough data points to plot
                    if data.empty or data.shape[0] < 2:
                        print(
                            f"The file {file_name} does not have enough data to plot.",
                        )
                        continue

                    # Plot the CDF
                    fig.add_trace(
                        go.Scatter(x=data.iloc[:, 0], y=data.iloc[:, 1], mode='lines', name=f'{legenda}',),
                    )
                except Exception as e:
                    print(f"Error processing the file {file_name}: {e}")

    # Graph configurations
    fig.update_layout(
        title=f'CDF Plot for {file_prefix}',
        xaxis_title=xaxis_title,
        yaxis_title='CDF',
        yaxis=dict(tickmode='array', tickvals=[0, 0.25, 0.5, 0.75, 1]),
        xaxis=dict(tickmode='linear', tick0=global_min, dtick=passo_xticks),
        legend_title="Labels",
    )

    # Show the figure if requested
    if show_plot:
        fig.show()

    # Save the figure if requested
    if save_file:
        fig_file_path = os.path.join(
            figs_folder, f"CDF_Plot_{file_prefix}.png",
        )
        fig.write_image(fig_file_path)
        print(f"Figure saved: {fig_file_path}")

# Specific functions for different metrics


def plot_bs_antenna_gain_towards_the_ue(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'IMT_CDF_of_BS_antenna_gain_towards_the_UE', passo_xticks, xaxis_title='Antenna Gain (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_coupling_loss(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'IMT_CDF_of_coupling_loss', passo_xticks, xaxis_title='Coupling Loss (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_dl_sinr(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'IMT_CDF_of_DL_SINR', passo_xticks, xaxis_title='DL SINR (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_dl_snr(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'IMT_CDF_of_DL_SNR', passo_xticks, xaxis_title='DL SNR (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_dl_throughput(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'IMT_CDF_of_DL_throughput', passo_xticks, xaxis_title='DL Throughput (Mbps)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_dl_transmit_power(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'IMT_CDF_of_DL_transmit_power', passo_xticks, xaxis_title='DL Transmit Power (dBm)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_imt_station_antenna_gain_towards_system(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'IMT_CDF_of_IMT_station_antenna_gain_towards_system', passo_xticks,
        xaxis_title='Antenna Gain (dB)', legends=legends, subfolders=subfolders, save_file=save_file,
        show_plot=show_plot,
    )


def plot_path_loss(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'IMT_CDF_of_path_loss', passo_xticks, xaxis_title='Path Loss (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_ue_antenna_gain_towards_the_bs(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'IMT_CDF_of_UE_antenna_gain_towards_the_BS', passo_xticks, xaxis_title='Antenna Gain (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_imt_to_system_path_loss(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'SYS_CDF_of_IMT_to_system_path_loss', passo_xticks, xaxis_title='Path Loss (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_system_antenna_towards_imt_stations(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'SYS_CDF_of_system_antenna_gain_towards_IMT_stations', passo_xticks,
        xaxis_title='Antenna Gain (dB)', legends=legends, subfolders=subfolders, save_file=save_file,
        show_plot=show_plot,
    )


def plot_system_inr(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'SYS_CDF_of_system_INR', passo_xticks, xaxis_title='INR (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_system_interference_power_from_imt_dl(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'SYS_CDF_of_system_interference_power_from_IMT_DL', passo_xticks,
        xaxis_title='Interference Power (dBm/MHz)', legends=legends, subfolders=subfolders, save_file=save_file,
        show_plot=show_plot,
    )


def plot_system_interference_power_from_imt_ul(
    base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True,
    show_plot=True,
):
    plot_cdf(
        base_dir, 'SYS_CDF_of_system_interference_power_from_IMT_UL', passo_xticks,
        xaxis_title='Interference Power (dBm/MHz)', legends=legends, subfolders=subfolders, save_file=save_file,
        show_plot=show_plot,
    )


def plot_system_pfd(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'SYS_CDF_of_system_PFD', passo_xticks, xaxis_title='PFD (dBW/m²)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )


def plot_inr_samples(base_dir, passo_xticks=5, legends=None, subfolders=None, save_file=True, show_plot=True):
    plot_cdf(
        base_dir, 'INR_samples', passo_xticks, xaxis_title='INR Samples (dB)',
        legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )

# Main function to identify labels and call the appropriate functions


def all_plots(base_dir, legends=None, subfolders=None, save_file=True, show_plot=False):
    plot_bs_antenna_gain_towards_the_ue(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_coupling_loss(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )
    plot_dl_sinr(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )
    plot_dl_snr(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )
    plot_dl_throughput(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )
    plot_dl_transmit_power(
        base_dir, legends=legends,
        subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_imt_station_antenna_gain_towards_system(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_path_loss(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )
    plot_ue_antenna_gain_towards_the_bs(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_imt_to_system_path_loss(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_system_antenna_towards_imt_stations(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_system_inr(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )
    plot_system_interference_power_from_imt_dl(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_system_interference_power_from_imt_ul(
        base_dir, legends=legends, subfolders=subfolders, save_file=save_file, show_plot=show_plot,
    )
    plot_system_pfd(
        base_dir, legends=legends, subfolders=subfolders,
        save_file=save_file, show_plot=show_plot,
    )


if __name__ == "__main__":
    # Example usage
    name = "imt_hibs_ras_2600_MHz"
    legends = None  # Replace with a list of legends if needed
    subfolders = None  # Replace with a list of specific subfolders if needed
    all_plots(
        name, legends=legends, subfolders=subfolders,
        save_file=True, show_plot=False,
    )

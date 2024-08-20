import matplotlib.pyplot as plt
import pandas as pd
import argparse

from plot_helpers import plotATLASCaloBoundary, plotEtaGuidingLines, plotATLASCaloLayers


def main(args):
    """
    @brief Plot the transport data for a single event in the z-r plane.
    """

    # Read data
    data = pd.read_json(args.input)

    # Create new figure
    plt.figure()
    # Loop over the tracks
    for track_id in data.track_id.unique():
        track_data = data[data.track_id == track_id]
        # Plot the transport in the z-r plane
        plt.plot(
            track_data.z, track_data.r, "o--", markersize=3, color="black", zorder=10
        )

    # Remember the track zoom region
    track_x_lim = plt.xlim()
    track_y_lim = plt.ylim()

    # Plot the ATLAS calorimeter boundary
    plotATLASCaloBoundary()
    # Plot the eta guiding lines
    plotEtaGuidingLines()
    # Plot the ATLAS calorimeter layers
    plotATLASCaloLayers(args.calo_layer_csv)

    if args.track_zoom:
        # Legend not supported in this mode
        plt.xlim(track_x_lim)
        plt.ylim(track_y_lim)
        plt.legend().remove()
    else:
        plt.ylim([0, 4])
        plt.xlim([-6.5, 6.5])

    plt.xlabel("Z [m]")
    plt.ylabel("R [m]")

    plt.text(0.01, 0.95, rf"{args.label}", transform=plt.gca().transAxes)

    plt.savefig(args.output, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot particle transport data.")
    parser.add_argument("--input", type=str, required=True, help="Input file path")
    parser.add_argument("--output", type=str, required=True, help="Output file path")
    parser.add_argument(
        "--calo_layer_csv",
        type=str,
        required=True,
        help="CSV file containing the ATLAS calorimeter layer data",
    )
    parser.add_argument(
        "--track_zoom", action="store_true", help="Zoom in to the track region"
    )
    parser.add_argument("--label", type=str, required=False, help="Label for the plot")

    main(parser.parse_args())

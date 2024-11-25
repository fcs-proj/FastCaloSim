# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

import matplotlib.pyplot as plt
import pandas as pd
import argparse

from plot_helpers import plotATLASCaloBoundary, plotEtaGuidingLines, atlas_calo_layers


def main(args):
    """
    @brief Plot the extrapolation data for a single event in the z-r plane.
    """
    # Read data
    data = pd.read_json(args.input)

    # Create new figure
    plt.figure()

    colorList = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    plotted_layers = set()

    # Loop over extrapolation states (tracks)
    for state_id in data.state_id.unique():
        state_data = data[data.state_id == state_id]
        # Loop over the layers
        markerList = ["o", "s", "p"]
        for layer_id in state_data.layer_id.unique():
            # Plot the extrapolation to the ID-calo boundary
            if layer_id == -1:
                sub_pos_data = state_data[state_data.subpos == "IDCaloBoundary"]
                if layer_id not in plotted_layers:
                    plt.plot(
                        sub_pos_data.z,
                        sub_pos_data.r,
                        "o",
                        markersize=3,
                        color="red",
                        zorder=10,
                        label="ID-Calo Boundary",
                        linestyle="",
                    )
                    plotted_layers.add(layer_id)

                else:
                    plt.plot(
                        sub_pos_data.z,
                        sub_pos_data.r,
                        "o",
                        markersize=3,
                        color="red",
                        zorder=10,
                        linestyle="",
                    )

                continue

            # Get data for the current layer and requested subposition
            layer_data = state_data[
                (state_data.layer_id == layer_id) & (state_data.subpos == args.subpos)
            ]

            # change marker after 10 layers (max number of matplotlib colors)
            marker = markerList[layer_id // 10]

            # Check if the layer has already been plotted and labeled
            if layer_id not in plotted_layers:
                plt.plot(
                    layer_data.z,
                    layer_data.r,
                    color=colorList[layer_id % 10],
                    marker=marker,
                    markersize=3,
                    label=atlas_calo_layers[layer_id].name,
                    linestyle="",
                )
                plotted_layers.add(layer_id)
            else:
                plt.plot(
                    layer_data.z,
                    layer_data.r,
                    color=colorList[layer_id % 10],
                    marker=marker,
                    markersize=3,
                    linestyle="",
                )

    # Plot the ATLAS calorimeter boundary
    plotATLASCaloBoundary()
    # Plot the eta guiding lines
    plotEtaGuidingLines()

    plt.ylim([0, 4])
    plt.xlim([-8.5, 8.5])

    plt.xlabel("Z [m]")
    plt.ylabel("R [m]")

    plt.text(0.01, 0.95, rf"{args.label}", transform=plt.gca().transAxes)

    plt.legend(frameon=False, loc=(1.05, 0.1), fontsize=12, ncol=2, columnspacing=0.5)

    plt.savefig(args.output, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot particle transport data.")
    parser.add_argument("--input", type=str, required=True, help="Input file path")
    parser.add_argument("--output", type=str, required=True, help="Output file path")
    parser.add_argument("--label", type=str, required=False, help="Label for the plot")
    parser.add_argument(
        "--subpos", type=str, default="ENT", help="Subposition to plot (ENT, MID, EXT)"
    )

    main(parser.parse_args())

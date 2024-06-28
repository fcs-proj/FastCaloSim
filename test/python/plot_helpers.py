import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


from dataclasses import dataclass


@dataclass
class Layer:
    id: int
    name: str


layer_dict = {
    0: "PreSamplerB",
    1: "EMB1",
    2: "EMB2",
    3: "EMB3",
    4: "PreSamplerE",
    5: "EME1",
    6: "EME2",
    7: "EME3",
    8: "HEC0",
    9: "HEC1",
    10: "HEC2",
    11: "HEC3",
    12: "TileBar0",
    13: "TileBar1",
    14: "TileBar2",
    15: "TileGap1",
    16: "TileGap2",
    17: "TileGap3",
    18: "TileExt0",
    19: "TileExt1",
    20: "TileExt2",
    21: "FCAL0",
    22: "FCAL1",
    23: "FCAL2",
}

atlas_calo_layers = [Layer(id=key, name=value) for key, value in layer_dict.items()]


def plotATLASCaloBoundary(linewidth=1, color="black", alpha=0.5):
    barrelR1 = 1148.0 / 1000
    barrelR2 = 120.0 / 1000
    barrelR3 = 41.0 / 1000
    endcapZ1 = 3550.0 / 1000
    endcapZ2 = 4587.0 / 1000
    endcapZ3 = 6783 / 1000
    # offset to close lines better
    # Note: this depends on the fig size
    # TODO: find a robust way to make this work for all fig sizes
    offset = 0.020
    roffset = 0

    # Horizontal lines
    plt.hlines(
        y=barrelR1,
        xmin=-endcapZ1 - offset,
        xmax=+endcapZ1 + offset,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.hlines(
        y=barrelR2,
        xmin=endcapZ1 - offset,
        xmax=+endcapZ2 + offset,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.hlines(
        y=barrelR2,
        xmin=-endcapZ2 - offset,
        xmax=-endcapZ1 + offset,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.hlines(
        y=barrelR3,
        xmin=endcapZ2 - offset,
        xmax=endcapZ3 + offset,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.hlines(
        y=barrelR3,
        xmin=-endcapZ3 - offset,
        xmax=-endcapZ2 + offset,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )

    plt.vlines(
        x=endcapZ1,
        ymin=barrelR2 - roffset,
        ymax=barrelR1,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.vlines(
        x=-endcapZ1,
        ymin=barrelR2 - roffset,
        ymax=barrelR1,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.vlines(
        x=endcapZ2,
        ymin=barrelR3 - roffset,
        ymax=barrelR2,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )
    plt.vlines(
        x=-endcapZ2,
        ymin=barrelR3 - roffset,
        ymax=barrelR2,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
    )


def plotATLASCaloLayers(input_layer_csv: str):
    df = pd.read_csv(
        input_layer_csv,
        header=None,
        index_col=False,
        names=[
            "layerID",
            "isBarrel",
            "eta",
            "rEnt",
            "rMid",
            "rExt",
            "zEnt",
            "zMid",
            "zExt",
        ],
    )
    markerList = ["o", "s", "p"]

    for layer in range(len(atlas_calo_layers)):
        rValues = np.array(list(df[df["layerID"] == layer]["rEnt"]))
        zValues = np.array(list(df[df["layerID"] == layer]["zEnt"]))

        # change marker after 10 layers (max number of matplotlib colors)
        marker = markerList[layer // 10]
        plt.plot(
            zValues / 1000,
            rValues / 1000,
            marker,
            markersize=2,
            label=atlas_calo_layers[layer].name,
        )


def plotEtaGuidingLines(annotate_labels=False):
    start_point_list = [
        [1.07915876e-16, 1.76240000e00],
        [-1.81705442, 1.54616455],
        [-4.14229855, 1.14211689],
        [-4.16136164, 0.41539365],
        [1.81705442, 1.54616455],
        [4.14229855, 1.14211689],
        [4.16136164, 0.41539365],
    ]
    position_list = [
        [7.02947263e-17, 1.14800000e00],
        [-1.34913097, 1.148],
        [-3.55, 0.978808],
        [-3.55, 0.35436657],
        [1.34913097, 1.148],
        [3.55, 0.978808],
        [3.55, 0.35436657],
    ]
    for idx, eta in enumerate([0, -1, -2, -3, 1, 2, 3]):
        plt.plot(
            [5 * start_point_list[idx][0], position_list[idx][0]],
            [5 * start_point_list[idx][1], position_list[idx][1]],
            linestyle="dotted",
            color="purple",
            alpha=0.4,
        )

        if annotate_labels:
            if eta >= 0:
                annotatePointZ = 1.1 * start_point_list[idx][0] + 0.3
                annotatePointR = 1.1 * start_point_list[idx][1] - 0.06
                plt.annotate(rf"$\eta$={eta}", xy=(annotatePointZ, annotatePointR))

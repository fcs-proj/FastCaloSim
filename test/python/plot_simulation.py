import pygeosimplify as pgs
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import os
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import tempfile


def main(args):
    # Create a figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Create a temporary directory
    tmp_dir = tempfile.TemporaryDirectory()

    # Read the input transport json data
    trans_df = pd.read_json(args.transport_input)

    # Read the input simulation json data
    sim_df = pd.read_json(args.sim_input)

    # Set the coordinate system branches for pygeosimplify
    pgs.set_coordinate_branch("XYZ", "isCartesian")
    pgs.set_coordinate_branch("EtaPhiR", "isCylindrical")
    pgs.set_coordinate_branch("EtaPhiZ", "isECCylindrical")

    # Plot the cells colored according to their enegy
    ax = pgs.plot_geometry(
        sim_df,
        cell_energy_col="cell_energy",
        unit_scale=1e-3,
        unit_scale_energy=1e-3,
        axis_labels=["x [m]", "y [m]", "z [m]"],
        energy_label="Cell energy [GeV]",
        color_map="gist_heat_r",
    )

    # Hide the cells first to only plot transport first
    set_cell_visibliity(ax, False)

    # Loop over tracks
    for track_id in trans_df.track_id.unique():
        track_df = trans_df[trans_df.track_id == track_id]

        # Initialize the first point
        previous_step = np.array([track_df.x[0], track_df.y[0], track_df.z[0]])
        ax.plot(
            *previous_step,
            color="black",
            markerfacecolor="black",
            linewidth=2,
            marker="o",
            markersize=5,
            alpha=0.5,
        )
        plt.gca().text2D(0.01, 0.95, rf"{args.label}", transform=plt.gca().transAxes)
        plt.savefig(f"{tmp_dir.name}/figure_0.png", dpi=300, bbox_inches="tight")

        # Iterate over the DataFrame from the second point onwards
        for i in range(1, len(track_df)):
            current_step = np.array([track_df.x[i], track_df.y[i], track_df.z[i]])

            # Plot only if the current step is sufficiently different from the previous step
            if np.linalg.norm(current_step - previous_step) >= 1e-3:
                ax.plot(
                    track_df.x[: i + 1],
                    track_df.y[: i + 1],
                    track_df.z[: i + 1],
                    color="black",
                    markerfacecolor="black",
                    linewidth=2,
                    marker="o",
                    markersize=5,
                    alpha=0.5,
                )
                plt.savefig(
                    f"{tmp_dir.name}/figure_{i}.png", dpi=300, bbox_inches="tight"
                )

                # Update previous step
                previous_step = current_step

    # Finally show the cells
    set_cell_visibliity(ax, True)

    # Save for the animation
    plt.savefig(
        f"{tmp_dir.name}/figure_{len(trans_df)}.png", dpi=300, bbox_inches="tight"
    )
    # Save for final resultt
    plt.savefig(args.output_png, dpi=300, bbox_inches="tight")

    # Make the gif
    make_gif(tmp_dir.name, args.output_gif, duration=100, fade_steps=30)

    # Close the temporary directory
    tmp_dir.cleanup()


def set_cell_visibliity(ax, visible):
    for collection in ax.collections:
        if isinstance(collection, Poly3DCollection):
            collection.set_visible(visible)


def make_gif(dir, output_path, duration=100, fade_steps=30):
    # Get all the PNG files in the specified folder
    images = [img for img in os.listdir(dir) if img.endswith(".png")]
    # Sort the images
    images = sorted(images, key=lambda x: int(x.split("_")[1].split(".")[0]))

    # Load images into a list
    frames = [Image.open(os.path.join(dir, img)) for img in images]

    # Create fade-in frames for the last frame
    last_frame = frames[-1]
    fade_frames = create_fade_frames(frames[-2], last_frame, fade_steps)

    # Combine all frames, adding fade-in frames at the end
    all_frames = frames[:-1] + fade_frames

    # Pause the last frame for a bit
    for _ in range(10):
        all_frames.append(last_frame)

    # Save the frames as a GIF
    all_frames[0].save(
        output_path,
        save_all=True,
        append_images=all_frames[1:],
        duration=duration,
        loop=0,
    )


def create_fade_frames(prev_frame, last_frame, steps):
    fade_frames = []
    for i in range(steps):
        alpha = i / float(steps)
        blended = Image.blend(prev_frame, last_frame, alpha)
        fade_frames.append(blended)
    return fade_frames


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot particle transport data.")
    parser.add_argument(
        "--sim_input",
        type=str,
        required=True,
        help="Input file path to simulation state data",
    )
    parser.add_argument(
        "--transport_input",
        type=str,
        required=True,
        help="Input file path to transport state data",
    )
    parser.add_argument(
        "--output_png", type=str, required=True, help="Output of the final png plot"
    )
    parser.add_argument(
        "--output_gif",
        type=str,
        required=True,
        help="Output of the final gif animation",
    )
    parser.add_argument("--label", type=str, required=False, help="Label for the plot")

    main(parser.parse_args())

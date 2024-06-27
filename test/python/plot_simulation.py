import pygeosimplify as pgs
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(args):
    
    # Read the input json data
    df = pd.read_json(args.input)
    
    pgs.set_coordinate_branch("XYZ", "isCartesian")
    pgs.set_coordinate_branch("EtaPhiR", "isCylindrical")
    pgs.set_coordinate_branch("EtaPhiZ", "isECCylindrical")
    
    # Plot the cells colored according to their enegy
    pgs.plot_geometry(
        df,
        cell_energy_col="cell_energy",
        unit_scale=1e-3,
        unit_scale_energy=1e-3,
        axis_labels=["x [m]", "y [m]", "z [m]"],
        energy_label="Cell energy [GeV]",
        color_map="gist_heat_r",
    )
    
    plt.gca().text2D(0.01, 0.95, fr"{args.label}", transform=plt.gca().transAxes)

    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot particle transport data.')
    parser.add_argument('--input', type=str, required=True, help='Input file path')
    parser.add_argument('--output', type=str, required=True, help='Output file path')
    parser.add_argument('--label', type=str, required=False, help='Label for the plot')

    
    main(parser.parse_args())
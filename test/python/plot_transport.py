import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

def plotATLASCaloBoundary(linewidth = 1, color='black', alpha=0.5):
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
    plt.hlines(y = barrelR1,  xmin = -endcapZ1 - offset, xmax= +endcapZ1 + offset, color=color,  alpha=alpha, linewidth=linewidth)
    plt.hlines(y = barrelR2,  xmin = endcapZ1 - offset, xmax= +endcapZ2 + offset, color=color,  alpha=alpha, linewidth=linewidth)
    plt.hlines(y = barrelR2,  xmin = -endcapZ2 - offset, xmax= -endcapZ1 + offset, color=color,  alpha=alpha, linewidth=linewidth)
    plt.hlines(y = barrelR3,  xmin = endcapZ2 - offset, xmax= endcapZ3 + offset, color=color,  alpha=alpha, linewidth=linewidth)
    plt.hlines(y = barrelR3,  xmin = -endcapZ3 - offset, xmax= -endcapZ2 + offset, color=color,  alpha=alpha, linewidth=linewidth)

    plt.vlines(x = endcapZ1,  ymin = barrelR2-roffset, ymax= barrelR1, color=color, alpha=alpha, linewidth=linewidth)
    plt.vlines(x = -endcapZ1, ymin = barrelR2-roffset, ymax= barrelR1, color=color, alpha=alpha, linewidth=linewidth)
    plt.vlines(x = endcapZ2,  ymin = barrelR3-roffset, ymax= barrelR2, color=color,  alpha=alpha, linewidth=linewidth)
    plt.vlines(x = -endcapZ2, ymin = barrelR3-roffset, ymax= barrelR2, color=color, alpha=alpha, linewidth=linewidth)
    
    
def plotATLASCaloLayers(input_layer_csv:str):
    
    layerID_to_name = {0: "PreSamplerB", 1: "EMB1", 2: "EMB2", 3: "EMB3", 4: "PreSamplerE", 5: "EME1", 6: "EME2", 7: "EME3", 8: "HEC0", 9: "HEC1", 10: "HEC2", 11: "HEC3", 12: "TileBar0", 13: "TileBar1", 14: "TileBar2", 15: "TileGap1", 16: "TileGap2", 17: "TileGap3", 18: "TileExt0", 19: "TileExt1", 20: "TileExt2", 21: "FCAL0", 22: "FCAL1", 23: "FCAL2"}

    df = pd.read_csv(input_layer_csv, header=None, index_col=False, names = ['layerID', 'isBarrel', 'eta', 'rEnt', 'rMid', 'rExt', 'zEnt', 'zMid', 'zExt'])
    markerList = ['o', 's', 'p']
    for sample in range(0,23):
        rValues = np.array(list(df[df['layerID'] == sample]['rEnt']))
        zValues = np.array(list(df[df['layerID'] == sample]['zEnt']))
        
        # change marker after 10 samples (max number of matplotlib colors)
        marker = markerList[sample // 10]        
        plt.plot(zValues/1000, rValues/1000, marker, markersize=2, label=layerID_to_name[sample])



def plotEtaGuidingLines(annotate_labels = False):
    start_point_list = [[1.07915876e-16, 1.76240000e+00],
                        [-1.81705442, 1.54616455],
                        [-4.14229855, 1.14211689],
                        [-4.16136164, 0.41539365],
                        [1.81705442, 1.54616455],
                        [4.14229855, 1.14211689],
                        [4.16136164, 0.41539365]]
    position_list = [[7.02947263e-17, 1.14800000e+00], 
                     [-1.34913097, 1.148],
                     [-3.55, 0.978808],
                     [-3.55, 0.35436657],
                     [1.34913097, 1.148],
                     [3.55, 0.978808],
                     [3.55, 0.35436657]]
    for idx, eta in enumerate([0, -1, -2, -3, 1, 2, 3]):

        plt.plot([5*start_point_list[idx][0], position_list[idx][0]], [5*start_point_list[idx][1], position_list[idx][1]], linestyle="dotted", color='purple', alpha=0.4)
            
        if annotate_labels:
            if eta >= 0:
                annotatePointZ =  1.1*start_point_list[idx][0] + 0.3
                annotatePointR =  1.1*start_point_list[idx][1] - 0.06
                plt.annotate(fr'$\eta$={eta}', xy=(annotatePointZ, annotatePointR))
 

def get_track_data_of_event(df, event_id, track_id) -> pd.DataFrame:
    """
    Retrieves the track data of a specific event and track ID from a DataFrame.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing the track data.
    - event_id (int): The ID of the event.
    - track_id (int): The ID of the track.

    Returns:
    - pd.DataFrame: The track data of the specified event and track ID.
    """
    df_event = df[df.event == event_id]
    return df_event[df_event.track == track_id]
                
def main(args):
              
    # Read the data
    data = pd.read_csv(args.input)

    # Convert to meters
    data.x = data.x / 1000
    data.y = data.y / 1000
    data.z = data.z / 1000
    data.r = data.r / 1000
        
    # Get list of events
    event_list = data.event.unique()

    # Loop over the events
    for event_id in event_list:
        # Create new figure
        plt.figure()
        # Loop over the tracks
        for track_id in data.track.unique():
            track_data = get_track_data_of_event(data, event_id, track_id)
            # Plot the transport in the z-r plane
            plt.plot(track_data.z, track_data.r, 'o--', markersize=3, color='black', zorder=10)
            
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
            plt.xlim([-6.5 ,6.5])
            
        plt.xlabel('Z [m]')
        plt.ylabel('R [m]')
        
        plt.text(0.01, 0.95, fr"{args.label}", transform=plt.gca().transAxes)

        # Save the plot
        if len(event_list) == 1: 
            output_path = args.output
        else:
            # Add a _event_i prior to file extension if more than one event
            output_path = args.output.split(".")[0] + f"_event_{event_id}." + args.output.split(".")[1]
            
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Plot particle transport data.')
    parser.add_argument('--input', type=str, required=True, help='Input file path')
    parser.add_argument('--output', type=str, required=True, help='Output file path')
    parser.add_argument('--calo_layer_csv', type=str, required=True, help='CSV file containing the ATLAS calorimeter layer data')
    parser.add_argument('--track_zoom', action='store_true', help='Zoom in to the track region')
    parser.add_argument('--label', type=str, required=False, help='Label for the plot')
    

    main(parser.parse_args())
    
    
        
        

import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

def plotATLASCaloBoundary():
    barrelR1 = 1148.0 / 1000
    barrelR2 = 120.0 / 1000
    barrelR3 = 41.0 / 1000
    endcapZ1 = 3550.0 / 1000
    endcapZ2 = 4587.0 / 1000
    endcapZ3 = 6783 / 1000
    offset = 0.0245 # offset to close lines better
    roffset = 0.005
    plt.hlines(y = barrelR1,  xmin = -endcapZ1 - offset, xmax= +endcapZ1 + offset, colors='black')
    plt.hlines(y = barrelR2,  xmin = endcapZ1 - offset, xmax= +endcapZ2 + offset, colors='black')
    plt.hlines(y = barrelR2,  xmin = -endcapZ2 - offset, xmax= -endcapZ1 + offset, colors='black')
    plt.hlines(y = barrelR3,  xmin = endcapZ2 - offset, xmax= endcapZ3 + offset, colors='black')
    plt.hlines(y = barrelR3,  xmin = -endcapZ3 - offset, xmax= -endcapZ2 + offset, colors='black')

    plt.vlines(x = endcapZ1,  ymin = barrelR2-roffset, ymax= barrelR1, colors='black')
    plt.vlines(x = -endcapZ1, ymin = barrelR2-roffset, ymax= barrelR1, colors='black')
    plt.vlines(x = endcapZ2,  ymin = barrelR3-roffset, ymax= barrelR2, colors='black')
    plt.vlines(x = -endcapZ2, ymin = barrelR3-roffset, ymax= barrelR2, colors='black')


def plotEtaGuidingLines():
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
            plt.plot(track_data.z, track_data.r, 'o--', markersize=3)
        
        # Plot the ATLAS calorimeter boundary
        plotATLASCaloBoundary()
        # Plot the eta guiding lins
        plotEtaGuidingLines()

        plt.xlabel('Z [m]')
        plt.ylabel('R [m]')

        plt.ylim([0, 4])
        plt.xlim([-6.5 ,6.5])

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

    main(parser.parse_args())
    
    
        
        

#!/usr/bin/env python
"""
The program plots a number of timeseries together for comparison
"""
# Import Python modules
from __future__ import division, print_function
import os
import sys
import argparse
from math import radians, cos, sin, asin, sqrt
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disables use of Tk/X11
import numpy as np
import matplotlib.pyplot as plt

# Import seismtools functions
from file_utilities import read_files

def plot_comparison(args, filenames, stations,
                    output_file, plot_title=None):
    """
    Plotting a comparison of multiple timeseries
    """
    all_styles = ['k', 'r', 'b', 'm', 'g', 'c', 'y', 'brown',
                  'gold', 'blueviolet', 'grey', 'pink']
    orientation = ['N/S', 'E/W', 'Up/Down']

    # Check number of input timeseries
    if len(stations) > len(all_styles):
        print("[ERROR]: Too many timeseries to plot!")
        sys.exit(-1)

    delta_ts = [station[0].dt for station in stations]
    xtmin = args.xmin
    xtmax = args.xmax
    min_is = [int(xtmin/delta_t) for delta_t in delta_ts]
    max_is = [int(xtmax/delta_t) for delta_t in delta_ts]

    # Create plot
    f, axarr = plt.subplots(nrows=3, ncols=3, figsize=(14, 9))

    # For each component: N/S, E/W, U/D
    for i in range(0, 3):
        
        title_acc = "Acceleration - %s" % (orientation[i])
        title_vel = "Velocity - %s" % (orientation[i])
        title_dis = "Displacement - %s" % (orientation[i])

        signals = [station[i] for station in stations]
        samples = [signal.samples for signal in signals]
        displs = [signal.displ for signal in signals]
        vels = [signal.velo for signal in signals]
        accs = [signal.accel for signal in signals]
        
        # cutting signal by bounds
        c_displs = [dis[min_i:max_i] for dis, min_i, max_i in zip(displs,
                                                                  min_is,
                                                                  max_is)]
        c_vels = [vel[min_i:max_i] for vel, min_i, max_i in zip(vels,
                                                                min_is,
                                                                max_is)]
        c_accs = [acc[min_i:max_i] for acc, min_i, max_i in zip(accs,
                                                                min_is,
                                                                max_is)]
        times = [np.arange(xtmin,
                           min(xtmax, (delta_t * sample)),
                           delta_t) for delta_t, sample in zip(delta_ts,
                                                               samples)]

        axarr[i][0] = plt.subplot2grid((3, 3), (i, 0))
        axarr[i][0].set_title(title_dis)
        axarr[i][0].grid(True)
        styles = all_styles[0:len(times)]
        for timeseries, c_dis, style in zip(times, c_displs, styles):
            axarr[i][0].plot(timeseries, c_dis, style)
        plt.xlim(xtmin, xtmax)

        axarr[i][1] = plt.subplot2grid((3, 3), (i, 1))
        axarr[i][1].set_title(title_vel)
        axarr[i][1].grid(True)
        styles = all_styles[0:len(times)]
        for timeseries, c_vel, style in zip(times, c_vels, styles):
            axarr[i][1].plot(timeseries, c_vel, style)
        plt.xlim(xtmin, xtmax)

        axarr[i][2] = plt.subplot2grid((3, 3), (i, 2))
        axarr[i][2].set_title(title_acc)
        axarr[i][2].grid(True)
        styles = all_styles[0:len(times)]
        for timeseries, c_acc, style in zip(times, c_accs, styles):
            axarr[i][2].plot(timeseries, c_acc, style)
        # Add labels to first plot
        if i == 0:
            plt.legend(filenames, prop={'size':6})
        plt.xlim(xtmin, xtmax)

    # Make nice plots with tight_layout
    f.tight_layout()

    # Add overall title if provided
    if plot_title is not None:
        st = plt.suptitle(plot_title, fontsize=16)
        # shift subplots down:
        #st.set_y(0.95)
        f.subplots_adjust(top=0.92)

    # All done, save plot
    if output_file.lower().endswith(".png"):
        fmt='png'
    elif output_file.lower().endswith(".pdf"):
        fmt='pdf'
    else:
        print("ERROR: Unknown format!")
        sys.exit(-1)

    plt.savefig(output_file, format=fmt,
                transparent=False, dpi=300)

def calculate_distance(epicenter, st_loc):
    """
    Calculates the distance between two pairs of lat, long coordinates
    using the Haversine formula
    """
    lat1 = radians(abs(epicenter[0]))
    lon1 = radians(abs(epicenter[1]))
    lat2 = radians(abs(st_loc[0]))
    lon2 = radians(abs(st_loc[1]))

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))

    # Radius of earth in kilometers
    r = 6371

    return c * r

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Creates comparison plots of "
                                     " a number of timeseries files.")
    parser.add_argument("-o", "--output", dest="outfile", required=True,
                        help="output png file")
    parser.add_argument("--epicenter-lat", dest="epicenter_lat", type=float,
                        help="earthquake epicenter latitude")
    parser.add_argument("--epicenter-lon", dest="epicenter_lon", type=float,
                        help="earthquake epicenter longitude")
    parser.add_argument("--st-lat", "--station-latitude", dest="st_lat",
                        type=float, help="station latitude")
    parser.add_argument("--st-lon", "--station-longitude", dest="st_lon",
                        type=float, help="station longitude")
    parser.add_argument("-s", "--station-name", "--station", dest="station",
                        help="station name")
    parser.add_argument("--station-list", dest="station_list",
                        help="station list with latitude and longitude")
    parser.add_argument("--xmin", dest="xmin", type=float,
                        help="xmin to plot")
    parser.add_argument("--xmax", dest="xmax", type=float,
                        help="xmax to plot")
    parser.add_argument('input_files', nargs='*')
    args = parser.parse_args()

    if args.st_lat is not None and args.st_lon is not None:
        args.st_loc = [args.st_lat, args.st_lon]
    else:
        args.st_loc = None
    if args.epicenter_lat is not None and args.epicenter_lon is not None:
        args.epicenter = [args.epicenter_lat, args.epicenter_lon]
    else:
        args.epicenter = None
    if args.xmin is None:
        args.xmin = 0.0
    if args.xmax is None:
        args.xmax = 30.0

    return args

def compare_timeseries_main():
    """
    Main function for compare_timeseries
    """
    # Parse command-line options
    args = parse_arguments()
    # Copy inputs
    output_file = args.outfile
    filenames = args.input_files

    # Set plot title
    plot_title = None
    if args.station is not None:
        plot_title = "%s" % (args.station)
    # Set title if station name provided and epicenter are provided
    if args.station is not None and args.epicenter is not None:
        # Calculate distance if locations are provided
        if args.st_loc is None and args.station_list is not None:
            # Find station coordinates from station list
            st_file = open(args.station_list, 'r')
            for line in st_file:
                line = line.strip()
                if not line:
                    # skip blank lines
                    continue
                if line.startswith("#") or line.startswith("%"):
                    # Skip comments
                    continue
                pieces = line.split()
                if len(pieces) < 3:
                    # Skip line with insufficient tokens
                    continue
                if pieces[2].lower() != args.station.lower():
                    # Not a match
                    continue
                # Match!
                args.st_loc = [float(pieces[1]), float(pieces[0])]
                break
            # All done processing station file
            st_file.close()

        if args.st_loc is not None:
            # Calculate distance here
            distance = calculate_distance(args.epicenter, args.st_loc)
            # Set plot title
            plot_title = "%s, Dist: ~%dkm" % (args.station,
                                              distance)

    # Read data
    _, stations = read_files(None, filenames)
    filenames = [os.path.basename(filename) for filename in filenames]

    # Create plot
    plot_comparison(args, filenames, stations,
                    output_file, plot_title=plot_title)

# ============================ MAIN ==============================
if __name__ == "__main__":
    compare_timeseries_main()
# end of main program

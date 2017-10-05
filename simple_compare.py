#!/usr/bin/env python
"""
# ==============================================================================
# The program is to plot a simple comparison among timeseries files
# ==============================================================================
"""
from __future__ import division, print_function
import os
import sys
import argparse
from math import radians, cos, sin, asin, sqrt

from ptools import read_file
from compare_signals import simple_plot, set_parameter

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
    parser.add_argument("-p", "--parameters", dest="params", required=True,
                        help="filter and plot parameters")
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

    return args

def simple_compare_main():
    """
    Main function for simple_compare
    """
    # Parse command-line options
    args = parse_arguments()
    # Copy inputs
    output_file = args.outfile
    filenames = args.input_files
    params = args.params.split()

    # Set all other parameters
    parameter = set_parameter(params)

    # Set plot title
    plot_title = None
    if args.station is not None:
        plot_title = "%s, Freq: %1.1f-%1.1fHz" % (args.station,
                                                  parameter[5],
                                                  parameter[6])

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
            plot_title = "%s, Dist: ~%dkm, Freq: %1.1f-%1.1fHz" % (args.station,
                                                                   distance,
                                                                   parameter[5],
                                                                   parameter[6])

    # Read data
    stations = [read_file(filename) for filename in filenames]
    filenames = [os.path.basename(filename) for filename in filenames]

    # Create plot
    simple_plot(parameter, filenames, stations,
                output_file=output_file, plot_title=plot_title)

# ============================ MAIN ==============================
if __name__ == "__main__":
    simple_compare_main()
# end of main program

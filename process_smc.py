#!/usr/bin/env python
"""
# ==========================================================================
# The program is to get filename/directory in various ways, and then processes
# found V1/V2 files --> generate column acceleration .txt files and .her files.
# ==========================================================================
"""
from __future__ import print_function
import os
import sys
from smc import load_smc_v1, load_smc_v2, print_smc, print_her, print_bbp
from seism import seism_station

def get_parameters():
    """
    This function gets the list of files to process and creates a list.
    Additionally, it asks the user to enter a destination directory where
    output files will be written.
    """
    file_list = []
    output_format = ''
    destination = ''

    # Check if user provided any parameters
    if len(sys.argv) == 4:
        output_format = sys.argv[1]
        destination = sys.argv[2]
        file_list = sys.argv[3:]
    elif len(sys.argv) == 1:
        # Get the output format the user wants
        while output_format != 'bbp' and output_format != 'her':
            output_format = raw_input('== Enter output format (bbp/her): ')
            output_format = output_format.lower()

        while not file_list:
            # ask user if filename is not provided in command-line
            file_list = raw_input('== Enter the file / directory name: ')
            file_list = file_list.split()

        while not destination:
            destination = raw_input('== Enter name of the directory to store outputs: ')
    else:
        print("Usage: %s output_format destination_directory input_file/directory" %
              (os.path.basename(sys.argv[0])))
        sys.exit(-1)
    # check existence of target directory
    if not os.path.exists(destination):
        os.makedirs(destination)
    else:
        # clear content of unprocessed and warning files if they exist
        if os.path.exists(os.path.join(destination, 'unprocessed.txt')):
            clear(os.path.join(destination, 'unprocessed.txt'))
        if os.path.exists(os.path.join(destination, 'warning.txt')):
            clear(os.path.join(destination, 'warning.txt'))

    # All done!
    return destination, file_list, output_format
# end of get_filename

def read_list(destination, file_list, output_format):
    """
    The function is to read a list of files/directory and check their
    types to call corresponding functions.
    """

    for cur_file in file_list:
        # if is a directory; read all files\directories in the directory
        # and append them to file_list
        if os.path.isdir(cur_file):
            for item in os.listdir(cur_file):
                filename = os.path.join(cur_file, item)
                if not filename in file_list:
                    file_list.append(filename)

        # if is an non-empty file
        elif os.path.isfile(cur_file) and os.stat(cur_file).st_size != 0:
            # if the file is V1/raw data file: generate text file
            # for acceleration, and .her file
            if (cur_file.upper().endswith(".V1")
                or cur_file.upper().endswith(".RAW")):

                # Read V1 file
                station = load_smc_v1(cur_file)

                # Check we have valid results
                if station and station.list:
                    if station.process_v1():
                        print_smc(destination, station)
                        if output_format == 'her':
                            print_her(destination, station)
                        elif output_format == 'bbp':
                            print_bbp(destination, station)
                        else:
                            print("Error: Unknown output format: %s!" %
                                  (output_format))
                        check_station(destination, station)
                    else:
                        print_message(destination, cur_file, 'unprocessed')

            # if the file is V2/processed data file; generate text file
            # for acceleration, and .her file
            elif cur_file.upper().endswith(".V2"):

                # Read V2 file
                station = load_smc_v2(cur_file)

                # Check we have valid results
                if station and station.list:
                    if station.process_v2():
                        print_smc(destination, station)
                        if output_format == 'her':
                            print_her(destination, station)
                        elif output_format == 'bbp':
                            print_bbp(destination, station)
                        else:
                            print("Error: Unknown output format: %s!" %
                                  (output_format))
                        check_station(destination, station)
                    else:
                        print_message(destination, cur_file, 'unprocessed')
            else:
                # Maybe we have a filelist, try to process it...
                try:
                    input_file = open(cur_file)
                except IOError as err:
                    print(err)
                    continue
                lines = input_file.read().split()
                input_file.close()
                if not '#filelist' in lines[0]:
                    # unrecognized file type, skip it!
                    continue
                else:
                    for item in lines[1:]:
                        if not item in file_list:
                            file_list.append(item)
        else:
            print("[ERROR]: no such file or directory: %s" % (cur_file))
# end of read_list

def print_message(destination, message, ftype):
    """
    The function is to generate a files containing warning/unprocessed
    messages for input files.
    """
    out_file = open(os.path.join(destination, '%s.txt' % (ftype)), 'a')
    out_file.write(message + "\n")
    out_file.close()
# end of print_message

def check_station(destination, station):
    """
    The function is to check the station name of each record,
    if it's in the location should be discarded, print warning.
    """
    # check instance
    if not isinstance(station, seism_station):
        return
    if not station.list:
        return

    discard = {'dam': 'Dam', 'Fire Sta': 'Fire Station',
               'Acosta Res': 'Acosta Res', 'Bldg': 'Building',
               'Br': 'Interchange Bridge'}
    name = station.name
    for key in discard:
        if key in name:
            filename = station.network + station.id + '.' + station.type
            msg = filename + " was processed, but it's from " + discard[key]
            print_message(destination, msg, 'warning')
            break
# end of check_station

def clear(filename):
    """
    This function clears the content of a file if it exists.
    """
    try:
        open(filename, 'w').close()
    except IOError as err:
        print(err)
# end of clear

def process_main():
    """
    Main function for the process_smc program
    """
    # Main function
    destination, file_list, output_format = get_parameters()
    read_list(destination, file_list, output_format)

# ============================ MAIN ==============================
if __name__ == "__main__":
    process_main()
# end of main program

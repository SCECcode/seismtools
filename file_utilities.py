#!/usr/bin/env python
"""
The program contains several input/output utility
functions used by other modules.
"""
# Import Python modules
from __future__ import division, print_function, absolute_import
import os
import sys
import numpy as np

# Import seismtools needed classes
from seism import seism_psignal

def reverse_up_down(station):
    """
    reverse up down component
    """
    # station has 3 components [ns, ew, ud]
    # only need to flip the 3rd one
    station[2].accel *= -1
    station[2].velo *= -1
    station[2].displ *= -1

    return station
# end of reverse_up_down

def scale_from_m_to_cm(station):
    # scales timeseries from meters to centimeters
    for i in range(0, len(station)):
        station[i].accel *= 100
        station[i].velo *= 100
        station[i].displ *= 100

    return station
# end of scale_from_m_to_cm

def read_files(obs_file, input_files):
    """
    Reads all input files
    """
    # read obs data
    obs_data = None
    if obs_file is not None:
        obs_data = read_file(obs_file)
        # Make sure we got it
        if not obs_data:
            print("[ERROR]: Reading obs file: %s!" % (obs_file))
            sys.exit(-1)
        # Fix units if needed
        if obs_file.lower().endswith(".bbp"):
            units = read_unit_bbp(obs_file)
            # If in meters, scale to cm
            if units == "m":
                obs_data = scale_from_m_to_cm(obs_data)
        else:
            # In Hercule files, for observation files data is already in
            # cm, so nothing to do here!
            #obs_data = reverse_up_down(obs_data)
            pass

    # reads signals
    stations = []
    for input_file in input_files:
        station = read_file(input_file)
        # Make sure we got it
        if not station:
            print("[ERROR]: Reading input file: %s!" % (input_file))
            sys.exit(-1)
        # Fix units if needed
        if input_file.lower().endswith(".bbp"):
            units = read_unit_bbp(input_file)
            # If in meters, scale to cm
            if units == "m":
                station = scale_from_m_to_cm(station)
        else:
            # Hercule file, need to scale and flip up/down component
            station = scale_from_m_to_cm(station)
            station = reverse_up_down(station)
        stations.append(station)

    # all done
    return obs_data, stations

def read_filelist(filelist):
    """
    This function reads the filelist provided by the user
    """
    station_list = []
    coor_x = []
    coor_y = []

    try:
        input_file = open(filelist, 'r')
    except IOError:
        print("[ERROR]: error loading filelist.")
        sys.exit(-1)

    for line in input_file:
        if not '#' in line:
            line = line.split()
            # Get station name and make substitution
            station_name = line[0]
            station_name = station_name.replace(".", "_")

            if len(line) == 1:
                # not containing coordinates
                station_list.append(station_name)
                coor_x.append(0.0)
                coor_y.append(0.0)
            elif len(line) == 3:
                # containing coordinates
                station_list.append(station_name)
                try:
                    coor_x.append(float(line[1]))
                    coor_y.append(float(line[2]))
                except ValueError:
                    coor_x.append(0.0)
                    coor_y.append(0.0)

    # Close the input file
    input_file.close()

    return station_list, coor_x, coor_y
# end of read_filelist

# ================================ READING ================================
def read_file(filename):
    """
    This function reads a timeseries file(s) in either bbp
    format (ends with .bbp) or hercules (ends otherwise)
    """
    if filename.lower().endswith(".bbp"):
        # Filename in bbp format
        return read_file_bbp(filename)
    # Otherwise use hercules format
    return read_file_her(filename)
# end of read_file

def read_file_bbp2(filename):
    """
    This function reads a bbp file and returns the timeseries in the
    format time, n/s, e/w, u/d tuple
    """
    time = np.array([])
    ns_comp = np.array([])
    ew_comp = np.array([])
    ud_comp = np.array([])

    try:
        input_file = open(filename, 'r')
        for line in input_file:
            line = line.strip()
            if line.startswith('#') or line.startswith('%'):
                # Skip comments
                continue
            # Trim in-line comments
            if line.find('#') > 0:
                line = line[:line.find('#')]
            if line.find('%') > 0:
                line = line[:line.find('%')]
            # Make them float
            pieces = line.split()
            pieces = [float(piece) for piece in pieces]
            time = np.append(time, pieces[0])
            ns_comp = np.append(ns_comp, pieces[1])
            ew_comp = np.append(ew_comp, pieces[2])
            ud_comp = np.append(ud_comp, pieces[3])
    except IOError:
        print("[ERROR]: error reading bbp file: %s" % (filename))
        sys.exit(1)

    # All done!
    return time, ns_comp, ew_comp, ud_comp
# end of read_file_bbp2

def read_file_bbp(filename):
    """
    This function reads timeseries data from a set of BBP files
    """
    # Get filenames for displacement, velocity and acceleration bbp files
    work_dir = os.path.dirname(filename)
    base_file = os.path.basename(filename)

    base_tokens = base_file.split('.')[0:-2]
    if not base_tokens:
        print("[ERROR]: Invalid BBP filename: %s" % (filename))
        sys.exit(1)
    dis_tokens = list(base_tokens)
    vel_tokens = list(base_tokens)
    acc_tokens = list(base_tokens)

    dis_tokens.append('dis')
    vel_tokens.append('vel')
    acc_tokens.append('acc')

    dis_tokens.append('bbp')
    vel_tokens.append('bbp')
    acc_tokens.append('bbp')

    dis_file = os.path.join(work_dir, '.'.join(dis_tokens))
    vel_file = os.path.join(work_dir, '.'.join(vel_tokens))
    acc_file = os.path.join(work_dir, '.'.join(acc_tokens))

    # Read 3 bbp files
    [time, dis_ns, dis_ew, dis_up] = read_file_bbp2(dis_file)
    [_, vel_ns, vel_ew, vel_up] = read_file_bbp2(vel_file)
    [_, acc_ns, acc_ew, acc_up] = read_file_bbp2(acc_file)

    samples = dis_ns.size
    delta_t = time[1]

    # samples, dt, data, acceleration, velocity, displacement
    psignal_ns = seism_psignal(samples, delta_t, np.c_[dis_ns, vel_ns, acc_ns],
                               'c', acc_ns, vel_ns, dis_ns)
    psignal_ew = seism_psignal(samples, delta_t, np.c_[dis_ew, vel_ew, acc_ew],
                               'c', acc_ew, vel_ew, dis_ew)
    psignal_up = seism_psignal(samples, delta_t, np.c_[dis_up, vel_up, acc_up],
                               'c', acc_up, vel_up, dis_up)

    station = [psignal_ns, psignal_ew, psignal_up]
    return station
# end of read_file_bbp

def read_file_her(filename):
    """
    The function is to read 10-column .her files.
    Return a list of psignals for each orientation.
    """
    time, dis_ns, dis_ew, dis_up = [np.array([], float) for _ in xrange(4)]
    vel_ns, vel_ew, vel_up = [np.array([], float) for _ in xrange(3)]
    acc_ns, acc_ew, acc_up = [np.array([], float) for _ in xrange(3)]

    try:
        (time, dis_ns, dis_ew, dis_up, vel_ns, vel_ew,
         vel_up, acc_ns, acc_ew, acc_up) = np.loadtxt(filename,
                                                      comments='#',
                                                      unpack=True)
    except IOError:
        print("[ERROR]: error loading her file.")
        return False

    samples = dis_ns.size
    delta_t = time[1]

    # samples, dt, data, acceleration, velocity, displacement
    psignal_ns = seism_psignal(samples, delta_t, np.c_[dis_ns, vel_ns, acc_ns],
                               'c', acc_ns, vel_ns, dis_ns)
    psignal_ew = seism_psignal(samples, delta_t, np.c_[dis_ew, vel_ew, acc_ew],
                               'c', acc_ew, vel_ew, dis_ew)
    psignal_up = seism_psignal(samples, delta_t, np.c_[dis_up, vel_up, acc_up],
                               'c', acc_up, vel_up, dis_up)

    station = [psignal_ns, psignal_ew, psignal_up]
    return station
# end of read_file_her

def read_unit_bbp(filename):
    """
    Get the units from the file's header
    Returns either "m" or "cm"
    """
    units = None

    try:
        input_file = open(filename, 'r')
        for line in input_file:
            if line.find("units=") > 0:
                units = line.split()[2]
                break
        input_file.close()
    except IOError:
        print("[ERROR]: No such file.")
        sys.exit(-1)

    # Make sure we got something
    if units is None:
        print("[ERROR]: Cannot find units in bbp file!")
        sys.exit(-1)

    # Figure out if we have meters or centimeters
    if units == "cm" or units == "cm/s" or units == "cm/s^2":
        return "cm"
    elif units == "m" or units == "m/s" or units == "m/s^2":
        return "m"

    # Invalid units in this file
    print("[ERROR]: Cannot parse units in bbp file!")
    sys.exit(-1)
# end of read_unit_bbp

def read_stamp(filename):
    """
    Get the time stamp from file's header
    """
    if filename.endswith(".bbp"):
        # File in bbp format
        return read_stamp_bbp(filename)
    # Otherwise use hercules format
    return read_stamp_her(filename)
# end of read_stamp

def read_stamp_bbp(filename):
    """
    Get the time stamp from the bbp file's header
    """
    try:
        input_file = open(filename, 'r')
        for line in input_file:
            if line.find("time=") > 0:
                stamp = line.split()[2].split(',')[-1].split(':')
                break
        input_file.close()
    except IOError:
        print("[ERROR]: No such file.")
        return []

    # Converting time stamps to floats
    stamp = [float(i) for i in stamp]
    return stamp
# end of read_stamp_bbp

def read_stamp_her(filename):
    """
    Get the time stamp from the her file's header
    """
    try:
        with open(filename) as input_file:
            try:
                header = input_file.readline().split()
                stamp = header[4].split(',')[-1].split(':')
                input_file.close()
            except IndexError:
                print("[ERROR]: missing time stamp.")
                return []
    except IOError:
        print("[ERROR]: No such file.")
        return []

    # converting time stamps to floats
    for i in range(0, len(stamp)):
        stamp[i] = float(stamp[i])
    return stamp
# end of read_stamp_her

# ================================ WRITING ==================================
def print_her(filename, station):
    # filename = 'processed-' + filename.split('/')[-1]
    try:
        out_f = open(filename, 'w')
    except IOError as e:
        print(e)
    dis_ns = station[0].displ.tolist()
    vel_ns = station[0].velo.tolist()
    acc_ns = station[0].accel.tolist()
    dis_ew = station[1].displ.tolist()
    vel_ew = station[1].velo.tolist()
    acc_ew = station[1].accel.tolist()
    dis_up = station[2].displ.tolist()
    vel_up = station[2].velo.tolist()
    acc_up = station[2].accel.tolist()

    # get a list of time incremented by dt
    time = [0.000]
    samples = station[0].samples
    dt = station[0].dt
    tmp = samples

    while tmp > 1:
        time.append(time[len(time)-1] + dt)
        tmp -= 1

    out_f.write('# missing header \n')

    descriptor = '{:>12}' + '  {:>12}'*9 + '\n'
    out_f.write(descriptor.format("# time",
                                  "dis_ns", "dis_ew", "dis_up",
                                  "vel_ns", "vel_ew", "vel_up",
                                  "acc_ns", "acc_ew", "acc_up")) # header

    descriptor = '{:>12.3f}' + '  {:>12.7f}'*9 + '\n'
    for c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 in zip(time,
                                                      dis_ns, dis_ew, dis_up,
                                                      vel_ns, vel_ew, vel_up,
                                                      acc_ns, acc_ew, acc_up):
        out_f.write(descriptor.format(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9))
    out_f.close()
# end of print_her

def print_bbp(input_file, output_file, station):
    """
    This function generates processed .bbp files for
    each of velocity/acceleration/displacement
    and copies the header of the input bbp file
    """
    output_dir = os.path.dirname(output_file)
    output_basename = os.path.basename(output_file)

    # Prepare data for output
    acc_ns = station[0].accel.tolist()
    vel_ns = station[0].velo.tolist()
    dis_ns = station[0].displ.tolist()
    acc_ew = station[1].accel.tolist()
    vel_ew = station[1].velo.tolist()
    dis_ew = station[1].displ.tolist()
    acc_up = station[2].accel.tolist()
    vel_up = station[2].velo.tolist()
    dis_up = station[2].displ.tolist()

    # Start with time = 0.0
    time = [0.000]
    samples = station[0].samples
    while samples > 1:
        time.append(time[len(time)-1] + station[0].dt)
        samples -= 1

    # Prepare to output
    out_data = [['dis', dis_ns, dis_ew, dis_up, 'displacement', 'cm'],
                ['vel', vel_ns, vel_ew, vel_up, 'velocity', 'cm/s'],
                ['acc', acc_ns, acc_ew, acc_up, 'acceleration', 'cm/s^2']]

    for data in out_data:
        if not output_basename.endswith('.bbp'):
            # Remove extension
            bbp_output_basename = os.path.splitext(output_basename)[0]
            bbp_output_filename = os.path.join(output_dir,
                                               "%s.%s.bbp" %
                                               (bbp_output_basename,
                                                data[0]))
            output_header = ["# Station: NoName",
                             "#    time= 00/00/00,00:00:00.00 UTC",
                             "#     lon= 0.00",
                             "#     lat= 0.00",
                             "#   units= %s" % (data[5]),
                             "#",
                             "# Data fields are TAB-separated",
                             "# Column 1: Time (s)",
                             "# Column 2: N/S component ground "
                             "%s (+ is 000)" % (data[4]),
                             "# Column 3: E/W component ground "
                             "%s (+ is 090)" % (data[4]),
                             "# Column 4: U/D component ground "
                             "%s (+ is upward)" % (data[4]),
                             "#"]
        else:
            # Read header of input file
            input_dirname = os.path.dirname(input_file)
            input_basename = os.path.basename(input_file)
            pieces = input_basename.split('.')
            pieces = pieces[0:-2]
            bbp_input_file = os.path.join(input_dirname,
                                          "%s.%s.bbp" %
                                          ('.'.join(pieces),
                                           data[0]))
            input_header = []
            in_fp = open(bbp_input_file, 'r')
            for line in in_fp:
                line = line.strip()
                if line.startswith("#"):
                    input_header.append(line)
            in_fp.close()

            # Compose new header
            output_header = []
            for item in input_header:
                if item.find("units=") > 0:
                    output_header.append("#   units= %s" % (data[5]))
                else:
                    output_header.append(item)

            pieces = output_basename.split('.')
            pieces = pieces[0:-2]
            bbp_output_filename = os.path.join(output_dir,
                                               "%s.%s.bbp" %
                                               ('.'.join(pieces),
                                                data[0]))
        # Write output file
        try:
            out_fp = open(bbp_output_filename, 'w')
        except IOError as e:
            print(e)
            continue

        # Write header
        for item in output_header:
            out_fp.write("%s\n" % (item))

        # Write timeseries
        for val_time, val_ns, val_ew, val_ud in zip(time, data[1],
                                                    data[2], data[3]):
            out_fp.write("%5.7f   %5.9e   %5.9e    %5.9e\n" %
                         (val_time, val_ns, val_ew, val_ud))

        # All done, close file
        out_fp.close()
        print("*Writing file: %s " % (bbp_output_filename))
# end of print_bbp

#!/usr/bin/env python
"""
# ===================================================================================
# The program contains general functions what may be used by other programs.
# ===================================================================================
"""
from __future__ import division, print_function, absolute_import
import numpy as np
from seism import s_filter, seism_psignal
from stools import seism_cutting, seism_appendzeros

def synchronize_all_stations(obs_data, stations, stamp, eqtimestamp, leading):
    """
    synchronize the stating time and ending time of data arrays
    obs_data = recorded data (optional); stations = simulation signal(s)
    """
    # If we have a recorded data time stamp
    if stamp is not None and obs_data is not None:
        start = stamp[0]*3600 + stamp[1]*60 + stamp[2]
        eq_time = eqtimestamp[0]*3600 + eqtimestamp[1]*60 + eqtimestamp[2]
        sim_start = eq_time - leading

        for i in range(0, 3):
            # synchronize the start time
            if start < sim_start:
                # data time < sim time < earthquake time; cutting data array
                obs_data[i] = seism_cutting('front', (sim_start - start),
                                            20, obs_data[i], False)
            elif start > eq_time:
                # sim time < earthquake time < data time; adding zeros in front
                obs_data[i] = seism_appendzeros('front', (start - eq_time),
                                                20, obs_data[i])
                for station in stations:
                    station[i] = seism_cutting('front', (eq_time - sim_start),
                                               20, station[i], False)
            else:
                # sim time < data time < earthquake time; adding zeros
                obs_data[i] = seism_appendzeros('front', (start - sim_start),
                                                20, obs_data[i])

    # synchronize the ending time
    if obs_data is not None:
        obs_dt = obs_data[0].dt
        obs_samples = obs_data[0].samples
        obs_time = obs_dt * obs_samples
    else:
        obs_time = None

    # Find target timeseries duration
    target_time = None
    if obs_time is not None:
        target_time = obs_time
    for station in stations:
        station_dt = station[0].dt
        station_samples = station[0].samples
        station_time = station_dt * station_samples
        if target_time is None:
            target_time = station_time
            continue
        target_time = min(target_time, station_time)

    # Work on obs_data
    if obs_data is not None:
        for i in range(0, 3):
            if obs_time > target_time:
                obs_data[i] = seism_cutting('end', (obs_time - target_time),
                                            20, obs_data[i], False)
        obs_samples = obs_data[0].samples
        obs_time = obs_dt * obs_samples

    # Work on simulated data
    for station in stations:
        for i in range(0, 3):
            sim_dt = station[i].dt
            sim_samples = station[i].samples
            sim_time = sim_dt * sim_samples
            if sim_time > target_time:
                station[i] = seism_cutting('end', (sim_time - target_time),
                                           20, station[i], False)

    # scale the data if they have one sample in difference after synchronizing
    total_samples = None
    if obs_data is not None:
        total_samples = obs_samples
    for station in stations:
        sim_samples = station[0].samples
        if total_samples is None:
            total_samples = sim_samples
            continue
        total_samples = max(sim_samples, total_samples)

    # For obs_data
    if obs_data is not None:
        for i in range(0, 3):
            if obs_data[i].samples == total_samples - 1:
                obs_data[i] = seism_appendzeros('end', obs_data[i].dt,
                                                20, obs_data[i])
    # For simulated data
    for station in stations:
        for i in range(0, 3):
            if station[i].samples == total_samples - 1:
                station[i] = seism_appendzeros('end', station[i].dt,
                                               20, station[i])

    return obs_data, stations
# end of synchronize_all_stations

def filter_data(psignal, fmin, fmax):
    """
    This function is used to filter with a bandpass filter between fmin/fmax
    """
    if not isinstance(psignal, seism_psignal):
        print("[ERROR]: found error filtering psignal.")
        return False
    delta_t = psignal.dt
    psignal.accel = s_filter(psignal.accel, delta_t, type='bandpass',
                             family='butter', fmin=fmin, fmax=fmax,
                             N=4, rp=0.1, rs=100)
    psignal.velo = s_filter(psignal.velo, delta_t, type='bandpass',
                            family='butter', fmin=fmin, fmax=fmax,
                            N=4, rp=0.1, rs=100)
    psignal.displ = s_filter(psignal.displ, delta_t, type='bandpass',
                             family='butter', fmin=fmin, fmax=fmax,
                             N=4, rp=0.1, rs=100)

    psignal.data = np.c_[psignal.displ, psignal.velo, psignal.accel]

    return psignal
# end of filter_data

def get_bands():
    """
    The function is to allow user specify sample rates.
    Without user input, sample rates are setting to default values.
    """
    freq_0 = 0.05
    freq_1 = 0.1
    freq_2 = 0.25
    freq_3 = 0.5
    freq_4 = 1
    freq_5 = 2
    freq_6 = 4
    bands = [freq_0, freq_1, freq_2, freq_3, freq_4, freq_5, freq_6]
    freqs = []
    flag = True

    while flag:
        flag = False
        freqs = raw_input('== Enter the sequence of '
                          'sample rates: ').replace(',', ' ').split()
        if not freqs:
            #setting to default values
            return bands

        if len(freqs) == 1:
            print("[ERROR]: invalid sample rates")
            flag = True
        else:
            bands = []
            for freq in freqs:
                try:
                    bands.append(float(freq))
                except ValueError:
                    print("[ERROR]: invalid sample rates")
                    flag = True
                    break

            for i in range(0, len(bands)-1):
                if bands[i] >= bands[i+1]:
                    print("[ERROR]: invalid sequence of sample rates")
                    flag = True
                    break
    return bands
# end of get_bands

def check_data(station):
    """
    Checks the data after rotation, process_dt, and synchronization
    to avoid encountering errors in gof_engine
    """
    for i in range(0, len(station)):
        signal = station[i]

        if signal.accel.size == 0:
            print("[ERROR]: Empty array after processing signals.")
            return False
        if signal.velo.size == 0:
            print("[ERROR]: Empty array after processing signals.")
            return False
        if signal.displ.size == 0:
            print("[ERROR]: Empty array after processing signals.")
            return False
        if np.isnan(np.sum(signal.accel)):
            print("[ERROR]: NaN data after processing signals.")
            return False
        if np.isnan(np.sum(signal.velo)):
            print("[ERROR]: NaN data after processing signals.")
            return False
        if np.isnan(np.sum(signal.displ)):
            print("[ERROR]: NaN data after processing signals.")
            return False
    return station
# end of check_data

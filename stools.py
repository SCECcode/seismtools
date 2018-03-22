#!/usr/bin/env python
"""
# =============================================================================
# The program contains general functions what may be used by other programs.
# Including: filter; integral; derivative; FAS;
# =============================================================================
"""
from __future__ import division, print_function
import sys
import numpy as np
import math
from scipy.signal import filtfilt, ellip, butter, kaiser
from scipy.integrate import cumtrapz

def integrate(data, dt):
    """
    compute derivative of a numpy array
    initial condition assumed 0
    result has same size as input
    """
    newdata = cumtrapz(data, dx=dt, initial=0) + data[0]*dt/2.0
    return newdata
    # data = np.cumsum(data*dt)
    # return data

def derivative(data, dt):
    """compute derivative of an numpy."""
    newdata = np.insert(data, 0, 0)
    newdata = np.diff(newdata)/dt
    return newdata

def s_filter(*args, **kwargs):
    """
    correct order for unlabeled arguments is data, dt;
    """
    data = np.array([], float)
    dt = 0.0
    fami = {'ellip': ellip, 'butter': butter}

    if len(args) == 2:
        data = args[0]
        dt = args[1]
    else:
        print("[ERROR]: filter missing data and dt.")
        return data

    if not isinstance(data, np.ndarray):
        print("[ERROR]: data input for filter is not an numpy array.")
        return data

    # default values
    N = 5
    rp = 0.1
    rs = 100
    Wn = 0.05/((1.0/dt)/2.0)
    fmin = 0.0
    fmax = 0.0
    a = np.array([], float)
    b = np.array([], float)

    if len(kwargs) > 0:
        if 'type' in kwargs:
            btype = kwargs['type']
        if 'N' in kwargs:
            N = kwargs['N']
        if 'rp' in kwargs:
            rp = kwargs['rp']
        if 'rs' in kwargs:
            rs = kwargs['rs']
        if 'Wn' in kwargs:
            Wn = kwargs['Wn']
        if 'fmin' in kwargs:
            fmin = kwargs['fmin']
            w_min = fmin/((1.0/dt)/2.0)
        if 'fmax' in kwargs:
            fmax = kwargs['fmax']
            w_max = fmax/((1.0/dt)/2.0)

        if fmin and fmax and btype == 'bandpass':
            Wn = [w_min, w_max]
        elif fmax and btype == 'lowpass':
            Wn = w_max
        elif fmin and btype == 'highpass':
            Wn = w_min

        # calling filter
        if kwargs['family'] == 'ellip':
            b, a = ellip(N=N, rp=rp, rs=rs, Wn=Wn, btype=btype, analog=False)
        elif kwargs['family'] == 'butter':
            b, a = butter(N=N, Wn=Wn, btype=btype, analog=False)

        data = filtfilt(b, a, data)
        return data
# end of s_filter

def smooth(data, factor):
    # factor = 3; c = 0.5, 0.25, 0.25
    # TODO: fix coefficients for factors other than 3
    c = 0.5/(factor-1)
    for i in range(1, data.size-1):
        data[i] = 0.5*data[i] + c*data[i-1] + c*data[i+1]
    return data

def FAS(data, dt, points, fmin, fmax, s_factor):
    afs = abs(np.fft.fft(data, points))*dt
    # freq = (1/signal.dt)*range(points)/points
    freq = (1/dt)*np.array(range(points))/points

    deltaf = (1/dt)/points

    inif = int(fmin/deltaf)
    endf = int(fmax/deltaf) + 1

    afs = afs[inif:endf]
    afs = smooth(afs, s_factor)
    freq = freq[inif:endf]
    return freq, afs

def get_points(samples):
    # points is the least base-2 number that is greater than max samples
    power = int(math.log(max(samples), 2)) + 1
    return 2**power
# end of get_points

def get_period(tmin, tmax):
    """ Return an array of period T """
    # tmin = 1/fmax
    # tmax = 1/fmin
    a = np.log10(tmin)
    b = np.log10(tmax)

    period = np.linspace(a, b, 20)
    period = np.power(10, period)
    return period

def max_osc_response(acc, dt, csi, period, ini_disp, ini_vel):
    signal_size = acc.size

    # initialize numpy arrays
    d = np.empty((signal_size))
    v = np.empty((signal_size))
    aa = np.empty((signal_size))

    d[0] = ini_disp
    v[0] = ini_vel

    w = 2*math.pi/period
    ww = w**2
    csicsi = csi**2
    dcsiw = 2*csi*w

    rcsi = math.sqrt(1-csicsi)
    csircs = csi/rcsi
    wd = w*rcsi
    ueskdt = -1/(ww*dt)
    dcsiew = 2*csi/w
    um2csi = (1-2*csicsi)/wd
    e = math.exp(-w*dt*csi)
    s = math.sin(wd*dt)
    c0 = math.cos(wd*dt)
    aa[0] = -ww*d[0]-dcsiw*v[0]

    ca = e*(csircs*s+c0)
    cb = e*s/wd
    cc = (e*((um2csi-csircs*dt)*s-(dcsiew+dt)*c0)+dcsiew)*ueskdt
    cd = (e*(-um2csi*s+dcsiew*c0)+dt-dcsiew)*ueskdt
    cap = -cb*ww
    cbp = e*(c0-csircs*s)
    ccp = (e*((w*dt/rcsi+csircs)*s+c0)-1)*ueskdt
    cdp = (1-ca)*ueskdt

    for i in range(1, signal_size):
        d[i] = ca*d[i-1]+cb*v[i-1]+cc*acc[i-1]+cd*acc[i]
        v[i] = cap*d[i-1]+cbp*v[i-1]+ccp*acc[i-1]+cdp*acc[i]
        aa[i] = -ww*d[i]-dcsiw*v[i]

    maxdisp = np.amax(np.absolute(d))
    maxvel = np.amax(np.absolute(v))
    maxacc = np.amax(np.absolute(aa))

    return maxdisp, maxvel, maxacc

def cal_acc_response(period, data, delta_ts):
    """
    # return the response for acceleration only
    """
    rsps = [[] for _ in delta_ts]
    for p in period:
        for rsp, timeseries, delta_t in zip(rsps, data, delta_ts):
            rsp.append(max_osc_response(timeseries, delta_t, 0.05,
                                        p, 0, 0)[-1])
    return rsps
# end of cal_acc_response

def taper(flag, m, samples):
    # m = samples for taper
    # samples = total samples
    window = kaiser(2*m+1, beta=14)

    if flag == 'front':
        # cut and replace the second half of window with 1s
        ones = np.ones(samples-m-1)
        window = window[0:(m+1)]
        window = np.concatenate([window, ones])

    elif flag == 'end':
        # cut and replace the first half of window with 1s
        ones = np.ones(samples-m-1)
        window = window[(m+1):]
        window = np.concatenate([ones, window])

    elif flag == 'all':
        ones = np.ones(samples-2*m-1)
        window = np.concatenate([window[0:(m+1)], ones, window[(m+1):]])

    # avoid concatenate error
    if window.size < samples:
        window = np.append(window, 1)

    if window.size != samples:
        print(window.size)
        print(samples)
        print("[ERROR]: taper and data do not have the same number of samples.")
        window = np.ones(samples)

    return window

def seism_appendzeros(flag, t_diff, m, signal):
    """adds zeros in the front and/or at the end of an numpy array
    apply taper before adding
    """
    # if not isinstance(signal, seism_psignal):
    #     return signal

    num = int(t_diff/signal.dt)
    zeros = np.zeros(num)

    if flag == 'front':
        # applying taper in the front
        if m != 0:
            window = taper('front', m, signal.samples)
            signal.accel = signal.accel*window
            signal.velo = signal.velo*window
            signal.displ = signal.displ*window

        # adding zeros in front of data
        signal.accel = np.append(zeros, signal.accel)
        signal.velo = np.append(zeros, signal.velo)
        signal.displ = np.append(zeros, signal.displ)

    elif flag == 'end':
        if m != 0:
            # applying taper in the front
            window = taper('end', m, signal.samples)
            signal.accel = signal.accel*window
            signal.velo = signal.velo*window
            signal.displ = signal.displ*window

        signal.accel = np.append(signal.accel, zeros)
        signal.velo = np.append(signal.velo, zeros)
        signal.displ = np.append(signal.displ, zeros)

    signal.samples += num
    return signal
# end of seism_appendzeros

def seism_cutting(flag, t_diff, m, signal, signal_flag):
    """cut data in the front or at the end of an numpy array
    apply taper after cutting
    """
    # if not isinstance(signal, seism_psignal):
    #     return signal

    num = int(t_diff/signal.dt)
    if num >= signal.samples:
        print("[ERROR]: fail to cut signal.")
        return signal

    if flag == 'front' and num != 0:
        # cutting signal
        if signal_flag == True:
            signal.data = signal.data[num:]
            signal.samples -= num
            window = taper('front', m, signal.samples)
            signal.data = signal.data*window
            return signal


        # cutting psignal
        signal.accel = signal.accel[num:]
        signal.velo = signal.velo[num:]
        signal.displ = signal.displ[num:]
        signal.samples -= num

        # applying taper at the front
        window = taper('front', m, signal.samples)
        signal.accel = signal.accel*window
        signal.velo = signal.velo*window
        signal.displ = signal.displ*window


    elif flag == 'end' and num != 0:
        num *= -1
         # cutting signal
        if signal_flag == True:
            signal.data = signal.data[:num]
            signal.samples += num
            window = taper('end', m, signal.samples)
            signal.data = signal.data*window
            return signal

        # cutting psignal
        signal.accel = signal.accel[:num]
        signal.velo = signal.velo[:num]
        signal.displ = signal.displ[:num]
        signal.samples += num

        # applying taper at the end
        window = taper('end', m, signal.samples)
        signal.accel = signal.accel*window
        signal.velo = signal.velo*window

    return signal
# end of seism_cutting

def scale_signal(signal, factor):
    """scale the data of given signal"""
    if not isinstance(signal.data, np.ndarray):
        print("[ERROR]: error in scale_signal; data is not an numpy array.")
        return signal
    signal.data = factor*signal.data
    return signal
# end of scale_signal

def correct_baseline(signal):
    if not isinstance(signal.data, np.ndarray):
        print("[ERROR]: error in correct_baseline; data is not an numpy array.")
        return signal

    # make average on first 10% of samples; minus average
    signal.data = signal.data - np.average(signal.data[0:int(signal.samples*0.1)])
    return signal
# end of correct_baseline

def polimod(x, y, n, m):
    """
    Polymod Fit polynomial to data - by J. Stewart 5/25/98

    polymod(x,y,n,m) finds the coefficients of a polynomial P(x) of
    degree n that fits the data. P(X(I))~=Y(I), in a least-squares sense
    but the first m terms of the polynomial are neglected in forming
    the fit (e.g. to use the squared and cubic terms in a third order
    polynomial, which has 4 coefficients, use n=3, m=1)

    The regression problem is formulated in matrix format as:

    y = G*m or

          3  2
    Y = [x  x  x  1] [p3
                      p2
                      p1
                      p0]

    where the vector p contains the coefficients to be found. For a
    3th order polynomial, matrix G would be:

    G = [x^3 x^2 X^1 1]

    where the number of rows in G equals the number of rows in x.
    """
    # Make sure the 2 vectors have the same size
    if len(x) != len(y):
        print("ERROR: X and Y vectors must be of same size!")
        sys.exit(-1)

    G = np.zeros((len(x), (n-m)))
    for i in range(0, len(x)):
        for j in range(0, (n-m)):
            G[i][j] = x[i] ** (j+1+m)

    # Transpose G
    GT = G.transpose()
    # Form solution see Section 2.2.2 of Geophysics 104 notes
    p = np.dot(np.dot(np.linalg.inv(np.dot(GT, G)), GT), y)
    # Polynomial coefficients are row vectors by convention
    return p

def baseline_function(acc, dt, gscale, ordern):
    """
    Integrates accelaration record and baseline corrects velocity and
    displacement time series using 5th order polynomial without the constant
    and linear term (only square, cubic, 4th and 5th order terms, so that
    the leading constants are applied to disp, vel, and acc)
    """
    # Use gscale to convert to cm/sec2
    acc = acc * gscale
    times = np.linspace(0, (len(acc) - 1) * dt, len(acc))

    # Integrate to get velocity and displacement
    vel = np.zeros(len(acc))
    dis = np.zeros(len(acc))

    vel[0] = (acc[0]/2.0) * dt
    for i in range(1, len(acc)):
        vel[i] = vel[i-1] + (((acc[i-1] + acc[i]) / 2.0) * dt)
    dis[0] = (vel[0]/2.0) * dt
    for i in range(1, len(vel)):
        dis[i] = dis[i-1] + (((vel[i-1] + vel[i]) / 2.0) * dt)

    if ordern == 10:
        p = polimod(times, dis, 10, 1)
        pd = [p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1], p[0], 0.0, 0.0]
        pv = [10*p[8], 9*p[7], 8*p[6], 7*p[5], 6*p[4], 5*p[3], 4*p[2],
              3*p[1], 2*p[0], 0.0]
        pa = [10*9*p[8], 9*8*p[7], 8*7*p[6], 7*6*p[5], 6*5*p[4],
              5*4*p[3], 4*3*p[2], 3*2*p[1], 2*1*p[0]]
    elif ordern == 5:
        p = polimod(times, dis, 5, 1)
        pd = [p[3], p[2], p[1], p[0], 0.0, 0.0]
        pv = [5*p[3], 4*p[2], 3*p[1], 2*p[0], 0.0]
        pa = [5*4*p[3], 4*3*p[2], 3*2*p[1], 2*1*p[0]]
    elif ordern == 3:
        p = polimod(times, dis, 3, 1)
        pd = [p[1], p[0], 0.0, 0.0]
        pv = [3*p[1], 2*p[0], 0.0]
        pa = [3*2*p[1], 2*1*p[0]]
    else:
        print("ERROR: Baseline function use order 3, 5, or 10!")
        sys.exit(-1)

    # Evalutate polynomial correction at each time step
    dcor = np.polyval(pd, times)
    vcor = np.polyval(pv, times)
    acor = np.polyval(pa, times)

    # Calculate corrected timeseries
    dmod = dis - dcor
    vmod = vel - vcor
    amod = acc - acor

    amod = amod / gscale

    return times, amod, vmod, dmod
# end of baseline_function

import time

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy.interpolate import splev, splrep
from numba import jit


gmx = '/usr/local/gromacs/bin/gmx'

@jit(nopython=True)
def autocorr_func(time_arr, data_arr):
    rng = np.arange(0, 200, time_arr[1] - time_arr[0])
    val = np.zeros(len(rng))
    N = len(time_arr)//2

    for dt in range(len(rng)):
        c = 0
        for t in range(N):
            c += data_arr[t]*data_arr[t+dt]/N
        val[dt] = c
    return rng, val


def read_xvg(fname):
    if not '.xvg' in fname:
        fname += '.xvg'
    with open(fname, 'r') as f:
        lines = f.readlines()
    lines = [line.split() for line in lines if not (line[0] in ['#', '@', '&'])]
    return np.array(lines).astype('float')

def clear():
    for fname in ['enecorr.xvg', 'energy.xvg', 'evisco.xvg', 'eviscoi.xvg', '#enecorr.xvg.1#', 'vis.xvg']:
        if os.path.exists(fname):
            os.system(f'rm "{fname}"')

def edr(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    os.chdir('/home/student/bmim/analysis/__tmp')
    clear()
    os.system(f'echo 216|{gmx} energy -f /home/student/bmim/anneal/{T}.{tmc}/nvt.edr -vis /home/student/bmim/analysis/__tmp/vis.xvg')
    if pars == 'stress_tensor':
        # time(ps), Pxx(bar), Pxy, Pxz, Pyx, Pyy, Pyz, Pzx, Pzy, Pzz, Temperature(K), Bond(kJ/mol), Pressure(bar)
        lines = read_xvg(f'energy.xvg')
    if pars == 'shear_viscosity':
        lines = read_xvg('#enecorr.xvg.1#')
    if pars == 'bulk_viscosity':
        lines = read_xvg('enecorr.xvg')
    clear()
    return lines

def average(analysis_type):
    # analysis_type = [stress_tensor, bulk_viscosity, shear_viscosity]
    if 'viscosity' in analysis_type:
        for T in [160, 180]:
            cnt = 0
            rng, val = [], []
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/{analysis_type}/{T}.{tmc}.npy'):
                    cnt += 1
                    data = np.load(f'/home/student/bmim/analysis/{analysis_type}/{T}.{tmc}.npy')[0].transpose()
                    if len(rng) == 0:
                        rng, val = data[0], data[1:][0]
                    else:
                        val += data[1:][0]
            val = val/cnt
            print(f'{analysis_type} = {round(sum(val[100:500]/400))} П,   {T} K')
            # rng: time in ps
            # val: viscosity
    if 'stress_tensor' in analysis_type:
        for T in [160, 180]:
            cnt = 0
            rng, P_ij = [], []
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/{analysis_type}/{T}.{tmc}.npy'):
                    cnt += 1
                    data = np.load(f'/home/student/bmim/analysis/{analysis_type}/{T}.{tmc}.npy')[0].transpose()[0:20]
                    if P_ij == []:
                        P_ij = data[1]
                    else:
                        P_ij = np.append(P_ij, data[1])
            rng = np.arange(0, len(P_ij+10)*0.1, 0.1)[0:len(P_ij)]
            spl = splrep(rng, P_ij)
            rng = np.arange(min(rng), max(rng), 0.01)
            P_ij = splev(rng, spl)

            #rng, P_ij: time(ps) and ij-component of the stress tensor(bar)
            #plt.plot(rng, P_ij)
            autocorr_rng, autocorr_val = autocorr_func(time_arr=rng, data_arr=P_ij)
            print(autocorr_rng)
            print(autocorr_val)
            plt.plot(autocorr_rng, autocorr_val, label='Stress tensor correlation function')
            plt.show()




import os, sys, time
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.fft import fft, fftfreq
from __additional import smooth, interp
from scipy.interpolate import splev, splrep
from scipy import linalg as LA
from numba import jit

gmx = '/usr/local/gromacs/bin/gmx'

def read_xvg(fname):
    if not '.xvg' in fname:
        fname += '.xvg'
    with open(fname, 'r') as f:
        lines = f.readlines()
    lines = [line.split() for line in lines if not (line[0] in ['#', '@', '&'])]
    return np.array(lines).astype('float')

def TMP_NO_vector(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    N_ndx = traj.topology.select('resname TMP and name N09')[0]
    O_ndx = traj.topology.select('resname TMP and name O0A')[0]
    return traj.xyz[:, [N_ndx, O_ndx]]

def TMP_rotacf(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    with open('/home/student/bmim/analysis/rotacf.ndx', 'w') as f:
        f.write('[ rotacf ]\n 20710 20711\n')
    os.system(f'{gmx} rotacf -f /home/student/bmim/anneal/{T}.{tmc}/nvt.xtc -s /home/student/bmim/anneal/{T}.{tmc}/nvt.tpr '
              f'-n /home/student/bmim/analysis/rotacf.ndx  -o /home/student/bmim/analysis/TMP_L/NO_rotacf/{T}.{tmc}.xvg -d -b 1000')
    return [read_xvg(f'/home/student/bmim/analysis/TMP_L/NO_rotacf/{T}.{tmc}.xvg')]


def van_hoff(vecs, dts):
    N = np.shape(vecs)[0]
    # нормировка
    for i in range(N):
        vecs[i] = vecs[i] / LA.norm(vecs[i])

    @jit(nopython=True)
    def delta_ang(dt):
        hist = np.zeros(181)
        for t in range(N//2):
            prod = vecs[t].dot(vecs[t+dt])
            if prod > 1:
                prod = 1
            if prod < -1:
                prod = -1
            delta = round(np.arccos(prod)/3.14*180)
            hist[delta] = hist[delta] + 1
        return hist

    # изменение угла
    hists = []
    for dt in dts:
        hists.append(delta_ang(dt))
    return hists


def analysis():
    for T in [160, 180, 190, 200, 220]:
        for tmc in range(5, 400, 5):
            if os.path.exists(f'/home/student/bmim/analysis/TMP_L/NO_vector/{T}.{tmc}.npy') and \
                    (not os.path.exists(f'/home/student/bmim/analysis/TMP_L/NO&fit/{T}.{tmc}.npy')):
                print(f'TMP_L analysis {T}.{tmc}')
                ### фитирование NO, fit: исходные данные и фит. Оба размерности (100001, 3)
                NO_arr = np.load(f'/home/student/bmim/analysis/TMP_L/NO_vector/{T}.{tmc}.npy')
                NO = (NO_arr[:, 1] - NO_arr[:, 0]) # [0:50000] ################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                fit = []
                for coord_i in range(3):
                    spl = splrep(np.arange(len(NO)), NO[:, coord_i], s=13.5, k=5)
                    fit.append(splev(np.arange(len(NO)), spl))
                fit = np.array(fit).transpose()

                ### van_hoff
                DL_dts = [1, 10, 100, 1000, 10000, 50000]
                SL_dts = [1, 10, 100, 1000, 10000, 50000]
                DL_van_hoff = van_hoff(NO, DL_dts)
                SL_van_hoff = van_hoff(fit, SL_dts)

                DL, SL = 0, 0
                FFT = []
                for coord_i in range(3):
                    x, y = np.arange(len(NO)), NO[:, coord_i]
                    y_fit = fit[:, coord_i]

                    ### DL
                    DL += sum((y - y_fit)**2)
                    """plt.plot(x, y, label='orginal', linewidth=1)
                    plt.plot(x, y_fit, label='interpolated', linewidth=1)
                    plt.legend(), plt.show(), time.sleep(100000)"""

                    ### SL
                    def f(x):
                        return splev(x, spl)
                    from scipy.misc import derivative
                    der = np.array([derivative(f, i, dx=1) for i in x])
                    SL += sum([(y_fit[i] - y_fit[i-1])**2 for i in range(len(y_fit))])

                    ### FFT
                    yf = np.absolute(fft(y))
                    xf = fftfreq(len(yf), 10000 / len(y))[:len(yf)//2]
                    yf = 1/len(yf)*yf[:len(yf)//2]
                    if FFT == []:
                        FFT = yf
                    else:
                        FFT += yf

                np.save(f'/home/student/bmim/analysis/TMP_L/NO_DL/{T}.{tmc}', [DL])
                np.save(f'/home/student/bmim/analysis/TMP_L/NO_DL_van_hoff/{T}.{tmc}', [DL_dts, DL_van_hoff])
                np.save(f'/home/student/bmim/analysis/TMP_L/NO_SL/{T}.{tmc}', [SL])
                np.save(f'/home/student/bmim/analysis/TMP_L/NO_SL_van_hoff/{T}.{tmc}', [SL_dts, SL_van_hoff])
                np.save(f'/home/student/bmim/analysis/TMP_L/NO_FFT/{T}.{tmc}.npy', [xf, FFT])
                np.save(f'/home/student/bmim/analysis/TMP_L/NO&fit/{T}.{tmc}.npy', [NO, fit])


def average(analysis_type):
    if analysis_type == 'DL_SL':
        for T in [160, 180, 190, 200, 220]:
            DL, SL, cnt = 0, 0, 0
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/TMP_L/NO_DL/{T}.{tmc}.npy'):
                    DL += np.load(f'/home/student/bmim/analysis/TMP_L/NO_DL/{T}.{tmc}.npy')[0]
                    SL += np.load(f'/home/student/bmim/analysis/TMP_L/NO_SL/{T}.{tmc}.npy')[0]
                    cnt += 1
            DL = DL / cnt
            SL = SL / cnt
            print(f'DL = {DL}, SL = {SL} || {T}K')

    if analysis_type == 'NO_FFT':
        for T in [160, 180]:
            cnt = 0
            rng, FFT = [], []
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/TMP_L/NO_FFT/{T}.{tmc}.npy'):
                    if FFT == []:
                        rng, FFT = np.load(f'/home/student/bmim/analysis/TMP_L/NO_FFT/{T}.{tmc}.npy')
                    else:
                        FFT = FFT + np.load(f'/home/student/bmim/analysis/TMP_L/NO_FFT/{T}.{tmc}.npy')[1]
                    cnt += 1
            FFT = FFT/cnt

            plt.plot(rng, FFT, label=str(T))
        plt.ylim([0, 0.001])
        plt.legend(), plt.show()

    if 'van_hoff' in analysis_type: ### NO_DL_van_hoff or NO_SL_van_hoff
        for T in [160, 180]:
            cnt = 0
            van_hoff = []
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/TMP_L/{analysis_type}/{T}.{tmc}.npy'):
                    dts = np.load(f'/home/student/bmim/analysis/TMP_L/{analysis_type}/{T}.{tmc}.npy', allow_pickle=True)[0]
                    val = np.load(f'/home/student/bmim/analysis/TMP_L/{analysis_type}/{T}.{tmc}.npy', allow_pickle=True)[1]
                    val = np.array([np.array(val[i]) for i in range(len(val))])
                    if van_hoff == []:
                        van_hoff = val
                    else:
                        van_hoff += val
                    cnt += 1
            van_hoff = van_hoff/cnt

            for i in range(len(dts)):
                plt.plot(np.arange(181), van_hoff[i], label=str(dts[i]))
            plt.title(f'{analysis_type} {T}')
            plt.legend(), plt.show()

    if analysis_type == 'NO_rotacf':
        tay = {160:100000, 180:10000}  # подгоняем время корреляции
        for T in [160, 180]:
            cnt = 0
            rng, rotacf = [], []
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/TMP_L/NO_rotacf/{T}.{tmc}.npy'):
                    cnt += 1
                    data = np.load(f'/home/student/bmim/analysis/TMP_L/NO_rotacf/{T}.{tmc}.npy')[0].transpose()
                    if rng == []:
                        rng, rotacf = data
                    else:
                        rotacf += data[1]
            rotacf = rotacf/cnt

            spl = splrep(rng, rotacf, s=0.5, k=5)
            rotacf = splev(rng, spl)

            plt.plot(rng, rotacf, label=str(T))
            plt.plot(rng, rotacf[0]*np.exp(-rng/tay[T]), label=f' exp, tau={tay[T]}ps, T={T}K')
        plt.xlabel('t, ps'), plt.ylabel('rotacf')
        plt.legend(), plt.show()
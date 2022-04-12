import os, sys, time
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.fft import fft, fftfreq
from __additional import smooth, interp
from scipy.interpolate import splev, splrep
from TMP_L import read_xvg, van_hoff
import mdtraj as md
from scipy.misc import derivative
from numba import jit


def clothest_BMI(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    if pars == 'clothest_BMI_chains':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C0N')
        distances = md.compute_distances(traj, pairs)
        clothest = np.argsort(sum(distances)/len(traj))[0:5]
        BMI_chains_ndx = traj.topology.select('resname BMI and name C0N')
        return traj.xyz[:, BMI_chains_ndx[clothest]]
    if pars == 'clothest_BMI_rings':
        pairs = traj.topology.select_pairs('resname TMP and name N09', 'resname BMI and name C01')
        distances = md.compute_distances(traj, pairs)
        clothest = np.argsort(sum(distances)/len(traj))[0:5]
        BMI_rings_ndx = traj.topology.select('resname BMI and name C01')
        return traj.xyz[:, BMI_rings_ndx[clothest]]

def stb(y_val):
    if max(y_val) - min(y_val) > 3:
        mean = (max(y_val) + min(y_val)) /2
        y_val[np.where(y_val < mean)[0]] += np.ones(len(np.where(y_val < mean)[0]))*6
    return y_val

def der(y_val):
    y_der = np.array([y_val[i+1] - y_val[i] for i in range(len(y_val) - 1)])
    return np.append(y_der, [y_der[-1]])

@jit(nopython=True)
def corr_f(TMP_velocity, BMI_velocity, dt):
    val, cnt = 0, 0
    for i in range(np.shape(TMP_velocity)[0]):
        for j in range(np.shape(BMI_velocity)[0]):
            if (j + dt) < np.shape(BMI_velocity)[0] and (j + dt) >= 0:
                val += TMP_velocity[i]*BMI_velocity[j + dt]
                cnt += 1
    return val/(cnt + 0.001)

def analysis(analysis_type): ### analysis_type = 'clothest_BMI_chains' or 'clothest_BMI_rings'
    for tmc in range(0, 400, 5):
        for T in [160, 180]:
            if os.path.exists(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/coords/{T}.{tmc}.npy') \
                    and os.path.exists(f'/home/student/bmim/analysis/TMP_L/NO&fit/{T}.{tmc}.npy'):
                print(f'BMI_TMP_correlation analysis {T}.{tmc}')

                loaded = np.load(f'/home/student/bmim/analysis/TMP_L/NO&fit/{T}.{tmc}.npy')
                np.save(f'/home/student/bmim/analysis/BMI-TMP_correlation/NO&fit/{T}.{tmc}.npy', loaded)
                NO = np.load(f'/home/student/bmim/analysis/BMI-TMP_correlation/NO&fit/{T}.{tmc}.npy')[0]
                NO_fit = np.load(f'/home/student/bmim/analysis/BMI-TMP_correlation/NO&fit/{T}.{tmc}.npy')[1]                            ### фит движения вектора NO, np.shape(NO_fit) = (100001, 3)

                clothest_BMI = np.load(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/coords/{T}.{tmc}.npy')[0:len(NO_fit)] #[0:len(NO_fit)] ### координаты 5 ближайших концов цепей BMI / колец BMI, обрезанные по времени до длины NO_fit
                BMI_fit = np.zeros((len(clothest_BMI), 5, 3))
                for mol_i in range(5):
                    for coord_i in range(3):
                        spl = splrep(np.arange(len(clothest_BMI)), stb(clothest_BMI[:, mol_i, coord_i]), s=9, k=5)
                        BMI_fit[:, mol_i, coord_i] = splev(np.arange(len(clothest_BMI)), spl)
                # NO,               np.shape(NO) = (100001, 3)
                # NO_fit,           np.shape(NO_fit) = (100001, 3)
                # clothest_BMI      np.shape(clothest_BMI) = (100001, 5, 3) - координаты 5 ближайших концов цепей BMI / колец BMI
                # BMI_fit           np.shape(BMI_fit) = (100001, 5, 3)
                # по времени оба фита имеют одинаковую длину
                # нужно вычислить скорость вращения NO  и  скорость движения BMI_fit. Найти их корреляцию


                ### STOCHASTIC
                # нужно вычислить стохастическую скорость вращения NO  и  скорость движения BMI_fit. Найти их корреляцию
                TMP_velocity = np.sqrt(sum([der(NO_fit[:, coord_i])**2 for coord_i in range(3)]))
                BMI_velocity = sum([np.sqrt(sum([der(BMI_fit[:, mol, coord_i])**2 for coord_i in range(3)])) for mol in range(5)])/5

                dts = np.linspace(-100, 100, 100).astype(int)
                corr_fs = []
                for dt in dts:
                    corr_fs.append(corr_f(TMP_velocity, BMI_velocity, dt))
                np.save(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/S_corr_func/{T}.{tmc}.npy', [dts, corr_fs])

                ### DYNAMIC
                # нужно вычислить динамическую скорость вращения NO  и  скорость движения BMI_fit. Найти их корреляцию
                TMP_velocity = np.sqrt(sum([stb(NO[:, coord_i]) ** 2 for coord_i in range(3)]))
                BMI_velocity = sum([np.sqrt(sum([stb(clothest_BMI[:, mol, coord_i]) ** 2 for coord_i in range(3)])) for mol in range(5)]) / 5

                dts = np.linspace(-100, 100, 100).astype(int)
                corr_fs = []
                for dt in dts:
                    corr_fs.append(corr_f(TMP_velocity, BMI_velocity, dt))
                np.save(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/D_corr_func/{T}.{tmc}.npy', [dts, corr_fs])


def average():
    for analysis_type in ['clothest_BMI_rings', 'clothest_BMI_chains', ]:
        for corr_func_type in ['S_corr_func', 'D_corr_func']:
            plot_data = []
            for T in [160, 180]:
                cnt = 0
                dts, corr_fs = [], []
                for tmc in range(5, 400, 5):
                    if os.path.exists(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/{corr_func_type}/{T}.{tmc}.npy'):
                        cnt += 1
                        if corr_fs == []:
                            dts, corr_fs = np.load(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/{corr_func_type}/{T}.{tmc}.npy')
                        else:
                            corr_fs += np.load(f'/home/student/bmim/analysis/BMI-TMP_correlation/{analysis_type}/{corr_func_type}/{T}.{tmc}.npy')[1]
                corr_fs = np.array(corr_fs)/(0.0001 + cnt)

                plot_data.append([dts, corr_fs, f'{analysis_type}, {corr_func_type}, {T}'])
                plt.plot(dts, corr_fs, label=f'{analysis_type}, {corr_func_type}, {T}')
            np.save(f'/home/student/bmim/analysis/BMI-TMP_correlation/plot_data/{analysis_type}_{corr_func_type}.npy', plot_data)
            plt.legend(), plt.show()
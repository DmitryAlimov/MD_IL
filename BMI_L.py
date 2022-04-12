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

gmx = '/usr/local/gromacs/bin/gmx'
mols = np.array([11, 14, 34, 36, 58, 62, 113, 186, 209, 212, 246, 252, 283, 296, 344, 380, 442, 446, 488, 497])
mol_ids = np.sort(np.append(mols*2, mols*2+np.ones(len(mols)))).astype(int) ### 22, 23, 28, 29, ...


def BMI_vector(traj, tmc, T, chunk, stride, data_folder, pars, cnt):
    if pars == 'C0N_N02_vector':
        ndx = np.array(traj.topology.select('resname BMI and name C0N N02'))[mol_ids]
    if pars == 'C01_C04_vector':
        ndx = np.array(traj.topology.select('resname BMI and name C01 C04'))[mol_ids]
    return traj.xyz[:, ndx]


def BMI_rotacf(vector_type):
    for traj in md.iterload('/home/student/bmim/anneal/160.5/nvt.xtc', top='/home/student/bmim/anneal/160.5/nvt.gro', chunk=1):
        if vector_type == 'C01_C04':
            ndx = traj.topology.select('resname BMI and name C01 C04')[mol_ids]
        if vector_type == 'C0N_N02':
            ndx = traj.topology.select('resname BMI and name C0N N02')[mol_ids]
        break

    # количество молекул BMI, для которых вычисляется rotacf = 20
    for tmc in range(5, 400, 5):
        for T in [160, 180]:
            for mol in range(len(mols)):
                with open('/home/student/bmim/analysis/rotacf.ndx', 'w') as f:
                    f.write('[ rotacf ]\n')
                    f.write(f'{ndx[2*mol]+1} {ndx[2*mol+1]+1}\n')
                os.system(f'echo 0|{gmx}  rotacf -f /home/student/bmim/anneal/{T}.{tmc}/nvt.xtc -s /home/student/bmim/anneal/{T}.{tmc}/nvt.tpr -n /home/student/bmim/analysis/rotacf.ndx -o /home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_rotacfs/{T}.{tmc}.{mol}.xvg -d')
                np.save(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_rotacfs/{T}.{tmc}.{mol}.npy',
                        read_xvg(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_rotacfs/{T}.{tmc}.{mol}.xvg'))


def analysis(vector_type): ### vector_type = 'C01_C04' or 'C0N_N02'
    for tmc in range(5, 400, 5):
        for T in [160, 180]:
            if os.path.exists(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_vector/{T}.{tmc}.npy')\
                    and (not os.path.exists(f'/home/student/bmim/analysis/BMI_L/C01_C04/vecs&fit/{T}.{tmc}.npy')):
                arr = np.load(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_vector/{T}.{tmc}.npy')

                vecss, fits = [], []
                DL_dts = [1, 10, 100, 1000, 10000, 50000]
                SL_dts = [1, 10, 100, 1000, 10000, 50000]
                DL_van_hoffs, SL_van_hoffs = [], []
                DLs, SLs = [], []
                FFTs = []

                for mol in range(len(mols)):
                    ### фитирование vecs, fit: исходные данные и фит. Оба размерности (100001, 3)
                    vecs = (arr[:, 2*mol+1] - arr[:, 2*mol]) #[0:5000] ################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    fit = []
                    for coord_i in range(3):
                        spl = splrep(np.arange(len(vecs)), vecs[:, coord_i], s=13.5, k=5)
                        fit.append(splev(np.arange(len(vecs)), spl))
                    fit = np.array(fit).transpose()
                    vecss.append(vecs), fits.append(fit)

                    ### van_hoff
                    DL_van_hoffs.append(van_hoff(vecs, DL_dts))
                    SL_van_hoffs.append(van_hoff(fit, SL_dts))

                    ### DL
                    DLs.append(sum([sum((vecs[:, i] - fit[:, i]) ** 2) for i in range(3)]))
                    """plt.plot(np.arange(len(vecs)), vecs[:, 0], label='orginal', linewidth=1)
                    plt.plot(np.arange(len(vecs)), fit[:, 0], label='interpolated', linewidth=1)
                    plt.legend(), plt.show()"""

                    ### SL
                    SLs.append(sum([sum((fit[i] - fit[i-1])**2) for i in np.arange(1, len(vecs))]))

                np.save(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_DL/{T}.{tmc}', [DLs])
                np.save(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_DL_van_hoff/{T}.{tmc}', [DL_dts, DL_van_hoffs])
                np.save(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_SL/{T}.{tmc}', [SLs])
                np.save(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_SL_van_hoff/{T}.{tmc}', [SL_dts, SL_van_hoffs])
                np.save(f'/home/student/bmim/analysis/BMI_L/{vector_type}/vecs&fit/{T}.{tmc}.npy', [vecss, fits])


def average(vector_type, analysis_type):
    if analysis_type == 'DL_SL':
        for T in [160, 180]:
            DL, SL, cnt = 0, 0, 0
            for tmc in range(0, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_DL/{T}.{tmc}.npy'):
                    loaded = np.load(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_DL/{T}.{tmc}.npy')[0]  ### лист из 20 значений DL для 20 выбранных молекул BMI для данной траектории T/tmc
                    DL += sum(loaded)/len(loaded)
                    loaded = np.load(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_SL/{T}.{tmc}.npy')[0]  ### лист из 20 значений SL для 20 выбранных молекул BMI для данной траектории T/tmc
                    SL += sum(loaded) / len(loaded)
                    cnt += 1
            DL, Sl = DL/cnt, SL/cnt
            print(f'DL = {DL}, SL = {SL} || {T}K')
    if 'van_hoff' in analysis_type: ### DL_van_hoff or SL_van_hoff
        for T in [160, 180]:
            cnt = 0
            van_hoff = []
            for tmc in range(5, 400, 5):
                if os.path.exists(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_{analysis_type}/{T}.{tmc}.npy'):
                    dts = np.load(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_{analysis_type}/{T}.{tmc}.npy', allow_pickle=True)[0]
                    val = np.load(f'/home/student/bmim/analysis/BMI_L/{vector_type}/{vector_type}_{analysis_type}/{T}.{tmc}.npy', allow_pickle=True)[1]
                    val = np.array([[val[i][j] for j in range(np.shape(val)[1])] for i in range(np.shape(val)[0])])
                    if van_hoff == []:
                        van_hoff = val
                    else:
                        van_hoff += val
                    cnt += 1
            van_hoff = van_hoff/cnt
            # np.shape(val) =  (20, 6, 181)         20 молекул BMI, 6 dts, 181 градусов

            for i in range(len(dts)):
                plt.plot(np.arange(181), sum(van_hoff)[i]/len(van_hoff), label=str(dts[i]))
            plt.title(f'{analysis_type} {T}K')
            plt.legend(), plt.show()
    if analysis_type == 'rotacf':
        val, cnt = [], 0
        for T in [160, 180]:
            for tmc in range(5, 400, 5):
                for mol in range(20):
                    if os.path.exists(f'/home/student/bmim/analysis/BMI_L/C01_C04/C01_C04_rotacfs/{T}.{tmc}.{mol}.npy'):
                        loaded = np.load(f'/home/student/bmim/analysis/BMI_L/C01_C04/C01_C04_rotacfs/{T}.{tmc}.{mol}.npy').transpose()
                        rng = loaded[0]
                        if val == []:
                            val = loaded[1]
                        else:
                            val += loaded[1]
                        cnt += 1
            val = val/cnt
            plt.plot(rng, val, label=str(T))
        plt.show()
import os, sys, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Process, Manager
man = Manager()
run_bool = man.list([0 for _ in range(100)])


def nm(tmc, T, ext):
    return f'/home/student/bmim/anneal/{T}.{tmc}/nvt.{ext}'

def calc(func, tmc, T, chunk=1000, stride=1, limit=1000000, data_folder='', pars=[], cnt=0):   #принимает (tmc, T) конкретной траектории, считает для нее func(traj) и сохраняет в data_folder
    run_bool[cnt+1] = 1
    data = []
    iter = 1
    for traj in md.iterload(nm(tmc, T, 'xtc'), top=nm(tmc, T, 'gro'), chunk=chunk, stride=stride):
        val = func(traj, tmc, T, chunk, stride, data_folder, pars, cnt)
        if np.size(data) == 0:
            data = val
        else:
            data = np.append(data, val, axis=0)
        iter += 1
        if iter*chunk*stride > limit:
            break
    np.save(f'{data_folder}/{T}.{tmc}.npy', data)
    run_bool[cnt+1] = 0

def MultiThreading(func, tmcs, Ts=[160, 180], chunk=1000, stride=1, limit=10000000, data_folder='', pars=[], nt=4, recalc=False):   #перебирает значения из tmcs и Ts и запускает для них calc(), nt: number of threads
    # nt: количество процессов, recalc: перезапись посчитанных распределений в папке  /data_folder/T.tmc.npy
    p_arr, cnt = [], 0

    for tmc in tmcs:
        for T in Ts:
            if cnt >= nt:
                while sum(run_bool) > 0:
                    time.sleep(1)
                for p in p_arr:
                    p.join()
                p_arr, cnt = [], 0

            run_bool[cnt + 1] = 1
            if os.path.exists(nm(tmc, T, 'gro')) and \
                    (not ((not recalc) and os.path.exists(f'{data_folder}/{T}.{tmc}.npy'))):
                print(f'{tmc} ns    {T}K     | {datetime.fromtimestamp(time.time()).strftime("%m/%d/%Y, %H:%M:%S")}')
                p_arr.append(Process(target=calc, args=[func, tmc, T, chunk, stride, limit, data_folder, pars, cnt]))
                p_arr[cnt].start()
                cnt += 1
    while sum(run_bool) > 0:
        time.sleep(1)
    for p in p_arr:
        p.join()
    print('terminated successfully')

def OneThread(func, tmcs, Ts=[160, 180], chunk=1000, stride=1, limit=10000000, data_folder='', pars=[], nt=1, recalc=False):   #перебирает значения из tmcs и Ts и запускает для них calc(), nt: number of threads
    for tmc in tmcs:
        for T in Ts:
            if os.path.exists(nm(tmc, T, 'gro')) and \
                    (not ((not recalc) and os.path.exists(f'{data_folder}/{T}.{tmc}.npy'))):
                calc(func, tmc, T, chunk=chunk, stride=stride, limit=limit, data_folder=data_folder, pars=pars)







### ============================================================================================================== ###
"""from rdf import rdf
for pars in ['BMI.C0N-BF.F', 'BMI.N02-BF.B', 'BMI.N02-BF.F', 'BMI.C01-BF.B', 'BMI.C01-BF.F', 'BF.B-BF.B', 'BMI.C0N-BMI.C0N', 'BMI.C01-BMI.C01', 'BMI.C0N-BF.B',]:
    print(pars)
    OneThread(func=rdf, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190, 200, 220], data_folder=f'/home/student/bmim/analysis/rdf/{pars}/', pars=pars, nt=1, chunk=25, limit=10000)
    
for pars in ['TMP.N-BMI.C0N', 'TMP.N-BMI.C06', 'TMP.N-BMI.N00', 'TMP.N-BMI.C03_C04', 'TMP.N-BMI', 'TMP.N-BF', 'TMP.N-BF.B', 'TMP.N-BF.F', 'TMP.N-BMI.C0J', 'TMP.N-BMI.COG']:
    print(pars)
    OneThread(func=rdf, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190, 200, 220], data_folder=f'/home/student/bmim/analysis/rdf/{pars}/', pars=pars, nt=1, chunk=1000)"""

from TMP_L import TMP_NO_vector, TMP_rotacf
#OneThread(func=TMP_NO_vector, tmcs=np.arange(5, 400, 5), Ts=[190], data_folder='/home/student/bmim/analysis/TMP_L/NO_vector/', nt=1)
#OneThread(func=TMP_rotacf, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190], data_folder='/home/student/bmim/analysis/TMP_L/NO_rotacf/', chunk=10, limit=11, nt=1)

from BMI_L import BMI_vector, BMI_rotacf
#OneThread(func=BMI_vector, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190], data_folder='/home/student/bmim/analysis/BMI_L/C01_C04/C01_C04_vector', pars='C01_C04_vector')
#OneThread(func=BMI_vector, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190], data_folder='/home/student/bmim/analysis/BMI_L/C0N_N02/C0N_N02_vector', pars='C0N_N02_vector')

from edr import edr
#MultiThreading(func=edr, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190], chunk=10, limit=11, pars='shear_viscosity', data_folder='/home/student/bmim/analysis/shear_viscosity/', nt=1)
# time(ps), Pxx(bar), Pxy, Pxz, Pyx, Pyy, Pyz, Pzx, Pzy, Pzz, Temperature(K), Bond(kJ/mol), Pressure(bar)
#MultiThreading(func=edr, tmcs=np.arange(20, 400, 5), Ts=[160, 180], chunk=10, limit=11, pars='stress_tensor', data_folder='/home/student/bmim/analysis/stress_tensor/', nt=1)

from clusters import BMI_chains
#MultiThreading(func=BMI_chains, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190], data_folder='/home/student/bmim/analysis/clusters/BMI_chains/', nt=4)

from BMI_TMP_correlation import clothest_BMI
OneThread(func=clothest_BMI, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190, 200, 220], pars='clothest_BMI_chains', chunk=5000, data_folder='/home/student/bmim/analysis/BMI-TMP_correlation/clothest_BMI_chains/coords/', nt=1)
#OneThread(func=clothest_BMI, tmcs=np.arange(5, 400, 5), Ts=[160, 180, 190, 200, 220], pars='clothest_BMI_rings', chunk=5000, data_folder='/home/student/bmim/analysis/BMI-TMP_correlation/clothest_BMI_rings/coords/', nt=1)


